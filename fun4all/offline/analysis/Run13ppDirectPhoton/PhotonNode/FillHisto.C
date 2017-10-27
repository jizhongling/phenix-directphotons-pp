#include "FillHisto.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsPhoton.h"

#include "EmcLocalRecalibrator.h"
#include "EmcLocalRecalibratorSasha.h"
#include "PhotonContainerClone.h"

#include <PhotonContainer.h>
#include <Photon.h>
#include <PhotonERT.h>
#include <SpinPattern.h>
#include <RunHeader.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TGraphErrors.h>

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <iterator>

using namespace std;

const double PI = TMath::Pi();

const float eMin = 0.3;
const float AsymCut = 0.8;

FillHisto::FillHisto(const string &name, const char *filename) :
  SubsysReco(name),
  hm(NULL),
  emcrecalib(NULL),
  emcrecalib_sasha(NULL),
  h_events(NULL),
  h3_tof(NULL),
  h3_tof_raw(NULL),
  h3_inv_mass_pi0calib(NULL),
  h3_inv_mass_pi0calib_raw(NULL),
  h3_trig(NULL),
  h3_trig_pion(NULL),
  hn_1photon(NULL),
  hn_2photon(NULL),
  hn_pion(NULL),
  hn_asym(NULL),
  hn_minv(NULL)
{
  datatype = ERT;

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 9999;

  // initialize crossing shift and spin pattern
  crossing_shift = 0;
  for(int i=0; i<120; i++)
  {
    spinpattern_blue[i] = 0;
    spinpattern_yellow[i] = 0;
  }
}

FillHisto::~FillHisto()
{
}

int FillHisto::Init(PHCompositeNode *topNode)
{
  // book histograms and graphs
  BookHistograms();

  // EMCal recalibration class
  emcrecalib = new EmcLocalRecalibrator();
  if(!emcrecalib)
  {
    cerr << "No emcrecalib" << endl;
    exit(1);
  }

  emcrecalib_sasha = new EmcLocalRecalibratorSasha();
  if(!emcrecalib_sasha)
  {
    cerr << "No emcrecalib_sasha" << endl;
    exit(1);
  }

  // read EMCal recalibration file
  EMCRecalibSetup();

  // read warnmap
  //ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  return EVENT_OK;
}

int FillHisto::InitRun(PHCompositeNode *topNode)
{
  SpinPattern *spinpat = findNode::getClass<SpinPattern>(topNode, "SpinPattern");
  if(!spinpat)
  {
    cerr << "No spinpattern" << endl;
    return ABORTRUN;
  }

  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!runheader)
  {
    cerr << "No runheader" << endl;
    return ABORTRUN;
  }

  // run number, time and crossing shift for each run
  runnumber = runheader->get_RunNumber();
  runtime = difftime( runheader->get_TimeStop(), runheader->get_TimeStart() );
  crossing_shift = spinpat->get_crossing_shift();

  // number of MinBias and spin pattern for each run
  for(int i=0; i<120; i++)
  {
    spinpattern_blue[i] = spinpat->get_spinpattern_blue(i);
    spinpattern_yellow[i] = spinpat->get_spinpattern_yellow(i);
  }

  // load EMCal recalibrations for run and fill
  emcrecalib->ReadEnergyCorrection( spinpat->get_runnumber() );
  emcrecalib->ReadTofCorrection( spinpat->get_fillnumber() );

  return EVENT_OK;
}

int FillHisto::process_event(PHCompositeNode *topNode)
{
  PhotonContainer *photoncont = findNode::getClass<PhotonContainer>(topNode, "PhotonContainer");
  if(!photoncont)
  {
    cerr << "No photoncont" << endl;
    return DISCARDEVENT;
  }

  float bbc_z = photoncont->get_bbc_z();

  /* Check trigger */
  if( datatype == ERT )
  {
    if( photoncont->get_bbcnarrow_live() && abs(bbc_z) < 10. )
    {
      if( photoncont->get_ert_a_scaled() )
        h_events->Fill("ert_a", 1.);
      if( photoncont->get_ert_b_scaled() )
        h_events->Fill("ert_b", 1.);
      if( photoncont->get_ert_c_scaled() )
        h_events->Fill("ert_c", 1.);
    }
  }

  else if( datatype == MB )
  {
    if( photoncont->get_bbcnovtx_scaled() )
    {
      h_events->Fill("bbc_novtx", 1.);
      if( abs(bbc_z) < 40. )
        h_events->Fill("bbc_40cm", 1.);
      if( abs(bbc_z) < 30. )
        h_events->Fill("bbc_30cm", 1.);
      if( abs(bbc_z) < 15. )
        h_events->Fill("bbc_15cm", 1.);
      if( abs(bbc_z) < 10. )
        h_events->Fill("bbc_10cm", 1.);
    }

    h_events->Fill("bbc", 1.);
    if( photoncont->get_ert_a_live() )
      h_events->Fill("bbc_ert_a", 1.);
    if( photoncont->get_ert_b_live() )
      h_events->Fill("bbc_ert_b", 1.);
    if( photoncont->get_ert_c_live() )
      h_events->Fill("bbc_ert_c", 1.);
  }

  // Run local recalibration of EMCal cluster data
  PhotonContainerClone *photoncont_raw = new PhotonContainerClone(photoncont);
  //emcrecalib->ApplyClusterCorrection( photoncont );
  emcrecalib_sasha->ApplyClusterCorrection( runnumber, photoncont );

  // Store TOF information for cluster as calibration check
  FillClusterTofSpectrum( photoncont_raw, "raw" );
  FillClusterTofSpectrum( photoncont );

  // Analyze pi0s events for crosscheck
  FillPi0InvariantMass( photoncont_raw, "raw" );
  FillPi0InvariantMass( photoncont );

  // Count events to calculate trigger efficiency
  FillTriggerEfficiency( photoncont );

  // Analyze single photon event
  FillSinglePhotonSpectrum( photoncont );

  // Analyze photon pair for direct photon event
  FillTwoPhotonSpectrum( photoncont );

  // Analyze photon pair for pi0 event
  FillPi0Spectrum( photoncont );

  delete photoncont_raw;

  return EVENT_OK;
}

int FillHisto::FillClusterTofSpectrum( const PhotonContainer *photoncont, const string &quali )
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  for( unsigned i = 0; i < nphotons; i++ )
  {
    Photon *photon = photoncont->GetPhoton(i);

    if( GetStatus(photon) == 0 )
    {
      int sector = anatools::GetSector(photon);
      double tof = photon->get_tof() - bbc_t0;
      double pT = anatools::Get_pT(photon);

      if ( quali == "raw" )
        h3_tof_raw->Fill( sector, pT, tof );
      else
        h3_tof->Fill( sector, pT, tof );
    }
  }

  return EVENT_OK;
}

int FillHisto::FillPi0InvariantMass( const PhotonContainer *photoncont, const string &quali )
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<nphotons; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        Photon *photon1 = photoncont->GetPhoton(i);
        Photon *photon2 = photoncont->GetPhoton(j);
        if( GetStatus(photon1) == 0
            && GetStatus(photon2) == 0
            && TestPhoton(photon1, bbc_t0)
            && TestPhoton(photon2, bbc_t0)
            && anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
        {
          int sector1 = anatools::GetSector(photon1);
          int sector2 = anatools::GetSector(photon2);
          if( sector1 != sector2 ) continue;
          //if( !anatools::SectorCheck(sector1,sector2) ) continue;

          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);

          if ( quali == "raw" )
            h3_inv_mass_pi0calib_raw->Fill(sector1, tot_pT, minv);
          else
            h3_inv_mass_pi0calib->Fill(sector1, tot_pT, minv);
        }
      }
  }

  return EVENT_OK;
}

int FillHisto::FillTriggerEfficiency( const PhotonContainer *photoncont )
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  /* Fire ERT on arm 0 (west) or 1 (east) */
  bool FireERT[2] = {};

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon = photoncont->GetPhoton(i);
    int arm = anatools::GetSector(photon) / 4;
    if( photon->get_trg1() ||
        photon->get_trg2() ||
        photon->get_trg3() )
      FireERT[arm] = true;
  }

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon1 = photoncont->GetPhoton(i);
    int arm = anatools::GetSector(photon1) / 4;
    int sector = anatools::GetSector( photon1 );
    double photon_pT = anatools::Get_pT( photon1 );

    /* Require the other arm to be fired */
    int oarm = ( arm==0 ? 1 : 0 );
    if( datatype == ERT && !FireERT[oarm] ) continue;

    if( GetStatus(photon1) == 0 &&
        TestPhoton(photon1, bbc_t0) )
    {
      v_used.push_back(i);

      h3_trig->Fill(photon_pT, sector, "all", 1.);
      if( photon1->get_trg1() )
        h3_trig->Fill(photon_pT, sector, "ERT4x4a", 1.);
      if( photon1->get_trg2() )
        h3_trig->Fill(photon_pT, sector, "ERT4x4b", 1.);
      if( photon1->get_trg3() )
        h3_trig->Fill(photon_pT, sector, "ERT4x4c", 1.);

      for(unsigned j=0; j<nphotons; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          Photon *photon2 = photoncont->GetPhoton(j);
          if( GetStatus(photon2) == 0 &&
              TestPhoton(photon2, bbc_t0) )
          {
            int sector2 = anatools::GetSector( photon2 );
            if( !anatools::SectorCheck(sector,sector2) )
              continue;

            double minv = anatools::GetInvMass(photon1, photon2);
            if( minv < 0.112 || minv > 0.162 )
              continue;

            double tot_pT = anatools::GetTot_pT(photon1, photon2);

            h3_trig_pion->Fill(tot_pT, sector, "all", 1.);
            if( photon1->get_trg1() ||
                photon2->get_trg1() )
              h3_trig_pion->Fill(tot_pT, sector, "ERT4x4a", 1.);
            if( photon1->get_trg2() ||
                photon2->get_trg2() )
              h3_trig_pion->Fill(tot_pT, sector, "ERT4x4b", 1.);
            if( photon1->get_trg3() ||
                photon2->get_trg3() )
              h3_trig_pion->Fill(tot_pT, sector, "ERT4x4c", 1.);
          }
        }
    }
  }

  return EVENT_OK;
}

int FillHisto::FillSinglePhotonSpectrum( const PhotonContainer *photoncont )
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 10. ) return DISCARDEVENT;

  /* Check trigger */
  if( datatype == ERT && !photoncont->get_ert_c_scaled() )
    return DISCARDEVENT;
  else if( datatype == MB && !photoncont->get_bbcnarrow_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  int pattern = GetPattern(photoncont);

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon = photoncont->GetPhoton(i);
    if( GetStatus(photon) == 0
        && TestPhoton(photon, bbc_t0) )
    {
      int sector = anatools::GetSector(photon);

      TLorentzVector pE = anatools::Get_pE(photon);
      double px = pE.Px();
      double py = pE.Py();
      double pz = pE.Pz();
      double pT = pE.Pt();
      double mom = pE.P(); 

      double eta = mom > 0. ? atan(pz/mom) : 9999.;
      double phi = px > 0. ? atan(py/px) : PI+atan(py/px);

      double fill_hn_1photon[] = {sector, pT, pattern, eta, phi};
      hn_1photon->Fill(fill_hn_1photon);
    }
  }

  return EVENT_OK;
}

int FillHisto::FillTwoPhotonSpectrum(const PhotonContainer *photoncont)
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  /* Check trigger */
  if( datatype == ERT && !photoncont->get_ert_c_scaled() )
    return DISCARDEVENT;
  else if( datatype == MB && !photoncont->get_bbcnarrow_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  int pattern = GetPattern(photoncont);

  for(unsigned i=0; i<nphotons; i++)
    for(unsigned j=0; j<nphotons; j++)
      if(j != i)
      {
        Photon *photon1 = photoncont->GetPhoton(i);
        Photon *photon2 = photoncont->GetPhoton(j);

        if( GetStatus(photon1) == 0
            && GetStatus(photon2) == 0
            && TestPhoton(photon1, bbc_t0)
            && TestPhoton(photon2, bbc_t0)
            && anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
        {
          int sector = anatools::GetSector(photon1);
          double pT = anatools::Get_pT(photon1);
          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);
          //double theta_cv = photon1->get_theta_cv();

          double fill_hn_2photon[] = {sector, pT, tot_pT, minv, pattern};
          hn_2photon->Fill(fill_hn_2photon);
        }
      }

  return EVENT_OK;
}

int FillHisto::FillPi0Spectrum(const PhotonContainer *photoncont)
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( abs(bbc_z) > 10. ) return DISCARDEVENT;

  /* Check trigger */
  if( datatype == ERT )
  {
    if( !photoncont->get_ert_c_scaled() || !photoncont->get_bbcnarrow_live() )
      return DISCARDEVENT;
  }
  else if( datatype == MB )
  {
    if( !photoncont->get_bbcnovtx_scaled() )
      return DISCARDEVENT;
  }

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<nphotons; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        Photon *photon1 = photoncont->GetPhoton(i);
        Photon *photon2 = photoncont->GetPhoton(j);
        double asym = anatools::GetAsymmetry_E(photon1, photon2);
        if( GetStatus(photon1) == 0 &&
            GetStatus(photon2) == 0 &&
            asym < AsymCut )
        {
          int sector1 = anatools::GetSector(photon1);
          int sector2 = anatools::GetSector(photon2);
          //if( sector1 != sector2 ) continue;
          if( !anatools::SectorCheck(sector1,sector2) ) continue;

          int sector = sector1;
          bool trig = photon1->get_trg3();
          if( photon2->get_E() > photon1->get_E() )
          {
            sector = sector2;
            trig = photon2->get_trg3();
          }
          if( datatype == ERT && !trig) continue;

          TLorentzVector pE1 = anatools::Get_pE(photon1);
          TLorentzVector pE2 = anatools::Get_pE(photon2);
          TLorentzVector tot_pE =  pE1 + pE2;
          double tot_pT = tot_pE.Pt();
          double minv = tot_pE.M();

          double fill_hn_pion_0[] = {sector, tot_pT, minv, 0.};
          hn_pion->Fill(fill_hn_pion_0);
          if( TestPhoton(photon1, bbc_t0) &&
              TestPhoton(photon2, bbc_t0) )
          {
            double fill_hn_pion_1[] = {sector, tot_pT, minv, 1.};
            hn_pion->Fill(fill_hn_pion_1);

            double dx = photon1->get_x() - photon2->get_x();
            double dy = photon1->get_y() - photon2->get_y();
            double dz = photon1->get_z() - photon2->get_z();
            double dR = sqrt(dx*dx + dy*dy + dz*dz);
            double fill_hn_asym[] = {sector, tot_pT, asym, dR, minv};
            hn_asym->Fill(fill_hn_asym);

            double angle = sqrt( 1. - cos( pE1.Angle(pE2.Vect()) ) );
            double fill_hn_minv[] = {sector, tot_pT, pE1.E(), pE2.E(), angle};
            hn_minv->Fill(fill_hn_minv);
          } // TestPhoton
        } // GetStatus and asymmetry cut
      } // second cluster loop
  } // first cluster loop

  return EVENT_OK;
}

int FillHisto::EndRun(const int runnumber)
{
  // Write histograms to file for this run
  // and reset histograms for the next run
  char filename[200];
  if( datatype == ERT )
    sprintf(filename, "histos-ERT/PhotonNode-%d.root", runnumber);
  else if( datatype == MB )
    sprintf(filename, "histos-MB/PhotonNode-%d.root", runnumber);
  hm->dumpHistos(filename);
  hm->Reset();

  return EVENT_OK;
}

int FillHisto::End(PHCompositeNode *topNode)
{
  delete hm;
  delete emcrecalib;
  delete emcrecalib_sasha;

  return EVENT_OK;
}

void FillHisto::SelectMB()
{
  datatype = MB;
  return;
}

void FillHisto::SelectERT()
{
  datatype = ERT;
  return;
}

void FillHisto::BookHistograms()
{
  /* Create HistogramManager */
  hm = new Fun4AllHistoManager("HistogramManager");

  /*
   * pT bins
   */
  const int n_pTbins = 31;
  const double pTbins[n_pTbins+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
    100.0};

  /*
   * phi sector
   */
  const double phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8,
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  /*
   * phi bins
   */
  const int n_phibins = 163;
  double phi_twr[n_phibins+1];
  for(int is=0; is<6; is++)
    for(int it=0; it<19; it++)
      phi_twr[19*is+it] = phi_sec[is] + 0.02 * ( it - 9 );
  for(int is=6; is<8; is++)
    for(int it=0; it<25; it++)
      phi_twr[114+25*(is-6)+it] = phi_sec[is] + 0.016 * ( it - 12 );

  /*
   * Events counter
   */
  if( datatype == ERT )
  {
    h_events = new TH1F("h_events", "Events counter", 3,0.5,3.5);
    h_events->GetXaxis()->SetBinLabel(1, "ert_a");
    h_events->GetXaxis()->SetBinLabel(2, "ert_b");
    h_events->GetXaxis()->SetBinLabel(3, "ert_c");
  }
  else if( datatype == MB )
  {
    h_events = new TH1F("h_events", "Events counter", 9,0.5,9.5);
    h_events->GetXaxis()->SetBinLabel(1, "bbc_novtx");
    h_events->GetXaxis()->SetBinLabel(2, "bbc_40cm");
    h_events->GetXaxis()->SetBinLabel(3, "bbc_30cm");
    h_events->GetXaxis()->SetBinLabel(4, "bbc_15cm");
    h_events->GetXaxis()->SetBinLabel(5, "bbc_10cm");
    h_events->GetXaxis()->SetBinLabel(6, "bbc");
    h_events->GetXaxis()->SetBinLabel(7, "bbc_ert_a");
    h_events->GetXaxis()->SetBinLabel(8, "bbc_ert_b");
    h_events->GetXaxis()->SetBinLabel(9, "bbc_ert_c");
  }
  hm->registerHisto( h_events, 1 );

  /*
   * 3D histogram of sector, pT and TOF to check TOF calibration
   */
  h3_tof = new TH3F("h3_tof", "TOF;EMCal sector;p_{T} [GeV];TOF [ns];", 8,-0.5,7.5, n_pTbins,0.,0., 1001,-100.05,100.05);
  h3_tof->GetYaxis()->Set(n_pTbins, pTbins);
  hm->registerHisto( h3_tof , 1 );

  /*
   * Using TOF from DST file without local recalibration.
   */
  h3_tof_raw = static_cast<TH3*>(h3_tof->Clone("h3_tof_raw"));
  hm->registerHisto( h3_tof_raw , 1 );

  /*
   * 3D histogram storing invariant mass of photon
   * candidate pairs in different sectors and pT bins; photon pair
   * selection particularly strict for asymmetry.
   * Require both photons to be in same sector. Used to check sector-by-
   * sector EMCal energy calibration.
   */
  h3_inv_mass_pi0calib = new TH3F("h3_inv_mass_pi0calib", "Photon pair invariant mass;EMCal sector;p_{T} [GeV];m_{inv} [GeV];", 8,-0.5,7.5, n_pTbins,0.,0., 300,0.,0.3);
  h3_inv_mass_pi0calib->GetYaxis()->Set(n_pTbins, pTbins);
  hm->registerHisto( h3_inv_mass_pi0calib , 1 );

  /*
   * Using energies from DST file without local recalibration.
   */
  h3_inv_mass_pi0calib_raw = static_cast<TH3*>(h3_inv_mass_pi0calib->Clone("h3_inv_mass_pi0calib_raw"));
  hm->registerHisto( h3_inv_mass_pi0calib_raw , 1 );

  /*
   * 3D histogram to count for trigger efficiency
   */
  h3_trig = new TH3F("h3_trig", "number of clusters;p_{T} [GeV];sector;trigger;", n_pTbins,0.,0., 8,-0.5,7.5, 4,0.5,4.5);
  h3_trig->GetXaxis()->Set(n_pTbins, pTbins);
  h3_trig->GetZaxis()->SetBinLabel(1, "all");
  h3_trig->GetZaxis()->SetBinLabel(2, "ERT4x4a");
  h3_trig->GetZaxis()->SetBinLabel(3, "ERT4x4b");
  h3_trig->GetZaxis()->SetBinLabel(4, "ERT4x4c");
  hm->registerHisto(h3_trig, 1);

  /*
   * Trigger efficiency for pion
   */
  h3_trig_pion = static_cast<TH3*>( h3_trig->Clone("h3_trig_pion") );
  hm->registerHisto(h3_trig_pion, 1);

  /* store single photon information
   *
   * - sector
   * - pT
   * - spin pattern
   * - eta
   * - phi
   *
   */
  const int nbins_hn_1photon[] = {8, n_pTbins, 3, 70, n_phibins};
  const double xmin_hn_1photon[] = {-0.5, 0., -1.5, -0.35, 0.};
  const double xmax_hn_1photon[] = {7.5, 0., 1.5, 0.35, 0.};
  hn_1photon = new THnSparseF("hn_1photon", "Single photon spectrum;sector;p_{T} [GeV];selection;#eta;#phi [rad];",
      5, nbins_hn_1photon, xmin_hn_1photon, xmax_hn_1photon);
  hn_1photon->SetBinEdges(1, pTbins);
  hn_1photon->SetBinEdges(4, phi_twr);
  hn_1photon->GetAxis(2)->SetBinLabel(1, "opposite");
  hn_1photon->GetAxis(2)->SetBinLabel(2, "other");
  hn_1photon->GetAxis(2)->SetBinLabel(3, "same");
  hm->registerHisto(hn_1photon, 1);

  /* store two photons information
   *
   * - sector
   * - single photon pT
   * - photon pair pT
   * - invariant mass
   * - spin pattern
   *
   */
  const int nbins_hn_2photon[] = {8, n_pTbins, n_pTbins, 700, 3};
  const double xmin_hn_2photon[] = {-0.5, 0., 0., 0., -1.5};
  const double xmax_hn_2photon[] = {7.5, 0., 0., 0.7, 1.5};
  hn_2photon = new THnSparseF("hn_2photon", "Two photons spectrum;sector;p^{photon}_{T} [GeV];p^{#pi^0}_{T};m_{inv} [GeV];selection;",
      5, nbins_hn_2photon, xmin_hn_2photon, xmax_hn_2photon);
  hn_2photon->SetBinEdges(1, pTbins);
  hn_2photon->SetBinEdges(2, pTbins);
  hn_2photon->GetAxis(4)->SetBinLabel(1, "opposite");
  hn_2photon->GetAxis(4)->SetBinLabel(2, "other");
  hn_2photon->GetAxis(4)->SetBinLabel(3, "same");
  hm->registerHisto(hn_2photon, 1);

  /* store pion information
   *
   * - sector
   * - pion pT
   * - invariant mass
   * - w/o ToF 
   *
   */
  const int nbins_hn_pion[] = {8, n_pTbins, 300, 2};
  const double xmin_hn_pion[] = {-0.5, 0., 0., -0.5};
  const double xmax_hn_pion[] = {7.5, 0., 0.3, 1.5};
  hn_pion = new THnSparseF("hn_pion", "#pi^{0} spectrum;sector;p^{#pi^0}_{T};m_{inv} [GeV];condition;",
      4, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(1, pTbins);
  hm->registerHisto(hn_pion, 1);

  /*
   * histogram for asymmetry and merging
   *
   * - sector
   * - pion pT
   * - asymmetry
   * - deltaR
   *
   */
  int nbins_hn_asym[] = {8, 60, 50, 50, 300};
  double xmin_hn_asym[] = {-0.5, 0., 0., 0., 0.};
  double xmax_hn_asym[] = {7.5, 30., 1., 50., 0.3};
  hn_asym = new THnSparseF("hn_asym", "Asymmetry and merging; sector; Reco p_{T} [GeV]; Asymmetry; deltaR [cm]; m_{inv} [GeV];",
      5, nbins_hn_asym, xmin_hn_asym, xmax_hn_asym);
  hm->registerHisto(hn_asym, 1);

  /* sotore pion invariant mass information
   *
   * - sector
   * - pion pT
   * - E1
   * - E2
   * - angle
   *
   */
  const int nbins_hn_minv[] = {8, n_pTbins, 300, 300, 100};
  const double xmin_hn_minv[] = {-0.5, 0., 0., 0., 0.};
  const double xmax_hn_minv[] = {7.5, 0., 30., 30., 1.};
  hn_minv = new THnSparseF("hn_minv", "Contribution to m_{inv};sector;p^{#pi^0}_{T} [GeV];E_{1} [GeV];E_{2} [GeV];#sqrt{1-cos(#theta)};",
      5, nbins_hn_minv, xmin_hn_minv, xmax_hn_minv);
  hn_minv->SetBinEdges(1, pTbins);
  //hm->registerHisto(hn_minv, 1);

  return;
}

void FillHisto::EMCRecalibSetup()
{
  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(1);

  string file_ecal_run = toad_loader->location("Run13pp_RunbyRun_Calib.dat");
  string file_tofmap = toad_loader->location("Run13pp510_EMC_TOF_Correction.root");

  if( datatype == ERT )
    emcrecalib->SelectERT();
  else if( datatype == MB )
    emcrecalib->SelectMB();

  emcrecalib->SetEnergyCorrectionFile( file_ecal_run );
  emcrecalib->SetTofCorrectionFile( file_tofmap );

  string _file_tcal = toad_loader->location("tcorr_run13pp500gev.txt");

  emcrecalib_sasha->anaGetCorrTof( _file_tcal.c_str() );

  delete toad_loader;
  return;
}

void FillHisto::ReadTowerStatus(const string &filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> biny >> binz >> status )
  {
    // count tower with bad status for PbSc and PbGl
    if ( status > 10 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status[sector][biny][binz] = status;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void FillHisto::ReadSashaWarnmap(const string &filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int ich = 0;
  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> ich >> status )
  {
    // Attention!! I use my indexing for warn map in this program!!!
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }

    // get tower location
    anatools::TowerLocation(ich, sector, biny, binz);

    // count tower with bad status for PbSc and PbGl
    if ( status > 0 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status[sector][biny][binz] = status;

    // mark edge towers
    if( anatools::Edge_cg(sector, biny, binz) )
      tower_status[sector][biny][binz] = 20;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

bool FillHisto::TestPhoton( const Photon *photon, double bbc_t0 )
{
  if( photon->get_E() > eMin &&
      abs( photon->get_tof() - bbc_t0 ) < 10. )
    return true;
  else
    return false;
}

int FillHisto::GetStatus(const Photon *photon)
{
  int sector, iypos, izpos;
  anatools::TowerLocation(photon, sector, iypos, izpos);

  return tower_status[sector][iypos][izpos];
}

int FillHisto::GetPattern(const PhotonContainer *photoncont)
{
  //int crossing = photoncont->get_crossing();
  int crossing = 0;
  crossing = (crossing + crossing_shift) % 120;
  int pattern = spinpattern_blue[crossing] * spinpattern_yellow[crossing];  
  if( abs(pattern) != 1 ) pattern = 0;

  return pattern;
}

bool FillHisto::TestTrackVeto(const PhotonERT *photon)
{
  // angle between EMC cluster and PC3 track
  double theta_cv = photon->get_theta_cv();

  // charge veto in medium dtheta
  double emc_e = photon->get_E();
  int sector = anatools::GetSector(photon);
  if(sector<6)
  {
    if( theta_cv > 4.22*pow(10,-4) - 1.16*pow(10,-2)*emc_e - 4.53*pow(10,-3)*pow(emc_e,2) &&
        theta_cv < 1.01*pow(10,-1) - 2.02*pow(10,-1)*emc_e + 1.51*pow(10,-1)*pow(emc_e,2) - 3.66*pow(10,-2)*pow(emc_e,3) )
      return false;
  }
  else
  {
    if( theta_cv > 1.27*pow(10,-2) - 2.41*pow(10,-2)*emc_e + 2.26*pow(10,-2)*pow(emc_e,2) &&
        theta_cv < 1.64*pow(10,-2) - 7.38*pow(10,-3)*emc_e + 1.45*pow(10,-1)*exp(-4*emc_e) )
      return false;
  }

  return true;
}
