#include "FillHistoMB.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsPhoton.h"

#include "EmcLocalRecalibratorMB.h"
#include "PhotonContainerCloneMB.h"

#include <PhotonContainerMB.h>
#include <PhotonMB.h>
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
#include <TH1.h>
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

FillHistoMB::FillHistoMB(const string &name, const char *filename) :
  SubsysReco(name),
  hm(NULL),
  emcrecalib(NULL),
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
  hn_minv(NULL),
  g_pileup_PbSc(NULL),
  g_pileup_PbGl(NULL),
  g_pileup_PbSc_notof(NULL),
  g_pileup_PbGl_notof(NULL)
{
  // set output file name
  outFileName = "histos/PhotonNode-";
  outFileName.append(filename);

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;

  irun = 0;  // label the point for TGraph
  runtime = -9999.;
  nmb = 0;  // number of MinBias data

  // initialize the CLOCK Live Counts and BBC narrow Live counts
  for(int i=0; i<6; i++)
    for(int j=0; j<1020; j++)
      n_clock_bbc[i][j] = 0;

  // initialize number of pion signal and background
  for(int i=0; i<2; i++)
  {
    npions_sig[i] = 0;
    npions_bg[i] = 0;
    npions_sig_notof[i] = 0;
    npions_bg_notof[i] = 0;
  }

  // initialize crossing shift and spin pattern
  crossing_shift = 0;
  for(int i=0; i<120; i++)
  {
    spinpattern_blue[i] = 0;
    spinpattern_yellow[i] = 0;
  }
}

FillHistoMB::~FillHistoMB()
{
}

int FillHistoMB::Init(PHCompositeNode *topNode)
{
  // book histograms and graphs
  BookHistograms();

  // EMCal recalibration class
  emcrecalib = new EmcLocalRecalibratorMB();
  if(!emcrecalib)
  {
    cerr << "No emcrecalib" << endl;
    exit(1);
  }

  // read EMCal recalibration file
  EMCRecalibSetup();

  // read CLOCK counts file
  ReadClockCounts("clock-counts.root");

  // read warnmap
  //ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  return EVENT_OK;
}

int FillHistoMB::InitRun(PHCompositeNode *topNode)
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

  // run time and crossing shift for each run
  int runnumber = runheader->get_RunNumber();
  runtime = difftime( runheader->get_TimeStop(), runheader->get_TimeStart() );
  crossing_shift = spinpat->get_crossing_shift();

  // number of MinBias and spin pattern for each run
  //nmb = 0;
  nmb = GetBBCNarrowLive(runnumber);
  for(int i=0; i<120; i++)
  {
    //nmb += spinpat->get_bbc_wide(i);
    spinpattern_blue[i] = spinpat->get_spinpattern_blue(i);
    spinpattern_yellow[i] = spinpat->get_spinpattern_yellow(i);
  }

  // number of pion signal and background in each run
  for(int i=0; i<2; i++)
  {
    npions_sig[i] = 0;
    npions_bg[i] = 0;
    npions_sig_notof[i] = 0;
    npions_bg_notof[i] = 0;
  }

  // load EMCal recalibrations for run and fill
  emcrecalib->ReadEnergyCorrection( spinpat->get_runnumber() );
  emcrecalib->ReadTofCorrection( spinpat->get_fillnumber() );

  return EVENT_OK;
}

int FillHistoMB::process_event(PHCompositeNode *topNode)
{
  PhotonContainerMB *photoncont = findNode::getClass<PhotonContainerMB>(topNode, "PhotonContainerMB");
  if(!photoncont)
  {
    cerr << "No photoncont" << endl;
    return DISCARDEVENT;
  }

  // Run local recalibration of EMCal cluster data
  PhotonContainerCloneMB *photoncont_raw = new PhotonContainerCloneMB(photoncont);
  emcrecalib->ApplyClusterCorrection( photoncont );

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

  // Count events to calculate pile up
  FillPileup( photoncont );

  delete photoncont_raw;

  return EVENT_OK;
}

int FillHistoMB::FillClusterTofSpectrum( const PhotonContainerMB *photoncont, const string &quali )
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  unsigned nphotons = photoncont->Size();

  for( unsigned i = 0; i < nphotons; i++ )
  {
    PhotonMB *photon = photoncont->GetPhoton(i);

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

int FillHistoMB::FillPi0InvariantMass( const PhotonContainerMB *photoncont, const string &quali )
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<nphotons; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        PhotonMB *photon1 = photoncont->GetPhoton(i);
        PhotonMB *photon2 = photoncont->GetPhoton(j);
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

int FillHistoMB::FillTriggerEfficiency( const PhotonContainerMB *photoncont )
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    PhotonMB *photon1 = photoncont->GetPhoton(i);
    int sector = anatools::GetSector( photon1 );
    double photon_pT = anatools::Get_pT( photon1 );

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
          PhotonMB *photon2 = photoncont->GetPhoton(j);
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

int FillHistoMB::FillSinglePhotonSpectrum( const PhotonContainerMB *photoncont )
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  /* Check trigger */
  if( !photoncont->get_bbcnarrow_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  int pattern = GetPattern(photoncont);

  for(unsigned i=0; i<nphotons; i++)
  {
    PhotonMB *photon = photoncont->GetPhoton(i);
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
      double phi = px > 0. ? atan(py/px) : 3.1416+atan(py/px);

      double fill_hn_1photon[] = {sector, pT, pattern, eta, phi};
      hn_1photon->Fill(fill_hn_1photon);
    }
  }

  return EVENT_OK;
}

int FillHistoMB::FillTwoPhotonSpectrum(const PhotonContainerMB *photoncont)
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  /* Check trigger */
  if( !photoncont->get_bbcnarrow_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  int pattern = GetPattern(photoncont);

  for(unsigned i=0; i<nphotons; i++)
    for(unsigned j=0; j<nphotons; j++)
      if(j != i)
      {
        PhotonMB *photon1 = photoncont->GetPhoton(i);
        PhotonMB *photon2 = photoncont->GetPhoton(j);
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

int FillHistoMB::FillPi0Spectrum(const PhotonContainerMB *photoncont)
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  h_events->Fill(1.);

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<nphotons; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        PhotonMB *photon1 = photoncont->GetPhoton(i);
        PhotonMB *photon2 = photoncont->GetPhoton(j);
        if( GetStatus(photon1) == 0 &&
            GetStatus(photon2) == 0 && 
            anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
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
          //if(!trig) continue;

          //TLorentzVector pE1 = anatools::Get_pE(photon1);
          //TLorentzVector pE2 = anatools::Get_pE(photon2);
          //TLorentzVector tot_pE =  pE1 + pE2;
          //double tot_px = tot_pE.Px();
          //double tot_py = tot_pE.Py();
          //double tot_pz = tot_pE.Pz();
          //double tot_pT = tot_pE.Pt();
          //double tot_mom = tot_pE.P();
          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);

          //double eta = tot_mom > 0. ? atan(tot_pz/tot_mom) : 9999.;
          //double phi = tot_px > 0. ? atan(tot_py/tot_px) : 3.1416+atan(tot_py/tot_px);

          /* Check trigger */
          if( photoncont->get_bbcnovtx_scaled() )
          {
            double fill_hn_pion_0[] = {sector, tot_pT, minv, 0.};
            hn_pion->Fill(fill_hn_pion_0);
            if( TestPhoton(photon1, bbc_t0) &&
                TestPhoton(photon2, bbc_t0) )
            {
              double fill_hn_pion_1[] = {sector, tot_pT, minv, 1.};
              hn_pion->Fill(fill_hn_pion_1);
            }
          }

          if( photoncont->get_bbcnarrow_scaled() )
          {
            double fill_hn_pion_2[] = {sector, tot_pT, minv, 2.};
            hn_pion->Fill(fill_hn_pion_2);
            if( TestPhoton(photon1, bbc_t0) &&
                TestPhoton(photon2, bbc_t0) )
            {
              double fill_hn_pion_3[] = {sector, tot_pT, minv, 3.};
              hn_pion->Fill(fill_hn_pion_3);
            }
          }

          //double angle = sqrt( 1. - cos( pE1.Angle(pE2.Vect()) ) );
          //double fill_hn_minv[] = {sector, tot_pT, pE1.E(), pE2.E(), angle};
          //hn_minv->Fill(fill_hn_minv);
        }
      }
  }

  return EVENT_OK;
}

int FillHistoMB::FillPileup(const PhotonContainerMB *photoncont)
{
  /* Get event global parameters */
  double bbc_t0 = photoncont->get_bbc_t0();

  if( !photoncont->get_bbcnarrow_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<nphotons; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        PhotonMB *photon1 = photoncont->GetPhoton(i);
        PhotonMB *photon2 = photoncont->GetPhoton(j);
        if( GetStatus(photon1) == 0 &&
            GetStatus(photon2) == 0 && 
            anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
        {
          int sector1 = anatools::GetSector(photon1);
          int sector2 = anatools::GetSector(photon2);
          //if( sector1 != sector2 ) continue;
          if( !anatools::SectorCheck(sector1,sector2) ) continue;

          int sector = sector1;
          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);

          if( TestPhoton(photon1, bbc_t0) &&
              TestPhoton(photon2, bbc_t0) )
          {
            if(tot_pT>2.)
            {
              if(minv>0.112 && minv<0.162)
                npions_sig[sector/6]++;
              else if( (minv>0.047 && minv<0.097) || (minv>0.177 && minv<0.227) )
                npions_bg[sector/6]++;
            }
          }

          if(tot_pT>2.)
          {
            if(minv>0.112 && minv<0.162)
              npions_sig_notof[sector/6]++;
            else if( (minv>0.047 && minv<0.097) || (minv>0.177 && minv<0.227) )
              npions_bg_notof[sector/6]++;
          }
        }
      }
  }

  return EVENT_OK;
}

int FillHistoMB::EndRun(const int runnumber)
{
  unsigned long long nclock = GetClockLive(runnumber);
  unsigned long scaledown = GetBBCNarrowScaledown(runnumber) + 1;
  //double nclock = 9.6e6 * runtime;
  if(nclock == 0 || nmb == 0)
    return DISCARDEVENT;

  double npions_PbSc = npions_sig[0] - npions_bg[0]/2;
  double npions_PbGl = npions_sig[1] - npions_bg[1]/2;
  double enpions_PbSc = sqrt( npions_sig[0] + npions_bg[0]/2 );
  double enpions_PbGl = sqrt( npions_sig[1] + npions_bg[1]/2 );

  npions_PbSc *= (double)scaledown;
  npions_PbGl *= (double)scaledown;
  enpions_PbSc *= (double)scaledown;
  enpions_PbGl *= (double)scaledown;

  g_pileup_PbSc->SetPoint(irun, (double)nmb/(double)nclock, npions_PbSc/(double)nmb);
  g_pileup_PbSc->SetPointError(irun, 0., enpions_PbSc/(double)nmb);
  g_pileup_PbGl->SetPoint(irun, (double)nmb/(double)nclock, npions_PbGl/(double)nmb);
  g_pileup_PbGl->SetPointError(irun, 0., enpions_PbGl/(double)nmb);

  double npions_PbSc_notof = npions_sig_notof[0] - npions_bg_notof[0]/2;
  double npions_PbGl_notof = npions_sig_notof[1] - npions_bg_notof[1]/2;
  double enpions_PbSc_notof = sqrt( npions_sig_notof[0] + npions_bg_notof[0]/2 );
  double enpions_PbGl_notof = sqrt( npions_sig_notof[1] + npions_bg_notof[1]/2 );

  npions_PbSc_notof *= (double)scaledown;
  npions_PbGl_notof *= (double)scaledown;
  enpions_PbSc_notof  *= (double)scaledown;
  enpions_PbGl_notof  *= (double)scaledown;

  g_pileup_PbSc_notof->SetPoint(irun, (double)nmb/(double)nclock, npions_PbSc_notof/(double)nmb);
  g_pileup_PbSc_notof->SetPointError(irun, 0., enpions_PbSc_notof/(double)nmb);
  g_pileup_PbGl_notof->SetPoint(irun, (double)nmb/(double)nclock, npions_PbGl_notof/(double)nmb);
  g_pileup_PbGl_notof->SetPointError(irun, 0., enpions_PbGl_notof/(double)nmb);

  irun++;

  char name[100];
  sprintf(name, "histos/PhotonNode-%d.root", runnumber);
  hm->dumpHistos(name);
  hn_pion->Reset();

  return EVENT_OK;
}

int FillHistoMB::End(PHCompositeNode *topNode)
{
  //hm->dumpHistos(outFileName);
  delete hm;
  delete emcrecalib;

  return EVENT_OK;
}

void FillHistoMB::BookHistograms()
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
  h_events = new TH1D("h_events", "Events counter", 1,0.5,1.5);
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
  h3_inv_mass_pi0calib = new TH3F("h3_inv_mass_pi0calib", "PhotonMB pair invariant mass;EMCal sector;p_{T} [GeV];m_{inv} [GeV];", 8,-0.5,7.5, n_pTbins,0.,0., 300,0.,0.3);
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
   * - eta
   * - phi
   *
   */
  //const int nbins_hn_pion[] = {8, n_pTbins, 300, 70, n_phibins};
  //const double xmin_hn_pion[] = {-0.5, 0., 0., -0.35, 0.};
  //const double xmax_hn_pion[] = {7.5, 0., 0.3, 0.35, 0.};
  const int nbins_hn_pion[] = {8, n_pTbins, 300, 4};
  const double xmin_hn_pion[] = {-0.5, 0., 0., -0.5};
  const double xmax_hn_pion[] = {7.5, 0., 0.3, 3.5};
  hn_pion = new THnSparseF("hn_pion", "#pi^{0} spectrum;sector;p^{#pi^0}_{T};m_{inv} [GeV];condition;",
      4, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(1, pTbins);
  //hn_pion->SetBinEdges(4, phi_twr);
  hm->registerHisto(hn_pion, 1);
  
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

  /*
   * graph to study pile up in PbSc with ToF cut
   */
  g_pileup_PbSc = new TGraphErrors(10);
  g_pileup_PbSc->SetNameTitle("g_pileup_PbSc", "PbSc;Nmb/Nclock;Npi0/Nmb;");
  hm->registerHisto(g_pileup_PbSc, 1);

  /*
   * graph to study pile up in PbGl with ToF cut
   */
  g_pileup_PbGl = (TGraphErrors*)g_pileup_PbSc->Clone("g_pileup_PbGl");
  g_pileup_PbGl->SetTitle("PbGl;Nmb/Nclock;Npi0/Nmb;");
  hm->registerHisto(g_pileup_PbGl, 1);

  /*
   * graph to study pile up in PbSc without ToF cut
   */
  g_pileup_PbSc_notof = (TGraphErrors*)g_pileup_PbSc->Clone("g_pileup_PbSc_notof");
  g_pileup_PbSc_notof->SetTitle("PbGl;Nmb/Nclock;Npi0/Nmb;");
  hm->registerHisto(g_pileup_PbSc_notof, 1);

  /*
   * graph to study pile up in PbGl without ToF cut
   */
  g_pileup_PbGl_notof = (TGraphErrors*)g_pileup_PbSc->Clone("g_pileup_PbGl_notof");
  g_pileup_PbGl_notof->SetTitle("PbGl;Nmb/Nclock;Npi0/Nmb;");
  hm->registerHisto(g_pileup_PbGl_notof, 1);

  return;
}

void FillHistoMB::EMCRecalibSetup()
{
  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(1);
  string file_ecal_run = toad_loader->location("Run13pp_RunbyRun_Calib.dat");
  string file_tofmap = toad_loader->location("Run13pp510_EMC_TOF_Correction.root");

  emcrecalib->SetEnergyCorrectionFile( file_ecal_run );
  emcrecalib->SetTofCorrectionFile( file_tofmap );

  delete toad_loader;
  return;
}

void FillHistoMB::ReadClockCounts(const string &filename)
{
  TOAD *toad_loader = new TOAD("PhotonNode");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;

  TFile *fin = new TFile( file_location.c_str() );
  TTree *t1 = (TTree*)fin->Get("t1");
  long long runno, clock, bbcnovtx_live, bbcnarrow_live;
  int bbcnovtx_scaledown, bbcnarrow_scaledown;
  t1->SetBranchAddress("runnumber", &runno);
  t1->SetBranchAddress("clock_live", &clock);
  t1->SetBranchAddress("bbcnovtx_live", &bbcnovtx_live);
  t1->SetBranchAddress("bbcnarrow_live", &bbcnarrow_live);
  t1->SetBranchAddress("bbcnovtx_scaledown", &bbcnovtx_scaledown);
  t1->SetBranchAddress("bbcnarrow_scaledown", &bbcnarrow_scaledown);
  int nentries = t1->GetEntries();
  for(int i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    n_clock_bbc[0][i] = runno;
    n_clock_bbc[1][i] = clock;
    n_clock_bbc[2][i] = bbcnovtx_live;
    n_clock_bbc[3][i] = bbcnarrow_live;
    n_clock_bbc[4][i] = bbcnovtx_scaledown;
    n_clock_bbc[5][i] = bbcnarrow_scaledown;
  }
  delete fin;

  delete toad_loader;
  return;
}

void FillHistoMB::ReadTowerStatus(const string &filename)
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

void FillHistoMB::ReadSashaWarnmap(const string &filename)
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

unsigned long long FillHistoMB::GetClockLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[1][i];

  return 0;
}

unsigned long long FillHistoMB::GetBBCNovtxLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[2][i];

  return 0;
}

unsigned long long FillHistoMB::GetBBCNarrowLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[3][i];

  return 0;
}

unsigned long FillHistoMB::GetBBCNovtxScaledown(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[4][i];

  return 0;
}

unsigned long FillHistoMB::GetBBCNarrowScaledown(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[5][i];

  return 0;
}

int FillHistoMB::GetPattern(const PhotonContainerMB *photoncont)
{
  int crossing = photoncont->get_crossing();
  crossing = (crossing + crossing_shift) % 120;
  int pattern = spinpattern_blue[crossing] * spinpattern_yellow[crossing];  
  if( abs(pattern) != 1 ) pattern = 0;

  return pattern;
}

int FillHistoMB::GetStatus(const PhotonMB *photon)
{
  int sector, iypos, izpos;
  anatools::TowerLocation(photon, sector, iypos, izpos);

  return tower_status[sector][iypos][izpos];
}

bool FillHistoMB::TestPhoton( const PhotonMB *photon, double bbc_t0 )
{
  if( photon->get_E() > eMin &&
      photon->get_tof() - bbc_t0 > -10. &&
      photon->get_tof() - bbc_t0 < 10. )
    return true;
  else
    return false;
}

//bool FillHistoMB::TestTrackVeto(const PhotonMB *photon)
//{
//  // angle between EMC cluster and PC3 track
//  double theta_cv = photon->get_theta_cv();
//
//  // charge veto in medium dtheta
//  double emc_e = photon->get_E();
//  int sector = anatools::GetSector(photon);
//  if(sector<6)
//  {
//    if( theta_cv > 4.22*pow(10,-4) - 1.16*pow(10,-2)*emc_e - 4.53*pow(10,-3)*pow(emc_e,2) &&
//        theta_cv < 1.01*pow(10,-1) - 2.02*pow(10,-1)*emc_e + 1.51*pow(10,-1)*pow(emc_e,2) - 3.66*pow(10,-2)*pow(emc_e,3) )
//      return false;
//  }
//  else
//  {
//    if( theta_cv > 1.27*pow(10,-2) - 2.41*pow(10,-2)*emc_e + 2.26*pow(10,-2)*pow(emc_e,2) &&
//        theta_cv < 1.64*pow(10,-2) - 7.38*pow(10,-3)*emc_e + 1.45*pow(10,-1)*exp(-4*emc_e) )
//      return false;
//  }
//
//  return true;
//}
