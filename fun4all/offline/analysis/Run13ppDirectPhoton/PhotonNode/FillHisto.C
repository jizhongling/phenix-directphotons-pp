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

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <RunHeader.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>

#include <ctime>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

const double PI = TMath::Pi();

const double eMin = 0.3;
const double AsymCut = 0.8;

FillHisto::FillHisto(const string &name) :
  SubsysReco(name),
  hm(NULL),
  emcrecalib(NULL),
  emcrecalib_sasha(NULL),
  h_events(NULL),
  h3_tof(NULL),
  h3_tof_raw(NULL),
  h3_minv(NULL),
  h3_minv_raw(NULL),
  h3_bbc(NULL),
  hn_bbc_pion(NULL),
  h3_ert(NULL),
  hn_ert_pion(NULL),
  hn_pion(NULL),
  hn_1photon(NULL),
  hn_2photon(NULL),
  hn_photonbg(NULL),
  runnumber(0),
  runtime(0.)
{
  datatype = ERT;

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;

  // initialize crossing shift and spin pattern
  crossing_shift = 0;
  for(int i=0; i<120; i++)
  {
    spinpattern_blue[i] = 0;
    spinpattern_yellow[i] = 0;
  }

  for(Int_t part=0; part<3; part++)
  {
    h2_photon_eta_phi[part] = NULL;
    h2_cluster_eta_phi[part] = NULL;
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

  double bbc_z = photoncont->get_bbc_z();

  /* Get BBC and ERT trigger counts */
  if( datatype == ERT )
  {
    if( photoncont->get_bbcnarrow_live() && fabs(bbc_z) < 10. )
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
    if( photoncont->get_bbcnarrow_scaled() )
    {
      h_events->Fill("bbc_narrow", 1.);
      if( fabs(bbc_z) < 10. )
      {
        h_events->Fill("bbc_narrow_10cm", 1.);
        if( photoncont->get_ert_c_live() )
          h_events->Fill("bbc_narrow_10cm_ert_c", 1.);
      }
    }

    if( photoncont->get_bbcnovtx_scaled() )
    {
      h_events->Fill("bbc_novtx", 1.);
      if( fabs(bbc_z) < 10. )
        h_events->Fill("bbc_novtx_10cm", 1.);
      if( photoncont->get_bbcnarrow_live() && fabs(bbc_z) < 10. )
      {
        h_events->Fill("bbc_novtx_narrow_10cm", 1.);
        if( photoncont->get_ert_c_live() )
          h_events->Fill("bbc_novtx_narrow_10cm_ert_c", 1.);
      }
    }
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

  // Count events to calculate BBC efficiency
  FillBBCEfficiency( photoncont );

  for( int itype = 0; itype < 3; itype++ )
  {
    // Count events to calculate ERT efficiency
    FillERTEfficiency( photoncont, itype );

    // Analyze photon for pi0 event
    FillPi0Spectrum( photoncont, itype );

    // Analyze photon for direct photon event
    FillPhotonSpectrum( photoncont, itype );
  }

  delete photoncont_raw;

  return EVENT_OK;
}

int FillHisto::FillClusterTofSpectrum(const PhotonContainer *photoncont, const string &quali)
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  for( unsigned i = 0; i < nphotons; i++ )
  {
    Photon *photon = photoncont->GetPhoton(i);
    if( GetStatus(photon) == 0 &&
        photon->get_prob() )
    {
      int sector = anatools::GetSector(photon);
      double tof = photon->get_tof() - bbc_t0;
      double pT = anatools::Get_pT(photon);

      if( quali == "raw" )
        h3_tof_raw->Fill((double)sector, pT, tof);
      else
        h3_tof->Fill((double)sector, pT, tof);
    }
  }

  return EVENT_OK;
}

int FillHisto::FillPi0InvariantMass(const PhotonContainer *photoncont, const string &quali)
{
  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;

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
        if( GetStatus(photon1) == 0 &&
            GetStatus(photon2) == 0 &&
            TestPhoton(photon1, bbc_t0) &&
            TestPhoton(photon2, bbc_t0) &&
            anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
        {
          int sector1 = anatools::GetSector(photon1);
          int sector2 = anatools::GetSector(photon2);
          if( sector1 != sector2 ) continue;

          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);

          if( quali == "raw" )
            h3_minv_raw->Fill((double)sector1, tot_pT, minv);
          else
            h3_minv->Fill((double)sector1, tot_pT, minv);
        }
      }
  }

  return EVENT_OK;
}

int FillHisto::FillBBCEfficiency(const PhotonContainer *photoncont)
{
  /* Check trigger */
  if( !photoncont->get_ert_b_scaled() )
    return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();
  vector<unsigned> v_used;

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon1 = photoncont->GetPhoton(i);
    if( GetStatus(photon1) == 0 &&
        photon1->get_prob() )
    {
      v_used.push_back(i);

      int sector = anatools::GetSector( photon1 );
      TLorentzVector pE = anatools::Get_pE(photon1);
      double photon_pT = pE.Pt();
      double ecore = pE.E();
      double eta = pE.Eta();
      double bbc_trig = 0.;

      h3_bbc->Fill((double)sector, photon_pT, 0.);
      if( photoncont->get_bbcnovtx_live() )
      {
        bbc_trig = 1.;
        h3_bbc->Fill((double)sector, photon_pT, 1.);
      }

      double fill_hn_photonbg[] = {(double)sector, photon_pT, ecore, eta, bbc_trig};
      hn_photonbg->Fill(fill_hn_photonbg);

      for(unsigned j=0; j<nphotons; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          Photon *photon2 = photoncont->GetPhoton(j);
          if( GetStatus(photon2) == 0 &&
              photon2->get_prob() )
          {
            int sector2 = anatools::GetSector( photon2 );
            if( !anatools::SectorCheck(sector,sector2) ) continue;

            double tot_pT = anatools::GetTot_pT(photon1, photon2);
            double minv = anatools::GetInvMass(photon1, photon2);

            double fill_hn_bbc_pion[] = {(double)sector, tot_pT, minv, 0.};
            hn_bbc_pion->Fill(fill_hn_bbc_pion);
            if( photoncont->get_bbcnovtx_live() )
            {
              fill_hn_bbc_pion[3] = 1.;
              hn_bbc_pion->Fill(fill_hn_bbc_pion);
            }
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int FillHisto::FillERTEfficiency(const PhotonContainer *photoncont, const int evtype)
{
  /* Check trigger */
  if( !IsEventType(evtype, photoncont) )
    return DISCARDEVENT;

  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( fabs(bbc_z) > 10. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();
  vector<unsigned> v_used;

  /* Fire ERT on arm 0 (west) or 1 (east) */
  bool FireERT[2] = {};
  if( datatype == ERT )
  {
    for(unsigned i=0; i<nphotons; i++)
    {
      Photon *photon = photoncont->GetPhoton(i);
      int arm = anatools::GetSector(photon) / 4;
      bool trig1 = photon->get_trg1();
      bool trig2 = photon->get_trg2();
      bool trig3 = photon->get_trg3();
      if( CheckGammaTrigger(evtype, trig1, trig2, trig3) )
        FireERT[arm] = true;
    }
  }

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon1 = photoncont->GetPhoton(i);
    int arm = anatools::GetSector(photon1) / 4;

    /* Require the other arm to be fired for ERT sample */
    int oarm = ( arm==0 ? 1 : 0 );
    if( datatype == ERT && !FireERT[oarm] ) continue;

    if( GetStatus(photon1) == 0 &&
        TestPhoton(photon1, bbc_t0) )
    {
      v_used.push_back(i);

      int sector = anatools::GetSector( photon1 );
      double photon_pT = anatools::Get_pT( photon1 );

      bool trig1 = photon1->get_trg1();
      bool trig2 = photon1->get_trg2();
      bool trig3 = photon1->get_trg3();

      h3_ert->Fill((double)sector, photon_pT, (double)evtype);
      if( CheckGammaTrigger(evtype, trig1, trig2, trig3) )
        h3_ert->Fill((double)sector, photon_pT, evtype+3.);

      for(unsigned j=0; j<nphotons; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          Photon *photon2 = photoncont->GetPhoton(j);
          if( GetStatus(photon2) == 0 &&
              TestPhoton(photon2, bbc_t0) )
          {
            int sector2 = anatools::GetSector( photon2 );
            if( !anatools::SectorCheck(sector,sector2) ) continue;

            if( photon2->get_E() > photon1->get_E() )
            {
              sector = sector2;
              trig1 = photon2->get_trg1();
              trig2 = photon2->get_trg2();
              trig3 = photon2->get_trg3();
            }

            double tot_pT = anatools::GetTot_pT(photon1, photon2);
            double minv = anatools::GetInvMass(photon1, photon2);

            double fill_hn_ert_pion[] = {(double)sector, tot_pT, minv, (double)evtype};
            hn_ert_pion->Fill(fill_hn_ert_pion);
            if( CheckGammaTrigger(evtype, trig1, trig2, trig3) )
            {
              fill_hn_ert_pion[3] = evtype + 3.;
              hn_ert_pion->Fill(fill_hn_ert_pion);
            }
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int FillHisto::FillPi0Spectrum(const PhotonContainer *photoncont, const int evtype)
{
  if( datatype == MB && evtype != 0 )
    return DISCARDEVENT;

  /* Check trigger */
  if( !IsEventType(evtype, photoncont) )
    return DISCARDEVENT;

  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  if( fabs(bbc_z) > 10. ) return DISCARDEVENT;

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
        if( GetStatus(photon1) == 0 &&
            GetStatus(photon2) == 0 &&
            anatools::GetAsymmetry_E(photon1, photon2) < AsymCut )
        {
          int sector1 = anatools::GetSector(photon1);
          int sector2 = anatools::GetSector(photon2);
          if( !anatools::SectorCheck(sector1,sector2) ) continue;

          int sector = sector1;
          bool trig1 = photon1->get_trg1();
          bool trig2 = photon1->get_trg2();
          bool trig3 = photon1->get_trg3();
          if( photon2->get_E() > photon1->get_E() )
          {
            sector = sector2;
            trig1 = photon2->get_trg1();
            trig2 = photon2->get_trg2();
            trig3 = photon2->get_trg3();
          }

          if( datatype == ERT && !CheckGammaTrigger(evtype, trig1, trig2, trig3) )
            continue;

          double tot_pT = anatools::GetTot_pT(photon1, photon2);
          double minv = anatools::GetInvMass(photon1, photon2);

          double fill_hn_pion[] = {(double)sector, tot_pT, minv, 0., (double)evtype};
          hn_pion->Fill(fill_hn_pion);
          if( fabs( photon1->get_tof() - bbc_t0 ) < 10. &&
              fabs( photon2->get_tof() - bbc_t0 ) < 10. )
          {
            fill_hn_pion[3] = 1.;
            hn_pion->Fill(fill_hn_pion);
          }
          if( photon1->get_prob() &&
              photon2->get_prob() )
          {
            fill_hn_pion[3] = 2.;
            hn_pion->Fill(fill_hn_pion);
          }
          if( TestPhoton(photon1, bbc_t0) &&
              TestPhoton(photon2, bbc_t0) )
          {
            fill_hn_pion[3] = 3.;
            hn_pion->Fill(fill_hn_pion);
          }
        } // GetStatus and asymmetry cut
      } // j loop
  } // i loop

  return EVENT_OK;
}

int FillHisto::FillPhotonSpectrum(const PhotonContainer *photoncont, const int evtype)
{
  if( datatype == MB && evtype != 0 )
    return DISCARDEVENT;

  /* Check trigger */
  if( !IsEventType(evtype, photoncont) )
    return DISCARDEVENT;

  /* Get event global parameters */
  double bbc_z = photoncont->get_bbc_z();
  double bbc_t0 = photoncont->get_bbc_t0();
  int pattern = GetPattern(photoncont);
  if( fabs(bbc_z) > 10. ) return DISCARDEVENT;

  unsigned nphotons = photoncont->Size();

  for(unsigned i=0; i<nphotons; i++)
  {
    Photon *photon1 = photoncont->GetPhoton(i);
    if( GetStatus(photon1) == 0 )
    {
      int sector = anatools::GetSector(photon1);
      int part = -1;
      if(sector < 0) part = -1;
      else if(sector < 4) part = 0;
      else if(sector < 6) part = 1;
      else if(sector < 8) part = 2;

      TLorentzVector pE = anatools::Get_pE(photon1);
      double pT = pE.Pt();
      double eta = pE.Eta();
      double phi = pE.Phi();
      if(sector >= 4)
      {
        pE.RotateZ(-PI);
        phi = pE.Phi() + PI;
      }

      bool trig1 = photon1->get_trg1();
      bool trig2 = photon1->get_trg2();
      bool trig3 = photon1->get_trg3();

      if( datatype == ERT && !CheckGammaTrigger(evtype, trig1, trig2, trig3) )
        continue;

      if( evtype == 2 && part >= 0 &&
          pT > 5. && pT < 10. )
        h2_cluster_eta_phi[part]->Fill(eta, phi);

      double fill_hn_1photon[] = {(double)sector, pT, (double)pattern, 0., (double)evtype};
      hn_1photon->Fill(fill_hn_1photon);
      if( fabs( photon1->get_tof() - bbc_t0 ) < 10. )
      {
        fill_hn_1photon[3] = 1.;
        hn_1photon->Fill(fill_hn_1photon);
      }
      if( photon1->get_prob() )
      {
        fill_hn_1photon[3] = 2.;
        hn_1photon->Fill(fill_hn_1photon);
      }
      if( TestPhoton(photon1, bbc_t0) )
      {
        if( evtype == 2 && part >= 0 &&
            pT > 5. && pT < 10. )
          h2_photon_eta_phi[part]->Fill(eta, phi);

        fill_hn_1photon[3] = 3.;
        hn_1photon->Fill(fill_hn_1photon);
      }

      for(unsigned j=0; j<nphotons; j++)
        if(j != i)
        {
          Photon *photon2 = photoncont->GetPhoton(j);
          if( GetStatus(photon2) != 0 ) continue;
          double minv = anatools::GetInvMass(photon1, photon2);

          double fill_hn_2photon[] = {(double)sector, pT, minv, (double)pattern, 0., (double)evtype};
          hn_2photon->Fill(fill_hn_2photon);
          if( fabs( photon1->get_tof() - bbc_t0 ) < 10. &&
              fabs( photon2->get_tof() - bbc_t0 ) < 10. )
          {
            fill_hn_2photon[4] = 1.;
            hn_2photon->Fill(fill_hn_2photon);
          }
          if( photon1->get_prob() && 
              photon2->get_prob() )
          {
            fill_hn_2photon[4] = 2.;
            hn_2photon->Fill(fill_hn_2photon);
          }
          if( TestPhoton(photon1, bbc_t0) &&
              TestPhoton(photon2, bbc_t0) )
          {
            fill_hn_2photon[4] = 3.;
            hn_2photon->Fill(fill_hn_2photon);
          }
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int FillHisto::EndRun(const int runnumber)
{
  // Write histograms to file for this run
  char fname[200];
  strcpy(fname, "");
  if( strcmp(Name(), "FillHisto_TAXI") == 0 )
  {
    if( datatype == ERT )
      sprintf(fname, "histos-ERT/PhotonNode-%d.root", runnumber);
    else if( datatype == MB )
      sprintf(fname, "histos-MB/PhotonNode-%d.root", runnumber);
  }
  else if( strcmp(Name(), "FillHisto_Sample") == 0 )
  {
    sprintf(fname, "histos-sample/PhotonNode-%d.root", runnumber);
  }
  if( strlen(fname) > 0 )
    hm->dumpHistos(fname);

  // Reset histograms for the next run
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

  /* pT bins */
  const int npT = 30;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  /* eta and phi bins step size */
  const double step[2] = {0.011, 0.008};

  /* eta bins */
  const int neta = 100;
  double etabin[2][neta+1];
  for(int part=0; part<2; part++)
    for(int it=0; it<=neta; it++)
      etabin[part][it] = step[part] * ( it - neta/2 );

  /* phi sector */
  const double phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8,
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  /* phi bins */
  const int nphi_sec[2] = {36, 48};
  const int nphi = (nphi_sec[0]/2+1)*6 + (nphi_sec[1]/2+1)*2 + 4 - 1;
  double phibin[nphi+1];
  int iphi = 0;
  for(int is=0; is<8; is++)
    for(int it=0; it<=nphi_sec[is/6]/2; it++)
      phibin[iphi++] = phi_sec[is] + step[is/6] * 2 * ( it - nphi_sec[is/6]/4 );
  phibin[iphi++] = -PI*3/16 - 0.02;
  phibin[iphi++] = PI*5/16 + 0.02;
  phibin[iphi++] = PI*11/16 - 0.02;
  phibin[iphi++] = PI*19/16 + 0.02;
  sort(phibin, phibin+nphi);

  /* Events counter */
  if( datatype == ERT )
  {
    h_events = new TH1F("h_events", "Events counter", 3,-0.5,2.5);
    h_events->GetXaxis()->SetBinLabel(1, "ert_a");
    h_events->GetXaxis()->SetBinLabel(2, "ert_b");
    h_events->GetXaxis()->SetBinLabel(3, "ert_c");
  }
  else if( datatype == MB )
  {
    h_events = new TH1F("h_events", "Events counter", 7,-0.5,6.5);
    h_events->GetXaxis()->SetBinLabel(1, "bbc_narrow");
    h_events->GetXaxis()->SetBinLabel(2, "bbc_narrow_10cm");
    h_events->GetXaxis()->SetBinLabel(3, "bbc_narrow_10cm_ert_c");
    h_events->GetXaxis()->SetBinLabel(4, "bbc_novtx");
    h_events->GetXaxis()->SetBinLabel(5, "bbc_novtx_10cm");
    h_events->GetXaxis()->SetBinLabel(6, "bbc_novtx_narrow_10cm");
    h_events->GetXaxis()->SetBinLabel(7, "bbc_novtx_narrow_10cm_ert_c");
  }
  hm->registerHisto(h_events);

  /* ToF calibration */
  h3_tof = new TH3F("h3_tof", "ToF;sector;p_{T} [GeV];tof [ns];", 8,-0.5,7.5, npT,0.,0., 1001,-100.05,100.05);
  h3_tof->GetYaxis()->Set(npT, pTbin);
  hm->registerHisto(h3_tof);

  h3_tof_raw = static_cast<TH3*>( h3_tof->Clone("h3_tof_raw") );
  hm->registerHisto(h3_tof_raw);

  /* Storing invariant mass of photon pairs in different sectors and pT bins.
   * Require both photons to be in same sector.
   * Used to check sector-by-sector EMCal energy calibration.
   */
  h3_minv = new TH3F("h3_minv", "Photon pair invariant mass;sector;p_{T} [GeV];m_{inv} [GeV];", 8,-0.5,7.5, npT,0.,0., 300,0.,0.3);
  h3_minv->GetYaxis()->Set(npT, pTbin);
  hm->registerHisto(h3_minv);

  h3_minv_raw = static_cast<TH3*>( h3_minv->Clone("h3_minv_raw") );
  hm->registerHisto(h3_minv_raw);

  /* BBC trigger efficiency for photon */
  h3_bbc = new TH3F("h3_bbc", "BBC efficiency;sector;p_{T} [GeV];w/o BBC;", 8,-0.5,7.5, npT,0.,0., 2,-0.5,1.5);
  h3_bbc->GetYaxis()->Set(npT, pTbin);
  h3_bbc->GetZaxis()->SetBinLabel(1, "all");
  h3_bbc->GetZaxis()->SetBinLabel(2, "bbc");
  hm->registerHisto(h3_bbc);

  /* BBC Trigger efficiency for pion */
  const int nbins_hn_bbc_pion[] = {8, npT, 300, 2};
  const double xmin_hn_bbc_pion[] = {-0.5, 0., 0., -0.5};
  const double xmax_hn_bbc_pion[] = {7.5, 0., 0.3, 1.5};
  hn_bbc_pion = new THnSparseF("hn_bbc_pion", "BBC efficiency;sector;p_{T} [GeV];m_{inv} [GeV];w/o BBC;",
      4, nbins_hn_bbc_pion, xmin_hn_bbc_pion, xmax_hn_bbc_pion);
  hn_bbc_pion->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_bbc_pion);

  /* ERT trigger efficiency for photon */
  h3_ert = new TH3F("h3_ert", "ERT efficiency;sector;p_{T} [GeV];w/o ERT;", 8,-0.5,7.5, npT,0.,0., 6,-0.5,5.5);
  h3_ert->GetYaxis()->Set(npT, pTbin);
  hm->registerHisto(h3_ert);

  /* ERT Trigger efficiency for pion */
  hn_ert_pion = static_cast<THnSparse*>( hn_bbc_pion->Clone("hn_ert_pion") );
  hn_ert_pion->SetTitle("ERT efficiency");
  hn_ert_pion->GetAxis(3)->SetTitle("w/o ERT");
  hn_ert_pion->GetAxis(3)->Set(6, -0.5, 5.5);
  hm->registerHisto(hn_ert_pion);

  for(Int_t part=0; part<3; part++)
  {
    h2_photon_eta_phi[part] = new TH2F(Form("h2_photon_eta_phi_part%d",part), "Photon #eta and #phi distribution;#eta;#phi;", neta,etabin[part/2], nphi,phibin);
    h2_cluster_eta_phi[part] = (TH2*)h2_photon_eta_phi[part]->Clone( Form("h2_cluster_eta_phi_part%d",part) );
    hm->registerHisto(h2_photon_eta_phi[part]);
    hm->registerHisto(h2_cluster_eta_phi[part]);
  }

  /* Store pion information */
  const int nbins_hn_pion[] = {8, npT, 300, 4, 3};
  const double xmin_hn_pion[] = {-0.5, 0., 0., -0.5, -0.5};
  const double xmax_hn_pion[] = {7.5, 0., 0.3, 3.5, 2.5};
  hn_pion = new THnSparseF("hn_pion", "#pi^{0} spectrum;sector;p_{T} [GeV];m_{inv} [GeV];cut;evtype;",
      5, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_pion);

  /* Store single photon information */
  const int nbins_hn_1photon[] = {8, npT, 3, 4, 3};
  const double xmin_hn_1photon[] = {-0.5, 0., -1.5, -0.5, -0.5};
  const double xmax_hn_1photon[] = {7.5, 0., 1.5, 3.5, 2.5};
  hn_1photon = new THnSparseF("hn_1photon", "Single photon spectrum;sector;p_{T} [GeV];pattern;cut;evtype;",
      5, nbins_hn_1photon, xmin_hn_1photon, xmax_hn_1photon);
  hn_1photon->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_1photon);

  /* Store two photons information */
  const int nbins_hn_2photon[] = {8, npT, 700, 3, 4, 3};
  const double xmin_hn_2photon[] = {-0.5, 0., 0., -1.5, -0.5, -0.5};
  const double xmax_hn_2photon[] = {7.5, 0., 0.7, 1.5, 3.5, 2.5};
  hn_2photon = new THnSparseF("hn_2photon", "Two photons spectrum;sector;p_{T} [GeV];m_{inv} [GeV];pattern;cut;evtype;",
      6, nbins_hn_2photon, xmin_hn_2photon, xmax_hn_2photon);
  hn_2photon->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_2photon);

  /* Store photon bg information */
  const int nbins_hn_photonbg[] = {8, npT, 100, 40, 2};
  const double xmin_hn_photonbg[] = {-0.5, 0., 0., -0.4, -0.5};
  const double xmax_hn_photonbg[] = {7.5, 0., 35., 0.4, 1.5};
  hn_photonbg = new THnSparseF("hn_photonbg", "Photon bg spectrum;sector;p_{T} [GeV];E [GeV];#eta;BBC;",
      5, nbins_hn_photonbg, xmin_hn_photonbg, xmax_hn_photonbg);
  hn_photonbg->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_photonbg);

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

bool FillHisto::IsEventType(const int evtype, const PhotonContainer *photoncont)
{
  if( datatype == ERT )
  {
    if( evtype == 0 && photoncont->get_ert_a_scaled() && photoncont->get_bbcnarrow_live() )
      return true;
    else if( evtype == 1 && photoncont->get_ert_b_scaled() && photoncont->get_bbcnarrow_live() )
      return true;
    else if( evtype == 2 && photoncont->get_ert_c_scaled() && photoncont->get_bbcnarrow_live() )
      return true;
  }
  else if( datatype == MB )
  {
    if( photoncont->get_bbcnarrow_scaled() )
      return true;
  }

  return false;
}

bool FillHisto::CheckGammaTrigger(const int evtype, const bool trig1, const bool trig2, const bool trig3)
{
  if( evtype == 0 )
    return trig1;
  else if( evtype == 1 )
    return trig2;
  else if( evtype == 2 )
    return trig3;

  return false;
}

bool FillHisto::TestPhoton(const Photon *photon, double bbc_t0)
{
  if( photon->get_E() > eMin &&
      photon->get_prob() &&
      fabs( photon->get_tof() - bbc_t0 ) < 10. )
    return true;

  return false;
}

bool FillHisto::TestTrackVeto(const PhotonERT *photon)
{
  // Angle between EMC cluster and PC3 track
  double theta_cv = photon->get_theta_cv();

  // Charge veto in medium dtheta
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
