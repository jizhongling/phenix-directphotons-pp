#include "PhotonHistos.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"
#include "AnaToolsTrigger.h"

#include "EmcLocalRecalibrator.h"
#include "EmcLocalRecalibratorSasha.h"
#include "SpinPattern.h"

#include <RunHeader.h>
#include <SpinDBOutput.hh>
#include <SpinDBContent.hh>

#include <PHGlobal.h>
#include <TrigLvl1.h>
#include <ErtOut.h>
#include <EmcIndexer.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <emcTowerContainer.h>
#include <emcTowerContent.h>
#include <PHCentralTrack.h>

#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <PHNodeIterator.h>

#include <TOAD.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <iostream>
#include <fstream>

using namespace std;

/* Some constants */
const double epsilon = TMath::Limits<float>::Epsilon();
const double PI = TMath::Pi();
const TVector2 v2_2PI(0.,2.*PI);

/* Some cuts for photon identification */
const double eMin = 0.3;
const double probMin = 0.02;
const double tofMax = 10.;
const double tofMaxIso = 50.;
const double AsymCut = 0.8;

/* Some cuts for isolation cut */
const double eClusMin = 0.15;
const double pTrkMin = 0.2;
const double pTrkMax = 15.;

/* Isolation cut cone angle and energy fraction */
const double cone_angle = 0.5;
const double eratio = 0.1;

/* Triger bit */
const unsigned bit_ppg = 0x70000000;
const unsigned bit_bbcnarrow = 0x00000010;
const unsigned bit_bbcnovtx = 0x00000002;
const unsigned bit_ert4x4[4] = {0x00000080, 0x00000040, 0x00000100, 0x000001C0};  // ert4x4a/b/c/or

PhotonHistos::PhotonHistos(const string &name, const char *filename) :
  SubsysReco(name),
  emcrecalib(NULL),
  emcrecalib_sasha(NULL),
  spinpattern(NULL),
  runnumber(0),
  fillnumber(0),
  hm(NULL)
{
  datatype = ERT;

  outFile = "PhotonHistos-";
  outFile.append(filename);

  /* Initialize array for tower status */
  for(int isector=0; isector<NSEC; isector++)
    for(int ibiny=0; ibiny<NY; ibiny++)
      for(int ibinz=0; ibinz<NZ; ibinz++)
      {
        tower_status_nils[isector][ibiny][ibinz] = 0;
        tower_status_sasha[isector][ibiny][ibinz] = 0;
      }

  h_events = NULL;
  h_prod = NULL;
  for(int ih=0; ih<nh_calib; ih++)
  {
    h2_tof[ih] = NULL;
    h2_minv[ih] = NULL;
  }
  for(int ih=0; ih<nh_bbc; ih++)
  {
    h_bbc[ih] = NULL;
    h2_bbc_pion[ih] = NULL;
  }
  for(int ih=0; ih<nh_ert; ih++)
  {
    h_ert[ih] = NULL;
    h2_ert_pion[ih] = NULL;
  }
  for(int ih=0; ih<nh_etwr; ih++)
    h3_etwr[ih] = NULL;
  for(int ih=0; ih<nh_dcpartqual; ih++)
  {
    h3_dcdphiz[ih] = NULL;
    h2_alphaboard[ih] = NULL;
  }
  for(int ih=0; ih<nh_dcquality; ih++)
    h3_dclive[ih] = NULL;
  for(int ih=0; ih<nh_dcpart; ih++)
    h2_emcdphiz[ih] = NULL;
  for(int ih=0; ih<nh_eta_phi; ih++)
    h2_eta_phi[ih] = NULL;
  for(int ih=0; ih<nh_pion; ih++)
    h2_pion[ih] = NULL;
  for(int ih=0; ih<nh_1photon; ih++)
    h_1photon[ih] = NULL;
  for(int ih=0; ih<nh_2photon; ih++)
  {
    h2_2photon[ih] = NULL;
    h2_2photon2pt[ih] = NULL;
  }
}

PhotonHistos::~PhotonHistos()
{
}

int PhotonHistos::Init(PHCompositeNode *topNode)
{
  /* Create a class to store spin information */
  spinpattern = new SpinPattern();
  if(!spinpattern)
  {
    cerr << "Failure to create spin pattern" << endl;
    exit(1);
  }

  /* EMCal recalibration class */
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

  /* Read EMCal recalibration file */
  EMCRecalibSetup();

  /* Read warnmap */
  ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  /* Book histograms and graphs */
  BookHistograms();

  return EVENT_OK;
}

int PhotonHistos::InitRun(PHCompositeNode *topNode)
{
  /* Get run number */
  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!runheader)
  {
    cerr << "No runheader" << endl;
    return ABORTRUN;
  }
  runnumber = runheader->get_RunNumber();

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  /* Initialize object to access spin DB */
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  /* Retrieve entry from Spin DB and get fill number */
  int qa_level = spin_out.GetDefaultQA(runnumber);
  spin_out.StoreDBContent(runnumber, runnumber, qa_level);
  spin_out.GetDBContentStore(spin_cont, runnumber);
  fillnumber = spin_cont.GetFillNumber();

  /* Load EMCal recalibrations for run and fill */
  emcrecalib->ReadEnergyCorrection( runnumber );
  emcrecalib->ReadTofCorrection( fillnumber );

  /* Update spinpattern */
  if( spin_out.CheckRunRow(runnumber,qa_level) == 1 &&
      spin_cont.GetRunNumber() == runnumber )
    UpdateSpinPattern(spin_cont);
  else
    spinpattern->Reset();

  return EVENT_OK;
}

int PhotonHistos::process_event(PHCompositeNode *topNode)
{
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  TrigLvl1 *data_triggerlvl1 = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  ErtOut *data_ert = findNode::getClass<ErtOut>(topNode, "ErtOut");
  emcClusterContainer *data_emccontainer_raw = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
  emcTowerContainer *data_emctwrcontainer = findNode::getClass<emcTowerContainer>(topNode, "emcHitContainer");
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");

  if(!data_global)
  {
    cerr << "No gbl" << endl;
    return DISCARDEVENT;
  }
  if(!data_triggerlvl1)
  {
    cerr << "No trg" << endl;
    return DISCARDEVENT;
  }
  if(!data_ert)
  {
    cerr << "No ert" << endl;
    return DISCARDEVENT;
  }
  if(!data_emccontainer_raw)
  {
    cerr << "No emcont" << endl;
    return DISCARDEVENT;
  }
  if(!data_emctwrcontainer)
  {
    cerr << "No emtwrcont" << endl;
    return DISCARDEVENT;
  }
  if(!data_tracks)
  {
    cerr << "\nNo tracker data" << endl;
    return DISCARDEVENT;
  }

  /* Check trigger */
  unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();
  if( (lvl1_live & bit_ppg) || (lvl1_scaled & bit_ppg) ) return DISCARDEVENT;

  /* Run local recalibration of EMCal cluster data */
  emcClusterContainer *data_emccontainer = data_emccontainer_raw->clone();
  //emcrecalib->ApplyClusterCorrection( data_emccontainer );
  emcrecalib_sasha->ApplyClusterCorrection( runnumber, data_emccontainer );

  /* Event counts */
  FillEventCounts(data_global, data_triggerlvl1);

  /* Store ToF information for cluster as calibration check */
  FillClusterTofSpectrum(data_emccontainer_raw, data_global, "raw");
  FillClusterTofSpectrum(data_emccontainer, data_global);

  /* Analyze pi0s events for crosscheck */
  FillPi0InvariantMass(data_emccontainer_raw, data_global, "raw");
  FillPi0InvariantMass(data_emccontainer, data_global);

  /* Count events to calculate BBC efficiency */
  FillBBCEfficiency(data_emccontainer, data_triggerlvl1);

  /* Fill EMCal cluster energy distribution on towers */
  //FillTowerEnergy(data_emccontainer_raw, data_emctwrcontainer, data_global, data_triggerlvl1);

  /* Fill DC track quality information */
  //FillTrackQuality(data_emccontainer, data_tracks, data_global, data_triggerlvl1);

  for(int itype=0; itype<4; itype++)
  {
    /* Count events to calculate ERT efficiency */
    FillERTEfficiency(data_emccontainer, data_tracks, data_global, data_triggerlvl1, data_ert, itype);

    /* Analyze photon for pi0 event */
    FillPi0Spectrum(data_emccontainer, data_tracks, data_global, data_triggerlvl1, data_ert, itype);

    /* Analyze photon for direct photon event */
    FillPhotonSpectrum(data_emccontainer, data_tracks, data_global, data_triggerlvl1, data_ert, itype);
  }

  /* Clean up */
  delete data_emccontainer;

  return EVENT_OK;
}

int PhotonHistos::FillEventCounts(const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1)
{
  /* Get trigger */
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  const unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();

  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();

  /* Fill ERT trigger counts */
  if( datatype == ERT )
  {
    if( (lvl1_live & bit_bbcnarrow) )
    {
      if( abs(bbc_z) < 10. )
      {
        if( lvl1_scaled & bit_ert4x4[0] )
          h_events->Fill("ert_a_10cm", 1.);
        if( lvl1_scaled & bit_ert4x4[1] )
          h_events->Fill("ert_b_10cm", 1.);
        if( lvl1_scaled & bit_ert4x4[2] )
          h_events->Fill("ert_c_10cm", 1.);
      }
      if( abs(bbc_z) < 30. )
      {
        if( lvl1_scaled & bit_ert4x4[0] )
          h_events->Fill("ert_a_30cm", 1.);
        if( lvl1_scaled & bit_ert4x4[1] )
          h_events->Fill("ert_b_30cm", 1.);
        if( lvl1_scaled & bit_ert4x4[2] )
          h_events->Fill("ert_c_30cm", 1.);
      }
    }
  }

  /* Fill BBC trigger counts */
  else if( datatype == MB )
  {
    if( lvl1_scaled & bit_bbcnarrow )
    {
      h_events->Fill("bbc_narrow", 1.);
      if( abs(bbc_z) < 10. )
      {
        h_events->Fill("bbc_narrow_10cm", 1.);
        if( lvl1_live & bit_ert4x4[2] )
          h_events->Fill("bbc_narrow_10cm_ert_c", 1.);
      }
      if( abs(bbc_z) < 30. )
      {
        h_events->Fill("bbc_narrow_30cm", 1.);
        if( lvl1_live & bit_ert4x4[2] )
          h_events->Fill("bbc_narrow_30cm_ert_c", 1.);
      }
    }
  }

  return EVENT_OK;
}

int PhotonHistos::FillClusterTofSpectrum(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const string &quali)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned ncluster = data_emccontainer->size();

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster = data_emccontainer->getCluster(i);
    if( IsGoodTower(cluster) &&
        cluster->ecore() > eMin &&
        cluster->prob_photon() > probMin )
    {
      int sector = anatools::GetSector(cluster);
      double tof = cluster->tofcorr() - bbc_t0;
      double pT = anatools::Get_pT(cluster);

      int raw = 0;
      if( quali == "raw" )
        raw = 1;

      if(sector >= 0 && sector <= 7)
      {
        int ih = sector + 8*raw;
        h2_tof[ih]->Fill(pT, tof);
      }
    }
  }

  return EVENT_OK;
}

int PhotonHistos::FillPi0InvariantMass(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const string &quali)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  unsigned ncluster = data_emccontainer->size();
  vector<unsigned> v_used;

  for(unsigned i=0; i<ncluster; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<ncluster; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
        emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
        if( IsGoodTower(cluster1) &&
            IsGoodTower(cluster2) &&
            TestPhoton(cluster1, bbc_t0) &&
            TestPhoton(cluster2, bbc_t0) &&
            anatools::GetAsymmetry_E(cluster1, cluster2) < AsymCut )
        {
          int sector1 = anatools::GetSector(cluster1);
          int sector2 = anatools::GetSector(cluster2);
          if( sector1 != sector2 ) continue;

          double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
          double minv = anatools::GetInvMass(cluster1, cluster2);

          int raw = 0;
          if( quali == "raw" )
            raw = 1;

          if(sector1 >= 0 && sector1 <= 7)
          {
            int ih = sector1 + 8*raw;
            h2_minv[ih]->Fill(tot_pT, minv);
          }
        }
      }
  }

  return EVENT_OK;
}

int PhotonHistos::FillBBCEfficiency(const emcClusterContainer *data_emccontainer, const TrigLvl1 *data_triggerlvl1)
{
  /* Check trigger */
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  const unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();
  if( !(lvl1_scaled & bit_ert4x4[1]) )
    return DISCARDEVENT;

  unsigned ncluster = data_emccontainer->size();
  vector<unsigned> v_used;

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
    if( IsGoodTower(cluster1) &&
        cluster1->ecore() > eMin &&
        cluster1->prob_photon() > probMin )
    {
      v_used.push_back(i);

      int sector = anatools::GetSector(cluster1);
      double pT = anatools::Get_pT(cluster1);

      int bbc_trig = 0;
      if( lvl1_live & bit_bbcnovtx )
        bbc_trig = 1;

      if(sector >= 0 && sector <= 7)
      {
        int ih = sector + 8*bbc_trig;
        h_bbc[ih]->Fill(pT);
      }

      for(unsigned j=0; j<ncluster; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
          if( IsGoodTower(cluster2) &&
              cluster2->ecore() > eMin &&
              cluster2->prob_photon() > probMin )
          {
            int sector2 = anatools::GetSector( cluster2 );
            if( !anatools::SectorCheck(sector,sector2) ) continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);

            int bbc_trig = 0;
            if( lvl1_live & bit_bbcnovtx )
              bbc_trig = 1;

            if(sector >= 0 && sector <= 7)
            {
              int ih = sector + 8*bbc_trig;
              h2_bbc_pion[ih]->Fill(tot_pT, minv);
            }
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::FillERTEfficiency(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype)
{
  if(evtype >= 3)
    return DISCARDEVENT;

  /* Check trigger */
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  const unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();
  if( (lvl1_live & bit_ppg) || (lvl1_scaled & bit_ppg) )
    return DISCARDEVENT;
  if( datatype == ERT && !( (lvl1_scaled & bit_ert4x4[evtype]) && (lvl1_live & bit_bbcnarrow) ) )
    return DISCARDEVENT;
  if( datatype == MB && !(lvl1_scaled & bit_bbcnarrow) )
    return DISCARDEVENT;

  anatools::TriggerMode triggermode = anatools::ERT_4x4or;
  if(evtype < 3)
    triggermode = (anatools::TriggerMode)evtype;

  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  /* bbc10cm = 1 for 10cm cut */
  int bbc10cm = 0;
  if( abs(bbc_z) < 10. )
    bbc10cm = 1;

  unsigned ncluster = data_emccontainer->size();
  vector<unsigned> v_used;

  /* Fire ERT on arm 0 (west) or 1 (east) */
  bool FireERT[2] = {};
  if( datatype == ERT )
  {
    int ert_n = data_ert->get_ERThit_N();
    for(int i=0; i<ert_n; i++)
    {
      int arm = data_ert->get_ERTarm(i);
      int trigmode = data_ert->get_ERTtrigmode(i);
      if( arm >= 0 && arm <= 1 && trigmode == evtype )
        FireERT[arm] = true;
    }
  }

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
    int arm = anatools::GetSector(cluster1) / 4;

    /* Require the other arm to be fired for ERT sample */
    int oarm = ( arm==0 ? 1 : 0 );
    if( datatype == ERT && !FireERT[oarm] ) continue;

    if( IsGoodTower(cluster1) &&
        TestPhoton(cluster1, bbc_t0) )
    {
      v_used.push_back(i);

      int sector = anatools::GetSector( cluster1 );
      double photon_pT = anatools::Get_pT( cluster1 );
      bool trig = anatools::PassERT(data_ert, cluster1, triggermode);

      int ert_trig = evtype;
      if( trig )
        ert_trig += 3;

      int isolated[4] = {};
      double econeTrk[4];
      double econeEM = SumEEmcal(cluster1, data_emccontainer, data_tracks, bbc_t0);
      SumPTrack(cluster1, data_tracks, econeTrk);
      for(int ival=0; ival<4; ival++)
      {
        double econe = econeEM + econeTrk[ival];
        if( econe < eratio * cluster1->ecore() )
          isolated[ival] = 1;
      }

      if(sector >= 0 && sector <= 7)
        for(int ival=0; ival<4; ival++)
        {
          int ih = sector + 8*ert_trig + 8*6*bbc10cm + 8*6*2*isolated[ival] + 8*6*2*2*ival;
          h_ert[ih]->Fill(photon_pT);
        }

      for(unsigned j=0; j<ncluster; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
          if( IsGoodTower(cluster2) &&
              TestPhoton(cluster2, bbc_t0) )
          {
            int sector2 = anatools::GetSector( cluster2 );
            if( !anatools::SectorCheck(sector,sector2) ) continue;

            if( cluster2->ecore() > cluster1->ecore() )
            {
              sector = sector2;
              trig = anatools::PassERT(data_ert, cluster2, triggermode);
            }

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);
            double tot_E = cluster1->ecore() + cluster2->ecore();

            int ert_trig = evtype;
            if( trig )
              ert_trig += 3;

            int isolated[4] = {};
            double econe[4];
            SumEPi0(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econe);
            for(int ival=0; ival<4; ival++)
              if( econe[ival] < eratio * tot_E )
                isolated[ival] = 1;

            if(sector >= 0 && sector <= 7)
              for(int ival=0; ival<4; ival++)
              {
                int ih = sector + 8*ert_trig + 8*6*bbc10cm + 8*6*2*isolated[ival] + 8*6*2*2*ival;
                h2_ert_pion[ih]->Fill(tot_pT, minv);
              }
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::FillTowerEnergy(const emcClusterContainer *data_emccontainer, const emcTowerContainer *data_emctwrcontainer,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();

  /* Check trigger */
  const unsigned bit_trig[4] = {0x00000002, 0x00000080, 0x00000040, 0x00000100};
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  int trig = 0;
  for(int i=0; i<4; i++)
    if( lvl1_live & bit_trig[i] )
      trig += pow(2,i);

  /* Get crossing number */
  int crossing = data_triggerlvl1->get_lvl1_clock_cross();
  int crossing_shift = spinpattern->get_crossing_shift();
  int bunch = (crossing + crossing_shift) % 120;
  if(bunch < 0) return DISCARDEVENT;

  int pattern_blue = spinpattern->get_spinpattern_blue(bunch);
  int pattern_yellow = spinpattern->get_spinpattern_yellow(bunch);

  /* Test whether at the beam abort gap or not */
  int collision = 1;
  if( pattern_blue == 10 || pattern_yellow == 10 )
    collision = 0;

  unsigned ncluster = data_emccontainer->size();
  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster = data_emccontainer->getCluster(i);
    if( IsGoodTower(cluster) &&
        cluster->ecore() > eMin )
    {
      double pT = anatools::Get_pT(cluster);
      double tof = cluster->tofcorr();
      int cid = cluster->towerid(0);
      int sector, iypos, izpos;
      EmcIndexer::decodeTowerId(cid, sector, izpos, iypos);  // notice izpos, iypos

      int cut = 0;
      if( cluster->prob_photon() > probMin )
        cut += 1;
      if( abs(tof - bbc_t0) < tofMax )
        cut += 2;
      if( abs(bbc_z) < 10. )
        cut += 4;

      for(int iz=0; iz<7; iz++)
        for(int iy=0; iy<7; iy++)
        {
          int twrid = EmcIndexer::getTowerId(sector, izpos+iz-3, iypos+iy-3);
          emcTowerContent *tower = data_emctwrcontainer->findTower(twrid);
          if(tower && sector >= 0 && sector <= 7)
          {
            double etwr = tower->Energy();
            int ih = sector + 8*collision + 8*2*cut + 8*2*8*trig;
            h3_etwr[ih]->Fill(pT, iy-3., iz-3., etwr);
          }
        }
    }
  }

  return EVENT_OK;
}

int PhotonHistos::FillTrackQuality(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1)
{
  /* Check trigger */
  int evtype = 2;  // ERT_4x4c
  if( !IsEventType(evtype, data_triggerlvl1) )
    return DISCARDEVENT;

  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  int npart = data_tracks->get_npart();
  for(int itrk=0; itrk<npart; itrk++)
  {
    unsigned quality = data_tracks->get_quality(itrk);
    double mom = data_tracks->get_mom(itrk);
    if( quality > 63 || mom != mom ||
        mom < pTrkMin || mom > pTrkMax )
      continue;

    /* phi and zed distributions */
    double dcphi = data_tracks->get_phi(itrk);
    double dczed = data_tracks->get_zed(itrk);
    double dcalpha = data_tracks->get_alpha(itrk);
    int dcns = dczed > 0. ? 0 : 1;
    int dcwe = dcphi < PI/2. ? 0 : 1;
    double dcboard = 0.;
    if( dcwe == 0 )
      dcboard = ( 0.573231 + dcphi - 0.0046 * cos( dcphi + 0.05721 ) ) / 0.01963496;
    else
      dcboard = ( 3.72402 - dcphi + 0.008047 * cos( dcphi + 0.87851 ) ) / 0.01963496;

    /* emcdphi and emcdz distributions */
    int ih = dcns + 2*dcwe + 2*2*quality;
    double dcdphi = data_tracks->get_emcdphi(itrk);
    double dcdz = data_tracks->get_emcdz(itrk);
    h3_dcdphiz[ih]->Fill(dcdphi, dcdz, mom);

    /* alpha-board */
    ih = dcns + 2*dcwe + 2*2*quality;
    h2_alphaboard[ih]->Fill(dcboard, dcalpha);

    /* DC+PC1 live area */
    ih = quality;
    h3_dclive[ih]->Fill(dczed, dcphi, mom);

    int charge = data_tracks->get_charge(itrk);
    double px = data_tracks->get_px(itrk);
    double py = data_tracks->get_py(itrk);
    double pt = sqrt( px*px + py*py );
    double prod = -charge * pt * dcalpha;
    h_prod->Fill(prod);
  }

  unsigned ncluster = data_emccontainer->size();
  for(unsigned iclus=0; iclus<ncluster; iclus++)
  {
    emcClusterContent *cluster = data_emccontainer->getCluster(iclus);
    if( cluster->ecore() > 5. )
    {
      int itrk_match = GetEmcMatchTrack(cluster, data_tracks);
      if( itrk_match >= 0 )
      {
        double dczed = data_tracks->get_zed(itrk_match);
        int dcarm = data_tracks->get_dcarm(itrk_match);
        int dcns = dczed > 0. ? 0 : 1;
        int dcwe = dcarm == 0 ? 1 : 0;

        double pemcx = data_tracks->get_pemcx(itrk_match);
        double pemcy = data_tracks->get_pemcy(itrk_match);
        double pemcz = data_tracks->get_pemcz(itrk_match);
        double emcx = cluster->x();
        double emcy = cluster->y();
        double emcz = cluster->z();
        double dcdphi = atan2(emcy,emcx) - atan2(pemcy,pemcx);
        double dcdz = emcz - pemcz;

        int ih = dcns + 2*dcwe;
        h2_emcdphiz[ih]->Fill(dcdphi, dcdz);
      }
    }
  }

  return EVENT_OK;
}

int PhotonHistos::FillPi0Spectrum(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype)
{
  /* Check trigger */
  if( !IsEventType(evtype, data_triggerlvl1) )
    return DISCARDEVENT;

  anatools::TriggerMode triggermode = anatools::ERT_4x4or;
  if(evtype < 3)
    triggermode = (anatools::TriggerMode)evtype;

  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  /* bbc10cm = 1 for 10cm cut */
  int bbc10cm = 0;
  if( abs(bbc_z) < 10. )
    bbc10cm = 1;

  /* Get crossing number */
  int crossing = data_triggerlvl1->get_lvl1_clock_cross();
  int crossing_shift = spinpattern->get_crossing_shift();
  int evenodd = ( crossing + crossing_shift ) % 2;
  int pattern = GetPattern(crossing);

  unsigned ncluster = data_emccontainer->size();
  vector<unsigned> v_used;

  for(unsigned i=0; i<ncluster; i++)
  {
    v_used.push_back(i);
    for(unsigned j=0; j<ncluster; j++)
      if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
      {
        emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
        emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
        if( IsGoodTower(cluster1) &&
            IsGoodTower(cluster2) &&
            cluster1->ecore() > eMin &&
            cluster2->ecore() > eMin &&
            anatools::GetAsymmetry_E(cluster1, cluster2) < AsymCut )
        {
          int sector1 = anatools::GetSector(cluster1);
          int sector2 = anatools::GetSector(cluster2);
          if( !anatools::SectorCheck(sector1,sector2) ) continue;

          int sector = sector1;
          bool trig = anatools::PassERT(data_ert, cluster1, triggermode);
          if( cluster2->ecore() > cluster1->ecore() )
          {
            sector = sector2;
            trig = anatools::PassERT(data_ert, cluster2, triggermode);
          }

          if( datatype == ERT && !trig )
            continue;

          double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
          double minv = anatools::GetInvMass(cluster1, cluster2);
          double tot_E = cluster1->ecore() + cluster2->ecore();

          int isolated[4] = {};
          double econe[4];
          SumEPi0(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econe);
          for(int ival=0; ival<4; ival++)
            if( econe[ival] < eratio * tot_E )
              isolated[ival] = 1;

          int tof = 0;
          int prob = 0;
          if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax &&
              abs( cluster2->tofcorr() - bbc_t0 ) < tofMax )
          {
            tof = 1;
            if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax/2. &&
                abs( cluster2->tofcorr() - bbc_t0 ) < tofMax/2. )
              tof = 2;
          }
          if( cluster1->prob_photon() > probMin &&
              cluster2->prob_photon() > probMin )
            prob = 1;

          if(sector >= 0 && sector <= 7)
            for(int ival=0; ival<4; ival++)
            {
              int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated[ival] + 8*2*3*2*tof + 8*2*3*2*3*prob + 8*2*3*2*3*2*evtype + 8*2*3*2*3*2*4*bbc10cm + 8*2*3*2*3*2*4*2*ival;
              h2_pion[ih]->Fill(tot_pT, minv);
            }
        } // IsGoodTower and asymmetry cut
      } // j loop
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::FillPhotonSpectrum(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype)
{
  /* Check trigger */
  if( !IsEventType(evtype, data_triggerlvl1) )
    return DISCARDEVENT;

  anatools::TriggerMode triggermode = anatools::ERT_4x4or;
  if(evtype < 3)
    triggermode = (anatools::TriggerMode)evtype;

  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  double bbc_t0 = data_global->getBbcTimeZero();
  if( abs(bbc_z) > 30. ) return DISCARDEVENT;

  /* bbc10cm = 1 for 10cm cut */
  int bbc10cm = 0;
  if( abs(bbc_z) < 10. )
    bbc10cm = 1;

  /* Get crossing number */
  int crossing = data_triggerlvl1->get_lvl1_clock_cross();
  int crossing_shift = spinpattern->get_crossing_shift();
  int evenodd = ( crossing + crossing_shift ) % 2;
  int pattern = GetPattern(crossing);

  unsigned ncluster = data_emccontainer->size();

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
    if( InFiducial(cluster1) &&
        cluster1->ecore() > eMin &&
        !DCChargeVeto(cluster1,data_tracks) )
    {
      int sector = anatools::GetSector(cluster1);
      int part = -1;
      if(sector < 0) part = -1;
      else if(sector < 4) part = 0;
      else if(sector < 6) part = 1;
      else if(sector < 8) part = 2;

      TLorentzVector pE = anatools::Get_pE(cluster1);
      double pT = pE.Pt();
      double eta = 9999.;
      if( pT > epsilon )
        eta = pE.Eta();
      double phi = pE.Phi();
      if(sector >= 4)
      {
        pE.RotateZ(-PI);
        phi = pE.Phi() + PI;
      }

      bool trig = anatools::PassERT(data_ert, cluster1, triggermode);
      if( datatype == ERT && !trig )
        continue;

      int isolated[4] = {};
      double econeTrk[4];
      double econeEM = SumEEmcal(cluster1, data_emccontainer, data_tracks, bbc_t0);
      SumPTrack(cluster1, data_tracks, econeTrk);
      for(int ival=0; ival<4; ival++)
      {
        double econe = econeEM + econeTrk[ival];
        if( econe < eratio * cluster1->ecore() )
          isolated[ival] = 1;
      }

      int tof = 0;
      int prob = 0;
      if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax )
      {
        tof = 1;
        if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax/2. )
          tof = 2;
      }
      if( cluster1->prob_photon() > probMin )
        prob = 1;
      if(tof == 0 || prob != 1)
        continue;

      if( evtype == 2 && part >= 0 &&
          pT > 5. && pT < 10. &&
          TestPhoton(cluster1, bbc_t0) )
        for(int ival=0; ival<4; ival++)
        {
          int ih = part + 3*isolated[ival] + 3*2*ival;
          h2_eta_phi[ih]->Fill(eta, phi);
        }

      if(sector >= 0 && sector <= 7)
        for(int ival=0; ival<4; ival++)
        {
          int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated[ival] + 8*2*3*2*evtype + 8*2*3*2*4*bbc10cm + 8*2*3*2*4*2*ival + 8*2*3*2*4*2*4*(tof-1);
          h_1photon[ih]->Fill(pT);
        }

      for(unsigned j=0; j<ncluster; j++)
        if(j != i)
        {
          emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
          if( !IsGoodTower(cluster2) ||
              cluster2->ecore() < eMin ||
              DCChargeVeto(cluster2,data_tracks) )
            continue;
          double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
          double minv = anatools::GetInvMass(cluster1, cluster2);

          int isoboth[4] = {};
          double econeTrk2[4];
          double econeEM2 = SumEEmcal(cluster2, data_emccontainer, data_tracks, bbc_t0);
          SumPTrack(cluster2, data_tracks, econeTrk2);
          for(int ival=0; ival<4; ival++)
          {
            double econe2 = econeEM2 + econeTrk2[ival];
            if( isolated[ival] && econe2 < eratio * cluster2->ecore() )
              isoboth[ival] = 1;
          }

          int isopair[5] = {};
          double econeEMPair1, econeEMPair2;
          SumEEmcal(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econeEMPair1, econeEMPair2);
          for(int ival=0; ival<4; ival++)
          {
            double econePair1 = econeEMPair1 + econeTrk[ival];
            double econePair2 = econeEMPair2 + econeTrk2[ival];
            if( econePair1 < eratio * cluster1->ecore() &&
                econePair2 < eratio * cluster2->ecore() )
              isopair[ival] = 1;
          }

          int tof = 0;
          int prob = 0;
          if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax &&
              abs( cluster2->tofcorr() - bbc_t0 ) < tofMax )
          {
            tof = 1;
            if( abs( cluster1->tofcorr() - bbc_t0 ) < tofMax/2. &&
                abs( cluster2->tofcorr() - bbc_t0 ) < tofMax/2. )
              tof = 2;
          }
          if( cluster1->prob_photon() > probMin && 
              cluster2->prob_photon() > probMin )
            prob = 1;
          if(tof == 0 || prob != 1)
            continue;

          if(sector >= 0 && sector <= 7)
            for(int ival=0; ival<4; ival++)
            {
              int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isoboth[ival] + 8*2*3*2*isopair[ival] + 8*2*3*2*2*evtype + 8*2*3*2*2*4*bbc10cm + 8*2*3*2*2*4*2*ival + 8*2*3*2*2*4*2*4*(tof-1);
              h2_2photon[ih]->Fill(pT, minv);
              h2_2photon2pt[ih]->Fill(tot_pT, minv);
            }
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::End(PHCompositeNode *topNode)
{
  hm->dumpHistos(outFile);
  delete spinpattern;
  delete emcrecalib;
  delete emcrecalib_sasha;
  delete hm;

  return EVENT_OK;
}

void PhotonHistos::SelectMB()
{
  datatype = MB;
  return;
}

void PhotonHistos::SelectERT()
{
  datatype = ERT;
  return;
}

void PhotonHistos::BookHistograms()
{
  /* Initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");

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
    h_events = new TH1F("h_events", "Events counter", 6,-0.5,5.5);
    h_events->GetXaxis()->SetBinLabel(1, "ert_a_10cm");
    h_events->GetXaxis()->SetBinLabel(2, "ert_b_10cm");
    h_events->GetXaxis()->SetBinLabel(3, "ert_c_10cm");
    h_events->GetXaxis()->SetBinLabel(4, "ert_a_30cm");
    h_events->GetXaxis()->SetBinLabel(5, "ert_b_30cm");
    h_events->GetXaxis()->SetBinLabel(6, "ert_c_30cm");
  }
  else if( datatype == MB )
  {
    h_events = new TH1F("h_events", "Events counter", 5,-0.5,4.5);
    h_events->GetXaxis()->SetBinLabel(1, "bbc_narrow");
    h_events->GetXaxis()->SetBinLabel(2, "bbc_narrow_10cm");
    h_events->GetXaxis()->SetBinLabel(3, "bbc_narrow_10cm_ert_c");
    h_events->GetXaxis()->SetBinLabel(4, "bbc_narrow_30cm");
    h_events->GetXaxis()->SetBinLabel(5, "bbc_narrow_30cm_ert_c");
  }
  hm->registerHisto(h_events);

  // ih = sector + 8*raw < 8*2
  for(int ih=0; ih<nh_calib; ih++)
  {
    /* ToF calibration */
    h2_tof[ih] = new TH2F(Form("h2_tof_%d",ih), "ToF;p_{T} [GeV];tof [ns];", npT,pTbin, 500,-50.,50.);
    hm->registerHisto(h2_tof[ih]);
    /* Storing invariant mass of photon pairs in different sectors and pT bins.
     * Require both photons to be in same sector.
     * Used to check sector-by-sector EMCal energy calibration.
     */
    h2_minv[ih] = new TH2F(Form("h2_minv_%d",ih), "Photon pair invariant mass;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_minv[ih]);
  }

  // ih = sector + 8*bbc_trig < 8*2
  for(int ih=0; ih<nh_bbc; ih++)
  {
    /* BBC trigger efficiency for photon */
    h_bbc[ih] = new TH1F(Form("h_bbc_%d",ih), "BBC efficiency;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_bbc[ih]);
    /* BBC Trigger efficiency for pion */
    h2_bbc_pion[ih] = new TH2F(Form("h2_bbc_pion_%d",ih), "BBC efficiency;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_bbc_pion[ih]);
    /* Store photon bg information */
  }

  // ih = sector + 8*ert_trig + 8*6*bbc10cm + 8*6*2*isolated[ival] + 8*6*2*2*ival < 8*6*2*2*4
  for(int ih=0; ih<nh_ert; ih++)
  {
    /* ERT trigger efficiency for photon */
    h_ert[ih] = new TH1F(Form("h_ert_%d",ih), "ERT efficiency;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_ert[ih]);
    /* ERT Trigger efficiency for pion */
    h2_ert_pion[ih] = new TH2F(Form("h2_ert_pion_%d",ih), "ERT efficiency;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_ert_pion[ih]);
  }

  /* Tower energy distribution */
  // ih = sector + 8*collision + 8*2*cut + 8*2*8*trig < 8*2*8*16
  //for(int ih=0; ih<nh_etwr; ih++)
  //{
  //  h3_etwr[ih] = new TH3F(Form("h3_etwr_%d",ih), "Tower energy;p_{T} [GeV];iy;iz;",
  //      npT,0.,0., 7,-3.5,3.5, 7,-3.5,3.5);
  //  h3_etwr[ih]->GetXaxis()->Set(npT, pTbin);
  //  h3_etwr[ih]->Sumw2();
  //  hm->registerHisto(h3_etwr[ih]);
  //}

  /* Track quality study */
  // ih = dcns + 2*dcwe + 2*2*quality < 2*2*64
  //for(int ih=0; ih<nh_dcpartqual; ih++)
  //{
  //  h3_dcdphiz[ih] = new TH3F(Form("h3_dcdphiz_%d",ih), "EMCal and DC phi and z deviation;dphi [rad];dz [cm];mom [GeV];",
  //      100,-0.2,0.2, 120,-60.,60., 30,0.,15.);
  //  hm->registerHisto(h3_dcdphiz[ih]);
  //  h2_alphaboard[ih] = new TH2F(Form("h2_alphaboard_%d",ih), "DC #alpha-board;board;#alpha;", 80,0.,80., 120,-0.6,0.6);
  //  hm->registerHisto(h2_alphaboard[ih]);
  //}

  // ih = quality < 64
  //for(int ih=0; ih<nh_dcquality; ih++)
  //{
  //  h3_dclive[ih] = new TH3F(Form("h3_dclive_%d",ih), "DC zed and phi distribution;zed [cm];phi [rad];mom [GeV]",
  //      200,-100.,100., 50,-1.,4., 30,0.,15.);
  //  hm->registerHisto(h3_dclive[ih]);
  //}

  //h_prod = new TH1F("h_prod", "Minus charge times pT times alpha;-charge*pT*#alpha [GeV]", 200,0.,0.2);
  //hm->registerHisto(h_prod);

  // ih = dcns + 2*dcwe < 2*2
  //for(int ih=0; ih<nh_dcpart; ih++)
  //{
  //  h2_emcdphiz[ih] = new TH2F(Form("h2_emcdphiz_%d",ih), "EMCal and DC phi and z deviation;dphi [rad];dz [cm];",
  //      100,-0.2,0.2, 120,-60.,60.);
  //  hm->registerHisto(h2_emcdphiz[ih]);
  //}

  /* Eta and phi distribution */
  // ih = part + 3*isolated[ival] + 3*2*ival < 3*2*4
  for(int ih=0; ih<nh_eta_phi; ih++)
  {
    h2_eta_phi[ih] = new TH2F(Form("h2_eta_phi_%d",ih), "#eta and #phi distribution;#eta;#phi;", neta,etabin[ih%3/2], nphi,phibin);
    hm->registerHisto(h2_eta_phi[ih]);
  }

  /* Store pi0 information */
  // ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated[ival] + 8*2*3*2*tof + 8*2*3*2*3*prob + 8*2*3*2*3*2*evtype + 8*2*3*2*3*2*4*bbc10cm + 8*2*3*2*3*2*4*2*ival < 8*2*3*2*3*2*4*2*4
  for(int ih=0; ih<nh_pion; ih++)
  {
    h2_pion[ih] = new TH2F(Form("h2_pion_%d",ih), "#pi^{0} spectrum;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_pion[ih]);
  }

  /* Store single photon information */
  // ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated[ival] + 8*2*3*2*evtype + 8*2*3*2*4*bbc10cm + 8*2*3*2*4*2*ival + 8*2*3*2*4*2*4*(tof-1) < 8*2*3*2*4*2*4*2
  for(int ih=0; ih<nh_1photon; ih++)
  {
    h_1photon[ih] = new TH1F(Form("h_1photon_%d",ih), "Single photon spectrum;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_1photon[ih]);
  }

  /* Store two photons information */
  // ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isoboth[ival] + 8*2*3*2*isopair[ival] + 8*2*3*2*2*evtype + 8*2*3*2*2*4*bbc10cm + 8*2*3*2*2*4*2*ival + 8*2*3*2*2*4*2*4*(tof-1) < 8*2*3*2*2*4*2*4*2
  for(int ih=0; ih<nh_2photon; ih++)
  {
    h2_2photon[ih] = new TH2F(Form("h2_2photon_%d",ih), "Two photons spectrum;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    h2_2photon2pt[ih] = (TH2*)h2_2photon[ih]->Clone(Form("h2_2photon2pt_%d",ih));
    hm->registerHisto(h2_2photon[ih]);
    hm->registerHisto(h2_2photon2pt[ih]);
  }

  return;
}

double PhotonHistos::SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
    const PHCentralTrack *data_tracks, double bbc_t0)
{ 
  /* Sum up all energy in cone around particle without one cluster */
  double econe = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus2 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to 'reference' particle
     * or on bad towers or lower than energy threshold */
    if( clus2->id() == cluster->id() ||
        IsBadTower(clus2) ||
        abs( clus2->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus2->ecore() < eClusMin )
      continue;

    /* 3 sigma charge veto */
    if( DCChargeVeto(clus2,data_tracks) )
      continue;

    /* Get cluster vector */
    TLorentzVector pE_part2 = anatools::Get_pE(clus2);
    if( pE_part2.Pt() < epsilon ) continue;
    TVector2 v2_part2 = pE_part2.EtaPhiVector();

    /* Check if cluster within cone */
    TVector2 v2_diff(v2_part2 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
      econe += clus2->ecore();
  }

  return econe;
}

void PhotonHistos::SumEEmcal(const emcClusterContent *cluster1, const emcClusterContent *cluster2,
    const emcClusterContainer *cluscont, const PHCentralTrack *data_tracks, double bbc_t0, double &econe1, double &econe2)
{ 
  /* Sum up all energy in cone around particle without two clusters */
  econe1 = 0.;
  econe2 = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref1 = anatools::Get_pE(cluster1);
  TLorentzVector pE_pref2 = anatools::Get_pE(cluster2);
  if( pE_pref1.Pt() < epsilon || pE_pref2.Pt() < epsilon ) return;
  TVector2 v2_pref1 = pE_pref1.EtaPhiVector();
  TVector2 v2_pref2 = pE_pref2.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus3 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to any of the two 'reference' particles
     * or on bad towers or lower than energy threshold */
    if( clus3->id() == cluster1->id() ||
        clus3->id() == cluster2->id() ||
        IsBadTower(clus3) ||
        abs( clus3->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus3->ecore() < eClusMin )
      continue;

    /* 3 sigma charge veto */
    if( DCChargeVeto(clus3,data_tracks) )
      continue;

    /* Get cluster vector */
    TLorentzVector pE_part3 = anatools::Get_pE(clus3);
    if( pE_part3.Pt() < epsilon ) continue;
    TVector2 v2_part3 = pE_part3.EtaPhiVector();

    /* Check if cluster within cone */
    TVector2 v2_diff(v2_part3 - v2_pref1);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
      econe1 += clus3->ecore();
    v2_diff = v2_part3 - v2_pref2;
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
      econe2 += clus3->ecore();
  }

  return;
}

void PhotonHistos::SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks, double econe[])
{ 
  /* Sum up all energy in cone around particle */
  for(int ival=0; ival<4; ival++)
    econe[ival] = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int npart = data_tracks->get_npart();

  for(int i=0; i<npart; i++)
  {
    int quality = data_tracks->get_quality(i);
    double px = data_tracks->get_px(i);
    double py = data_tracks->get_py(i);
    double pz = data_tracks->get_pz(i);
    double mom = data_tracks->get_mom(i);
    if( !TMath::Finite(px+py+pz+mom) )
      continue;

    /* Test if track passes the momentum cuts */
    if( mom < pTrkMin || mom > pTrkMax )
      continue;

    /* Get track vector */
    TVector3 v3_track(px, py, pz);
    if( v3_track.Pt() < epsilon ) continue;
    TVector2 v2_track = v3_track.EtaPhiVector();

    /* Add track energy from clusters within cone range */
    TVector2 v2_diff(v2_track - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
    {
      econe[0] += mom;
      if( quality > 3 )
        econe[1] += mom;
      if( quality == 63 || quality == 31 || quality == 51 )
        econe[2] += mom;
    }
  }

  return;
}

void PhotonHistos::SumEPi0(const emcClusterContent *cluster1, const emcClusterContent *cluster2,
    const emcClusterContainer *cluscont, const PHCentralTrack *data_tracks, double bbc_t0, double econe[])
{ 
  /* Sum up all energy in cone around pi0 */
  double econeEM = 0.;
  for(int ival=0; ival<4; ival++)
    econe[ival] = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster1) + anatools::Get_pE(cluster2);
  if( pE_pref.Pt() < epsilon ) return;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus3 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to any of the two 'reference' particles
     * or on bad towers or lower than energy threshold */
    if( clus3->id() == cluster1->id() ||
        clus3->id() == cluster2->id() ||
        IsBadTower(clus3) ||
        abs( clus3->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus3->ecore() < eClusMin )
      continue;

    /* 3 sigma charge veto */
    if( DCChargeVeto(clus3,data_tracks) )
      continue;

    /* Get cluster vector */
    TLorentzVector pE_part3 = anatools::Get_pE(clus3);
    if( pE_part3.Pt() < epsilon ) continue;
    TVector2 v2_part3 = pE_part3.EtaPhiVector();

    /* Check if cluster within cone */
    TVector2 v2_diff(v2_part3 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
      econeEM += clus3->ecore();
  }

  int npart = data_tracks->get_npart();

  for(int i=0; i<npart; i++)
  {
    int quality = data_tracks->get_quality(i);
    double px = data_tracks->get_px(i);
    double py = data_tracks->get_py(i);
    double pz = data_tracks->get_pz(i);
    double mom = data_tracks->get_mom(i);
    if( !TMath::Finite(px+py+pz+mom) )
      continue;

    /* Test if track passes the momentum cuts */
    if( mom < pTrkMin || mom > pTrkMax )
      continue;

    /* Get track vector */
    TVector3 v3_track(px, py, pz);
    if( v3_track.Pt() < epsilon ) continue;
    TVector2 v2_track = v3_track.EtaPhiVector();

    /* Add track energy from clusters within cone range */
    TVector2 v2_diff(v2_track - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
    {
      econe[0] += mom;
      if( quality > 3 )
        econe[1] += mom;
      if( quality == 63 || quality == 31 || quality == 51 )
        econe[2] += mom;
    }
  }

  for(int ival=0; ival<4; ival++)
    econe[ival] += econeEM;

  return;
}

bool PhotonHistos::IsEventType(const int evtype, const TrigLvl1 *data_triggerlvl1)
{
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  const unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();

  if( datatype == ERT )
  {
    if( evtype < 3 && (lvl1_scaled & bit_ert4x4[evtype]) && (lvl1_live & bit_bbcnarrow) )
      return true;
    else if( evtype == 3 && (lvl1_live & bit_ert4x4[evtype]) )
      return true;
  }
  else if( datatype == MB )
  {
    if( evtype == 0 && (lvl1_scaled & bit_bbcnarrow) )
      return true;
  }

  return false;
}

bool PhotonHistos::TestPhoton(const emcClusterContent *cluster, double bbc_t0)
{
  if( cluster->ecore() > eMin &&
      abs( cluster->tofcorr() - bbc_t0 ) < tofMax &&
      cluster->prob_photon() > probMin )
    return true;
  else
    return false;
}

bool PhotonHistos::DCChargeVeto(const emcClusterContent *cluster, const PHCentralTrack *data_tracks)
{
  /* 3 sigma charge veto */
  int itrk_match = GetEmcMatchTrack(cluster, data_tracks);
  if( itrk_match >= 0 )
    return true;
  else 
    return false;
}

bool PhotonHistos::InFiducial(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_sasha[sector][iypos][izpos] == 0 )
    return true;
  else
    return false;
}

bool PhotonHistos::IsGoodTower(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_sasha[sector][iypos][izpos] == 0 ||
      tower_status_sasha[sector][iypos][izpos] == 30 )
    return true;
  else
    return false;
}

bool PhotonHistos::IsBadTower(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  //if( tower_status_nils[sector][iypos][izpos] >= 50 )
  if( tower_status_sasha[sector][iypos][izpos] > 0 &&
      tower_status_sasha[sector][iypos][izpos] < 20 )
    return true;
  else
    return false;
}

int PhotonHistos::GetPattern(int crossing)
{
  int pattern = 0;
  int crossing_shift = spinpattern->get_crossing_shift();
  int bunch = (crossing + crossing_shift) % 120;
  if(bunch >= 0)
  {
    int pattern_blue = spinpattern->get_spinpattern_blue(bunch);
    int pattern_yellow = spinpattern->get_spinpattern_yellow(bunch);
    pattern = pattern_blue * pattern_yellow;  
  }
  if( abs(pattern) > 1 ) pattern = 0;
  pattern += 1;

  return pattern;
}

int PhotonHistos::GetEmcMatchTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks)
{
  int itrk_match = -1;
  double dzmin = 9999.;
  double mommax = 0.;

  TVector3 v3_cluster(cluster->x(), cluster->y(), cluster->z());

  int npart = data_tracks->get_npart();
  for(int itrk=0; itrk<npart; itrk++)
  {
    TVector3 v3_track(data_tracks->get_pemcx(itrk), data_tracks->get_pemcy(itrk), data_tracks->get_pemcz(itrk));
    double dphi = abs((v3_track-v3_cluster).Phi());
    double dz = abs((v3_track-v3_cluster).Z());
    double mom = data_tracks->get_mom(itrk);
    if( dphi > 0.06 || dz > 13. || !TMath::Finite(mom) )
      continue;

    if( itrk_match != -1 )
    {
      if( dz < 8. && dz < dzmin )
      {
        itrk_match = itrk;
        dzmin = dz;
        mommax = mom;
      }
      else if( dzmin >= 8. && mom > mommax )
      {
        itrk_match = itrk;
        dzmin = dz;
        mommax = mom;
      }
    }
    else
    {
      itrk_match = itrk;
      dzmin = dz;
      mommax = mom;
    }
  }

  return itrk_match;
}

void PhotonHistos::EMCRecalibSetup()
{
  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(0);

  string file_ecal_run = toad_loader->location("Run13pp_RunbyRun_Calib.dat");
  string file_tofmap = toad_loader->location("Run13pp510_EMC_TOF_Correction.root");

  emcrecalib->SetEnergyCorrectionFile( file_ecal_run );
  emcrecalib->SetTofCorrectionFile( file_tofmap );

  string _file_ecal = toad_loader->location("ecorr_run13pp500gev.txt");
  string _file_ecal_run = toad_loader->location("ecorr_run_run13pp500gev.txt");
  string _file_tcal = toad_loader->location("tcorr_run13pp500gev.txt");

  emcrecalib_sasha->anaGetCorrCal( _file_ecal.c_str() );
  emcrecalib_sasha->anaGetCorrCal_run( _file_ecal_run.c_str() );
  emcrecalib_sasha->anaGetCorrTof( _file_tcal.c_str() );

  delete toad_loader;
  return;
}

void PhotonHistos::ReadTowerStatus(const string &filename)
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
    tower_status_nils[sector][biny][binz] = status;
  }

  //cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void PhotonHistos::ReadSashaWarnmap(const string &filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int ich = 0;
  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> ich >> status )
  {
    /* Attention!! I use my indexing for warn map in this program!!! */
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }

    /* Get tower location */
    anatools::TowerLocation(ich, sector, biny, binz);

    /* Count tower with bad status for PbSc and PbGl */
    if ( status > 0 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status_sasha[sector][biny][binz] = status;

    /* Mark edge towers */
    if( anatools::Edge_cg(sector, biny, binz) &&
        tower_status_sasha[sector][biny][binz] == 0 )
      tower_status_sasha[sector][biny][binz] = 20;
    /* Mark fiducial arm */
    if( anatools::ArmEdge_cg(sector, biny, binz) &&
        tower_status_sasha[sector][biny][binz] == 0 )
      tower_status_sasha[sector][biny][binz] = 30;
  }

  //cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void PhotonHistos::UpdateSpinPattern(SpinDBContent &spin_cont)
{
  /* Get spin info */
  int runnumber = spin_cont.GetRunNumber();
  int qa_level = spin_cont.GetQALevel();
  int fillnumber = spin_cont.GetFillNumber();
  int badrunqa = spin_cont.GetBadRunFlag();
  int crossing_shift = spin_cont.GetCrossingShift();

  double pb, pbstat, pbsyst;
  double py, pystat, pysyst;

  spin_cont.GetPolarizationBlue(1, pb, pbstat, pbsyst);
  spin_cont.GetPolarizationYellow(1, py, pystat, pysyst);

  int badbunch[120];
  int spinpattern_blue[120];
  int spinpattern_yellow[120];

  long long bbc_narrow[120];
  long long bbc_wide[120];
  long long zdc_narrow[120];
  long long zdc_wide[120];

  for(int i=0; i<120; i++)
  {
    badbunch[i] = spin_cont.GetBadBunchFlag(i);
    spinpattern_blue[i] = spin_cont.GetSpinPatternBlue(i);
    spinpattern_yellow[i] = spin_cont.GetSpinPatternYellow(i);

    bbc_narrow[i] = spin_cont.GetScalerBbcVertexCut(i);
    bbc_wide[i] = spin_cont.GetScalerBbcNoCut(i);
    zdc_narrow[i] = spin_cont.GetScalerZdcNarrow(i);
    zdc_wide[i] = spin_cont.GetScalerZdcWide(i);
  }

  /* Update spinpattern */
  int current_qa_level = spinpattern->get_qa_level();
  if(current_qa_level < qa_level)
  {
    spinpattern->Reset();

    spinpattern->set_runnumber(runnumber);
    spinpattern->set_qa_level(qa_level);
    spinpattern->set_fillnumber(fillnumber);
    spinpattern->set_badrunqa(badrunqa);
    spinpattern->set_crossing_shift(crossing_shift);

    spinpattern->set_pb(pb);
    spinpattern->set_pbstat(pbstat);
    spinpattern->set_pbsyst(pbsyst);
    spinpattern->set_py(py);
    spinpattern->set_pystat(pystat);
    spinpattern->set_pysyst(pysyst);

    for(int i=0; i<120; i++)
    {
      spinpattern->set_badbunch(i, badbunch[i]);
      spinpattern->set_spinpattern_blue(i, spinpattern_blue[i]);
      spinpattern->set_spinpattern_yellow(i, spinpattern_yellow[i]);

      spinpattern->set_bbc_narrow(i, bbc_narrow[i]);
      spinpattern->set_bbc_wide(i, bbc_wide[i]);
      spinpattern->set_zdc_narrow(i, zdc_narrow[i]);
      spinpattern->set_zdc_wide(i, zdc_wide[i]);
    }
  }

  return;
}
