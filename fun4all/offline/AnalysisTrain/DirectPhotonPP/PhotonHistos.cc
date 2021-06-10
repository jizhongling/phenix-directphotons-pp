#include "PhotonHistos.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"
#include "AnaToolsTrigger.h"

#include "EmcLocalRecalibrator.h"
#include "EmcLocalRecalibratorSasha.h"
#include "EMCWarnmapChecker.h"
#include "DCDeadmapChecker.h"
#include "SpinPattern.h"

#include <RunHeader.h>
#include <SpinDBOutput.hh>
#include <SpinDBContent.hh>

#include <PHGlobal.h>
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
const double tofMaxIso = 10.;
const double tofMaxALL = 15.;
const double AsymCut = 0.8;

/* Some cuts for isolation cut */
const double eClusMin = 0.3;
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

/* pT bins for ALL */
double PhotonHistos::pTbin_pol[] = { 2.0,
  2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
  10.0, 12.0, 15.0, 20.0, 30.0 };

PhotonHistos::PhotonHistos(const string &name, const char *filename) :
  SubsysReco(name),
  emcrecalib(nullptr),
  emcrecalib_sasha(nullptr),
  emcwarnmap(nullptr),
  dcdeadmap(nullptr),
  spinpattern(nullptr),
  runnumber(0),
  fillnumber(0),
  hm(nullptr)
{
  datatype = ERT;

  outFile = "PhotonHistos-";
  outFile.append(filename);

  h_events = nullptr;
  h_prod = nullptr;
  for(int ih=0; ih<nh_calib; ih++)
  {
    h2_tof[ih] = nullptr;
    h2_minv[ih] = nullptr;
  }
  for(int ih=0; ih<nh_bbc; ih++)
  {
    h_bbc[ih] = nullptr;
    h2_bbc_pion[ih] = nullptr;
  }
  for(int ih=0; ih<nh_ertsm; ih++)
    h2_ertsm[ih] = nullptr;
  for(int ih=0; ih<nh_ert; ih++)
  {
    h_ert[ih] = nullptr;
    h2_ert_pion[ih] = nullptr;
  }
  for(int ih=0; ih<nh_dcpartqual; ih++)
  {
    h3_dcdphiz[ih] = nullptr;
    h2_alphaboard[ih] = nullptr;
  }
  for(int ih=0; ih<nh_dcgood; ih++)
    h3_dclive[ih] = nullptr;
  for(int ih=0; ih<nh_pion; ih++)
    h2_pion[ih] = nullptr;
  for(int ih=0; ih<nh_pion_pol; ih++)
    h2_pion_pol[ih] = nullptr;
  for(int ih=0; ih<nh_eta_phi; ih++)
    h2_eta_phi[ih] = nullptr;
  for(int ih=0; ih<nh_1photon; ih++)
    h_1photon[ih] = nullptr;
  for(int ih=0; ih<nh_1photon_pol; ih++)
    h_1photon_pol[ih] = nullptr;
  for(int ih=0; ih<nh_photon_bunch; ih++)
    h_photon_bunch[ih] = nullptr;
  for(int ih=0; ih<nh_2photon; ih++)
  {
    h2_2photon[ih] = nullptr;
    h2_2photon2pt[ih] = nullptr;
  }
  for(int ih=0; ih<nh_2photon_pol; ih++)
  {
    h2_2photon_pol[ih] = nullptr;
    h2_2photon2pt_pol[ih] = nullptr;
  }
  for(int ih=0; ih<nh_mul_pion; ih++)
  {
    h2_mul_pion_sig[ih] = nullptr;
    h2_mul_pion_bg[ih] = nullptr;
  }
  for(int ih=0; ih<nh_mul_photon; ih++)
    h2_mul_photon[ih] = nullptr;
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

  /* Initialize warnmap checker */
  emcwarnmap = new EMCWarnmapChecker();
  if(!emcwarnmap)
  {
    cerr << "No emcwarnmap" << endl;
    exit(1);
  }

  /* Initialize DC deadmap checker */
  dcdeadmap = new DCDeadmapChecker();
  if(!dcdeadmap)
  {
    cerr << "No dcdeadmap" << endl;
    exit(1);
  }

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
  emcrecalib->ReadEnergyCorrection(runnumber);
  emcrecalib->ReadTofCorrection(fillnumber);

  /* Set DC deadmap index by runnumber */
  dcdeadmap->SetMapByRunnumber(runnumber);

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
  emcClusterContainer *data_emccontainer[2];
  for(int i=0; i<2; i++)
    data_emccontainer[i] = data_emccontainer_raw->clone();
  emcrecalib_sasha->ApplyClusterCorrection( runnumber, data_emccontainer[0] );
  emcrecalib->ApplyClusterCorrection( data_emccontainer[1] );

  /* Event counts */
  FillEventCounts(data_global, data_triggerlvl1);

  /* Store ToF information for cluster as calibration check */
  FillClusterTofSpectrum(data_emccontainer_raw, data_global, "raw");
  FillClusterTofSpectrum(data_emccontainer[0], data_global);

  /* Analyze pi0s events for crosscheck */
  FillPi0InvariantMass(data_emccontainer_raw, data_global, "raw");
  FillPi0InvariantMass(data_emccontainer[0], data_global);

  /* Count events to calculate BBC efficiency */
  FillBBCEfficiency(data_emccontainer[0], data_triggerlvl1);

  /* Count events to calculate ERT efficiency */
  FillERTEfficiency(data_emccontainer[0], data_tracks, data_global, data_triggerlvl1, data_ert);

  /* Fill DC track quality information */
  FillTrackQuality(data_emccontainer[0], data_tracks, data_global, data_triggerlvl1, data_ert);

  /* Analyze photon for pi0 event */
  for(int i=0; i<2; i++)
    FillPi0Spectrum(i, data_emccontainer[i], data_tracks, data_global, data_triggerlvl1, data_ert);

  /* Analyze photon for direct photon event */
  for(int i=0; i<2; i++)
    FillPhotonSpectrum(i, data_emccontainer[i], data_tracks, data_global, data_triggerlvl1, data_ert);

  /* Clean up */
  for(int i=0; i<2; i++)
    delete data_emccontainer[i];

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
    for(int iert=0; iert<3; iert++)
      if( lvl1_scaled & bit_ert4x4[iert] )
      {
        if( BBC10cm(data_global, data_triggerlvl1) )
          h_events->Fill(Form("ert_%c_10cm", 97+iert), 1.);
        if( BBC30cm(data_global, data_triggerlvl1) )
          h_events->Fill(Form("ert_%c_30cm", 97+iert), 1.);
      }
  }

  /* Fill BBC trigger counts */
  else if( datatype == MB )
  {
    if( lvl1_scaled & bit_bbcnarrow )
    {
      h_events->Fill("bbc_narrow", 1.);
      if( fabs(bbc_z) < 10. )
      {
        h_events->Fill("bbc_narrow_10cm", 1.);
        if( lvl1_live & bit_ert4x4[2] )
          h_events->Fill("bbc_narrow_10cm_ert_c", 1.);
      }
    if( lvl1_scaled & bit_bbcnovtx )
    {
      h_events->Fill("bbc_novtx", 1.);
      if( fabs(bbc_z) < 30. )
      {
        h_events->Fill("bbc_novtx_30cm", 1.);
        if( lvl1_live & bit_ert4x4[2] )
          h_events->Fill("bbc_novtx_30cm_ert_c", 1.);
      }
    }
    }
  }

  return EVENT_OK;
}

int PhotonHistos::FillClusterTofSpectrum(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const string &quali)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;
  double bbc_t0 = data_global->getBbcTimeZero();

  unsigned ncluster = data_emccontainer->size();

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster = data_emccontainer->getCluster(i);
    if( emcwarnmap->IsGoodTower(cluster) &&
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
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;
  double bbc_t0 = data_global->getBbcTimeZero();

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
        if( emcwarnmap->IsGoodTower(cluster1) &&
            emcwarnmap->IsGoodTower(cluster2) &&
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
    if( emcwarnmap->IsGoodTower(cluster1) &&
        cluster1->ecore() > eMin &&
        cluster1->prob_photon() > probMin )
    {
      v_used.push_back(i);

      int part = anatools::GetPart(cluster1);
      double pT = anatools::Get_pT(cluster1);

      int bbc_trig = (lvl1_live & bit_bbcnovtx) ? 1 : 0;

      if(part >= 0)
      {
        int ih = part + 3*bbc_trig;
        h_bbc[ih]->Fill(pT);
      }

      for(unsigned j=0; j<ncluster; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
          if( emcwarnmap->IsGoodTower(cluster2) &&
              cluster2->ecore() > eMin &&
              cluster2->prob_photon() > probMin &&
              anatools::GetAsymmetry_E(cluster1, cluster2) < AsymCut )
          {
            int part2 = anatools::GetPart( cluster2 );
            if( part2 != part ) continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);

            if(part >= 0)
            {
              int ih = part + 3*bbc_trig;
              h2_bbc_pion[ih]->Fill(tot_pT, minv);
            }
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::FillERTEfficiency(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert)
{
  /* Get event global parameters */
  if( !BBC10cm(data_global, data_triggerlvl1) )
    return DISCARDEVENT;
  double bbc_t0 = data_global->getBbcTimeZero();

  unsigned ncluster = data_emccontainer->size();
  vector<unsigned> v_used;

  /* Fire ERT on arm 0 (west) or 1 (east) */
  bool FireERT[3][2] = {};  // FireERT[evtype][arm]
  if( datatype == ERT )
  {
    int ert_n = data_ert->get_ERThit_N();
    for(int i=0; i<ert_n; i++)
    {
      int arm = data_ert->get_ERTarm(i);
      int trigmode = data_ert->get_ERTtrigmode(i);
      if( trigmode >= 0 && trigmode <= 2 &&
          arm >= 0 && arm <= 1 )
        FireERT[trigmode][arm] = true;
    }
  }

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
    int sector = anatools::GetSector(cluster1);
    int part = anatools::GetPart(cluster1);
    int arm = sector / 4;
    double photon_pT = anatools::Get_pT(cluster1);

    int ert_trig[3] = {};  // ert_trig[evtype]
    for(int evtype=0; evtype<3; evtype++)
      ert_trig[evtype] = anatools::PassERT(data_ert, cluster1, (anatools::TriggerMode)evtype) ? 1 : 0;

    /* Require the other arm to be fired for ERT sample */
    int oarm = ( arm==0 ? 1 : 0 );

    /* ERT efficiency for each of super modules */
    for(int evtype=0; evtype<3; evtype++)
      if( datatype == MB || FireERT[evtype][oarm] )
      {
        int ih = sector + 8*ert_trig[evtype] + 8*2*evtype;
        h2_ertsm[ih]->Fill(photon_pT, (double)anatools::GetSM(cluster1));
      }

    /* For pi0 consider good towers */
    if( emcwarnmap->IsGoodTower(cluster1) &&
        TestPhoton(cluster1, bbc_t0) )
    {
      v_used.push_back(i);

      /* For photon consider fiducial towers */
      if( emcwarnmap->InFiducial(cluster1) )
      {
        for(int checkmap=0; checkmap<2; checkmap++)
        {
          dcdeadmap->Checkmap(checkmap);

          int isolated[3] = {};
          double econeEM[2], econeTrk[3];
          SumEEmcal(cluster1, data_emccontainer, data_tracks, bbc_t0, econeEM);
          SumPTrack(cluster1, data_tracks, econeTrk);
          for(int ival=0; ival<3; ival++)
          {
            double econe = (ival == 0 ? econeEM[0] : econeEM[1]) + econeTrk[ival];
            if( econe < eratio * cluster1->ecore() )
              isolated[ival] = 1;
          }

          for(int evtype=0; evtype<3; evtype++)
            if( part >= 0 && (datatype == MB || FireERT[evtype][oarm]) )
              for(int ival=0; ival<3; ival++)
                if( ival == 0 || !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
                {
                  int ih = part + 3*ert_trig[evtype] + 3*2*evtype + 3*2*3*checkmap + 3*2*3*2*isolated[ival] + 3*2*3*2*2*ival;
                  h_ert[ih]->Fill(photon_pT);
                }
        } // checkmap
      } // InFiducial

      for(unsigned j=0; j<ncluster; j++)
        if( j != i && find(v_used.begin(), v_used.end(), j) == v_used.end() )
        {
          emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
          if( emcwarnmap->IsGoodTower(cluster2) &&
              TestPhoton(cluster2, bbc_t0) &&
              anatools::GetAsymmetry_E(cluster1, cluster2) < AsymCut )
          {
            int part2 = anatools::GetPart(cluster2);
            if( part2 != part ) continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);
            double tot_E = cluster1->ecore() + cluster2->ecore();

            if( cluster2->ecore() > cluster1->ecore() )
              for(int evtype=0; evtype<3; evtype++)
                ert_trig[evtype] = anatools::PassERT(data_ert, cluster2, (anatools::TriggerMode)evtype) ? 1 : 0;

            for(int checkmap=0; checkmap<2; checkmap++)
            {
              dcdeadmap->Checkmap(checkmap);

              int isolated[3] = {};
              double econe[3];
              SumEPi0(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econe);
              for(int ival=0; ival<3; ival++)
                if( econe[ival] < eratio * tot_E )
                  isolated[ival] = 1;

              for(int evtype=0; evtype<3; evtype++)
                if( part >= 0 && (datatype == MB || FireERT[evtype][oarm]) )
                  for(int ival=0; ival<3; ival++)
                  {
                    int ih = part + 3*ert_trig[evtype] + 3*2*evtype + 3*2*3*checkmap + 3*2*3*2*isolated[ival] + 3*2*3*2*2*ival;
                    h2_ert_pion[ih]->Fill(tot_pT, minv);
                  }
            } // checkmap
          } // check photon2
        } // j loop
    } // check photon1
  } // i loop

  return EVENT_OK;
}

int PhotonHistos::FillTrackQuality(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert)
{
  /* Require ERT_4x4c and BBC 30cm vertex cut */
  double bbc_z = data_global->getBbcZVertex();
  if( !IsEventType(2, data_triggerlvl1) ||
      fabs(bbc_z) > 30. )
    return DISCARDEVENT;

  int npart = data_tracks->get_npart();
  for(int itrk=0; itrk<npart; itrk++)
  {
    unsigned quality = data_tracks->get_quality(itrk);
    int iqual = 0;
    if( quality == 63 || quality == 31 || quality == 51 )
      iqual = 2;
    else if(quality > 3)
      iqual = 1;
    double mom = data_tracks->get_mom(itrk);
    if( quality > 63 || !TMath::Finite(mom) ||
        mom < pTrkMin || mom > pTrkMax )
      continue;

    /* phi and zed distributions */
    double dcphi = data_tracks->get_phi(itrk);
    double dczed = data_tracks->get_zed(itrk);
    double dcalpha = data_tracks->get_alpha(itrk);
    int dcarm = data_tracks->get_dcarm(itrk);
    int dcns = dczed > 0. ? 0 : 1;
    int dcwe = dcarm == 1 ? 0 : 1;
    double dcboard = 0.;
    if( dcwe == 0 )
      dcboard = ( 0.573231 + dcphi - 0.0046 * cos( dcphi + 0.05721 ) ) / 0.01963496;
    else
      dcboard = ( 3.72402 - dcphi + 0.008047 * cos( dcphi + 0.87851 ) ) / 0.01963496;

    /* emcdphi and emcdz distributions */
    int ih = dcns + 2*dcwe + 2*2*iqual;
    double dcdphi = data_tracks->get_emcdphi(itrk);
    double dcdz = data_tracks->get_emcdz(itrk);
    h3_dcdphiz[ih]->Fill(dcdphi, dcdz, mom);

    /* alpha-board */
    ih = dcns + 2*dcwe + 2*2*iqual;
    h2_alphaboard[ih]->Fill(dcboard, dcalpha);

    /* DC+PC1 live area */
    if(quality > 3)
    {
      dcdeadmap->Checkmap(true);
      int ih = dcdeadmap->IsDead(data_tracks, itrk) ? 0 : 1;
      h3_dclive[ih]->Fill(dczed, dcphi, mom);
    }

    int charge = data_tracks->get_charge(itrk);
    double px = data_tracks->get_px(itrk);
    double py = data_tracks->get_py(itrk);
    double pt = sqrt( px*px + py*py );
    double prod = -charge * pt * dcalpha;
    h_prod->Fill(prod);
  }

  return EVENT_OK;
}

int PhotonHistos::FillPi0Spectrum(const int ical, const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;
  double bbc_t0 = data_global->getBbcTimeZero();

  /* Get crossing number */
  const int crossing = data_triggerlvl1->get_lvl1_clock_cross();
  const int crossing_shift = spinpattern->get_crossing_shift();
  const int bunch = ( crossing + crossing_shift ) % 120;
  const int evenodd = bunch % 2;
  const int pattern[3] = {
    spinpattern->get_spinpattern_blue(bunch),
    spinpattern->get_spinpattern_yellow(bunch),
    GetPattern(crossing)};

  /* Count event multiplicity */
  // mul[beam][ipT]
  int mul_sig[3][npT_pol+1] = {};
  int mul_bg[3][npT_pol+1] = {};

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
        if( cluster1->ecore() > eMin &&
            cluster2->ecore() > eMin )
        {
          bool trig[4] = {};
          for(int evtype=0; evtype<4; evtype++)
          {
            anatools::TriggerMode trigmode = anatools::ERT_4x4or;
            if(evtype < 3)
              trigmode = (anatools::TriggerMode)evtype;
            if( cluster1->ecore() > cluster2->ecore() )
              trig[evtype] = anatools::PassERT(data_ert, cluster1, trigmode);
            else
              trig[evtype] = anatools::PassERT(data_ert, cluster2, trigmode);
          }
          if(datatype == ERT && !trig[3]) continue;

          int part1 = anatools::GetPart(cluster1);
          int part2 = anatools::GetPart(cluster2);
          double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
          double minv = anatools::GetInvMass(cluster1, cluster2);
          double tot_E = cluster1->ecore() + cluster2->ecore();

          int tof = 0;
          int prob = 0;
          if( fabs( cluster1->tofcorr() - bbc_t0 ) < tofMax &&
              fabs( cluster2->tofcorr() - bbc_t0 ) < tofMax )
            tof = 1;
          if( cluster1->prob_photon() > probMin &&
              cluster2->prob_photon() > probMin )
            prob = 1;

          if(ical == 0)
            for(int checkmap=0; checkmap<2; checkmap++)
            {
              dcdeadmap->Checkmap(checkmap);

              int isolated[3] = {};
              double econe[3];
              SumEPi0(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econe);
              for(int ival=0; ival<3; ival++)
                if( econe[ival] < eratio * tot_E )
                  isolated[ival] = 1;

              for(int evtype=0; evtype<3; evtype++)
                if( part1 >= 0 && part1 == part2 &&
                    IsEventType(evtype, data_triggerlvl1) && trig[evtype] &&
                    BBC10cm(data_global, data_triggerlvl1) &&
                    emcwarnmap->IsGoodTower(cluster1) &&
                    emcwarnmap->IsGoodTower(cluster2) &&
                    anatools::GetAsymmetry_E(cluster1, cluster2) < AsymCut )
                  for(int ival=0; ival<3; ival++)
                  {
                    int ih = part1 + 3*evtype + 3*3*tof + 3*3*2*prob + 3*3*2*2*checkmap + 3*3*2*2*2*isolated[ival] + 3*3*2*2*2*2*ival;
                    h2_pion[ih]->Fill(tot_pT, minv);
                  } // evtype, ival
            } // checkmap

          if( ical == 1 && evenodd >= 0 && prob == 1 &&
              IsEventType(3, data_triggerlvl1) && trig[3] &&
              emcwarnmap->PassCut(cluster1) &&
              emcwarnmap->PassCut(cluster2) &&
              fabs( cluster1->tofcorr() - bbc_t0 ) < tofMaxALL &&
              fabs( cluster2->tofcorr() - bbc_t0 ) < tofMaxALL &&
              PassChargeVeto(cluster1) &&
              PassChargeVeto(cluster2) )
          {
            int ipt = 0;
            for(ipt=0; ipt<npT_pol; ipt++)
              if( tot_pT < pTbin_pol[ipt+1] )
                break;

            for(int beam=0; beam<3; beam++)
              if( abs(pattern[beam]) == 1 )
              {
                int pol = pattern[beam] > 0 ? 1 : 0;
                int ih = beam + 3*evenodd + 3*2*pol;
                h2_pion_pol[ih]->Fill(tot_pT, minv);

                if( tot_pT > pTbin_pol[0] )
                {
                  if( minv > 0.112 && minv < 0.162 )
                    mul_sig[beam][ipt]++;
                  else if( (minv > 0.047 && minv < 0.097) ||
                      (minv > 0.177 && minv < 0.227) )
                    mul_bg[beam][ipt]++;
                }
              } // fill pol
          } // fill pol, calc mul
        } // min energy cut
      } // j loop
  } // i loop

  if( ical == 1 && evenodd >= 0 )
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double xpt = ( pTbin_pol[ipt] + pTbin_pol[ipt+1] ) / 2.;
      for(int beam=0; beam<3; beam++)
      {
        int ih = beam + 3*evenodd;
        h2_mul_pion_sig[ih]->Fill(xpt, (double)mul_sig[beam][ipt]);
        h2_mul_pion_bg[ih]->Fill(xpt, (double)mul_bg[beam][ipt]);
      }
    }

  return EVENT_OK;
}

int PhotonHistos::FillPhotonSpectrum(const int ical, const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
    const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert)
{
  /* Get event global parameters */
  double bbc_z = data_global->getBbcZVertex();
  if( fabs(bbc_z) > 30. ) return DISCARDEVENT;
  double bbc_t0 = data_global->getBbcTimeZero();

  /* Get crossing number */
  const int crossing = data_triggerlvl1->get_lvl1_clock_cross();
  const int crossing_shift = spinpattern->get_crossing_shift();
  const int bunch = ( crossing + crossing_shift ) % 120;
  const int evenodd = bunch % 2;
  const int pattern[3] = {
    spinpattern->get_spinpattern_blue(bunch),
    spinpattern->get_spinpattern_yellow(bunch),
    GetPattern(crossing)};

  /* Count event multiplicity */
  // mul[imul][beam][checkmap][ipt]
  int mul_photon[6][3][2][npT_pol+1] = {};

  unsigned ncluster = data_emccontainer->size();

  for(unsigned i=0; i<ncluster; i++)
  {
    emcClusterContent *cluster1 = data_emccontainer->getCluster(i);
    if( emcwarnmap->InFiducial(cluster1) &&
        TestPhoton(cluster1, bbc_t0) )
    {
      bool trig[4] = {};
      for(int evtype=0; evtype<4; evtype++)
      {
        anatools::TriggerMode trigmode = anatools::ERT_4x4or;
        if(evtype < 3)
          trigmode = (anatools::TriggerMode)evtype;
        trig[evtype] = anatools::PassERT(data_ert, cluster1, trigmode);
      }
      if(datatype == ERT && !trig[3]) continue;

      int part = anatools::GetPart(cluster1);
      double pT = anatools::Get_pT(cluster1);

      for(int checkmap=0; checkmap<2; checkmap++)
      {
        dcdeadmap->Checkmap(checkmap);

        int isolated[3] = {};
        double econeEM[2], econeTrk[3];
        SumEEmcal(cluster1, data_emccontainer, data_tracks, bbc_t0, econeEM);
        SumPTrack(cluster1, data_tracks, econeTrk);
        for(int ival=0; ival<3; ival++)
        {
          double econe = (ival == 0 ? econeEM[0] : econeEM[1]) + econeTrk[ival];
          if( econe < eratio * cluster1->ecore() )
            isolated[ival] = 1;
        }

        if( ical == 0 && part >= 0 && pT > 5. && pT < 10. &&
            IsEventType(2, data_triggerlvl1) && trig[2] &&
            BBC10cm(data_global, data_triggerlvl1) )
          for(int ival=0; ival<3; ival++)
            if( ival == 0 || !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
            {
              TLorentzVector pE = anatools::Get_pE(cluster1);
              double eta = pE.Eta();
              double phi = pE.Phi();
              if(part >= 1)
              {
                pE.RotateZ(-PI);
                phi = pE.Phi() + PI;
                pE.RotateZ(PI);
              }

              int ih = part + 3*checkmap + 3*2*isolated[ival] + 3*2*2*ival;
              h2_eta_phi[ih]->Fill(eta, phi);
            }

        for(int evtype=0; evtype<3; evtype++)
          if( ical == 0 && part >= 0 &&
              IsEventType(evtype, data_triggerlvl1) && trig[evtype] &&
              BBC10cm(data_global, data_triggerlvl1) )
            for(int ival=0; ival<3; ival++)
              if( ival == 0 || !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
              {
                int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated[ival] + 3*3*2*2*ival;
                h_1photon[ih]->Fill(pT);
              }

        if( evenodd >= 0 &&
            IsEventType(3, data_triggerlvl1) && trig[3] &&
            !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
        {
          int ipt = 0;
          for(ipt=0; ipt<npT_pol; ipt++)
            if( pT < pTbin_pol[ipt+1] )
              break;

          for(int beam=0; beam<3; beam++)
            if( abs(pattern[beam]) == 1 )
            {
              int pol = pattern[beam] > 0 ? 1 : 0;
              int ih = beam + 3*evenodd + 3*2*pol + 3*2*2*checkmap + 3*2*2*2*isolated[1] + 3*2*2*2*2*ical;
              h_1photon_pol[ih]->Fill(pT);

              if(beam == 2)
              {
                int ih = isolated[1] + 6*evenodd + 6*2*(bunch/2) + 6*2*60*checkmap + 6*2*60*2*ical;
                h_photon_bunch[ih]->Fill(pT);
              }

              if( pT > pTbin_pol[0] )
                mul_photon[isolated[1]][beam][checkmap][ipt]++;
            } // beam
        } // fill pol, calc mul

        for(unsigned j=0; j<ncluster; j++)
          if(j != i)
          {
            emcClusterContent *cluster2 = data_emccontainer->getCluster(j);
            if( !emcwarnmap->IsGoodTower(cluster2) ||
                !TestPhoton(cluster2, bbc_t0) )
              continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);

            int isopair[3] = {};
            double econeEMPair[2];
            SumEEmcal(cluster1, cluster2, data_emccontainer, data_tracks, bbc_t0, econeEMPair);
            for(int ival=0; ival<3; ival++)
            {
              double econePair = (ival == 0 ? econeEMPair[0] : econeEMPair[1]) + econeTrk[ival];
              if( econePair < eratio * cluster1->ecore() )
                isopair[ival] = 1;
            }

            for(int evtype=0; evtype<3; evtype++)
              if( ical == 0 && part >= 0 &&
                  IsEventType(evtype, data_triggerlvl1) && trig[evtype] &&
                  BBC10cm(data_global, data_triggerlvl1) )
                for(int ival=0; ival<3; ival++)
                  if( ival == 0 || !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
                  {
                    int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated[ival] + 3*3*2*2*isopair[ival] + 3*3*2*2*2*ival;
                    h2_2photon[ih]->Fill(pT, minv);
                    h2_2photon2pt[ih]->Fill(tot_pT, minv);
                  }

            if( evenodd >= 0 &&
                IsEventType(3, data_triggerlvl1) && trig[3] &&
                !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
            {
              int ipt[2] = {};
              double pTcmp[2] = {pT, tot_pT};
              for(int pttype=0; pttype<2; pttype++)
                for(ipt[pttype]=0; ipt[pttype]<npT_pol; ipt[pttype]++)
                  if( pTcmp[pttype] < pTbin_pol[ipt[pttype]+1] )
                    break;

              for(int beam=0; beam<3; beam++)
                if( abs(pattern[beam]) == 1 )
                {
                  int pol = pattern[beam] > 0 ? 1 : 0;
                  int ih = beam + 3*evenodd + 3*2*pol + 3*2*2*checkmap + 3*2*2*2*isolated[1] + 3*2*2*2*2*isopair[1] + 3*2*2*2*2*2*ical;
                  h2_2photon_pol[ih]->Fill(pT, minv);
                  h2_2photon2pt_pol[ih]->Fill(tot_pT, minv);

                  for(int pttype=0; pttype<2; pttype++)
                  {
                    int ih = 6*evenodd + 6*2*(bunch/2) + 6*2*60*checkmap + 6*2*60*2*ical;
                    double fill_pT = pttype ? tot_pT : pT;
                    if( minv > 0.112 && minv < 0.162 )
                    {
                      if(beam == 2)
                        h_photon_bunch[2+2*pttype+ih]->Fill(fill_pT);
                      if( pTcmp[pttype] > pTbin_pol[0] )
                        mul_photon[2+2*pttype][beam][checkmap][ipt[pttype]]++;
                    }
                    else if( (minv > 0.047 && minv < 0.097) ||
                        (minv > 0.177 && minv < 0.227) )
                    {
                      if(beam == 2)
                        h_photon_bunch[3+2*pttype+ih]->Fill(fill_pT);
                      if( pTcmp[pttype] > pTbin_pol[0] )
                        mul_photon[3+2*pttype][beam][checkmap][ipt[pttype]]++;
                    }
                  } // pttype
                } // fill pol
            } // fill pol, calc mul
          } // j loop
      } // checkmap
    } // check photon1
  } // i loop

  if( evenodd >= 0 )
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double xpt = ( pTbin_pol[ipt] + pTbin_pol[ipt+1] ) / 2.;
      for(int imul=0; imul<6; imul++)
        for(int beam=0; beam<3; beam++)
          for(int checkmap=0; checkmap<2; checkmap++)
          {
            int ih = imul + 6*beam + 6*3*evenodd + 6*3*2*checkmap + 6*3*2*2*ical;
            h2_mul_photon[ih]->Fill(xpt, (double)mul_photon[imul][beam][checkmap][ipt]);
          }
    }

  return EVENT_OK;
}

int PhotonHistos::End(PHCompositeNode *topNode)
{
  hm->dumpHistos(outFile);
  delete hm;
  delete emcrecalib;
  delete emcrecalib_sasha;
  delete emcwarnmap;
  delete dcdeadmap;
  delete spinpattern;

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
    h_events = new TH1F("h_events", "Events counter", 6,-0.5,5.5);
    h_events->GetXaxis()->SetBinLabel(1, "bbc_narrow");
    h_events->GetXaxis()->SetBinLabel(2, "bbc_narrow_10cm");
    h_events->GetXaxis()->SetBinLabel(3, "bbc_narrow_10cm_ert_c");
    h_events->GetXaxis()->SetBinLabel(4, "bbc_novtx");
    h_events->GetXaxis()->SetBinLabel(5, "bbc_novtx_30cm");
    h_events->GetXaxis()->SetBinLabel(6, "bbc_novtx_30cm_ert_c");
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

  // ih = part + 3*bbc_trig < 3*2
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

  // ih = sector + 8*ert_trig[evtype] + 8*2*evtype < 8*2*3
  for(int ih=0; ih<nh_ertsm; ih++)
  {
    /* ERT trigger efficiency for each of super modules */
    h2_ertsm[ih] = new TH2F(Form("h2_ertsm_%d",ih), "ERT efficiency;p_{T} [GeV];SM;", npT,pTbin, 32,-0.5,31.5);
    hm->registerHisto(h2_ertsm[ih]);
  }

  // ih = part + 3*ert_trig[evtype] + 3*2*evtype + 3*2*3*checkmap + 3*2*3*2*isolated[ival] + 3*2*3*2*2*ival < 3*2*3*2*2*3
  for(int ih=0; ih<nh_ert; ih++)
  {
    /* ERT trigger efficiency for photon */
    h_ert[ih] = new TH1F(Form("h_ert_%d",ih), "ERT efficiency;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_ert[ih]);
    /* ERT Trigger efficiency for pion */
    h2_ert_pion[ih] = new TH2F(Form("h2_ert_pion_%d",ih), "ERT efficiency;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_ert_pion[ih]);
  }

  /* Track quality study */
  // ih = dcns + 2*dcwe + 2*2*iqual < 2*2*3
  for(int ih=0; ih<nh_dcpartqual; ih++)
  {
    h3_dcdphiz[ih] = new TH3F(Form("h3_dcdphiz_%d",ih), "EMCal and DC phi and z deviation;dphi [rad];dz [cm];mom [GeV];",
        100,-0.1,0.1, 100,-50.,50., 30,0.,15.);
    hm->registerHisto(h3_dcdphiz[ih]);
    h2_alphaboard[ih] = new TH2F(Form("h2_alphaboard_%d",ih), "DC #alpha-board;board;#alpha;", 82*4,-1.5,80.5, 120,-0.6,0.6);
    hm->registerHisto(h2_alphaboard[ih]);
  }

  // ih = isGoodDC < 2
  for(int ih=0; ih<nh_dcgood; ih++)
  {
    h3_dclive[ih] = new TH3F(Form("h3_dclive_%d",ih), "DC zed and phi distribution;zed [cm];phi [rad];mom [GeV];",
        200,-100.,100., 50,-1.,4., 30,0.,15.);
    hm->registerHisto(h3_dclive[ih]);
  }

  h_prod = new TH1F("h_prod", "Minus charge times pT times alpha;-charge*pT*#alpha [GeV]", 200,0.,0.2);
  hm->registerHisto(h_prod);

  /* Store pi0 information */
  // ih = part1 + 3*evtype + 3*3*tof + 3*3*2*prob + 3*3*2*2*checkmap + 3*3*2*2*2*isolated[ival] + 3*3*2*2*2*2*ival < 3*3*2*2*2*2*3
  for(int ih=0; ih<nh_pion; ih++)
  {
    h2_pion[ih] = new TH2F(Form("h2_pion_%d",ih), "#pi^{0} spectrum;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    hm->registerHisto(h2_pion[ih]);
  }

  /* Store polarized pi0 information */
  // ih = beam + 3*evenodd + 3*2*pol < 3*2*2
  for(int ih=0; ih<nh_pion_pol; ih++)
  {
    h2_pion_pol[ih] = new TH2F(Form("h2_pion_pol_%d",ih), "Polarized spectrum;p_{T} [GeV];m_{inv} [GeV];",
        300,0.,30., 300,0.,0.3);
    hm->registerHisto(h2_pion_pol[ih]);
  }

  /* Eta and phi distribution */
  // ih = part + 3*checkmap + 3*2*isolated[ival] + 3*2*2*ival < 3*2*2*3
  for(int ih=0; ih<nh_eta_phi; ih++)
  {
    h2_eta_phi[ih] = new TH2F(Form("h2_eta_phi_%d",ih), "#eta and #phi distribution;#eta;#phi;", neta,etabin[ih%3/2], nphi,phibin);
    hm->registerHisto(h2_eta_phi[ih]);
  }

  /* Store single photon information */
  // ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated[ival] + 3*3*2*2*ival < 3*3*2*2*3
  for(int ih=0; ih<nh_1photon; ih++)
  {
    h_1photon[ih] = new TH1F(Form("h_1photon_%d",ih), "Single photon spectrum;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_1photon[ih]);
  }

  /* Store polarized single photons information */
  // ih = beam + 3*evenodd + 3*2*pol + 3*2*2*checkmap + 3*2*2*2*isolated[1] + 3*2*2*2*2*ical < 3*2*2*2*2*2
  for(int ih=0; ih<nh_1photon_pol; ih++)
  {
    h_1photon_pol[ih] = new TH1F(Form("h_1photon_pol_%d",ih), "Polarized single photon spectrum;p_{T} [GeV];", 300,0.,30.);
    hm->registerHisto(h_1photon_pol[ih]);
  }

  /* Store two photons information */
  // ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated[ival] + 3*3*2*2*isopair[ival] + 3*3*2*2*2*ival < 3*3*2*2*2*3
  for(int ih=0; ih<nh_2photon; ih++)
  {
    h2_2photon[ih] = new TH2F(Form("h2_2photon_%d",ih), "Two photons spectrum;p_{T} [GeV];m_{inv} [GeV];", npT,pTbin, 300,0.,0.3);
    h2_2photon2pt[ih] = (TH2*)h2_2photon[0]->Clone(Form("h2_2photon2pt_%d",ih));
    hm->registerHisto(h2_2photon[ih]);
    hm->registerHisto(h2_2photon2pt[ih]);
  }

  /* Store polarized two photons information */
  // ih = beam + 3*evenodd + 3*2*pol + 3*2*2*checkmap + 3*2*2*2*isolated[1] + 3*2*2*2*2*isopair[1] + 3*2*2*2*2*2*ical < 3*2*2*2*2*2*2
  for(int ih=0; ih<nh_2photon_pol; ih++)
  {
    h2_2photon_pol[ih] = (TH2*)h2_pion_pol[0]->Clone(Form("h2_2photon_pol_%d",ih));
    h2_2photon2pt_pol[ih] = (TH2*)h2_pion_pol[0]->Clone(Form("h2_2photon2pt_pol_%d",ih));
    hm->registerHisto(h2_2photon_pol[ih]);
    hm->registerHisto(h2_2photon2pt_pol[ih]);
  }

  /* Store pion event multiplicity information */
  // ih = beam + 3*evenodd < 3*2
  for(int ih=0; ih<nh_mul_pion; ih++)
  {
    h2_mul_pion_sig[ih] = new TH2F(Form("h2_mul_pion_sig_%d",ih), "Event multiplity;p_{T} [GeV];Multiplicity;", npT_pol,pTbin_pol, 20,-0.5,19.5);
    h2_mul_pion_bg[ih] = (TH2*)h2_mul_pion_sig[0]->Clone(Form("h2_mul_pion_bg_%d",ih));
    hm->registerHisto(h2_mul_pion_sig[ih]);
    hm->registerHisto(h2_mul_pion_bg[ih]);
  }

  /* Store isolated photon event multiplicity information */
  // ih = imul + 6*beam + 6*3*evenodd + 6*3*2*checkmap + 6*3*2*2*ical < 6*3*2*2*2
  for(int ih=0; ih<nh_mul_photon; ih++)
  {
    h2_mul_photon[ih] = (TH2*)h2_mul_pion_sig[0]->Clone(Form("h2_mul_photon_%d",ih));
    hm->registerHisto(h2_mul_photon[ih]);
  }

  /* Store polarized photon information for bunch shuffling */
  // ih = imul + 6*evenodd + 6*2*(bunch/2) + 6*2*60*checkmap + 6*2*60*2*ical < 6*2*60*2*2
  for(int ih=0; ih<nh_photon_bunch; ih++)
  {
    h_photon_bunch[ih] = new TH1F(Form("h_photon_bunch_%d",ih), "Polarized photon spectrum;p_{T} [GeV];", npT_pol,pTbin_pol);
    hm->registerHisto(h_photon_bunch[ih]);
  }

  return;
}

void PhotonHistos::SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
    const PHCentralTrack *data_tracks, double bbc_t0, double econe[])
{ 
  /* Sum up all energy in cone around particle without one cluster */
  for(int ival=0; ival<2; ival++)
    econe[ival] = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus2 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to 'reference' particle
     * or on bad towers or lower than energy threshold */
    if( clus2->id() == cluster->id() ||
        emcwarnmap->IsBadTower(clus2) ||
        fabs( clus2->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus2->ecore() < eClusMin )
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
    {
      /* 3 sigma charge veto */
      econe[0] += clus2->ecore();
      if( !dcdeadmap->ChargeVeto(clus2, data_tracks) )
        econe[1] += clus2->ecore();
    }
  }

  return;
}

void PhotonHistos::SumEEmcal(const emcClusterContent *cluster, const emcClusterContent *cluster_part,
    const emcClusterContainer *cluscont, const PHCentralTrack *data_tracks, double bbc_t0, double econe[])
{ 
  /* Sum up all energy in cone around particle without the partner cluster */
  for(int ival=0; ival<2; ival++)
    econe[ival] = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus3 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to any of the two 'reference' particles
     * or on bad towers or lower than energy threshold */
    if( clus3->id() == cluster->id() ||
        clus3->id() == cluster_part->id() ||
        emcwarnmap->IsBadTower(clus3) ||
        fabs( clus3->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus3->ecore() < eClusMin )
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
    {
      econe[0] += clus3->ecore();
      /* 3 sigma charge veto */
      if( !dcdeadmap->ChargeVeto(clus3, data_tracks) )
        econe[1] += clus3->ecore();
    }
  }

  return;
}

void PhotonHistos::SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks, double econe[])
{ 
  /* Sum up all energy in cone around particle */
  for(int ival=0; ival<3; ival++)
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

    /* Test if track passes the momentum cuts and deadmap */
    if( !TMath::Finite(px+py+pz+mom) ||
        mom < pTrkMin || mom > pTrkMax ||
        dcdeadmap->IsDead(data_tracks, i) )
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
  double econeEM[2] = {};
  for(int ival=0; ival<3; ival++)
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
        emcwarnmap->IsBadTower(clus3) ||
        fabs( clus3->tofcorr() - bbc_t0 ) > tofMaxIso ||
        clus3->ecore() < eClusMin )
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
    {
      econeEM[0] += clus3->ecore();
      /* 3 sigma charge veto */
      if( !dcdeadmap->ChargeVeto(clus3, data_tracks) )
        econeEM[1] += clus3->ecore();
    }
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

    /* Test if track passes the momentum cuts and deadmap */
    if( mom < pTrkMin || mom > pTrkMax ||
        dcdeadmap->IsDead(data_tracks, i) )
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
      if( quality > 3 )
        econe[1] += mom;
      if( quality == 63 || quality == 31 || quality == 51 )
        econe[2] += mom;
    }
  }

  for(int ival=0; ival<3; ival++)
    econe[ival] += (ival == 0 ? econeEM[0] : econeEM[1]);

  return;
}

bool PhotonHistos::IsEventType(int evtype, const TrigLvl1 *data_triggerlvl1)
{
  const unsigned lvl1_scaled = data_triggerlvl1->get_lvl1_trigscaled();

  if( datatype == ERT && evtype < 4 &&
      (lvl1_scaled & bit_ert4x4[evtype]) )
    return true;

  if( datatype == MB && evtype == 0 &&
      (lvl1_scaled & bit_bbcnarrow) )
    return true;

  return false;
}

bool PhotonHistos::BBC10cm(const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1)
{
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  double bbc_z = data_global->getBbcZVertex();

  if( datatype == ERT &&
      (lvl1_live & bit_bbcnarrow) &&
      fabs(bbc_z) < 10. )
    return true;

  if( datatype == MB &&
      fabs(bbc_z) < 10. )
    return true;

  return false;
}

bool PhotonHistos::BBC30cm(const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1)
{
  const unsigned lvl1_live = data_triggerlvl1->get_lvl1_triglive();
  double bbc_z = data_global->getBbcZVertex();

  if( datatype == ERT &&
      (lvl1_live & bit_bbcnovtx) &&
      fabs(bbc_z) < 30. )
    return true;

  if( datatype == MB &&
      fabs(bbc_z) < 30. )
    return true;

  return false;
}

bool PhotonHistos::TestPhoton(const emcClusterContent *cluster, double bbc_t0)
{
  if( cluster->ecore() > eMin &&
      fabs( cluster->tofcorr() - bbc_t0 ) < tofMax &&
      cluster->prob_photon() > probMin )
    return true;
  else
    return false;
}

bool PhotonHistos::PassChargeVeto(const emcClusterContent *cluster)
{
  /* Angle between EMC cluster and PC3 track */
  double theta_cv = anatools::GetTheta_CV(cluster);

  /* Charge veto in medium dtheta */
  double emc_e = cluster->ecore();
  int sector = anatools::CorrectClusterSector(cluster->arm(), cluster->sector());
  if(sector < 6)
  {
    if( theta_cv > 4.22e-4 - 1.16e-2*emc_e - 4.53e-3*pow(emc_e,2) &&
        theta_cv < 1.01e-1 - 2.02e-1*emc_e + 1.51e-1*pow(emc_e,2) - 3.66e-2*pow(emc_e,3) )
      return false;
  }
  else
  {
    if( theta_cv > 1.27e-2 - 2.41e-2*emc_e + 2.26e-2*pow(emc_e,2) &&
        theta_cv < 1.64e-2 - 7.38e-3*emc_e + 1.45e-1*exp(-4.*emc_e) )
      return false;
  }

  return true;
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

  if( abs(pattern) != 1 )
    return 0;

  return pattern;
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

      /* For Run13, GL1p swap narrow and wide.
       * So I swap narrow and wide to correct the error. */
      spinpattern->set_bbc_narrow(i, bbc_wide[i]);
      spinpattern->set_bbc_wide(i, bbc_narrow[i]);
      spinpattern->set_zdc_narrow(i, zdc_wide[i]);
      spinpattern->set_zdc_wide(i, zdc_narrow[i]);
    }
  }

  return;
}
