#include "HadronResponse.h"

#include "AnaTrk.h"
#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"

#include <DCDeadmapChecker.h>

#include <PHGlobal.h>
#include <PHCentralTrack.h>
#include <McEvalSingleList.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <boost/foreach.hpp>

using namespace std;

double HadronResponse::vpT[] = { 0.0,
  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
  5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

/* global constants */
const int PHOTON_PID = 1;
const int POSITRON_PID = 2;
const int ELECTRON_PID = 3;
const int PIZERO_PID = 7;

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

HadronResponse::HadronResponse(const string &name,
    const char *filename,
    const double pt_init):
  SubsysReco(name),
  h_pt_weight(NULL),
  pt_start(pt_init),
  pt_count(0),
  weight_pythia(0.),
  dcdeadmap(NULL),
  hm(NULL),
  hn_dc(NULL),
  hn_emcal(NULL),
  hn_1photon(NULL),
  hn_2photon(NULL),
  hn_cluster(NULL)
{
  /* construct output file names */
  outFileName = "histos/HadronResponse-";
  outFileName.append(filename);

  /* initialize array for tower status */
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
      {
        tower_status_nils[isector][ibiny][ibinz] = 0;
        tower_status_sasha[isector][ibiny][ibinz] = 0;
      }

  for(int ih=0; ih<nh_dcpart; ih++)
    h2_alphaboard[ih] = NULL;
}

HadronResponse::~HadronResponse()
{
}

int HadronResponse::Init(PHCompositeNode *topNode)
{
  /* Get weight histogram */
  TFile *f_proc_pt = new TFile("../AnaPHPythiaDirectPhoton-macros/PureHard-histo.root");
  if( f_proc_pt && f_proc_pt->IsOpen() )
  {
    TH2 *h2_proc_pt = (TH2*)f_proc_pt->Get("h2_proc_pt");
    if(h2_proc_pt)
      h_pt_weight = h2_proc_pt->ProjectionX("h_pt_weight");
  }

  if(!h_pt_weight)
  {
    cout << "Cannot get weight histogram" << endl;
    exit(1);
  }

  /* initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  /* create and register histograms */
  BookHistograms();

  /* read warnmap */
  ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  /* initialize DC deadmap checker */
  dcdeadmap = new DCDeadmapChecker();

  return EVENT_OK;
}

void HadronResponse::InitBatch()
{
  /* Set Pythia weight */
  double pt_low = pt_start + pt_count/2 * 0.1;
  int ipt_cut = h_pt_weight->GetXaxis()->FindBin(3.);
  int ipt_low = h_pt_weight->GetXaxis()->FindBin(pt_low);
  int ipt_high = h_pt_weight->GetXaxis()->FindBin(50.);
  weight_pythia = h_pt_weight->Integral(ipt_low,ipt_high) /
    h_pt_weight->Integral(ipt_cut,ipt_high);
  pt_count++;

  return;
}

int HadronResponse::process_event(PHCompositeNode *topNode)
{
  /* central track reco info */
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!data_tracks)
  {
    cout << "Cannot find PHCentralTrack" << endl;
    return DISCARDEVENT;
  }

  /* central track sim info */
  McEvalSingleList *mctrk = findNode::getClass<McEvalSingleList>(topNode, "McSingle");
  if(!mctrk)
  {
    cout << "Cannot find McEvalSingleList" << endl;
    return DISCARDEVENT;
  }

  /* emc track info */
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  /* emc cluster info */
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

  /* Set DC deadmap */
  dcdeadmap->SetMapByRandom();

  int nPart = data_tracks->get_npart();
  int nMC = mctrk->get_McEvalSingleTrackN();
  //if(nPart <= 0 || nMC <= 0)
  //  return DISCARDEVENT;
  //int *cglMC = new int[nPart];

  // associate reco tracks with mc info
  //for(int iPart=0; iPart<nPart; iPart++)
  //{
  //  cglMC[iPart] = -1;
  //  double dpMin = 9999.;
  //  for(int iMC=0; iMC<nMC; iMC++)
  //  {
  //    double pxMC = mctrk->get_momentumx(iMC);
  //    double pyMC = mctrk->get_momentumy(iMC);
  //    double pzMC = mctrk->get_momentumz(iMC);
  //    int nReco = mctrk->get_Nreco(iMC);
  //    for(int iReco=0; iReco<nReco; iReco++)
  //    {
  //      int jPart = mctrk->get_recoid(iMC, iReco);
  //      if(jPart == iPart)
  //      {
  //        double pxReco = data_tracks->get_px(iPart);
  //        double pyReco = data_tracks->get_py(iPart);
  //        double pzReco = data_tracks->get_pz(iPart);
  //        double dp = sqrt( (pxMC-pxReco)*(pxMC-pxReco) + (pyMC-pyReco)*(pyMC-pyReco) + (pzMC-pzReco)*(pzMC-pzReco) );
  //        if(dp < dpMin)
  //        {
  //          cglMC[iPart] = iMC;
  //          dpMin = dp;
  //        }
  //      } // jPart
  //    } // iReco
  //  } // iMC
  //} // iPart

  /* loop over cgl reco tracks */
  for(int iPart=0; iPart<nPart; iPart++)
  {
    unsigned quality = data_tracks->get_quality(iPart);
    double dcphi = data_tracks->get_phi(iPart);
    double dczed = data_tracks->get_zed(iPart);
    double dcalpha = data_tracks->get_alpha(iPart);
    int dcarm = data_tracks->get_dcarm(iPart);
    int dcns = dczed > 0. ? 0 : 1;
    int dcwe = dcarm == 1 ? 0 : 1;
    double dcboard = 0.;
    if( dcwe == 0 )
      dcboard = ( 0.573231 + dcphi - 0.0046 * cos( dcphi + 0.05721 ) ) / 0.01963496;
    else
      dcboard = ( 3.72402 - dcphi + 0.008047 * cos( dcphi + 0.87851 ) ) / 0.01963496;

    double pxReco = data_tracks->get_px(iPart);
    double pyReco = data_tracks->get_py(iPart);
    double pzReco = data_tracks->get_pz(iPart);
    double ptReco = sqrt( pxReco*pxReco + pyReco*pyReco );
    double momReco = sqrt( pxReco*pxReco + pyReco*pyReco + pzReco*pzReco );

    if( quality <= 3 || !TMath::Finite(pxReco+pyReco+pzReco) )
      continue;

    /* Test if track passes the momentum cuts */
    if( momReco > pTrkMin && momReco < pTrkMax )
    {
      int ih = dcns + 2*dcwe;
      h2_alphaboard[ih]->Fill(dcboard, dcalpha, weight_pythia);
    }

    double fill_hn_dc[] = {ptReco, momReco, (double)dcns, (double)dcwe, 1.};
    hn_dc->Fill(fill_hn_dc, weight_pythia);
  }

  /* loop over cgl truth tracks */
  for(int iMC=0; iMC<nMC; iMC++)
  {
    double dczed = mctrk->get_zed(iMC);
    double dcphi = mctrk->get_phi(iMC);
    int dcns = dczed > 0. ? 0 : 1;
    int dcwe = dcphi < PI/2. ? 0 : 1;

    double pxMC = mctrk->get_momentumx(iMC);
    double pyMC = mctrk->get_momentumy(iMC);
    double pzMC = mctrk->get_momentumz(iMC);
    double ptMC = sqrt( pxMC*pxMC + pyMC*pyMC );
    double momMC = sqrt( pxMC*pxMC + pyMC*pyMC + pzMC*pzMC );

    double fill_hn_dc[] = {ptMC, momMC, (double)dcns, (double)dcwe, 0.};
    hn_dc->Fill(fill_hn_dc, weight_pythia);
  }

  /* associate cluster with track
   * map key is trkno */
  typedef map<int,AnaTrk*> map_Ana_t;
  map_Ana_t track_list;

  /* store photon */
  vector<AnaTrk*> photon;

  /* number of tracks and clusters */
  int nemctrk = emctrkcont->size();
  int nemcclus = emccluscont->size();
  if( nemctrk<=0 || nemcclus<=0 ) return DISCARDEVENT;

  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *track = new AnaTrk(emctrk, emccluscont, (int*)tower_status_sasha);
    if(track)
      track_list.insert( make_pair(track->trkno,track) );
  }

  /* analyze emc clusters */
  for(int iclus=0; iclus<nemcclus; iclus++)
  {
    emcClusterContent *cluster1 = emccluscont->getCluster(iclus);
    if( IsGoodTower(cluster1) &&
        cluster1->ecore() > eMin )
    {
      int sector = anatools::GetSector(cluster1); 
      double pT = anatools::Get_pT(cluster1);

      double fill_hn_emcal[] = {pT, cluster1->ecore(), (double)sector, 1.};
      hn_emcal->Fill(fill_hn_emcal, weight_pythia);

      if( DCChargeVeto(cluster1,data_tracks) )
      {
        double fill_hn_emcal[] = {pT, cluster1->ecore(), (double)sector, 2.};
        hn_emcal->Fill(fill_hn_emcal, weight_pythia);
      }

      if( InFiducial(cluster1) &&
          cluster1->prob_photon() > probMin &&
          !DCChargeVeto(cluster1,data_tracks) )
      {
        double econeEM = SumEEmcal(cluster1, emccluscont, data_tracks);
        double econeTrk = SumPTrack(cluster1, data_tracks);
        double econe = econeEM + econeTrk;

        int isolated = 0;
        if( econe < eratio * cluster1->ecore() )
          isolated = 1;

        double fill_hn_1photon[] = {pT, (double)sector, (double)isolated};
        hn_1photon->Fill(fill_hn_1photon, weight_pythia);

        for(int jclus=0; jclus<nemcclus; jclus++)
          if(jclus != iclus)
          {
            emcClusterContent *cluster2 = emccluscont->getCluster(jclus);
            if( !IsGoodTower(cluster2) ||
                cluster2->ecore() < eMin ||
                cluster2->prob_photon() < probMin ||
                DCChargeVeto(cluster2,data_tracks) )
              continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);

            int isoboth = 0;
            double econeEM2 = SumEEmcal(cluster2, emccluscont, data_tracks);
            double econeTrk2 = SumPTrack(cluster2, data_tracks);
            double econe2 = econeEM2 + econeTrk2;
            if( isolated && econe2 < eratio * cluster2->ecore() )
              isoboth = 1;

            int isopair = 0;
            double econeEMPair1, econeEMPair2;
            SumEEmcal(cluster1, cluster2, emccluscont, data_tracks, econeEMPair1, econeEMPair2);
            double econePair1 = econeEMPair1 + econeTrk;
            double econePair2 = econeEMPair2 + econeTrk2;
            if( econePair1 < eratio * cluster1->ecore() &&
                econePair2 < eratio * cluster2->ecore() )
              isopair = 1;

            double fill_hn_2photon[] = {pT, tot_pT, minv, (double)sector, (double)isoboth, (double)isopair};
            hn_2photon->Fill(fill_hn_2photon, weight_pythia);
          } // jclus
      } // fiducial & prob & charge veto
    } // good tower & ecore
  } // iclus

  /* analyze emc tracks */
  BOOST_FOREACH( map_Ana_t::value_type &trk_el, track_list )
  {
    AnaTrk *trk = trk_el.second;
    if( trk->cid < 0 ||
        trk->trkedep < eMin )
      continue;

    double fill_hn_emcal[] = {trk->trkpt, trk->emctrk->get_ekin(), (double)trk->sector, 0.};
    hn_emcal->Fill(fill_hn_emcal, weight_pythia);

    if( trk->pid != PHOTON_PID )
      continue;

    double econeEM = SumEEmcal(trk->emcclus, emccluscont, data_tracks);
    double econeTrk = SumPTrack(trk->emcclus, data_tracks);
    double econeReco = econeEM + econeTrk;

    /* loop over second emc tracks */
    double econeMC = 0.;
    BOOST_FOREACH( map_Ana_t::value_type &trk2_el, track_list )
    {
      AnaTrk *trk2 = trk2_el.second;
      if( trk2->trkno != trk->trkno &&
          trk->trkvp.Angle(trk2->trkvp) < cone_angle )
        econeMC += trk2->emctrk->get_ekin();
    }

    double fill_hn_cluster[] = {trk->trkpt, trk->cluspt, econeMC, econeReco, (double)trk->sector};
    hn_cluster->Fill(fill_hn_cluster, weight_pythia);
  }

  /* clear associated list */
  //delete[] cglMC;
  BOOST_FOREACH( map_Ana_t::value_type &trk_el, track_list )
    delete trk_el.second;

  return EVENT_OK;
}

int HadronResponse::End(PHCompositeNode *topNode)
{
  /* write histogram output to ROOT file */
  hm->dumpHistos();
  delete hm;
  delete dcdeadmap;

  return EVENT_OK;
}

void HadronResponse::BookHistograms()
{
  /* for dc alpha-board */
  // ih = dcns + 2*dcwe < 2*2
  for(int ih=0; ih<nh_dcpart; ih++)
  {
    h2_alphaboard[ih] = new TH2F(Form("h2_alphaboard_%d",ih), "DC #alpha-board;board;#alpha;", 82*4,-1.5,80.5, 120,-0.6,0.6);
    hm->registerHisto(h2_alphaboard[ih]);
  }

  /* for dc track */
  int nbins_hn_dc[] = {npT, 300, 2, 2, 2};
  double xmin_hn_dc[] = {0., 0., -0.5, -0.5, -0.5};
  double xmax_hn_dc[] = {0., 30.,1.5, 1.5, 1.5};
  hn_dc = new THnSparseF("hn_dc", "DC track info;p_{T} [GeV];E [GeV];NS;WE;Reco;",
      5, nbins_hn_dc, xmin_hn_dc, xmax_hn_dc);
  hn_dc->SetBinEdges(0, vpT);
  hn_dc->Sumw2();
  hm->registerHisto(hn_dc);

  /* for emcal cluster */
  int nbins_hn_emcal[] = {npT, 300, 8, 3};
  double xmin_hn_emcal[] = {0., 0., -0.5, -0.5};
  double xmax_hn_emcal[] = {0., 30., 7.5, 2.5};
  hn_emcal = new THnSparseF("hn_emcal", "EMCal info;p_{T} [GeV];E [GeV];Sector;Reco;",
      4, nbins_hn_emcal, xmin_hn_emcal, xmax_hn_emcal);
  hn_emcal->SetBinEdges(0, vpT);
  hn_emcal->Sumw2();
  hm->registerHisto(hn_emcal);

  /* for emcal reco one photon */
  int nbins_hn_1photon[] = {npT, 8, 2};
  double xmin_hn_1photon[] = {0., -0.5, -0.5};
  double xmax_hn_1photon[] = {0., 7.5, 1.5};
  hn_1photon = new THnSparseF("hn_1photon", "EMCal one photon;p_{T} [GeV];Sector;Isolated;",
      3, nbins_hn_1photon, xmin_hn_1photon, xmax_hn_1photon);
  hn_1photon->SetBinEdges(0, vpT);
  hn_1photon->Sumw2();
  hm->registerHisto(hn_1photon);

  /* for emcal reco two photons */
  int nbins_hn_2photon[] = {npT, npT, 300, 8, 2, 2};
  double xmin_hn_2photon[] = {0., 0., 0., -0.5, -0.5, -0.5};
  double xmax_hn_2photon[] = {0., 0., 0.3, 7.5, 1.5, 1.5};
  hn_2photon = new THnSparseF("hn_2photon", "EMCal one photon;p_{T} [GeV];Sector;Isolated;",
      6, nbins_hn_2photon, xmin_hn_2photon, xmax_hn_2photon);
  hn_2photon->SetBinEdges(0, vpT);
  hn_2photon->SetBinEdges(1, vpT);
  hn_2photon->Sumw2();
  hm->registerHisto(hn_2photon);

  /* for isolation cone  */
  int nbins_hn_cluster[] = {npT, npT, 400, 400, 8};
  double xmin_hn_cluster[] = {0., 0., 0., 0., -0.5};
  double xmax_hn_cluster[] = {0., 0., 40., 40., 7.5};
  hn_cluster = new THnSparseF("hn_cluster", "Isolation cone info;p^{truth}_{T} [GeV];p^{reco}_{T} [GeV];E_{truth} [GeV];E_{reco} [GeV];Sector;",
      5, nbins_hn_cluster, xmin_hn_cluster, xmax_hn_cluster);
  hn_cluster->SetBinEdges(0, vpT);
  hn_cluster->SetBinEdges(1, vpT);
  hn_cluster->Sumw2();
  hm->registerHisto(hn_cluster);

  return;
}

double HadronResponse::SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
    const PHCentralTrack *data_tracks)
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
     * or on bad towers or lower than energy threshold
     * and with 3 sigma charge veto */
    if( clus2->id() == cluster->id() ||
        IsBadTower(clus2) ||
        clus2->ecore() < eClusMin ||
        DCChargeVeto(clus2,data_tracks) )
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

void HadronResponse::SumEEmcal(const emcClusterContent *cluster1, const emcClusterContent *cluster2,
    const emcClusterContainer *cluscont, const PHCentralTrack *data_tracks, double &econe1, double &econe2)
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
     * or on bad towers or lower than energy threshold 
     * and with 3 sigma charge veto */
    if( clus3->id() == cluster1->id() ||
        clus3->id() == cluster2->id() ||
        IsBadTower(clus3) ||
        abs( clus3->tofcorr() ) > tofMaxIso ||
        clus3->ecore() < eClusMin ||
        DCChargeVeto(clus3,data_tracks) )
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

double HadronResponse::SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks)
{ 
  /* Sum up all energy in cone around particle */
  double econe = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int npart = data_tracks->get_npart();

  for(int i=0; i<npart; i++)
  {
    int quality = data_tracks->get_quality(i);
    double px = data_tracks->get_px(i);
    double py = data_tracks->get_py(i);
    double pz = data_tracks->get_pz(i);
    double mom = data_tracks->get_mom(i);

    /* Test if track passes quality and the momentum cuts 
     * and the DC deadmap */
    if( quality <= 3 || !TMath::Finite(px+py+pz+mom) ||
        mom < pTrkMin || mom > pTrkMax ||
        IsDCDead(data_tracks, i) )
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
      econe += mom;
  }

  return econe;
}

bool HadronResponse::DCChargeVeto(const emcClusterContent *cluster, const PHCentralTrack *data_tracks)
{
  /* 3 sigma charge veto */
  int itrk_match = GetEmcMatchTrack(cluster, data_tracks);
  if( itrk_match >= 0 )
    return true;
  else 
    return false;
}

bool HadronResponse::InFiducial(const emcClusterContent *cluster)
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

bool HadronResponse::IsGoodTower(const emcClusterContent *cluster)
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

bool HadronResponse::IsBadTower(const emcClusterContent *cluster)
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

bool HadronResponse::IsDCDead(const PHCentralTrack *data_tracks, int itrk)
{
  /* phi and zed distributions */
  double dcphi = data_tracks->get_phi(itrk);
  double dczed = data_tracks->get_zed(itrk);
  double dcalpha = data_tracks->get_alpha(itrk);
  int dcarm = data_tracks->get_dcarm(itrk);
  string dcns = dczed > 0. ? "N" : "S";
  string dcwe = dcarm == 1 ? "W" : "E";
  string nswe = dcns + dcwe;
  double dcboard = 0.;
  if( dcwe == "W" )
    dcboard = ( 0.573231 + dcphi - 0.0046 * cos( dcphi + 0.05721 ) ) / 0.01963496;
  else
    dcboard = ( 3.72402 - dcphi + 0.008047 * cos( dcphi + 0.87851 ) ) / 0.01963496;

  return dcdeadmap->IsDead(nswe, dcboard, dcalpha);
}

int HadronResponse::GetEmcMatchTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks)
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
    if( dphi > 0.015 ||
        !TMath::Finite(mom) ||
        IsDCDead(data_tracks,itrk) )
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

void HadronResponse::ReadTowerStatus(const string& filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  unsigned int sector = 0;
  unsigned int biny = 0;
  unsigned int binz = 0;
  unsigned int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> biny >> binz >> status )
  {
    /* count tower with bad status for PbSc and PbGl */
    if ( status > 10 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status_nils[sector][biny][binz] = status;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void HadronResponse::ReadSashaWarnmap(const string &filename)
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
    /* Attention!! I use my indexing for warn map in this program!!! */
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }

    /* get tower location */
    anatools::TowerLocation(ich, sector, biny, binz);

    /* count tower with bad status for PbSc and PbGl */
    if ( status > 0 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status_sasha[sector][biny][binz] = status;

    /* mark edge towers */
    if( anatools::Edge_cg(sector, biny, binz) )
      tower_status_sasha[sector][biny][binz] = 20;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}
