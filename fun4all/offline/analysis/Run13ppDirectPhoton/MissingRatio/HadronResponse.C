#include "HadronResponse.h"

#include <AnaToolsCluster.h>
#include "ERTSimTrigger.h"
#include <EMCWarnmapChecker.h>
#include <DCDeadmapChecker.h>
#include "AnaTrk.h"

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <PHCentralTrack.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <iostream>

using namespace std;

/* Global constants */
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

HadronResponse::HadronResponse(const string &name, const char *filename):
  SubsysReco(name),
  weight_pythia(1.),
  ertsim(nullptr),
  emcwarnmap(nullptr),
  dcdeadmap(nullptr),
  hm(nullptr),
  h_events(nullptr),
  hn_1photon(nullptr),
  hn_2photon(nullptr),
  hn_prob_photon(nullptr),
  hn_dcdphiz(nullptr),
  hn_alphaboard(nullptr),
  hn_dclive(nullptr),
  h_prod(nullptr)
{
  /* Construct output file names */
  outFileName = "histos/HadronResponse-";
  outFileName.append(filename);

  for(int ih=0; ih<nh_eta_phi; ih++)
    h2_eta_phi[ih] = nullptr;
}

HadronResponse::~HadronResponse()
{
}

int HadronResponse::Init(PHCompositeNode *topNode)
{
  /* Create and register histograms */
  BookHistograms();

  /* Initialize ERT sim trigger */
  ertsim = new ERTSimTrigger();
  if(!ertsim)
  {
    cerr << "No ertsim" << endl;
    exit(1);
  }

  /* Initialize EMC warnmap checker */
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

  return EVENT_OK;
}

int HadronResponse::process_event(PHCompositeNode *topNode)
{
  /* Count events */
  h_events->Fill(1.);

  /* EMC track truth info */
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  /* EMC cluster reco info */
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

  /* Central track reco info */
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!data_tracks)
  {
    cout << "Cannot find PHCentralTrack" << endl;
    return DISCARDEVENT;
  }

  /* ERT Lvl1 trigger */
  int ERTLvl1 = 0;

  /* Set DC deadmap */
  dcdeadmap->SetMapByEvent();

  /* Loop over EMC reco clusters */
  int nemcclus = emccluscont->size();
  for(int iclus=0; iclus<nemcclus; iclus++)
  {
    emcClusterContent *cluster1 = emccluscont->getCluster(iclus);
    int trig = ertsim->Fired(cluster1) ? 1 : 0;
    if(trig) ERTLvl1 = 1;

    for(int checkmap=0; checkmap<2; checkmap++)
    {
      if( emcwarnmap->InFiducial(cluster1) &&
          TestPhoton(cluster1) &&
          !dcdeadmap->ChargeVeto(cluster1, data_tracks) )
      {
        int sector = anatools::GetSector(cluster1);
        int part = anatools::GetPart(cluster1);
        double pT = anatools::Get_pT(cluster1);

        double econeEM = SumEEmcal(cluster1, emccluscont, data_tracks);
        double econeTrk = SumPTrack(cluster1, data_tracks);
        double econe = econeEM + econeTrk;

        int isolated = 0;
        if( econe < eratio * cluster1->ecore() )
          isolated = 1;

        if( pT > 5. && pT < 10. )
        {
          TLorentzVector pE = anatools::Get_pE(cluster1);
          double eta = pE.Eta();
          double phi = pE.Phi();
          if(part >= 1)
          {
            pE.RotateZ(-PI);
            phi = pE.Phi() + PI;
          }

          int ih = part + 3*isolated + 3*2*checkmap + 3*2*2*trig;
          h2_eta_phi[ih]->Fill(eta, phi);
        }

        double fill_hn_1photon[] = {pT, (double)sector, (double)isolated, (double)checkmap, (double)trig};
        hn_1photon->Fill(fill_hn_1photon, weight_pythia);

        for(int jclus=0; jclus<nemcclus; jclus++)
          if(jclus != iclus)
          {
            emcClusterContent *cluster2 = emccluscont->getCluster(jclus);
            if( !emcwarnmap->IsGoodTower(cluster2) ||
                !TestPhoton(cluster2) )
              continue;

            double tot_pT = anatools::GetTot_pT(cluster1, cluster2);
            double minv = anatools::GetInvMass(cluster1, cluster2);

            int isopair = 0;
            double econeEMPair;
            SumEEmcal(cluster1, cluster2, emccluscont, data_tracks, econeEMPair);
            double econePair = econeEMPair + econeTrk;
            if( econePair < eratio * cluster1->ecore() )
              isopair = 1;

            double fill_hn_2photon[] =
            {pT, tot_pT, minv, (double)sector, (double)isolated, (double)isopair, (double)checkmap, (double)trig};
            hn_2photon->Fill(fill_hn_2photon, weight_pythia);
          } // jclus
      } // checkmap
    } // check photon1
  } // iclus

  /* Loop over EMC truth tracks */
  int nemctrk = emctrkcont->size();
  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *anatrk = new AnaTrk(emctrk, emccluscont);
    if( anatrk && anatrk->cid >= 0 && !anatrk->decayed )
    {
      int isPhoton = anatrk->pid == PHOTON_PID ? 1 : 0;
      int passProb = anatrk->prob_photon > probMin ? 1 : 0;
      double fill_hn_prob_photon[] =
      {anatrk->cluspt, (double)anatrk->sector, (double)isPhoton, (double)passProb};
      hn_prob_photon->Fill(fill_hn_prob_photon, weight_pythia);
    }
  }

  /* Loop over cgl reco tracks */
  int npart = data_tracks->get_npart();
  for(int itrk=0; itrk<npart; itrk++)
  {
    unsigned quality = data_tracks->get_quality(itrk);
    double mom = data_tracks->get_mom(itrk);
    if( quality <= 3 || !TMath::Finite(mom) ||
        mom < pTrkMin || mom > pTrkMax )
      continue;

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

    double dcdphi = data_tracks->get_emcdphi(itrk);
    double dcdz = data_tracks->get_emcdz(itrk);
    double fill_hn_dcdphiz[] = {dcdphi, dcdz, mom, (double)dcns, (double)dcwe, (double)ERTLvl1};
    hn_dcdphiz->Fill(fill_hn_dcdphiz, weight_pythia);

    double fill_hn_alphaboard[] = {dcboard, dcalpha, (double)dcns, (double)dcwe, (double)ERTLvl1};
    hn_alphaboard->Fill(fill_hn_alphaboard, weight_pythia);

    dcdeadmap->Checkmap(true);
    int isDCGood = dcdeadmap->IsDead(data_tracks, itrk) ? 0 : 1;
    double fill_hn_dclive[] = {dczed, dcphi, mom, (double)isDCGood, (double)ERTLvl1};
    hn_dclive->Fill(fill_hn_dclive, weight_pythia);

    int charge = data_tracks->get_charge(itrk);
    double px = data_tracks->get_px(itrk);
    double py = data_tracks->get_py(itrk);
    double pt = sqrt( px*px + py*py );
    double prod = -charge * pt * dcalpha;
    h_prod->Fill(prod);
  }

  return EVENT_OK;
}

int HadronResponse::End(PHCompositeNode *topNode)
{
  /* Write histogram output to ROOT file */
  hm->dumpHistos();
  delete hm;
  delete ertsim;
  delete emcwarnmap;
  delete dcdeadmap;

  return EVENT_OK;
}

void HadronResponse::BookHistograms()
{
  /* Initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

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

  /* Count events */
  h_events = new TH1F("h_events", "Events count", 1, 0.5, 1.5);
  hm->registerHisto(h_events);

  /* Eta and phi distribution */
  // ih = part + 3*isolated + 3*2*checkmap + 3*2*2*trig < 3*2*2*2
  for(int ih=0; ih<nh_eta_phi; ih++)
  {
    h2_eta_phi[ih] = new TH2F(Form("h2_eta_phi_%d",ih), "#eta and #phi distribution;#eta;#phi;", neta,etabin[ih%3/2], nphi,phibin);
    hm->registerHisto(h2_eta_phi[ih]);
  }

  /* EMCal reco one photon */
  int nbins_hn_1photon[] = {npT, 8, 2, 2, 2};
  double xmin_hn_1photon[] = {0., -0.5, -0.5, -0.5, -0.5};
  double xmax_hn_1photon[] = {0., 7.5, 1.5, 1.5, 1.5};
  hn_1photon = new THnSparseF("hn_1photon", "EMCal one photon;p_{T} [GeV];Sector;Isolated;Checkmap;Trig;",
      5, nbins_hn_1photon, xmin_hn_1photon, xmax_hn_1photon);
  hn_1photon->SetBinEdges(0, pTbin);
  hn_1photon->Sumw2();
  hm->registerHisto(hn_1photon);

  /* EMCal reco two photons */
  int nbins_hn_2photon[] = {npT, npT, 300, 8, 2, 2, 2, 2};
  double xmin_hn_2photon[] = {0., 0., 0., -0.5, -0.5, -0.5, -0.5, -0.5};
  double xmax_hn_2photon[] = {0., 0., 0.3, 7.5, 1.5, 1.5, 1.5, 1.5};
  hn_2photon = new THnSparseF("hn_2photon", "EMCal two photon;p^{1photon}_{T} [GeV];p^{2photon}_{T} [GeV];m_{inv} [GeV];Sector;Isolated;Isopair;Checkmap;Trig;",
      8, nbins_hn_2photon, xmin_hn_2photon, xmax_hn_2photon);
  hn_2photon->SetBinEdges(0, pTbin);
  hn_2photon->SetBinEdges(1, pTbin);
  hn_2photon->Sumw2();
  hm->registerHisto(hn_2photon);

  /* EMCal prob_photon */
  int nbins_hn_prob_photon[] = {npT, 8, 2, 2};
  double xmin_hn_prob_photon[] = {0., -0.5, -0.5, -0.5};
  double xmax_hn_prob_photon[] = {0., 7.5, 1.5, 1.5};
  hn_prob_photon = new THnSparseF("hn_prob_photon", "EMCal prob_photon;p_{T} [GeV];Sector;IsPhoton;PassProb;",
      4, nbins_hn_prob_photon, xmin_hn_prob_photon, xmax_hn_prob_photon);
  hn_prob_photon->SetBinEdges(0, pTbin);
  hn_prob_photon->Sumw2();
  hm->registerHisto(hn_prob_photon);

  /* DC dphi and dz */
  int nbins_hn_dcdphiz[] = {100, 100, 30, 2, 2, 2};
  double xmin_hn_dcdphiz[] = {-0.1, -50., 0., -0.5, -0.5, -0.5};
  double xmax_hn_dcdphiz[] = {0.1, 50., 15, 1.5, 1.5, 1.5};
  hn_dcdphiz = new THnSparseF("hn_dcdphiz", "EMCal and DC phi and z deviation;dphi [rad];dz [cm];mom [GeV];NS;WE;ERT;",
      6, nbins_hn_dcdphiz, xmin_hn_dcdphiz, xmax_hn_dcdphiz);
  hn_dcdphiz->Sumw2();
  hm->registerHisto(hn_dcdphiz);

  /* DC alpha-board */
  int nbins_hn_alphaboard[] = {82*4, 120, 2, 2, 2};
  double xmin_hn_alphaboard[] = {-1.5, -0.6, -0.5, -0.5, -0.5};
  double xmax_hn_alphaboard[] = {80.5, 0.6, 1.5, 1.5, 1.5};
  hn_alphaboard = new THnSparseF("hn_alphaboard", "DC #alpha-board;board;#alpha;NS;WE;ERT;",
      5, nbins_hn_alphaboard, xmin_hn_alphaboard, xmax_hn_alphaboard);
  hn_alphaboard->Sumw2();
  hm->registerHisto(hn_alphaboard);

  /* DC live area */
  int nbins_hn_dclive[] = {200, 50, 30, 2, 2};
  double xmin_hn_dclive[] = {-100., -1., 0., -0.5, -0.5};
  double xmax_hn_dclive[] = {100., 4., 15., 1.5, 1.5};
  hn_dclive = new THnSparseF("hn_dclive", "DC zed and phi distribution;zed [cm];phi [rad];mom [GeV];isDCGood;ERT;",
      5, nbins_hn_dclive, xmin_hn_dclive, xmax_hn_dclive);
  hn_dclive->Sumw2();
  hm->registerHisto(hn_dclive);

  h_prod = new TH1F("h_prod", "Minus charge times pT times alpha;-charge*pT*#alpha [GeV]", 200,0.,0.2);
  hm->registerHisto(h_prod);

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
        emcwarnmap->IsBadTower(clus2) ||
        //fabs( clus2->tofcorr() ) > tofMaxIso ||
        clus2->ecore() < eClusMin ||
        dcdeadmap->ChargeVeto(clus2, data_tracks) )
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

void HadronResponse::SumEEmcal(const emcClusterContent *cluster, const emcClusterContent *cluster_part,
    const emcClusterContainer *cluscont, const PHCentralTrack *data_tracks, double &econe)
{ 
  /* Sum up all energy in cone around particle without the partner cluster */
  econe = 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(cluster);
  if( pE_pref.Pt() < epsilon ) return;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nclus = cluscont->size();

  for (int iclus=0; iclus<nclus; iclus++)
  {
    emcClusterContent *clus3 = cluscont->getCluster(iclus);

    /* Skip if pointer identical to any of the two 'reference' particles
     * or on bad towers or lower than energy threshold 
     * and with 3 sigma charge veto */
    if( clus3->id() == cluster->id() ||
        clus3->id() == cluster_part->id() ||
        emcwarnmap->IsBadTower(clus3) ||
        //fabs( clus3->tofcorr() ) > tofMaxIso ||
        clus3->ecore() < eClusMin ||
        dcdeadmap->ChargeVeto(clus3, data_tracks) )
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
      econe += clus3->ecore();
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
      econe += mom;
  }

  return econe;
}

bool HadronResponse::TestPhoton(const emcClusterContent *cluster)
{
  /* Do not use ToF cut. ToF is wrong for high pT clusters */
  if( cluster->ecore() > eMin &&
      //fabs( cluster->tofcorr() ) < tofMax &&
      cluster->prob_photon() > probMin )
    return true;
  else
    return false;
}
