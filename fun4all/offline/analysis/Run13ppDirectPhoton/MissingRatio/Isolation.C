#include "Isolation.h"

#include "AnaTrk.h"
#include <AnaToolsTowerID.h>
#include <AnaToolsCluster.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <PHGlobal.h>
#include <PHCentralTrack.h>
//#include <PHSnglCentralTrack.h>
//#include <McEvalSingleList.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TMath.h>
#include <THnSparse.h>

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>

using namespace std;

/* Some constants */
const double epsilon = TMath::Limits<float>::Epsilon();
const double PI = TMath::Pi();
const TVector2 v2_2PI(0., 2.*PI);

/* Global constants */
const int PHOTON_PID = 1;
const int POSITRON_PID = 2;
const int ELECTRON_PID = 3;
const int PIZERO_PID = 7;

const double eMin = 0.3;
const double AsymCut = 0.8;

Isolation::Isolation(const string &name, const char *filename):
  SubsysReco(name),
  hm(NULL),
  hn_photon(NULL)
{
  /* Construct output file names */
  outFileName = "histos/Isolation-";
  outFileName.append(filename);
}

Isolation::~Isolation()
{
}

int Isolation::Init(PHCompositeNode *topNode)
{
  /* Create and register histograms */
  BookHistograms();

  return EVENT_OK;
}

int Isolation::process_event(PHCompositeNode *topNode)
{
  /* Global info */
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
  {
    cout << "Cannot find PHGlobal" << endl;
    return DISCARDEVENT;
  }

  /* Central tracks */
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!data_tracks)
  {
    cerr << "Cannot find PHCentralTrack" << endl;
    return DISCARDEVENT;
  }

  /* Track info */
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  /* Cluster info */
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

  /* Number of tracks */
  int nemctrk = emctrkcont->size();

  /* Loop over all emcGeaTrack */
  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    /* Associate cluster to track */
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *anatrk = new AnaTrk(emctrk, emccluscont);
    if(!anatrk) continue;

    /* Get the associated cluster */
    emcGeaClusterContent *emcclus = anatrk->emcclus;

    /* Test if particle is photon and associated cluster
     * and on good towers */
    if( !emcclus ||
        anatrk->cid < 0 ||
        anatrk->pid != PHOTON_PID )
    {
      delete anatrk;
      continue;
    }

    /* Test if particle is prompt photon */
    int prompt = 0;
    if( anatrk->anclvl == 0 )
      prompt = 1;

    /* Fill for different cone angle and energy fraction */
    for(int icone=0; icone<11; icone++)
      for(int ie=0; ie<20; ie++)
      {
        double rcone = icone * 0.1;
        double econe = SumEEmcal(anatrk, rcone) + SumPTrack(anatrk, data_tracks, rcone);
        double re = (ie+1) * 0.02;
        int isolated = 0;
        if( econe < re * emcclus->ecore() )
          isolated = 1;

        double fill_hn_photon[] = {anatrk->cluspt, rcone, re, (double)isolated, (double)prompt};
        hn_photon->Fill(fill_hn_photon);
      } // icone, ie

    delete anatrk;
  } // anatrk

  return EVENT_OK;
}

int Isolation::End(PHCompositeNode *topNode)
{
  /* Write histogram output to ROOT file */
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void Isolation::BookHistograms()
{
  /* Initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  /* pT bins  */
  const int npT = 30;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  const int nbins_hn_photon[] = {npT, 11, 20, 2, 2};
  const double xmin_hn_photon[] = {0., -0.05, 0.01, -0.5, -0.5};
  const double xmax_hn_photon[] = {0., 1.05, 0.41, 1.5, 1.5};
  hn_photon = new THnSparseF("hn_photon", "Photon counts;p_{T} [GeV];rcone;renergy;isolated;prompt;",
      5, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, pTbin);
  hm->registerHisto(hn_photon);

  return;
}

double Isolation::SumEEmcal(const AnaTrk *anatrk, double rcone)
{
  /* Sum up all energy in cone around particle */
  double econe = 0;

  /* Get associated cluster and its container */
  emcGeaClusterContainer *emccluscont = anatrk->emccluscont;
  emcGeaClusterContent *emcclus_pref = anatrk->emcclus;
  if(!emccluscont || !emcclus_pref)
    return 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(emcclus_pref);
  if( pE_pref.Pt() < epsilon ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nemcclus = emccluscont->size();

  for (int iclus=0; iclus<nemcclus; iclus++)
  {
    emcGeaClusterContent *emcclus2 = emccluscont->getCluster(iclus);

    /* Skip if pointer identical to 'reference' particle
     * or on bad towers or lower than energy threshold */
    if( emcclus2->id() == anatrk->cid ||
        anatrk->cid < 0 ||
        emcclus2->ecore() < 0.15 )
      continue;

    /* Get cluster vector */
    TLorentzVector pE_part2 = anatools::Get_pE(emcclus2);
    if( pE_part2.Pt() < epsilon ) continue;
    TVector2 v2_part2 = pE_part2.EtaPhiVector();

    /* Check if cluster within cone */
    TVector2 v2_diff(v2_part2 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < rcone )
      econe += emcclus2->ecore();
  }

  return econe;
}

double Isolation::SumPTrack(const AnaTrk *anatrk, const PHCentralTrack *tracks, double rcone)
{
  /* Sum up all energy in cone around particle */
  double econe = 0;

  /* Get associated cluster */
  emcGeaClusterContent *emcclus_pref = anatrk->emcclus;
  if(!emcclus_pref)
    return 0.;

  /* Get reference vector */
  TLorentzVector pE_pref = anatools::Get_pE(emcclus_pref);
  if( pE_pref.Pt() < epsilon ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int ntrk = tracks->get_npart();

  for (int itrk=0; itrk<ntrk; itrk++)
  {
    double px = tracks->get_px(itrk);
    double py = tracks->get_py(itrk);
    double pz = tracks->get_pz(itrk);
    double mom = tracks->get_mom(itrk);
    //cout << px << "\t" << py << "\t" << pz << "\t" << mom << endl;

    /* Test if track passes momentum cuts */
    if( mom < 0.2 || mom > 15. )
      continue;

    /* Get track vector */
    TVector3 v3_part2(px, py, pz);
    if( v3_part2.Pt() < epsilon ) continue;
    TVector2 v2_part2 = v3_part2.EtaPhiVector();

    /* Check if particle within cone */
    TVector2 v2_diff(v2_part2 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < rcone )
      econe += mom;
  }

  return econe;
}
