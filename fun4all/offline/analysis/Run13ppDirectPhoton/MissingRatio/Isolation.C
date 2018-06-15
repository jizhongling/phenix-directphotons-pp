#include "Isolation.h"

#include "AnaTrk.h"
#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"

#include <PHGlobal.h>
#include <PHCentralTrack.h>
//#include <PHSnglCentralTrack.h>
//#include <McEvalSingleList.h>

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

#include <TMath.h>
#include <THnSparse.h>

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// global constants
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
  // construct output file names
  outFileName = "histos/Isolation-";
  outFileName.append(filename);

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;
}

Isolation::~Isolation()
{
}

int Isolation::Init(PHCompositeNode *topNode)
{
  // create and register histograms
  BookHistograms();

  // read warnmap
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  return EVENT_OK;
}

int Isolation::process_event(PHCompositeNode *topNode)
{
  // global info
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
  {
    cout << "Cannot find PHGlobal" << endl;
    return DISCARDEVENT;
  }

  // central tracks
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!data_tracks)
  {
    cerr << "Cannot find PHCentralTrack" << endl;
    return DISCARDEVENT;
  }

  // track info
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  // cluster info
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

  // number of tracks
  int nemctrk = emctrkcont->size();

  // Loop over all emcGeaTrack
  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    // Associate cluster to track
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *anatrk = new AnaTrk(emctrk, emccluscont, (int*)tower_status);
    if(!anatrk) continue;

    // Get the associated cluster
    emcGeaClusterContent *emcclus = anatrk->emcclus;

    // Test if particle is photon and associated cluster
    // and on good towers
    if( !emcclus ||
        anatrk->cid < 0 ||
        anatrk->pid != PHOTON_PID )
    {
      delete anatrk;
      continue;
    }

    // Test if particle is prompt photon
    int prompt = 0;
    if( anatrk->anclvl == 0 )
      prompt = 1;

    // Fill for different cone angle and energy fraction
    for(int icone=0; icone<10; icone++)
      for(int ie=0; ie<20; ie++)
      {
        double rcone = (icone+1) * 0.1;
        double econe = SumEEmcal(anatrk, rcone) + SumPTrack(anatrk, data_tracks, rcone);
        double re = (ie+1) * 0.02;
        int isolated = 0;
        if( econe < re * emcclus->ecore() )
          isolated = 1;

        double fill_hn_photon[] = {anatrk->cluspt, rcone, re, isolated, prompt};
        hn_photon->Fill(fill_hn_photon);
      } // icone, ie

    delete anatrk;
  } // anatrk

  return EVENT_OK;
}

int Isolation::End(PHCompositeNode *topNode)
{
  // write histogram output to ROOT file
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void Isolation::BookHistograms()
{
  // Initialize histogram manager
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  // pT bins 
  const int npT = 30;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  const int nbins_hn_photon[] = {npT, 10, 20, 2, 2};
  const double xmin_hn_photon[] = {0., 0.05, 0.01, -0.5, -0.5};
  const double xmax_hn_photon[] = {0., 1.05, 0.41, 1.5, 1.5};
  hn_photon = new THnSparseF("hn_photon", "Photon counts;p_{T} [GeV];rcone;renergy;isolated;prompt;",
      5, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, pTbin);
  hm->registerHisto(hn_photon);

  return;
}

void Isolation::ReadSashaWarnmap(const string &filename)
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

double Isolation::SumEEmcal(const AnaTrk *anatrk, double rcone)
{
  // Sum up all energy in cone around particle
  double econe = 0;

  // Get associated cluster and its container
  emcGeaClusterContainer *emccluscont = anatrk->emccluscont;
  emcGeaClusterContent *emcclus_pref = anatrk->emcclus;
  if(!emccluscont || !emcclus_pref)
    return 0.;

  // Get reference vector
  TLorentzVector pE_pref = anatools::Get_pE(emcclus_pref);
  if( pE_pref.Pt() < 0.01 ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int nemcclus = emccluscont->size();

  for (int iclus=0; iclus<nemcclus; iclus++)
  {
    emcGeaClusterContent *emcclus2 = emccluscont->getCluster(iclus);

    // Skip if pointer identical to 'reference' particle
    // or on bad towers or lower than energy threshold
    if( emcclus2->id() == anatrk->cid ||
        anatrk->cid < 0 ||
        emcclus2->ecore() < 0.15 )
      continue;

    // Get cluster vector
    TLorentzVector pE_part2 = anatools::Get_pE(emcclus2);
    if( pE_part2.Pt() < 0.01 ) continue;
    TVector2 v2_part2 = pE_part2.EtaPhiVector();

    // Check if cluster within cone
    if( (v2_part2-v2_pref).Mod() < rcone )
      econe += emcclus2->ecore();
  }

  return econe;
}

double Isolation::SumPTrack(const AnaTrk *anatrk, const PHCentralTrack *tracks, double rcone)
{
  // Sum up all energy in cone around particle
  double econe = 0;

  // Get associated cluster
  emcGeaClusterContent *emcclus_pref = anatrk->emcclus;
  if(!emcclus_pref)
    return 0.;

  // Get reference vector
  TLorentzVector pE_pref = anatools::Get_pE(emcclus_pref);
  if( pE_pref.Pt() < 0.01 ) return econe;
  TVector2 v2_pref = pE_pref.EtaPhiVector();

  int ntrk = tracks->get_npart();

  for (int itrk=0; itrk<ntrk; itrk++)
  {
    double px = tracks->get_px(itrk);
    double py = tracks->get_py(itrk);
    double pz = tracks->get_pz(itrk);
    double mom = tracks->get_mom(itrk);
    //cout << px << "\t" << py << "\t" << pz << "\t" << mom << endl;

    // Test if track passes momentum cuts
    if( mom < 0.2 || mom > 15. )
      continue;

    // Get track vector
    TVector3 v3_part2(px, py, pz);
    if( v3_part2.Pt() < 0.01 ) continue;
    TVector2 v2_part2 = v3_part2.EtaPhiVector();

    // Check if particle within cone
    if( (v2_part2-v2_pref).Mod() < rcone )
      econe += mom;
  }

  return econe;
}
