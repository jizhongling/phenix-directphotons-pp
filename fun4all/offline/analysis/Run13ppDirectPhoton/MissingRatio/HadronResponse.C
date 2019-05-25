#include "HadronResponse.h"

#include "AnaTrk.h"
#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"

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

// global constants
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
  hm(NULL),
  hn_dc(NULL),
  hn_emcal(NULL),
  hn_cluster(NULL)
{
  // construct output file names
  outFileName = "histos/HadronResponse-";
  outFileName.append(filename);

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
      {
        tower_status_nils[isector][ibiny][ibinz] = 0;
        tower_status_sasha[isector][ibiny][ibinz] = 0;
      }
}

HadronResponse::~HadronResponse()
{
}

int HadronResponse::Init(PHCompositeNode *topNode)
{
  // initialize histogram manager
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  // create and register histograms
  BookHistograms();

  // read warnmap
  ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  return EVENT_OK;
}

int HadronResponse::process_event(PHCompositeNode *topNode)
{
  // central track reco info
  PHCentralTrack *data_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!data_tracks)
  {
    cout << "Cannot find PHCentralTrack" << endl;
    return DISCARDEVENT;
  }

  // central track sim info
  McEvalSingleList *mctrk = findNode::getClass<McEvalSingleList>(topNode, "McSingle");
  if(!mctrk)
  {
    cout << "Cannot find McEvalSingleList" << endl;
    return DISCARDEVENT;
  }

  // emc track info
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  // emc cluster info
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

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

  // loop over cgl reco tracks
  for(int iPart=0; iPart<nPart; iPart++)
  {
    unsigned quality = data_tracks->get_quality(iPart);
    if(quality < 3) continue;

    double dczed = data_tracks->get_zed(iPart);
    int dcarm = data_tracks->get_dcarm(iPart);
    int dcns = dczed > 0. ? 0 : 1;
    int dcwe = dcarm == 0 ? 1 : 0;

    double pxReco = data_tracks->get_px(iPart);
    double pyReco = data_tracks->get_py(iPart);
    double pzReco = data_tracks->get_pz(iPart);
    double ptReco = sqrt( pxReco*pxReco + pyReco*pyReco );
    double momReco = sqrt( pxReco*pxReco + pyReco*pyReco + pzReco*pzReco );

    double fill_hn_dc[] = {ptReco, momReco, (double)dcns, (double)dcwe, 1.};
    hn_dc->Fill(fill_hn_dc);
  }

  // loop over cgl truth tracks
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
    hn_dc->Fill(fill_hn_dc);
  }

  // associate cluster with track
  // map key is trkno
  typedef map<int,AnaTrk*> map_Ana_t;
  map_Ana_t track_list;

  // store photon
  vector<AnaTrk*> photon;

  // number of tracks and clusters
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

  // analyze emc clusters
  for(int iclus=0; iclus<nemcclus; iclus++)
  {
    emcClusterContent *cluster = emccluscont->getCluster(iclus);
    if( IsGoodTower(cluster) &&
        cluster->ecore() > eMin )
    {
      int sector = anatools::GetSector(cluster); 
      double pT = anatools::Get_pT(cluster);

      double fill_hn_emcal[] = {pT, cluster->ecore(), (double)sector, 1.};
      hn_emcal->Fill(fill_hn_emcal);
    }
  }

  // analyze emc tracks
  BOOST_FOREACH( map_Ana_t::value_type &trk_el, track_list )
  {
    AnaTrk *trk = trk_el.second;
    if( trk->cid < 0 ||
        trk->trkedep < eMin )
      continue;

    double fill_hn_emcal[] = {trk->trkpt, trk->emctrk->get_ekin(), (double)trk->sector, 0.};
    hn_emcal->Fill(fill_hn_emcal);

    if( trk->pid != PHOTON_PID )
      continue;

    double econeEM = SumEEmcal(trk->emcclus, emccluscont, data_tracks);
    double econeTrk = SumPTrack(trk->emcclus, data_tracks);
    double econeReco = econeEM + econeTrk;

    // loop over second emc tracks
    double econeMC = 0.;
    BOOST_FOREACH( map_Ana_t::value_type &trk2_el, track_list )
    {
      AnaTrk *trk2 = trk2_el.second;
      if( trk2->trkno != trk->trkno &&
          trk->trkvp.Angle(trk2->trkvp) < cone_angle )
        econeMC += trk2->emctrk->get_ekin();
    }

    double fill_hn_cluster[] = {trk->trkpt, trk->cluspt, econeMC, econeReco, (double)trk->sector};
    hn_cluster->Fill(fill_hn_cluster);
  }

  // clear associated list
  //delete[] cglMC;
  BOOST_FOREACH( map_Ana_t::value_type &trk_el, track_list )
    delete trk_el.second;

  return EVENT_OK;
}

int HadronResponse::End(PHCompositeNode *topNode)
{
  // write histogram output to ROOT file
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void HadronResponse::BookHistograms()
{
  // for dc track
  int nbins_hn_dc[] = {npT, 300, 2, 2, 2};
  double xmin_hn_dc[] = {0., 0., -0.5, -0.5, -0.5};
  double xmax_hn_dc[] = {0., 30.,1.5, 1.5, 1.5};
  hn_dc = new THnSparseF("hn_dc", "DC track info;p_{T} [GeV];E [GeV];NS;WE;Reco;",
      5, nbins_hn_dc, xmin_hn_dc, xmax_hn_dc);
  hn_dc->SetBinEdges(0, vpT);
  hm->registerHisto(hn_dc);

  // for emcal cluster
  int nbins_hn_emcal[] = {npT, 300, 8, 2};
  double xmin_hn_emcal[] = {0., 0., -0.5, -0.5};
  double xmax_hn_emcal[] = {0., 30., 7.5, 1.5};
  hn_emcal = new THnSparseF("hn_emcal", "EMCal info;p_{T} [GeV];E [GeV];Sector;Reco;",
      4, nbins_hn_emcal, xmin_hn_emcal, xmax_hn_emcal);
  hn_emcal->SetBinEdges(0, vpT);
  hm->registerHisto(hn_emcal);

  // for isolation cone 
  int nbins_hn_cluster[] = {npT, npT, 400, 400, 8};
  double xmin_hn_cluster[] = {0., 0., 0., 0., -0.5};
  double xmax_hn_cluster[] = {0., 0., 40., 40., 7.5};
  hn_cluster = new THnSparseF("hn_cluster", "Isolation cone info;p^{truth}_{T} [GeV];p^{reco}_{T} [GeV];E_{truth} [GeV];E_{reco} [GeV];Sector;",
      5, nbins_hn_cluster, xmin_hn_cluster, xmax_hn_cluster);
  hn_cluster->SetBinEdges(0, vpT);
  hn_cluster->SetBinEdges(1, vpT);
  hm->registerHisto(hn_cluster);

  return;
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
    // count tower with bad status for PbSc and PbGl
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
    tower_status_sasha[sector][biny][binz] = status;

    // mark edge towers
    if( anatools::Edge_cg(sector, biny, binz) )
      tower_status_sasha[sector][biny][binz] = 20;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
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
     * or on bad towers or lower than energy threshold */
    if( clus2->id() == cluster->id() ||
        IsBadTower(clus2) ||
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
    if( quality <= 3 )
      continue;

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
      econe += mom;
  }

  return econe;
}
