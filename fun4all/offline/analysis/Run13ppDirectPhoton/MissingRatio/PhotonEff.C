#include "PhotonEff.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsCluster.h"

#include <PHGlobal.h>
//#include <PHCentralTrack.h>
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
#include <TF1.h>
#include <TH2.h>
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
const double PI = TMath::Pi();

PhotonEff::PhotonEff(const string &name, const char *filename):
  SubsysReco(name),
  hm(NULL),
  hn_1photon(NULL)
{
  // construct output file names
  outFileName = "histos/PhotonEff-";
  outFileName.append(filename);

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;

  // function for pT weight for direct photon
  cross = new TF1("cross", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]", 0, 30);
  cross->SetParameters(255., 5.98, 0.273, 14.43);

  // initialize histograms
  for(int part=0; part<3; part++)
    h2_photon_eta_phi[part] = NULL;
}

PhotonEff::~PhotonEff()
{
}

int PhotonEff::Init(PHCompositeNode *topNode)
{
  // create and register histograms
  BookHistograms();

  // read warnmap
  //ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("warn_all_run13pp500gev.dat");

  return EVENT_OK;
}

int PhotonEff::process_event(PHCompositeNode *topNode)
{
  // global info
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
  {
    cout << "Cannot find PHGlobal" << endl;
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

  /* Get event global parameters */
  float bbc_t0 = data_global->getBbcTimeZero();

  // initialize weight of this event
  // it will depend on photon pT
  double photonpt = 0.;
  double weight = 1.;

  // workaround to avoid "set but not used" compiler warning
  if (weight) {};

  unsigned nemctrk = emctrkcont->size();
  for(unsigned itrk=0; itrk<nemctrk; itrk++)
  {
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    if( emctrk->get_pid() == PHOTON_PID && emctrk->get_anclvl() == 0 )
    {
      photonpt = emctrk->get_pt();
      if( photonpt > 1. )
        weight = cross->Eval(photonpt);
      else
        weight = cross->Eval(1.);
      break;
    }
  }

  unsigned nemccluster = emccluscont->size();
  for(unsigned i=0; i<nemccluster; i++)
  {
    emcClusterContent *emccluster = emccluscont->getCluster(i);
    if( emccluster->ecore() > eMin &&
        GetStatus(emccluster) == 0 )
    {
      int sector = anatools::GetSector(emccluster);
      int part = -1;
      if(sector < 0) part = -1;
      else if(sector < 4) part = 0;
      else if(sector < 6) part = 1;
      else if(sector < 8) part = 2;

      TLorentzVector pE = anatools::Get_pE(emccluster);
      double recopt = pE.Pt();
      double eta = pE.Eta();
      double phi = pE.Phi();
      if(sector >= 4)
      {
        pE.RotateZ(-PI);
        phi = pE.Phi() + PI;
      }

      double tof = emccluster->tofcorr();
      double prob = emccluster->prob_photon();

      double fill_hn_1photon[] = {(double)sector, photonpt, recopt, 0., tof, tof-bbc_t0};
      hn_1photon->Fill(fill_hn_1photon);
      if( fabs(tof - bbc_t0) < 10. )
      {
        fill_hn_1photon[3] = 1.;
        hn_1photon->Fill(fill_hn_1photon);
      }
      if( prob > 0.02 )
      {
        fill_hn_1photon[3] = 2.;
        hn_1photon->Fill(fill_hn_1photon);
      }
      if( fabs(tof - bbc_t0) < 10. && prob > 0.02 )
      {
        if( part >= 0 && recopt > 5. && recopt < 10. )
          h2_photon_eta_phi[part]->Fill(eta, phi);

        fill_hn_1photon[3] = 3.;
        hn_1photon->Fill(fill_hn_1photon);
      }
    }
  }

  return EVENT_OK;
}

int PhotonEff::End(PHCompositeNode *topNode)
{
  // write histogram output to ROOT file
  hm->dumpHistos(outFileName);
  delete hm;

  return EVENT_OK;
}

void PhotonEff::BookHistograms()
{
  // initialize histogram manager
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

  /* Store single photon information */
  const int nbins_hn_1photon[] = {8, npT, npT, 4, 60, 60};
  const double xmin_hn_1photon[] = {-0.5, 0., 0., -0.5, -30., -30.};
  const double xmax_hn_1photon[] = {7.5, 0., 0., 3.5, 30., 30.};
  hn_1photon = new THnSparseF("hn_1photon", "Single photon spectrum;sector;p^{truth}_{T} [GeV];reco p^{reco}_{T} [GeV];cut;tof [ns];tofdiff [ns];",
      6, nbins_hn_1photon, xmin_hn_1photon, xmax_hn_1photon);
  hn_1photon->SetBinEdges(1, pTbin);
  hn_1photon->SetBinEdges(2, pTbin);
  hm->registerHisto(hn_1photon);

  for(Int_t part=0; part<3; part++)
  {
    h2_photon_eta_phi[part] = new TH2F(Form("h2_photon_eta_phi_part%d",part), "Photon #eta and #phi distribution;#eta;#phi;", neta,etabin[part/2], nphi,phibin);
    h2_photon_eta_phi[part]->Sumw2();
    hm->registerHisto(h2_photon_eta_phi[part]);
  }

  return;
}

void PhotonEff::ReadTowerStatus(const string& filename)
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
    tower_status[sector][biny][binz] = status;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void PhotonEff::ReadSashaWarnmap(const string &filename)
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

int PhotonEff::GetStatus(const emcClusterContent *emccluster)
{
  int arm = emccluster->arm();
  int rawsector = emccluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = emccluster->iypos();
  int izpos = emccluster->izpos();

  return tower_status[sector][iypos][izpos];
}
