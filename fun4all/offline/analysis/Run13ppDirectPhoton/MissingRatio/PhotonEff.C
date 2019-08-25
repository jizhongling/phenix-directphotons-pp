#include "PhotonEff.h"

#include <AnaToolsTowerID.h>
#include <AnaToolsCluster.h>
#include <EMCWarnmapChecker.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <PHGlobal.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

#include <TF1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <iostream>
#include <fstream>

using namespace std;

/* Global constants */
const int PHOTON_PID = 1;
const int POSITRON_PID = 2;
const int ELECTRON_PID = 3;
const int PIZERO_PID = 7;

const double eMin = 0.3;
const double PI = TMath::Pi();

PhotonEff::PhotonEff(const string &name, const char *filename):
  SubsysReco(name),
  emcwarnmap(NULL),
  hm(NULL),
  hn_1photon(NULL)
{
  /* Construct output file names */
  outFileName = "histos/PhotonEff-";
  outFileName.append(filename);

  /* Function for pT weight for direct photon */
  cross = new TF1("cross", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]", 0, 30);
  cross->SetParameters(255., 5.98, 0.273, 14.43);

  /* Initialize histograms */
  for(int part=0; part<3; part++)
    h2_photon_eta_phi[part] = NULL;
}

PhotonEff::~PhotonEff()
{
}

int PhotonEff::Init(PHCompositeNode *topNode)
{
  /* Create and register histograms */
  BookHistograms();

  /* Initialize EMC warnmap checker */
  emcwarnmap = new EMCWarnmapChecker();
  if(!emcwarnmap)
  {
    cerr << "No emcwarnmap" << endl;
    exit(1);
  }

  return EVENT_OK;
}

int PhotonEff::process_event(PHCompositeNode *topNode)
{
  /* Global info */
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
  {
    cout << "Cannot find PHGlobal" << endl;
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

  /* Get event global parameters */
  float bbc_t0 = data_global->getBbcTimeZero();

  /* Initialize weight of this event
   * It will depend on photon pT */
  double photonpt = 0.;
  double weight = 1.;

  /* Workaround to avoid "set but not used" compiler warning */
  if(weight);

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
        emcwarnmap->IsGoodTower(emccluster) )
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
  /* Write histogram output to ROOT file */
  hm->dumpHistos(outFileName);
  delete hm;
  delete emcwarnmap;

  return EVENT_OK;
}

void PhotonEff::BookHistograms()
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
