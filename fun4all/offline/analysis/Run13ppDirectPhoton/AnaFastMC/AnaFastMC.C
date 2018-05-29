#include "AnaFastMC.h"

#include "AnaToolsTowerID.h"

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <Fun4AllReturnCodes.h>
#include <getClass.h>
#include <TOAD.h>
#include <phool.h>
#include <Fun4AllHistoManager.h>

#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>

#include <TPythia6.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8) 
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TGenPhaseSpace.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

const double PI = TMath::Pi();
const double mPi0 = 0.1349770;

const double eMin = 0.3;
const double AsymCut = 0.8;

AnaFastMC::AnaFastMC(const string &name):
  SubsysReco(name),
  outFileName("histos/AnaFastMC-"),
  mcmethod(PHParticleGen),
  warnmap(Nils),
  phpythiaheader(NULL),
  phpythia(NULL),
  hm(NULL),
  h_pion(NULL),
  h_photon(NULL),
  hn_missing(NULL),
  hn_pion(NULL),
  hn_photon(NULL)
{
  // initialize histograms
  for(int part=0; part<3; part++)
  {
    h2_pion_eta_phi[part] = NULL;
    h2_photon_eta_phi[part] = NULL;
  }

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;

  // function for pT weight for pi0
  cross_pi0 = new TF1("cross_pi0", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross_pi0->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  // function for pT weight for direct photon
  cross_ph = new TF1("cross_ph", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]", 0, 30);
  cross_ph->SetParameters(255., 5.98, 0.273, 14.43);

  NPart = 0;
  NPeak = 0;
  for(int i=0; i<MAXPEAK; i++)
  {
    Vpart[i].SetPxPyPzE(0.,0.,0.,0.);
    itw_part[i] = -1;
  }
}

AnaFastMC::~AnaFastMC()
{
}

int AnaFastMC::Init(PHCompositeNode *topNode)
{
  // create and register histograms
  BookHistograms();

  // read warnmap
  if( warnmap == Nils )
    ReadTowerStatus("Warnmap_Run13pp510.txt");
  else if( warnmap == Sasha )
    ReadSashaWarnmap("warn_all_run13pp500gev.dat");
  else if( warnmap == Sim )
    ReadSimWarnmap("dead_eff_run13pp500gev.dat");

  return EVENT_OK;
}

int AnaFastMC::process_event(PHCompositeNode *topNode)
{
  /* initialize parameters for pi0 decay */
  static TLorentzVector beam_pi0;
  static TLorentzVector *pG1;
  static TLorentzVector *pG2;
  vector<TLorentzVector> incident;

  static TGenPhaseSpace event;
  static double masses[2] = {0., 0.};

  double pt_pi0 = 0.;
  double weight_pi0 = 1.;

  /* initialize parameters for direct photon */
  static TLorentzVector beam_ph;

  double pt_ph = 0.;
  double weight_ph = 1.;

  /* Initialize parameters for clusters' information */
  int sector[MAXPEAK] = {};
  int iy[MAXPEAK] = {};
  int iz[MAXPEAK] = {};

  /* Select method to generate input particle */
  /* Use PHParticleGen */
  if( mcmethod == PHParticleGen )
  {
    // Get PYTHIA Header
    phpythiaheader = findNode::getClass<PHPythiaHeader>(topNode,"PHPythiaHeader");
    if (!phpythiaheader)
    {
      cout << PHWHERE << "Unable to get PHPythiaHeader, is Node missing?" << endl;
      return ABORTEVENT;
    }

    // Get PYTHIA Particles
    phpythia = findNode::getClass<PHPythiaContainer>(topNode,"PHPythia");
    if (!phpythia)
    {
      cout << PHWHERE << "Unable to get PHPythia, is Node missing?" << endl;
      return ABORTEVENT;
    }

    int npart = phpythia->size();
    for (int ipart=0; ipart<npart; ipart++)
    {
      TMCParticle *part = phpythia->getParticle(ipart);
      double px = part->GetPx();
      double py = part->GetPy();
      double pz = part->GetPz();
      double energy = part->GetEnergy();

      if( part->GetKF() == 111 )
      {
        beam_pi0.SetPxPyPzE(px,py,pz,energy);
        pt_pi0 = anatools::GetPt(part);
      }
      else if( part->GetKF() == 22 )
      {
        beam_ph.SetPxPyPzE(px,py,pz,energy);
        pt_ph = anatools::GetPt(part);
      }
    }
  }

  /* OR: Use FastMC (i.e. random number generator for pi0 and direct photon eta, phi, pt) */
  else if( mcmethod == FastMC )
  {
    const double eta_max = 0.5;
    double phi = 2.0*PI*gRandom->Rndm();
    double eta = eta_max*2*gRandom->Rndm() - eta_max;

    const double emin = 0.; // min pi0 energy
    const double emax = 40.; // max pi0 energy
    double pt = gRandom->Rndm() * (emax-emin) + emin;
    double px = pt*cos(phi);
    double py = pt*sin(phi);

    // This is if we use uniform rapidity (eta here is y)
    double pz_pi0 = sqrt((pt*pt+mPi0*mPi0)*(exp(2*eta)-1)*(exp(2*eta)-1)/4./(exp(2*eta)));  // pz_pi0 = mt*gamma*beta = mt*sinh(eta)
    double pz_ph = sqrt(pt*pt*(exp(2*eta)-1)*(exp(2*eta)-1)/4./(exp(2*eta)));  // pz_ph = pt*gamma*beta = pt*sinh(eta)
    if( eta < 0 )
    {
      pz_pi0 = -pz_pi0;
      pz_ph = -pz_ph;
    }
    double energy_pi0 = sqrt(pz_pi0*pz_pi0+pt*pt+mPi0*mPi0);
    double energy_ph = sqrt(pz_ph*pz_ph+pt*pt);

    beam_pi0.SetPxPyPzE(px,py,pz_pi0,energy_pi0);
    pt_pi0 = pt;

    beam_ph.SetPxPyPzE(px,py,pz_ph,energy_ph);
    pt_ph = pt;
  }

  /* Let pi0 decay into two photons */
  event.SetDecay(beam_pi0, 2, masses);
  event.Generate();

  pG1 = event.GetDecay(0);
  pG2 = event.GetDecay(1);
  incident.push_back(*pG1);
  incident.push_back(*pG2);

  /* Get event weight for pi0 */
  if(pt_pi0 > 1.)
    weight_pi0 = cross_pi0->Eval(pt_pi0);
  else
    weight_pi0 = cross_pi0->Eval(1.);

  /* Fill histogram for all generated pi0 */
  h_pion->Fill( pt_pi0, weight_pi0 );

  /* Initialize parameters for pi0_sim */
  float ptsim_pi0 = 0.;
  float minv = 0.;
  float dist = 9999.;
  ResetTowerEnergy();

  /* determine if pi0 decay photons within a tower that is NOT flagged bad (set boolean 'acc_pi0')
   * and fill ptsim_pi0, minv, and dist values */
  bool acc_pi0 = pi0_sim( pG1, pG2, ptsim_pi0, minv, dist );

  /* BEGIN pi0 acceptance calculation */

  /* Loop over all clusters in calorimeter */
  for(int i=0; i<MAXPEAK; i++)
    if(itw_part[i] >= 0)
    {
      /* Get sector, iy and iz information */
      sector[i] = -1;
      iy[i] = -1;
      iz[i] = -1;
      anatools::TowerLocation(itw_part[i], sector[i], iy[i], iz[i]);
    }

  /* Missing one photon */
  if(NPart == 1)
  {
    /* parameters for clusters from pi0 decay photons */
    int sec1 = sector[0];
    double e1 = Vpart[0].E();

    if( e1 > eMin )
    {
      double pt_reco = Vpart[0].Pt();
      double fill_hn_missing[] = {pt_pi0, pt_reco, sec1, NPart, 1.};
      hn_missing->Fill(fill_hn_missing, weight_pi0);
    }
  }

  /* Two photons are in acceptance */
  if(acc_pi0)
  {
    /* Get NPeak */
    GetNpeak();

    /* parameters for two clusters from pi0 decay photons */
    int sec1 = sector[0];
    int sec2 = sector[1];
    double e1 = Vpart[0].E();
    double e2 = Vpart[1].E();

    /* for missing ratio */
    if(e1 > eMin && e2 > eMin)
      for(int iph=0; iph<2; iph++)
      {
        double pt_reco = Vpart[iph].Pt();
        double fill_hn_missing[] = {pt_pi0, pt_reco, sector[iph], NPart, NPeak};
        hn_missing->Fill(fill_hn_missing, weight_pi0);
      }

    /* Check in the same detector part
     * and pass the energy and asymmetry cuts */
    if( anatools::SectorCheck(sec1,sec2) &&
        e1 > eMin && e2 > eMin && abs(e1-e2)/(e1+e2) < AsymCut )
    {
      /* Fill eta and phi distribution */
      for(int iph=0; iph<2; iph++)
      {
        int part = -1;
        if(sec1 < 0) part = -1;
        else if(sec1 < 4) part = 0;
        else if(sec1 < 6) part = 1;
        else if(sec1 < 8) part = 2;
        
        double pt_reco = Vpart[iph].Pt();
        double eta = Vpart[iph].Eta();
        double phi = Vpart[iph].Phi();
        if(sec1 >= 4)
        {
          Vpart[iph].RotateZ(-PI);
          phi = Vpart[iph].Phi() + PI;
        }

        if( part >= 0 && pt_reco > 5. && pt_reco < 10. )
          h2_pion_eta_phi[part]->Fill(-eta, phi, weight_ph);
      }

      /* Fill pi0 histogram
       * NPeak=2: No merging for two photons
       * NPeak=1: Merging for two photons */
      double fill_hn_pion[] = {pt_pi0, ptsim_pi0, minv, sec1, NPeak};
      hn_pion->Fill(fill_hn_pion, weight_pi0);
    }
  } // acc_pi0

  /* Get event weight for direct photon */
  if(pt_ph > 1.)
    weight_ph = cross_ph->Eval(pt_ph);
  else
    weight_ph = cross_ph->Eval(1.);

  /* Fill histogram for all generated direct photons */
  h_photon->Fill( beam_ph.Pt(), weight_ph );

  /* Initialize parameters for photon_sim */
  float ptsim_ph = 0.;
  ResetTowerEnergy();

  /* determine if direct photons within a tower that is NOT flagged bad (set boolean 'acc_ph')
   * and fill ptsim */
  bool acc_ph = photon_sim( &beam_ph, ptsim_ph );

  /* BEGIN direct photon acceptance calculation */
  if(acc_ph)
  {
    /* Loop over all clusters in calorimeter */
    for(int i=0; i<MAXPEAK; i++)
      if(itw_part[i] >= 0)
      {
        /* Get sector, iy and iz information */
        sector[i] = -1;
        iy[i] = -1;
        iz[i] = -1;
        anatools::TowerLocation(itw_part[i], sector[i], iy[i], iz[i]);
      }

    /* parameters for direct photons */
    int sec1 = sector[0];
    double e1 = Vpart[0].E();

    /* Fill photon histogram */
    if(e1 > eMin)
    {
      int part = -1;
      if(sec1 < 0) part = -1;
      else if(sec1 < 4) part = 0;
      else if(sec1 < 6) part = 1;
      else if(sec1 < 8) part = 2;

      double pt_truth = beam_ph.Pt();
      double pt_reco = Vpart[0].Pt();
      double eta = Vpart[0].Eta();
      double phi = Vpart[0].Phi();
      if(sec1 >= 4)
      {
        Vpart[0].RotateZ(-PI);
        phi = Vpart[0].Phi() + PI;
      }

      if( part >= 0 && pt_reco > 5. && pt_reco < 10. )
        h2_photon_eta_phi[part]->Fill(-eta, phi, weight_ph);

      double fill_hn_photon[] = {pt_truth, pt_reco, sec1};
      hn_photon->Fill(fill_hn_photon, weight_ph);
    }
  } // acc_ph

  return EVENT_OK;
}

int AnaFastMC::End(PHCompositeNode *topNode)
{
  // write histogram output to ROOT file
  hm->dumpHistos(outFileName);
  delete hm;

  return EVENT_OK;
}


void AnaFastMC::BookHistograms()
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
  
  h_pion = new TH1F("h_pion", "Total pion count;p_{T} [GeV];", npT, pTbin);
  h_pion->Sumw2();
  hm->registerHisto(h_pion);

  h_photon = new TH1F("h_photon", "Total photon count;p_{T} [GeV];", npT, pTbin);
  h_photon->Sumw2();
  hm->registerHisto(h_photon);

  for(Int_t part=0; part<3; part++)
  {
    h2_pion_eta_phi[part] = new TH2F(Form("h2_pion_eta_phi_part%d",part), "Pion #eta and #phi distribution;#eta;#phi;", neta,etabin[part/2], nphi,phibin);
    h2_photon_eta_phi[part] = new TH2F(Form("h2_photon_eta_phi_part%d",part), "Photon #eta and #phi distribution;#eta;#phi;", neta,etabin[part/2], nphi,phibin);
    h2_pion_eta_phi[part]->Sumw2();
    h2_photon_eta_phi[part]->Sumw2();
    hm->registerHisto(h2_pion_eta_phi[part]);
    hm->registerHisto(h2_photon_eta_phi[part]);
  }

  int nbins_hn_missing[] = {npT, npT, 8, 3, 3};
  double xmin_hn_missing[] = {0., 0., -0.5, -0.5, -0.5};
  double xmax_hn_missing[] = {0., 0., 7.5, 2.5, 2.5};
  hn_missing = new THnSparseF("hn_missing", "#pi^{0} missing ratio;p^{#pi^{0}}_{T} [GeV];p^{#gamma}_{T} [GeV];sector;NPart;NPeak;",
      5, nbins_hn_missing, xmin_hn_missing, xmax_hn_missing);
  hn_missing->SetBinEdges(0, pTbin);
  hn_missing->SetBinEdges(1, pTbin);
  hn_missing->Sumw2();
  hm->registerHisto(hn_missing);

  int nbins_hn_pion[] = {npT, npT, 300, 8, 3};
  double xmin_hn_pion[] = {0., 0., 0., -0.5, -0.5};
  double xmax_hn_pion[] = {0., 0., 0.3, 7.5, 2.5};
  hn_pion = new THnSparseF("hn_pion", "EMCal pion count;Truth p_{T} truth [GeV];p_{T} reco [GeV];m_{inv} [GeV];sector;NPeak;",
      5, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(0, pTbin);
  hn_pion->SetBinEdges(1, pTbin);
  hn_pion->Sumw2();
  hm->registerHisto(hn_pion);

  int nbins_hn_photon[] = {npT, npT, 8};
  double xmin_hn_photon[] = {0., 0., -0.5};
  double xmax_hn_photon[] = {0., 0., 7.5};
  hn_photon = new THnSparseF("hn_photon", "EMCal photon count;p_{T} truth [GeV/c];p_{T} reco [GeV];sector;",
      3, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, pTbin);
  hn_photon->SetBinEdges(1, pTbin);
  hn_photon->Sumw2();
  hm->registerHisto(hn_photon);

  return;
}

void AnaFastMC::ReadTowerStatus(const string &filename)
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

void AnaFastMC::ReadSashaWarnmap(const string &filename)
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

void AnaFastMC::ReadSimWarnmap(const string &filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  unsigned int sector = 0;
  unsigned int biny = 0;
  unsigned int binz = 0;

  TOAD *toad_loader = new TOAD("AnaFastMC");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> binz >> biny )
  {
    // count tower with bad status for PbSc and PbGl
    if( sector < 6 ) nBadSc++;
    else nBadGl++;

    tower_status[sector][biny][binz] = 1;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

bool AnaFastMC::pi0_sim(TLorentzVector *pG1, TLorentzVector *pG2, float& ptsim, float& mm, float& dist)
  // returns false if one or two gammas out of acceptance
{
  const float mcorr = 0.993; // Hisa's effect of reducing pi0 mass due to conversion and background pi0's (from other decays)

  ptsim = 0.;
  mm = 0.;

  double e1 = pG1->E();
  double px1 = pG1->Px();
  double py1 = pG1->Py();
  double pz1 = pG1->Pz();
  double e2 = pG2->E();
  double px2 = pG2->Px();
  double py2 = pG2->Py();
  double pz2 = pG2->Pz();

  // Simulate EMCal postion resolution
  bool acc1 = Gamma_Pos(px1,py1,pz1);
  bool acc2 = Gamma_Pos(px2,py2,pz2);

  // Simulate EMCal energy resolution
  float e1t = 0.;
  float e2t = 0.;
  int itw01 = -1;
  int itw02 = -1;
  float ximp1,yimp1,zimp1;
  float ximp2,yimp2,zimp2;

  if( acc1 ) acc1 = Gamma_En(px1,py1,pz1,e1t,itw01,ximp1,yimp1,zimp1);
  if( acc1 ) {
    float k1 = e1t/e1;
    //    k1 = 1; // !!!
    pG1->SetPxPyPzE(px1*k1,py1*k1,pz1*k1,e1*k1);
  }

  if( acc2 ) acc2 = Gamma_En(px2,py2,pz2,e2t,itw02,ximp2,yimp2,zimp2);
  if( acc2 ) {
    float k2 = e2t/e2;
    //    k2 = 1; // !!!
    pG2->SetPxPyPzE(px2*k2,py2*k2,pz2*k2,e2*k2);
  }

  dist = sqrt((ximp1-ximp2)*(ximp1-ximp2)+
      (yimp1-yimp2)*(yimp1-yimp2)+
      (zimp1-zimp2)*(zimp1-zimp2));

  if( acc1 && acc2 ) { 
    TLorentzVector pPi0 = *pG1 + *pG2;
    mm = pPi0.M()*mcorr;
    ptsim = pPi0.Pt();
  }

  if( acc1 ) {
    Vpart[NPart] = *pG1;
    itw_part[NPart] = itw01;
    NPart++;
  }
  if( acc2 ) {
    Vpart[NPart] = *pG2;
    itw_part[NPart] = itw02;
    NPart++;
  }

  return (acc1 && acc2);
}

bool AnaFastMC::photon_sim(TLorentzVector *pG1, float& ptsim)
  // returns false if this photon out of acceptance
{
  ptsim = 0.;

  double e1 = pG1->E();
  double px1 = pG1->Px();
  double py1 = pG1->Py();
  double pz1 = pG1->Pz();

  // Simulate EMCal postion resolution
  bool acc1 = Gamma_Pos(px1,py1,pz1);

  // Simulate EMCal energy resolution
  float e1t = 0.;
  int itw01 = -1;
  float ximp1,yimp1,zimp1;

  if( acc1 ) acc1 = Gamma_En(px1,py1,pz1,e1t,itw01,ximp1,yimp1,zimp1);
  if( acc1 ) {
    float k1 = e1t/e1;
    //    k1 = 1; // !!!
    pG1->SetPxPyPzE(px1*k1,py1*k1,pz1*k1,e1*k1);
    ptsim = pG1->Pt();

    Vpart[NPart] = *pG1;
    itw_part[NPart] = itw01;
    NPart++;
  }

  return acc1;
}

void AnaFastMC::ResetTowerEnergy()
{
  NPart = 0;
  NPeak = 0;
  for(int i=0; i<MAXPEAK; i++)
  {
    Vpart[i].SetPxPyPzE(0.,0.,0.,0.);
    itw_part[i] = -1;
  }

  for( int is=0; is<NSEC; is++ )
    for( int iy=0; iy<NY; iy++ )
      for( int iz=0; iz<NZ; iz++ )
        eTwr[is][iy][iz] = 0.;
}

void AnaFastMC::FillTowerEnergy( int sec, int iy, int iz, float e )
{
  if( sec<0 || sec>NSEC-1 || iy<0 || iy>NY-1 || iz<0 || iz>NZ-1 ) return;
  eTwr[sec][iy][iz] += e;
  return;
}

float AnaFastMC::GetETwr( int sec, int iy, int iz )
{
  if( sec<0 || sec>NSEC-1 || iy<0 || iy>NY-1 || iz<0 || iz>NZ-1 ) return 0;
  return eTwr[sec][iy][iz];
}

int AnaFastMC::GetNpeak()
{
  const float eThresh = 0.02;

  int npeak=0;
  float e;
  for( int is=0; is<NSEC; is++ ) {
    for( int iy=0; iy<NY; iy++ ) {
      for( int iz=0; iz<NZ; iz++ ) {
        e = eTwr[is][iy][iz];
        if( e>eThresh ) {
          bool bpeak = true;
          // Check NO angles
          //    if( GetETwr(is,iy-1,iz)>e || GetETwr(is,iy+1,iz)>e || 
          //        GetETwr(is,iy,iz-1)>e || GetETwr(is,iy,iz+1)>e ) bpeak=false;
          // Check 3x3 matrix
          for( int dy=-1; dy<=1; dy++ ) {
            for( int dz=-1; dz<=1; dz++ ) {
              if( GetETwr( is, iy+dy, iz+dz ) > e ) bpeak=false;
            }
          }
          if( bpeak ) {
            npeak++;
          }
          //    if( e>0.01 ) printf("E=%f %d\n",e,bpeak);
        } // if( e> eTresh )
      }
    }
  }
  NPeak = npeak;
  return npeak;
}

bool AnaFastMC::GetImpactSectorTower(Double_t px, Double_t py, Double_t pz,  
    int& sec, int& iz, int& iy, float& zz, float& yy, 
    float& phi0, float& ximp, float& yimp, float& zimp )
//
// Returns: sector number (sec), tower coordinates in sector (iz,iy)
//          and hit position in tower (zz,yy) in tower units;
//          All that - for shower CG
//          and impact position on sector face
// 
// Edge cut is done here
//
{
  static float zvert = 0;

  // EMCal Geometry
  // X axis is towards W0

  const Double_t xsec_sc = 507;
  const Double_t zsize_sc = 5.58; // Tower size in Z
  const Double_t ysize_sc = 5.57; // Tower size in Y
  //  const Double_t zsize_sc = 6.5; // Tower size in Z
  //  const Double_t ysize_sc = 6.5; // Tower size in Y

  const Double_t xsec_gl = 540;
  const Double_t zsize_gl = 4.092; // Tower size in Z
  const Double_t ysize_gl = 4.106; // Tower size in Y

  // Sector phi position: 0<phi<2*PI; -phi=2*PI-phi
  const Double_t phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8, 
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  Double_t xsec, zsize, ysize;

  sec=-1; iz=-1; iy=-1; zz=0; yy=0;

  TVector3 vv(px,py,pz);
  Double_t phi = vv.Phi(); // -PI<phi<PI !!

  if( px > 0 ) { // West Arm
    for( int is=0; is<4; is++ ) {
      if( TMath::Abs(phi-phi_sec[is])<PI/16 ) sec = is;
    }
    if( sec < 0 ) return false;    
    // Make -PI/16 < phi < +PI/16
    vv.RotateZ(-phi_sec[sec]); phi = vv.Phi();
  }
  else { // East Arm
    // Rotate on -PI
    vv.RotateZ(-PI); phi = vv.Phi();
    for( int is=4; is<8; is++ ) {
      if( TMath::Abs(phi-(phi_sec[is]-PI))<PI/16 ) sec = is;
    }
    if( sec < 0 ) return false;    
    // Make -PI/16 < phi < +PI/16
    vv.RotateZ(-(phi_sec[sec]-PI)); phi = vv.Phi();
    // revert Phi and Theta, because numbering on East is reversed
    phi = -phi;
    Double_t theta = PI-vv.Theta();
    vv.SetPhi(phi);
    vv.SetTheta(theta);
  }
  phi0 = phi;

  float dl, x0;
  int nZ, nY;
  if( sec < 6 ) { // PbSc
    xsec = xsec_sc;
    zsize = zsize_sc;
    ysize = ysize_sc;
    dl = zsize * ( 1.93 + 0.383*log(vv.Mag()) );
    x0 = 2.;
    nZ = nZ_sc;
    nY = nY_sc;
  }
  else { // PbGl
    xsec = xsec_gl;
    zsize = zsize_gl;
    ysize = ysize_gl;
    dl = zsize * ( 1.93 + 0.383*log(vv.Mag()) ); // !!!!! Check it !!!!!
    x0 = 2.; // !!!!! Check it !!!!!
    nZ = nZ_gl;
    nY = nY_gl;
  }
  Double_t zsec_min = -zsize*nZ/2;
  Double_t ysec_min = -ysize*nY/2;

  // Center of gravity shift along Z
  Double_t sinZ = vv.CosTheta(); // I work with PI/2-Theta
  float dz = dl*sinZ;
  float dz_gamma = x0*sinZ; // Shift of gamma related to electron
  // Center of gravity shift along Y
  float dy = dl*phi; // aY ~ sin(aY)
  float dy_gamma = x0*phi; // Shift of gamma related to electron

  Double_t ysec = xsec*vv.Py()/vv.Px();
  Double_t zsec = xsec*vv.Pz()/vv.Px();
  zsec += zvert; // Zvert shift
  // Photons are deeper than electrons
  ysec += dy_gamma;
  zsec += dz_gamma;

  // For absolute (smeared by resolution!) impact position calculation
  TVector3 vv_imp(xsec,ysec,zsec);
  // Rotation on phi_sec[sec] should be made here !!!!!
  ximp = xsec;
  yimp = ysec;
  zimp = zsec;

  //  printf("E=%f: (%f,%f)  (%f,%f)\n",vv.Mag(),zsec,ysec,dz,dy);

  // Systematics in position measurements
  //  ysec *= 1.01;
  //  zsec *= 1.01;

  if( zsec <= zsec_min || ysec <= ysec_min ) {
    sec = -1;
    return false;
  }

  // Impact tower
  iz = int((zsec-zsec_min)/zsize);
  iy = int((ysec-ysec_min)/ysize);

  if( iz >= nZ || iy >= nY ) { 
    sec=-1; iz=-1; iy=-1;
    return false;
  }

  // Calculate max energy tower (center of gravity)
  Double_t ysec_cg = ysec + dy;
  Double_t zsec_cg = zsec + dz;
  int iz_cg, iy_cg;
  if( zsec_cg <= zsec_min ) {sec=-1; return false;}
  else {
    iz_cg = int((zsec_cg-zsec_min)/zsize);
    if( iz_cg >= nZ ) {sec=-1; return false;}
  }
  if( ysec_cg <= ysec_min ) {sec=-1; return false;}
  else {
    iy_cg = int((ysec_cg-ysec_min)/ysize);
    if( iy_cg >= nY ) {sec=-1; return false;}
  }

  //  if( CheckWarnMap(sec,iy,iz) ) return false;
  // Check CG tower (maximal energy)
  iz = iz_cg;
  iy = iy_cg;
  if( tower_status[sec][iy][iz] > 0 ) return false; 

  zz = (zsec_cg-zsec_min)/zsize-iz;
  yy = (ysec_cg-ysec_min)/ysize-iy;
  if( zz<0 || zz>1 || yy<0 || yy>1 ) printf("Error in GetImpactSectorTower: %f %f\n",zz,yy);

  return true;
}

bool AnaFastMC::GetShower(Double_t px, Double_t py, Double_t pz, float& eout, int& itw )
  // Return false if outside acceptance
  // Called by Gamma_En(...)
{
  //  const float bwidth = 0.15; // Slope parameter in shower shape
  //const float bwidth = 0.20; // Slope parameter in shower shape
  //  const float bwidth = 0.005; // for xcg=0
  //  const float bwidth = 1.; // for xcg~x
  const float corr = 0.933;
  //static TRandom* gen = new TRandom();

  const float thresh = 0.010;
  //  const float thresh = 0.030; //!!!!!

  eout = 0;
  itw = -1;
  if( TMath::Sqrt(TMath::Abs(px*px+py*py+pz*pz)) <= 0.01 ) return false;

  int sec, iz0, iy0;
  float zz, yy; // hit position in tower frame (in tower units)
  float phi;
  float ximp, yimp, zimp;
  // if outside acceptance - return
  if( !GetImpactSectorTower(px,py,pz, sec,iz0,iy0,zz,yy,phi,ximp,yimp,zimp) ) return false;

  if( sec<6 ) { 
    itw = sec*2592 + iy0*72 + iz0;
  }
  else {
    itw = 15552 + (sec-6)*4608 + iy0*96 + iz0;
  }

  float en = sqrt(px*px+py*py+pz*pz);

  // The rest is valid only for PbSc (tuning is needed for PbGl) !!!!!

  // Calculate angle and energy dependent shower parameters

  TVector3 vv(px,py,pz);
  Double_t sinZ = vv.CosTheta(); // I work with PI/2-Theta

  //  Double_t aY = TMath::Abs(vv.Phi());
  //  while( aY > PI/16 ) aY -= PI/16;
  //  Double_t sinY = aY; // aY ~ sin(aY)

  Double_t sinY = TMath::Sin(phi);
  float sin2a = sinZ*sinZ + sinY*sinY;
  if( sin2a > 0.5*0.5 ) printf("Something wrong in GetShower: too big angles %f %f\n",sinZ,sinY);
  float lgE = 0;
  if( en > 0.01 ) lgE = log(en);
  else sin2a = 0;

  float par1=0.59-(1.45+0.13*lgE)*sin2a;
  float par2=0.265+(0.80+0.32*lgE)*sin2a;
  float par3=0.25+(0.45-0.036*lgE)*sin2a;
  float par4=0.42;

  // !!!!! Modify profile
  //  par3 = 0.;
  //  par4 *= 2.;
  //  par2 = 0.6; // as K.Okada

  // Calculate center of gravity
  float dz0, dy0, dz, dy, r1, r2, r3;
  float tt = 0.98 + 0.98*sqrt(en);
  float bx = 0.20 + tt*sinZ*sinZ;
  float by = 0.20 + tt*sinY*sinY;
  dz0 = zz - 0.5;
  dy0 = yy - 0.5;
  dz0 = float(0.5*TMath::SinH(Double_t(dz0/bx))/TMath::SinH(0.5/bx));
  dy0 = float(0.5*TMath::SinH(Double_t(dy0/by))/TMath::SinH(0.5/by));

  // Calculate ecore

  float et;
  float esum = 0;
  float ecore = 0;
  float et5x5[5][5];
  for( int iz=-2; iz<=2; iz++ ) {
    for( int iy=-2; iy<=2; iy++ ) {
      //       dz = dz0 + iz;
      //       dy = dy0 + iy;
      dz = iz - dz0;
      dy = iy - dy0;
      r2=dz*dz+dy*dy;
      r1=sqrt(r2);
      r3=r2*r1;
      et = par1*exp(-r3/par2)+par3*exp(-r1/par4);
      //      et = GetEnergyFromShowerProfile(en,dz,dy); // If from loaded profile
      et5x5[iy+2][iz+2] = et;
      esum += et;
      //      if( !CheckTwrMap(sec,iy0+iy,iz0+iz) ) printf("Rejected: %d %d %d\n",sec,iy0,iz0);
      if( et > 0.02 && et*en > thresh ) 
        ecore += et;
    }
  }
  eout = ecore*en/esum; // normalization to real shower energy
  eout /= corr; // tune to data global calibration

  for( int iy=-2; iy<=2; iy++ ) {
    for( int iz=-2; iz<=2; iz++ ) {
      et = et5x5[iy+2][iz+2] * en/esum;
      FillTowerEnergy(sec,iy0+iy,iz0+iz,et);
    }
  }

  return true;
}

bool AnaFastMC::Gamma_En(Double_t px, Double_t py, Double_t pz, float& eout, int& itw,
    float& ximp, float& yimp, float& zimp)
// returns false if out of acceptance
{
  // PbSc
  static float a_sc = 0.078; // official
  //  const float b_sc = 0.030; // sigma/E = a/sqrt(E) + b, official
  //  static float b_sc = 0.055; // sigma/E = a/sqrt(E) + b
  static float b_sc = 0.040; // sigma/E = a/sqrt(E) + b, for eta
  static float bc_sc = 0.012; // this is artificial const. term introduced by ecore calculations in GetShower(...); should be subtracted when GetShower(...) used
  static float cnoise_sc = 0.015; // Noise term
  //  static float corr_sc = 1.000; // For Gamma
  static float corr_sc = 1.011; // for MB

  // PbGl
  static float a_gl = 0.085;
  //  static float b_gl = 0.065; // for pi0
  static float b_gl = 0.037; // for eta
  static float bc_gl = 0.;
  static float cnoise_gl = 0.030; // Noise term
  //  static float corr_gl = 1.007;
  static float corr_gl = 1.009;

  float a, b, bc, corr, cnoise;

  //const float Eback = 0.1; // mean energy per fired tower
  //const float Rload = 0.08; // 25% of towers fired
  //int Ntwr = 4; // number of towers in Core
  float et, et1; // et2;

  eout = 0;
  itw = -1;

  float ein = sqrt(px*px+py*py+pz*pz);
  if( ein <= 0 ) return false;

  int sec, iz0, iy0;
  float zz, yy; // hit position in tower frame (in tower units)
  float phi;
  GetImpactSectorTower(px,py,pz, sec,iz0,iy0,zz,yy,phi,ximp,yimp,zimp);
  if( sec<0 ) return false;

  if( sec<6 ) { // PbSc
    a = a_sc;
    b = b_sc;
    bc = bc_sc;
    cnoise = cnoise_sc;
    corr = corr_sc;
  }
  else { // PbGl
    a = a_gl;
    b = b_gl;
    bc = bc_gl;
    cnoise = cnoise_gl;
    corr = corr_gl;
  }

  //  float res = sqrt(a*a*ein + b*b*ein*ein + cnoise*cnoise);
  //  float res = sqrt(a*a*ein + b*b*ein*ein + 0.03*0.03);
  float res = sqrt(a*a*ein + b*b*ein*ein - bc*bc*ein*ein);
  //    res = sqrt(a*a*ein + b*b*ein*ein - bc*bc*ein*ein + 0.03*0.03);

  res = sqrt(res*res+0.01*0.01); // 10 MeV noise per cluster

  et = gRandom->Gaus(ein,res);
  if( et <= 0 ) et=0;
  //  et = ein; //!!!!!

  // Switching from High Gain to Low Gain
  if( et > 1 ) {
    et *= gRandom->Gaus(1,0.03);
  }

  float k = et/ein;
  if( !GetShower(px*k,py*k,pz*k,et1,itw) ) {
    return false;
  }

  //  if(sec>5) et1 *= (1-0.135644*exp(-1.10904*et1)); // Add'l attenuation for PbGl
  //  if(sec>5) et1 *= 1.; // Add'l attenuation for PbGl
  //  if(sec>5) et1 *= pow(et1/2,+2./500.); // Add'l attenuation for PbGl
  //  if(sec>5) et1 *= pow(et1/2,+2./150.); // Add'l attenuation for PbGl
  if(sec>5) et1 *= pow(et1/2,+2./120.); // Add'l attenuation for PbGl
  else et1 *= pow(et1,-2./800.); // Add'l attenuation for PbSc
  //  else et1 *= pow(et1/2,+2./9999999.); // Add'l attenuation for PbSc; normalization for 2 GeV

  // Non-linearity near pedestal/threshold (data is undercorrected for it)
  if( sec<=5 ) et1 *= (1-0.003/et1/et1); 
  else         et1 *= (1-0.005/et1/et1);

  // Shower overlap effect
  /*
     et2 = 0;
     Ntwr=4;
     if( ein < 0.5 ) Ntwr=3;
     if( ein < 0.2 ) Ntwr=2;
     if( ein < 0.1 ) Ntwr=1;
  //  Ntwr=4;
  for( int i=0; i<Ntwr; i++ ) {
  if( gRandom->Rndm() < Rload ) 
  et2 += gRandom->Exp(Eback);
  }
  //  eout = et1+et2;
  */
  eout = et1*corr;
  //  eout = et; 
  //  eout = ein;
  return true;
}

bool AnaFastMC::Gamma_Pos(Double_t& px, Double_t& py, Double_t& pz )
{
  // PbSc
  const float a_sc = 0.16;
  //  const float a_sc = 0.24; // Default
  //  const float a_sc = 0.14;
  const float b_sc = 0.90;
  //  const float b_sc = 0.59; // Default
  const float c_sc = 2.00; // Pos. res. = a+b/sqrt(E) (+) (c*sinT)
  const float xsec_sc = 507;

  // PbGl
  //  const float a_gl = 0.161;
  //  const float a_gl = 0.24; // Default
  const float a_gl = 0.14;
  //  const float b_gl = 0.673; // Default
  const float b_gl = 1.15;
  // "c_gl" should be checked !!!!!
  const float c_gl = 2.; // Pos. res. = a (+) b/sqrt(E) (+) (c*sinT)
  const float xsec_gl = 540;

  float a, b, c, xsec;

  int sec, iz0, iy0;
  float zz, yy; // hit position in tower frame (in tower units)
  float phi;
  float ximp, yimp, zimp;
  GetImpactSectorTower(px,py,pz, sec,iz0,iy0,zz,yy,phi,ximp,yimp,zimp);
  if( sec<0 ) return false;

  if( sec<6 ) { // PbSc
    a = a_sc;
    b = b_sc;
    c = c_sc;
    xsec = xsec_sc;
  }
  else { // PbGl
    a = a_gl;
    b = b_gl;
    c = c_gl;
    xsec = xsec_gl;
  }

  Double_t res = 0., dl;
  TVector3 vv(px,py,pz);
  Double_t En = vv.Mag();

  if( sec<6 ) res = a + b/TMath::Sqrt(En);
  else        res = sqrt(a*a + b*b/En);

  // Long. fluctuations
  dl = gRandom->Gaus(0,c);

  // Resolution along Z
  Double_t sinZ = vv.CosTheta(); // I work with PI/2-Theta
  //  Double_t resz = 0.65*res * ( 1.+ 2.*TMath::Gaus(zz-0.5,0.,0.18) ); // To account for worse resolution in tower center
  Double_t resz = res;
  Double_t dz = gRandom->Gaus(0,resz);
  dz += dl*sinZ;

  // Resolution along Y
  /*
     Double_t aY = TMath::Abs(vv.Phi());
     int nrot=0;
     while( aY > PI/16 ) { aY -= PI/16; nrot++; }
     Double_t dy = gRandom->Gaus(0,res);
     dy += dl*aY; // aY ~ sin(aY)

     TVector3 dv(0,dy,dz);
     vv.RotateZ(-nrot*PI/16);

     vv.SetMag(TMath::Abs(En/vv.Px()*xsec)); // Now in cm
     TVector3 VV = vv + dv;
     VV.SetMag(En); // Conserve energy keeping Theta and Phi
     VV.RotateZ(nrot*PI/16);
     */

  //  Double_t resy = 0.65*res * ( 1.+ 2.*TMath::Gaus(yy-0.5,0.,0.18) ); // To account for worse resolution in tower center
  Double_t resy = res;
  Double_t dy = gRandom->Gaus(0,resy);
  dy += dl*phi; // aY ~ sin(aY)

  // Rotate according to dy in sector frame
  vv.RotateZ(dy/xsec); // Approximate: correct within 3%

  // Shift according to dz
  TVector3 dv(0,0,dz);
  if(vv.Px()*vv.Px()+vv.Py()*vv.Py() > 0.)
    vv.SetMag(En/TMath::Sqrt(vv.Px()*vv.Px()+vv.Py()*vv.Py())*xsec); // Now in cm
  TVector3 VV = vv + dv;
  VV.SetMag(En); // Conserve energy keeping Theta and Phi

  px = VV.Px();
  py = VV.Py();
  pz = VV.Pz();

  return true;
}
