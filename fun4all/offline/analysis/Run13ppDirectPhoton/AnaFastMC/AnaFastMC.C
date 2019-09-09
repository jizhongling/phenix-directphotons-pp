#include "AnaFastMC.h"

#include "PtWeights.h"
#include <AnaToolsTowerID.h>
#include <EMCWarnmapChecker.h>
#include <DCDeadmapChecker.h>

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

#include <PHPythiaContainer.h>
#include <PHPyCommon.h>

#include <TPythia6.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8) 
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TGenPhaseSpace.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

/* Some constants */
const double PI = TMath::Pi();
const TVector2 v2_2PI(0., 2.*PI);
const double mPi0 = 0.1349770;
const double mEta = 0.547862;

/* Hisa's effect of reducing pi0 mass due to conversion and background pi0's (from other decays) */
const double mcorr = 0.993;

/* Some cuts for photon identification */
const double eMin = 0.3;
const double AsymCut = 0.8;

/* Some cuts for isolation cut */
const double eClusMin = 0.15;
const double pTrkMin = 0.2;
const double pTrkMax = 15.;

/* Isolation cut cone angle and energy fraction */
double cone_angle = 0.5;
double eratio = 0.1;

AnaFastMC::AnaFastMC(const string &name):
  SubsysReco(name),
  outFileName("histos/AnaFastMC-"),
  mcmethod(FastMC),
  phpythia(nullptr),
  weight_pythia(1.),
  ptweights(nullptr),
  pdg_db(nullptr),
  emcwarnmap(nullptr),
  dcdeadmap(nullptr),
  hm(nullptr),
  h_events(nullptr),
  h_pion(nullptr),
  h_photon(nullptr),
  h_photon_eta050(nullptr),
  h_photon_eta025(nullptr),
  h_isophoton_eta050(nullptr),
  h_isophoton_eta025(nullptr),
  h3_isopi0(nullptr),
  h3_isoeta(nullptr),
  hn_pion(nullptr),
  hn_missing(nullptr),
  hn_missing_eta(nullptr),
  hn_hadron(nullptr),
  hn_photon(nullptr),
  hn_geom(nullptr),
  hn_isolated(nullptr)
{
  /* Initialize histograms */
  for(int ih=0; ih<nh_eta_phi; ih++)
  {
    h2_pion_eta_phi[ih] = nullptr;
    h2_photon_eta_phi[ih] = nullptr;
  }

  /* Initialize array for tower status */
  for(int sec=0; sec<NSEC; sec++)
    for(int iy=0; iy<NY; iy++)
      for(int iz=0; iz<NZ; iz++)
      {
        tower_status_sim[sec][iy][iz] = 0;
        eTwr[sec][iy][iz] = 0.;
      }

  NPart = 0;
  NPeak = 0;
  for(int i=0; i<MAXPEAK; i++)
  {
    Vpart[i].SetPxPyPzE(0.,0.,0.,0.);
    itw_part[i] = -1;
    sec_part[i] = -1;
  }
}

AnaFastMC::~AnaFastMC()
{
}

int AnaFastMC::Init(PHCompositeNode *topNode)
{
  /* Create and register histograms */
  BookHistograms();

  /* pT weights calculator */
  PtWeights *ptweights = new PtWeights();
  if(!ptweights)
  {
    cerr << "No ptweights" << endl;
    exit(1);
  }

  /* Get PDGDatabase object */
  pdg_db = new TDatabasePDG();
  if(!pdg_db)
  {
    cerr << "No pdg_db" << endl;
    exit(1);
  }

  /* Initialize EMC warnmap checker */
  emcwarnmap = new EMCWarnmapChecker();
  if(!emcwarnmap)
  {
    cerr << "No emcwarnmap" << endl;
    exit(1);
  }

  /* Read sim warnmap */
  ReadSimWarnmap();

  /* Initialize DC deadmap checker */
  dcdeadmap = new DCDeadmapChecker();
  if(!dcdeadmap)
  {
    cerr << "No dcdeadmap" << endl;
    exit(1);
  }

  return EVENT_OK;
}

int AnaFastMC::process_event(PHCompositeNode *topNode)
{
  /* Count events */
  h_events->Fill(1.);

  /* Use FastMC input (i.e. random number generator for pi0, eta and direct photon eta, phi, pt) */
  if( mcmethod == FastMC )
    FastMCInput();

  /* Use PHParticleGen input */
  else if( mcmethod == PHParticleGen )
    PythiaInput(topNode);

  return EVENT_OK;
}

void AnaFastMC::FastMCInput()
{
  /* Initialize parameters for pi0 and eta decay and direct photon */
  TLorentzVector beam_pi0;
  TLorentzVector beam_eta;
  TLorentzVector beam_ph;
  double weight_pi0 = 1.;
  double weight_eta = 1.;
  double weight_ph = 1.;

  const double eta_max = 0.5;
  double phi = 2.0*PI*gRandom->Rndm();
  double eta = eta_max*2*gRandom->Rndm() - eta_max;

  const double emin = 0.; // min pi0 energy
  const double emax = 40.; // max pi0 energy
  double pt = gRandom->Rndm() * (emax-emin) + emin;
  double px = pt*cos(phi);
  double py = pt*sin(phi);

  /* This is if we use uniform rapidity (eta here is y) */
  double pz_pi0 = sqrt((pt*pt+mPi0*mPi0)*(exp(2*eta)-1)*(exp(2*eta)-1)/4./(exp(2*eta)));  // pz_pi0 = mt*gamma*beta = mt*sinh(eta)
  double pz_eta = sqrt((pt*pt+mEta*mEta)*(exp(2*eta)-1)*(exp(2*eta)-1)/4./(exp(2*eta)));  // pz_eta = mt*gamma*beta = mt*sinh(eta)
  double pz_ph = sqrt(pt*pt*(exp(2*eta)-1)*(exp(2*eta)-1)/4./(exp(2*eta)));  // pz_ph = pt*gamma*beta = pt*sinh(eta)
  if( eta < 0 )
  {
    pz_pi0 = -pz_pi0;
    pz_eta = -pz_eta;
    pz_ph = -pz_ph;
  }
  double energy_pi0 = sqrt(pz_pi0*pz_pi0+pt*pt+mPi0*mPi0);
  double energy_eta = sqrt(pz_eta*pz_eta+pt*pt+mEta*mEta);
  double energy_ph = sqrt(pz_ph*pz_ph+pt*pt);

  beam_pi0.SetPxPyPzE(px,py,pz_pi0,energy_pi0);
  beam_eta.SetPxPyPzE(px,py,pz_eta,energy_eta);
  beam_ph.SetPxPyPzE(px,py,pz_ph,energy_ph);

  /* Code for pi0 acceptance and merging rate */

  /* Fill Vpart[], itwr_part[] and sec_part[]
   * and set NPart and NPeak */
  pi0_sim(beam_pi0);

  /* Get event weight for pi0 */
  if(pt > 1.)
    weight_pi0 = ptweights->EvalPi0(pt);
  else
    weight_pi0 = ptweights->EvalPi0(1.);

  /* Fill histogram for all generated pi0 */
  h_pion->Fill(pt, weight_pi0);

  /* Two photons are in acceptance */
  if(NPart == 2)
  {
    /* Parameters for two clusters from pi0 decay photons */
    int sec1 = sec_part[0];
    int sec2 = sec_part[1];
    double e1 = Vpart[0].E();
    double e2 = Vpart[1].E();

    /* Check warnmap, in the same detector part
     * and pass the energy and asymmetry cuts */
    if( !CheckWarnMap(itw_part[0]) &&
        !CheckWarnMap(itw_part[1]) &&
        anatools::SectorCheck(sec1,sec2) &&
        e1 > eMin && e2 > eMin &&
        fabs(e1-e2)/(e1+e2) < AsymCut )
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
      } // iph

      /* Fill pi0 histogram
       * NPeak=2: No merging for two photons
       * NPeak=1: Merging for two photons */
      TLorentzVector pPi0 = Vpart[0] + Vpart[1];
      double ptsim = pPi0.Pt();
      double minv = pPi0.M()*mcorr;
      double fill_hn_pion[] = {pt, ptsim, minv, (double)sec1, (double)NPeak};
      hn_pion->Fill(fill_hn_pion, weight_pi0);
    } // !CheckWarnMap
  } // NPart

  /* Code for pi0 missing ratio and self veto power */

  for(int iph=0; iph<2; iph++)
    if( emcwarnmap->InFiducial(itw_part[iph]) &&
        Vpart[iph].E() > eMin )
    {
      int npart = 1;
      int npeak = 1;
      if( emcwarnmap->IsGoodTower(itw_part[1-iph]) &&
          Vpart[1-iph].E() > eMin )
      {
        npart = NPart;
        npeak = NPeak;
      }

      double ptsim = pt;
      if(NPart == 2)
        ptsim = (Vpart[0] + Vpart[1]).Pt();

      double fill_hn_missing[] = {ptsim, Vpart[iph].Pt(), (double)sec_part[iph], (double)npart, (double)npeak};
      hn_missing->Fill(fill_hn_missing, weight_pi0);

      if( npart == 1 )
      {
        int isolated = 1;
        h3_isopi0->Fill(Vpart[iph].Pt(), (double)sec_part[iph], (double)isolated);
      }
      else if( npart == 2 && npeak == 2 )
      {
        int isolated = 0;
        double angle = Vpart[0].Angle( Vpart[1].Vect() );
        double econe = angle < cone_angle ? Vpart[1].E() : 0.;
        if( econe < eratio * Vpart[0].E() )
          isolated = 1;
        h3_isopi0->Fill(Vpart[iph].Pt(), (double)sec_part[iph], (double)isolated);
      } // npart and npeak
    } // iph

  /* Code for eta self veto power */

  /* Fill Vpart[], itwr_part[] and sec_part[]
   * and set NPart and NPeak */
  pi0_sim(beam_eta);

  /* Get event weight for eta */
  double pteff = sqrt(pt*pt + mEta*mEta - mPi0*mPi0);
  if(pt > 1.)
    weight_eta = ptweights->EvalPi0(pteff);
  else
    weight_eta = ptweights->EvalPi0(1.);

  for(int iph=0; iph<2; iph++)
    if( emcwarnmap->InFiducial(itw_part[iph]) &&
        Vpart[iph].E() > eMin )
    {
      int npart = 1;
      int npeak = 1;
      if( emcwarnmap->IsGoodTower(itw_part[1-iph]) &&
          Vpart[1-iph].E() > eMin )
      {
        npart = NPart;
        npeak = NPeak;
      }

      double ptsim = pt;
      if(NPart == 2)
        ptsim = (Vpart[0] + Vpart[1]).Pt();

      double fill_hn_missing_eta[] = {ptsim, Vpart[iph].Pt(), (double)sec_part[iph], (double)npart, (double)npeak};
      hn_missing_eta->Fill(fill_hn_missing_eta, weight_eta);

      if( npart == 1 )
      {
        int isolated = 1;
        h3_isoeta->Fill(Vpart[iph].Pt(), (double)sec_part[iph], (double)isolated);
      }
      else if( npart == 2 && npeak == 2 )
      {
        int isolated = 0;
        double angle = Vpart[0].Angle( Vpart[1].Vect() );
        double econe = angle < cone_angle ? Vpart[1].E() : 0.;
        if( econe < eratio * Vpart[0].E() )
          isolated = 1;
        h3_isoeta->Fill(Vpart[iph].Pt(), (double)sec_part[iph], (double)isolated);
      } // npart and npeak
    } // iph

  /* Code for direct photon acceptance */

  /* Fill Vpart[], itwr_part[] and sec_part[]
   * and set NPart and NPeak */
  photon_sim(beam_ph);

  /* Get event weight for direct photon */
  if(pt > 1.)
    weight_ph = ptweights->EvalPhoton(pt);
  else
    weight_ph = ptweights->EvalPhoton(1.);

  /* Fill histogram for all generated direct photons */
  h_photon->Fill(pt, weight_ph);

  /* Direct photon is in fiducial */
  if( emcwarnmap->InFiducial(itw_part[0]) )
  {
    /* Parameters for direct photon */
    int sec1 = sec_part[0];
    double e1 = Vpart[0].E();
    double ptsim = Vpart[0].Pt();

    /* For eta and phi distribution and acceptance */
    if(e1 > eMin)
    {
      int part = -1;
      if(sec1 < 0) part = -1;
      else if(sec1 < 4) part = 0;
      else if(sec1 < 6) part = 1;
      else if(sec1 < 8) part = 2;

      double eta = Vpart[0].Eta();
      double phi = Vpart[0].Phi();
      if(sec1 >= 4)
      {
        Vpart[0].RotateZ(-PI);
        phi = Vpart[0].Phi() + PI;
      }

      if( part >= 0 && ptsim > 5. && ptsim < 10. )
        h2_photon_eta_phi[part]->Fill(-eta, phi, weight_ph);

      double fill_hn_photon[] = {pt, ptsim, (double)sec1};
      hn_photon->Fill(fill_hn_photon, weight_ph);
    } // e1
  } // InFiducial

  return;
}

void AnaFastMC::PythiaInput(PHCompositeNode *topNode)
{
  /* Get PYTHIA */
  phpythia = findNode::getClass<PHPythiaContainer>(topNode,"PHPythia");
  if(!phpythia)
  {
    cout << PHWHERE << "Unable to get PHPythia, is Node missing?" << endl;
    return;
  }

  /* Set DC deadmap */
  dcdeadmap->SetMapByEvent();

  /* Loop over all signal particles */
  int npart = phpythia->size();
  for(int ipart=0; ipart<npart; ipart++)
  {
    TMCParticle *particle = phpythia->getParticle(ipart);
    TMCParticle *parent = phpythia->getParent(particle);

    /* Get particle code */
    int id = particle->GetKF();

    /* Put particle's momentum and energy into TLorentzVector */
    TLorentzVector pE_part(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
    double pt = pE_part.Pt();
    double eta = pE_part.Eta();

    /* Only consider high-pT particles in central acceptance */
    if( pt < 2. || fabs(eta) > 1. )
      continue;

    int hadronid = -1;
    if( id == PY_PIZERO )
      hadronid = 0;
    else if( id == PY_ETA )
      hadronid = 1;
    else if( id == 223 )  // omega
      hadronid = 2;
    else if( id == 331 )  // eta prime
      hadronid = 3;

    if( hadronid >= 0 )
    {
      double fill_hn_hadron[] = {pt, (double)hadronid};
      hn_hadron->Fill(fill_hn_hadron, weight_pythia);
    }

    /* Test if particle is a stable prompt photon */
    if( id != PY_GAMMA ||
        particle->GetKS() != 1 ||
        ( parent && abs(parent->GetKF()) > 100 ) )
      continue;

    /* Fill Vpart[], itwr_part[] and sec_part[]
     * and set NPart and NPeak */
    photon_sim(pE_part);

    /* Copy parameters for reconstructed photon */
    bool InAcc = emcwarnmap->InFiducial(itw_part[0]);
    int sec1 = sec_part[0];
    TLorentzVector pE_reco(Vpart[0]);

    /* Get isolation cone energy */
    double econe_all, econe_emc, econe_trk[3];
    SumETruth(particle, InAcc, econe_all, econe_emc, econe_trk);

    /* Fill histogram for all prompt photons with |eta| < 0.5 */
    if( fabs(eta) < 0.5 )
    {
      h_photon_eta050->Fill(pt, weight_pythia);

      /* Fill histogram for all prompt photons with |eta| < 0.25 */
      if( fabs(eta) < 0.25 )
        h_photon_eta025->Fill(pt, weight_pythia);

      /* Fill histogram for all isolated prompt photons with |eta| < 0.5 */
      if( econe_all < eratio * pE_part.E() )
      {
        h_isophoton_eta050->Fill(pt, weight_pythia);

        /* Fill histogram for all isolated prompt photons with |eta| < 0.25 */
        if( fabs(eta) < 0.25 )
          h_isophoton_eta025->Fill(pt, weight_pythia);
      }
    }

    /* Fill histogram for prompt photons in accpetance */
    if( InAcc )
    {
      /* Parameters for prompt photon */
      double ptsim = pE_reco.Pt();

      /* All prompt photons in acceptance */
      double fill_hn_geom[] = {pt, ptsim, (double)sec1};
      hn_geom->Fill(fill_hn_geom, weight_pythia);

      /* Eta and phi distribution and acceptance for isolated prompt photons */
      int part = -1;
      if(sec1 < 0) part = -1;
      else if(sec1 < 4) part = 0;
      else if(sec1 < 6) part = 1;
      else if(sec1 < 8) part = 2;

      double eta = pE_reco.Eta();
      double phi = pE_reco.Phi();
      if(sec1 >= 4)
      {
        pE_reco.RotateZ(-PI);
        phi = pE_reco.Phi() + PI;
      }

      for(int ival=0; ival<3; ival++)
        if( econe_emc + econe_trk[ival] < eratio * pE_part.P() &&
            pE_reco.P() > eMin )
        {
          if( part >= 0 && ptsim > 5. && ptsim < 10. )
          {
            int ih = part + 3*ival;
            h2_photon_eta_phi[ih]->Fill(-eta, phi, weight_pythia);
          }

          double fill_hn_isolated[] = {pt, ptsim, (double)sec1, (double)ival};
          hn_isolated->Fill(fill_hn_isolated, weight_pythia);
        } // Iso & e1
    } // InFiducial
  } // ipart

  return;
}

int AnaFastMC::End(PHCompositeNode *topNode)
{
  /* Write histogram output to ROOT file */
  hm->dumpHistos(outFileName);
  delete hm;
  delete ptweights;
  delete pdg_db;
  delete emcwarnmap;
  delete dcdeadmap;

  return EVENT_OK;
}

void AnaFastMC::SumETruth(const TMCParticle *pref, bool prefInAcc,
    double &econe_all, double &econe_emc, double econe_trk[])
{
  /* Sum up all energy in cone around particle */
  econe_all = 0.;
  econe_emc = 0.;
  for(int ival=0; ival<3; ival++)
    econe_trk[ival] = 0.;

  /* Get reference vector */
  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());
  TVector2 v2_pref = v3_pref.EtaPhiVector();

  /* Loop over all particles */
  int npart = phpythia->size();
  for(int ipart2=0; ipart2<npart; ipart2++)
  {
    TMCParticle *part2 = phpythia->getParticle(ipart2);

    /* Only consider stable particles and
     * skip if pointer identical to 'reference' particle */
    if( part2 == pref || part2->GetKS() != 1 )
      continue;

    /* Get particle vector */
    TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());
    TVector2 v2_part2 = v3_part2.EtaPhiVector();

    /* Check if particle within cone */
    TVector2 v2_diff(v2_part2 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < cone_angle )
    {
      econe_all += part2->GetEnergy();

      /* Reference particle should in acceptance and
       * only consider interacting particles in cone */
      int id = abs(part2->GetKF());
      if( prefInAcc && id != PY_MU && id != PY_NU_E && id != PY_NU_MU )
      {
        TParticlePDG *pdg_part2 = pdg_db->GetParticle( part2->GetKF() );
        int charge = pdg_part2->Charge();
        double mom = v3_part2.Mag();

        /* For neutral particles sum kinetic energy in EMCal */
        if( charge == 0 && mom > eClusMin )
        {
          /* Fill itwr_part[] and sec_part[]
           * and set NPart */
          geom_sim(v3_part2);

          /* Test if particle is in acceptance and not on hot tower */
          if( NPart == 1 && !emcwarnmap->IsBadTower(itw_part[0]) )
            econe_emc += GetEMCResponse(id, mom);
        }

        /* For charged particles sum kinetic energy in EMCal
         * with hadron response */
        if( charge != 0 && mom > eClusMin )
          econe_trk[0] += GetEMCResponse(id, mom);

        /* For charged particles sum momentum in DC
         * if particle is in DC acceptance w/o deadmap */
        if( charge != 0 && mom > pTrkMin && mom < pTrkMax )
        {
          dcdeadmap->Checkmap(false);
          if( InDCAcceptance(v3_part2, charge) )
          {
            econe_trk[1] += mom;
            dcdeadmap->Checkmap(true);
            if( InDCAcceptance(v3_part2, charge) )
              econe_trk[2] += mom;
          }
        }
      } // prefInAcc and id
    } // cone_angle
  } // ipart2

  return;
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

  /* Eta and phi bins step size */
  const double step[2] = {0.011, 0.008};

  /* Eta bins */
  const int neta = 100;
  double etabin[2][neta+1];
  for(int part=0; part<2; part++)
    for(int it=0; it<=neta; it++)
      etabin[part][it] = step[part] * ( it - neta/2 );

  /* Phi sector */
  const double phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8,
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  /* Phi bins */
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

  /* Use FastMC and PHParticleGen input */
  h_events = new TH1F("h_events", "Events count", 1, 0.5, 1.5);
  hm->registerHisto(h_events);

  h_photon = new TH1F("h_photon", "Total photon count;p_{T} [GeV];", npT, pTbin);
  h_photon->Sumw2();
  hm->registerHisto(h_photon);

  // ih = part + 3*ival < 3*3
  for(int ih=0; ih<nh_eta_phi; ih++)
  {
    h2_photon_eta_phi[ih] = new TH2F(Form("h2_photon_eta_phi_part%d",ih), "Photon #eta and #phi distribution;#eta;#phi;", neta,etabin[ih%3/2], nphi,phibin);
    h2_photon_eta_phi[ih]->Sumw2();
    hm->registerHisto(h2_photon_eta_phi[ih]);
  }

  const int nbins_hn_photon[] = {npT, npT, 8};
  const double xmin_hn_photon[] = {0., 0., -0.5};
  const double xmax_hn_photon[] = {0., 0., 7.5};
  hn_photon = new THnSparseF("hn_photon", "EMCal photon count;p_{T} truth [GeV];p_{T} reco [GeV];sector;",
      3, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, pTbin);
  hn_photon->SetBinEdges(1, pTbin);
  hn_photon->Sumw2();
  hm->registerHisto(hn_photon);

  /* Use FastMC input */
  if( mcmethod == FastMC )
  {
    h_pion = new TH1F("h_pion", "Total pion count;p_{T} [GeV];", npT, pTbin);
    h_pion->Sumw2();
    hm->registerHisto(h_pion);

    for(int part=0; part<3; part++)
    {
      h2_pion_eta_phi[part] = new TH2F(Form("h2_pion_eta_phi_part%d",part), "Pion #eta and #phi distribution;#eta;#phi;", neta,etabin[part/2], nphi,phibin);
      h2_pion_eta_phi[part]->Sumw2();
      hm->registerHisto(h2_pion_eta_phi[part]);
    }

    h3_isopi0 = new TH3F("h3_isopi0", "Self veto for #pi^{0};p_{T} [GeV];sector;isolated;",
        npT,0.,0., 8,-0.5,7.5, 2,-0.5,1.5);
    h3_isopi0->GetXaxis()->Set(npT, pTbin);
    hm->registerHisto(h3_isopi0);

    h3_isoeta = static_cast<TH3*>( h3_isopi0->Clone("h3_isoeta") );
    h3_isoeta->SetTitle("Self veto for #eta");
    hm->registerHisto(h3_isoeta);

    const int nbins_hn_pion[] = {npT, npT, 300, 8, 3};
    const double xmin_hn_pion[] = {0., 0., 0., -0.5, -0.5};
    const double xmax_hn_pion[] = {0., 0., 0.3, 7.5, 2.5};
    hn_pion = new THnSparseF("hn_pion", "EMCal pion count;p_{T} truth [GeV];p_{T} reco [GeV];m_{inv} [GeV];sector;NPeak;",
        5, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
    hn_pion->SetBinEdges(0, pTbin);
    hn_pion->SetBinEdges(1, pTbin);
    hn_pion->Sumw2();
    hm->registerHisto(hn_pion);

    const int nbins_hn_missing[] = {npT, npT, 8, 3, 3};
    const double xmin_hn_missing[] = {0., 0., -0.5, -0.5, -0.5};
    const double xmax_hn_missing[] = {0., 0., 7.5, 2.5, 2.5};
    hn_missing = new THnSparseF("hn_missing", "#pi^{0} missing ratio;p^{#pi^{0}}_{T} [GeV];p^{#gamma}_{T} [GeV];sector;NPart;NPeak;",
        5, nbins_hn_missing, xmin_hn_missing, xmax_hn_missing);
    hn_missing->SetBinEdges(0, pTbin);
    hn_missing->SetBinEdges(1, pTbin);
    hn_missing->Sumw2();
    hm->registerHisto(hn_missing);

    hn_missing_eta = (THnSparse*)hn_missing->Clone("hn_missing_eta");
    hn_missing_eta->SetTitle("#eta missing ratio;p^{#eta}_{T} [GeV];p^{#gamma}_{T};sector;NPart;NPeak;");
    hn_missing_eta->Sumw2();
    hm->registerHisto(hn_missing_eta);
  }

  /* Use PHParticleGen input */
  if( mcmethod == PHParticleGen )
  {
    h_photon_eta050 = (TH1*)h_photon->Clone("h_photon_eta050");
    h_photon_eta050->Sumw2();
    hm->registerHisto(h_photon_eta050);

    h_photon_eta025 = (TH1*)h_photon->Clone("h_photon_eta025");
    h_photon_eta025->Sumw2();
    hm->registerHisto(h_photon_eta025);

    h_isophoton_eta050 = (TH1*)h_photon->Clone("h_isophoton_eta050");
    h_isophoton_eta050->Sumw2();
    hm->registerHisto(h_isophoton_eta050);

    h_isophoton_eta025 = (TH1*)h_photon->Clone("h_isophoton_eta025");
    h_isophoton_eta025->Sumw2();
    hm->registerHisto(h_isophoton_eta025);

    hn_geom = (THnSparse*)hn_photon->Clone("hn_geom");
    hn_geom->Sumw2();
    hm->registerHisto(hn_geom);

    const int nbins_hn_hadron[] = {npT, 4};
    const double xmin_hn_hadron[] = {0., -0.5};
    const double xmax_hn_hadron[] = {0., 3.5};
    hn_hadron = new THnSparseF("hn_hadron", "Hadron count;p_{T} [GeV];Hadron ID;",
        2, nbins_hn_hadron, xmin_hn_hadron, xmax_hn_hadron);
    hn_hadron->SetBinEdges(0, pTbin);
    hn_hadron->Sumw2();
    hm->registerHisto(hn_hadron);

    const int nbins_hn_isolated[] = {npT, npT, 8, 3};
    const double xmin_hn_isolated[] = {0., 0., -0.5, -0.5};
    const double xmax_hn_isolated[] = {0., 0., 7.5, 2.5};
    hn_isolated = new THnSparseF("hn_isolated", "EMCal photon count;p_{T} truth [GeV];p_{T} reco [GeV];sector;ival;",
        4, nbins_hn_isolated, xmin_hn_isolated, xmax_hn_isolated);
    hn_isolated->SetBinEdges(0, pTbin);
    hn_isolated->SetBinEdges(1, pTbin);
    hn_isolated->Sumw2();
    hm->registerHisto(hn_isolated);
  }

  return;
}

void AnaFastMC::ReadSimWarnmap()
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int sector = 0;
  int biny = 0;
  int binz = 0;

  TOAD *toad_loader = new TOAD("AnaFastMC");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("dead_eff_run13pp500gev.dat");
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> binz >> biny )
  {
    /* Count tower with bad status for PbSc and PbGl */
    if( sector < 6 ) nBadSc++;
    else nBadGl++;

    tower_status_sim[sector][biny][binz] = 1;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void AnaFastMC::pi0_sim(const TLorentzVector &beam_pi0)
{
  /* Initialize parameters */
  ResetClusters();
  ResetTowerEnergy();

  /* Decay to two photons */
  static TGenPhaseSpace event;
  const double masses[2] = {0., 0.};

  /* Let pi0 decay into two photons */
  TLorentzVector beam_copy(beam_pi0);
  event.SetDecay(beam_copy, 2, masses);
  event.Generate();

  /* Get pi0 decay photons */
  const TLorentzVector *pG1 = event.GetDecay(0);
  const TLorentzVector *pG2 = event.GetDecay(1);

  double e1 = pG1->E();
  double px1 = pG1->Px();
  double py1 = pG1->Py();
  double pz1 = pG1->Pz();
  double e2 = pG2->E();
  double px2 = pG2->Px();
  double py2 = pG2->Py();
  double pz2 = pG2->Pz();

  /* Simulate EMCal postion resolution */
  bool acc1 = Gamma_Pos(px1,py1,pz1);
  bool acc2 = Gamma_Pos(px2,py2,pz2);

  /* Simulate EMCal energy resolution */
  double e1t = 0.;
  double e2t = 0.;
  int itw01 = -1;
  int itw02 = -1;
  double ximp1,yimp1,zimp1;
  double ximp2,yimp2,zimp2;

  if( acc1 ) acc1 = Gamma_En(px1,py1,pz1,e1t,itw01,ximp1,yimp1,zimp1);
  if( acc1 ) {
    double k1 = e1t/e1;
    //    k1 = 1; // !!!
    Vpart[NPart].SetPxPyPzE(px1*k1,py1*k1,pz1*k1,e1*k1);
    itw_part[NPart] = itw01;
    NPart++;
  }

  if( acc2 ) acc2 = Gamma_En(px2,py2,pz2,e2t,itw02,ximp2,yimp2,zimp2);
  if( acc2 ) {
    double k2 = e2t/e2;
    //    k2 = 1; // !!!
    Vpart[NPart].SetPxPyPzE(px2*k2,py2*k2,pz2*k2,e2*k2);
    itw_part[NPart] = itw02;
    NPart++;
  }

  /* Loop over all clusters in calorimeter */
  for(int i=0; i<MAXPEAK; i++)
    if(itw_part[i] >= 0)
    {
      /* Get sector */
      int iy, iz;
      anatools::TowerLocation(itw_part[i], sec_part[i], iy, iz);
    }

  /* Get NPeak */
  GetNpeak();

  return;
}

void AnaFastMC::photon_sim(const TLorentzVector &beam_ph)
{
  /* Initialize parameters */
  ResetClusters();
  ResetTowerEnergy();

  double e1 = beam_ph.E();
  double px1 = beam_ph.Px();
  double py1 = beam_ph.Py();
  double pz1 = beam_ph.Pz();

  /* Simulate EMCal postion resolution */
  bool acc1 = Gamma_Pos(px1,py1,pz1);

  /* Simulate EMCal energy resolution */
  double e1t = 0.;
  int itw01 = -1;
  double ximp1,yimp1,zimp1;

  if( acc1 ) acc1 = Gamma_En(px1,py1,pz1,e1t,itw01,ximp1,yimp1,zimp1);
  if( acc1 ) {
    double k1 = e1t/e1;
    //    k1 = 1; // !!!
    Vpart[NPart].SetPxPyPzE(px1*k1,py1*k1,pz1*k1,e1*k1);
    itw_part[NPart] = itw01;
    NPart++;
  }

  /* Loop over all clusters in calorimeter */
  for(int i=0; i<MAXPEAK; i++)
    if(itw_part[i] >= 0)
    {
      /* Get sector */
      int iy, iz;
      anatools::TowerLocation(itw_part[i], sec_part[i], iy, iz);
    }

  /* Get NPeak */
  GetNpeak();

  return;
}

void AnaFastMC::geom_sim(const TVector3 &beam)
{
  /* Initialize parameters */
  ResetClusters();

  double px = beam.Px();
  double py = beam.Py();
  double pz = beam.Pz();

  int sec, iz0, iy0;
  double zz, yy; // hit position in tower frame (in tower units)
  double phi;
  double ximp, yimp, zimp;

  bool acc = GetImpactSectorTower(px,py,pz, sec,iz0,iy0,zz,yy,phi,ximp,yimp,zimp);
  if( acc ) {
    int itw = anatools::TowerID(sec, iy0, iz0);
    if( itw >= 0 ) {
      itw_part[NPart] = itw;
      sec_part[NPart] = sec;
      NPart++;
    }
  }

  return;
}

void AnaFastMC::ResetClusters()
{
  NPart = 0;
  NPeak = 0;
  for(int i=0; i<MAXPEAK; i++)
  {
    Vpart[i].SetPxPyPzE(0.,0.,0.,0.);
    itw_part[i] = -1;
    sec_part[i] = -1;
  }
  return;
}

void AnaFastMC::ResetTowerEnergy()
{
  for( int is=0; is<NSEC; is++ )
    for( int iy=0; iy<NY; iy++ )
      for( int iz=0; iz<NZ; iz++ )
        eTwr[is][iy][iz] = 0.;
  return;
}

void AnaFastMC::FillTowerEnergy( int sec, int iy, int iz, double e )
{
  if( sec<0 || sec>NSEC-1 || iy<0 || iy>NY-1 || iz<0 || iz>NZ-1 ) return;
  eTwr[sec][iy][iz] += e;
  return;
}

double AnaFastMC::GetETwr( int sec, int iy, int iz )
{
  if( sec<0 || sec>NSEC-1 || iy<0 || iy>NY-1 || iz<0 || iz>NZ-1 ) return 0;
  return eTwr[sec][iy][iz];
}

int AnaFastMC::GetNpeak()
{
  const double eThresh = 0.02;

  int npeak=0;
  double e;
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

bool AnaFastMC::CheckWarnMap( int itower )
{
  if( itower < 0 || itower >= n_twrs ) return false;
  int sec, iy, iz;
  anatools::TowerLocation(itower, sec, iy, iz);
  if( tower_status_sim[sec][iy][iz] != 0 )
    return true;
  return false;
}

bool AnaFastMC::InDCAcceptance( const TVector3 &v3_part, int charge )
{
  double px = v3_part.Px();
  double py = v3_part.Py();
  double pz = v3_part.Pz();
  double pt = sqrt(px*px + py*py);
  double alpha = -charge / pt;

  double eta = v3_part.Eta();
  double phi = v3_part.Phi();
  if(phi < -PI/2.) phi += PI*2.;
  string ns = pz > 0. ? "N" : "S";
  string we = phi < PI/2. ? "W" : "E";
  string nswe = ns + we;
  double board = 0.;
  if( we == "W" )
    board = ( 0.573231 + phi - 0.0046 * cos( phi + 0.05721 ) ) / 0.01963496;
  else
    board = ( 3.72402 - phi + 0.008047 * cos( phi + 0.87851 ) ) / 0.01963496;

  if( fabs(eta) < 0.35 &&
      ( ( phi > -PI*3./16. && phi < PI*5./16. ) ||
        ( phi > PI*11./16. && phi < PI*19./16. )
      ) &&
      !dcdeadmap->IsDead(nswe, board, alpha)
    )
    return true;

  return false;
}

double AnaFastMC::GetEMCResponse(int id, double mom)
  // EMCal hadron response
{
  if( id == PY_GAMMA ||
      id == PY_ELECTRON )
    return mom;

  double eout = 0.;

  if( gRandom->Rndm() < 0.43 )
  {
    eout = gRandom->Gaus(0.28,0.04);
    if( eout < 0. ) eout = 0.;
  }
  else
  {
    double emean = 0.26 + 0.44*mom - 0.0026*mom*mom;
    double ss = (0.41+0.007*mom)*emean;
    eout = gRandom->Gaus(emean,ss);
    if( eout < 0. ) eout = 0.;
  }

  return eout;
}

bool AnaFastMC::GetImpactSectorTower(double px, double py, double pz,  
    int& sec, int& iz, int& iy, double& zz, double& yy, 
    double& phi0, double& ximp, double& yimp, double& zimp )
  //
  // Returns: sector number (sec), tower coordinates in sector (iz,iy)
  //          and hit position in tower (zz,yy) in tower units;
  //          All that - for shower CG
  //          and impact position on sector face
  // 
  // Edge cut is NOT done here
  //
{
  static double zvert = 0;

  // EMCal Geometry
  // X axis is towards W0

  const double xsec_sc = 507;
  const double zsize_sc = 5.58; // Tower size in Z
  const double ysize_sc = 5.57; // Tower size in Y
  //  const double zsize_sc = 6.5; // Tower size in Z
  //  const double ysize_sc = 6.5; // Tower size in Y

  const double xsec_gl = 540;
  const double zsize_gl = 4.092; // Tower size in Z
  const double ysize_gl = 4.106; // Tower size in Y

  // Sector phi position: 0<phi<2*PI; -phi=2*PI-phi
  const double phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8, 
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  double xsec, zsize, ysize;

  sec=-1; iz=-1; iy=-1; zz=0; yy=0;

  TVector3 vv(px,py,pz);
  double phi = vv.Phi(); // -PI<phi<PI !!

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
    double theta = PI-vv.Theta();
    vv.SetPhi(phi);
    vv.SetTheta(theta);
  }
  phi0 = phi;

  double dl, x0;
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
  double zsec_min = -zsize*nZ/2;
  double ysec_min = -ysize*nY/2;

  // Center of gravity shift along Z
  double sinZ = vv.CosTheta(); // I work with PI/2-Theta
  double dz = dl*sinZ;
  double dz_gamma = x0*sinZ; // Shift of gamma related to electron
  // Center of gravity shift along Y
  double dy = dl*phi; // aY ~ sin(aY)
  double dy_gamma = x0*phi; // Shift of gamma related to electron

  double ysec = xsec*vv.Py()/vv.Px();
  double zsec = xsec*vv.Pz()/vv.Px();
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
  double ysec_cg = ysec + dy;
  double zsec_cg = zsec + dz;
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
  //if( CheckWarnMap(sec,iy,iz) ) return false;

  zz = (zsec_cg-zsec_min)/zsize-iz;
  yy = (ysec_cg-ysec_min)/ysize-iy;
  if( zz<0 || zz>1 || yy<0 || yy>1 ) printf("Error in GetImpactSectorTower: %f %f\n",zz,yy);

  return true;
}

bool AnaFastMC::GetShower(double px, double py, double pz, double& eout, int& itw )
  // Return false if outside acceptance
  // Called by Gamma_En(...)
{
  //  const double bwidth = 0.15; // Slope parameter in shower shape
  //const double bwidth = 0.20; // Slope parameter in shower shape
  //  const double bwidth = 0.005; // for xcg=0
  //  const double bwidth = 1.; // for xcg~x
  const double corr = 0.933;
  //static TRandom* gen = new TRandom();

  const double thresh = 0.010;
  //  const double thresh = 0.030; //!!!!!

  eout = 0;
  itw = -1;
  if( TMath::Sqrt(TMath::Abs(px*px+py*py+pz*pz)) <= 0.01 ) return false;

  int sec, iz0, iy0;
  double zz, yy; // hit position in tower frame (in tower units)
  double phi;
  double ximp, yimp, zimp;
  // if outside acceptance - return
  if( !GetImpactSectorTower(px,py,pz, sec,iz0,iy0,zz,yy,phi,ximp,yimp,zimp) ) return false;

  if( sec<6 ) { 
    itw = sec*2592 + iy0*72 + iz0;
  }
  else {
    itw = 15552 + (sec-6)*4608 + iy0*96 + iz0;
  }

  double en = sqrt(px*px+py*py+pz*pz);

  // The rest is valid only for PbSc (tuning is needed for PbGl) !!!!!

  // Calculate angle and energy dependent shower parameters

  TVector3 vv(px,py,pz);
  double sinZ = vv.CosTheta(); // I work with PI/2-Theta

  //  double aY = TMath::Abs(vv.Phi());
  //  while( aY > PI/16 ) aY -= PI/16;
  //  double sinY = aY; // aY ~ sin(aY)

  double sinY = TMath::Sin(phi);
  double sin2a = sinZ*sinZ + sinY*sinY;
  if( sin2a > 0.5*0.5 ) printf("Something wrong in GetShower: too big angles %f %f\n",sinZ,sinY);
  double lgE = 0;
  if( en > 0.01 ) lgE = log(en);
  else sin2a = 0;

  double par1=0.59-(1.45+0.13*lgE)*sin2a;
  double par2=0.265+(0.80+0.32*lgE)*sin2a;
  double par3=0.25+(0.45-0.036*lgE)*sin2a;
  double par4=0.42;

  // !!!!! Modify profile
  //  par3 = 0.;
  //  par4 *= 2.;
  //  par2 = 0.6; // as K.Okada

  // Calculate center of gravity
  double dz0, dy0, dz, dy, r1, r2, r3;
  double tt = 0.98 + 0.98*sqrt(en);
  double bx = 0.20 + tt*sinZ*sinZ;
  double by = 0.20 + tt*sinY*sinY;
  dz0 = zz - 0.5;
  dy0 = yy - 0.5;
  dz0 = double(0.5*TMath::SinH(double(dz0/bx))/TMath::SinH(0.5/bx));
  dy0 = double(0.5*TMath::SinH(double(dy0/by))/TMath::SinH(0.5/by));

  // Calculate ecore

  double et;
  double esum = 0;
  double ecore = 0;
  double et5x5[5][5];
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

bool AnaFastMC::Gamma_En(double px, double py, double pz, double& eout, int& itw,
    double& ximp, double& yimp, double& zimp)
  // returns false if out of acceptance
{
  // PbSc
  static double a_sc = 0.078; // official
  //  const double b_sc = 0.030; // sigma/E = a/sqrt(E) + b, official
  //  static double b_sc = 0.055; // sigma/E = a/sqrt(E) + b
  static double b_sc = 0.040; // sigma/E = a/sqrt(E) + b, for eta
  static double bc_sc = 0.012; // this is artificial const. term introduced by ecore calculations in GetShower(...); should be subtracted when GetShower(...) used
  static double cnoise_sc = 0.015; // Noise term
  //  static double corr_sc = 1.000; // For Gamma
  static double corr_sc = 1.011; // for MB

  // PbGl
  static double a_gl = 0.085;
  //  static double b_gl = 0.065; // for pi0
  static double b_gl = 0.037; // for eta
  static double bc_gl = 0.;
  static double cnoise_gl = 0.030; // Noise term
  //  static double corr_gl = 1.007;
  static double corr_gl = 1.009;

  double a, b, bc, corr, cnoise;

  //const double Eback = 0.1; // mean energy per fired tower
  //const double Rload = 0.08; // 25% of towers fired
  //int Ntwr = 4; // number of towers in Core
  double et, et1; // et2;

  eout = 0;
  itw = -1;

  double ein = sqrt(px*px+py*py+pz*pz);
  if( ein <= 0 ) return false;

  int sec, iz0, iy0;
  double zz, yy; // hit position in tower frame (in tower units)
  double phi;
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

  if(cnoise);
  //  double res = sqrt(a*a*ein + b*b*ein*ein + cnoise*cnoise);
  //  double res = sqrt(a*a*ein + b*b*ein*ein + 0.03*0.03);
  double res = sqrt(a*a*ein + b*b*ein*ein - bc*bc*ein*ein);
  //    res = sqrt(a*a*ein + b*b*ein*ein - bc*bc*ein*ein + 0.03*0.03);

  res = sqrt(res*res+0.01*0.01); // 10 MeV noise per cluster

  et = gRandom->Gaus(ein,res);
  if( et <= 0 ) et=0;
  //  et = ein; //!!!!!

  // Switching from High Gain to Low Gain
  if( et > 1 ) {
    et *= gRandom->Gaus(1,0.03);
  }

  double k = et/ein;
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

bool AnaFastMC::Gamma_Pos(double& px, double& py, double& pz )
{
  // PbSc
  const double a_sc = 0.16;
  //  const double a_sc = 0.24; // Default
  //  const double a_sc = 0.14;
  const double b_sc = 0.90;
  //  const double b_sc = 0.59; // Default
  const double c_sc = 2.00; // Pos. res. = a+b/sqrt(E) (+) (c*sinT)
  const double xsec_sc = 507;

  // PbGl
  //  const double a_gl = 0.161;
  //  const double a_gl = 0.24; // Default
  const double a_gl = 0.14;
  //  const double b_gl = 0.673; // Default
  const double b_gl = 1.15;
  // "c_gl" should be checked !!!!!
  const double c_gl = 2.; // Pos. res. = a (+) b/sqrt(E) (+) (c*sinT)
  const double xsec_gl = 540;

  double a, b, c, xsec;

  int sec, iz0, iy0;
  double zz, yy; // hit position in tower frame (in tower units)
  double phi;
  double ximp, yimp, zimp;
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

  double res = 0., dl;
  TVector3 vv(px,py,pz);
  double En = vv.Mag();

  if( sec<6 ) res = a + b/TMath::Sqrt(En);
  else        res = sqrt(a*a + b*b/En);

  // Long. fluctuations
  dl = gRandom->Gaus(0,c);

  // Resolution along Z
  double sinZ = vv.CosTheta(); // I work with PI/2-Theta
  //  double resz = 0.65*res * ( 1.+ 2.*TMath::Gaus(zz-0.5,0.,0.18) ); // To account for worse resolution in tower center
  double resz = res;
  double dz = gRandom->Gaus(0,resz);
  dz += dl*sinZ;

  // Resolution along Y
  /*
     double aY = TMath::Abs(vv.Phi());
     int nrot=0;
     while( aY > PI/16 ) { aY -= PI/16; nrot++; }
     double dy = gRandom->Gaus(0,res);
     dy += dl*aY; // aY ~ sin(aY)

     TVector3 dv(0,dy,dz);
     vv.RotateZ(-nrot*PI/16);

     vv.SetMag(TMath::Abs(En/vv.Px()*xsec)); // Now in cm
     TVector3 VV = vv + dv;
     VV.SetMag(En); // Conserve energy keeping Theta and Phi
     VV.RotateZ(nrot*PI/16);
     */

  //  double resy = 0.65*res * ( 1.+ 2.*TMath::Gaus(yy-0.5,0.,0.18) ); // To account for worse resolution in tower center
  double resy = res;
  double dy = gRandom->Gaus(0,resy);
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
