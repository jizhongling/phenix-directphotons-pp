#include "AnaPHPythiaHistos.h"

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <Fun4AllReturnCodes.h>
#include <Fun4AllHistoManager.h>
#include <getClass.h>

#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>
#include <PHPyCommon.h>

#include <TPythia6.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8) 
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <TMath.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <cstdlib>
#include <iostream>

using namespace std;

/* Some constants */
const double PI = TMath::Pi();
const TVector2 v2_2PI(0., 2.*PI);

/* Some cuts for photon identification */
const double ptMin = 3.;
const double eMin = 0.3;

/* Some cuts for isolation cut */
const double eClusMin = 0.15;
const double eTrkMin = 0.2;
const double eTrkMax = 15.;

AnaPHPythiaHistos::AnaPHPythiaHistos(const string &name, const char *filename):
  SubsysReco(name),
  outFileName(filename),
  phpythiaheader(nullptr),
  phpythia(nullptr),
  hm(nullptr),
  h2_proc_pt(nullptr),
  hn_photon(nullptr),
  hn_corr(nullptr)
{
}

AnaPHPythiaHistos::~AnaPHPythiaHistos()
{
}

int AnaPHPythiaHistos::Init(PHCompositeNode *topNode)
{
  /* Create histograms */
  BookHistograms();

  return EVENT_OK;
}

int AnaPHPythiaHistos::process_event(PHCompositeNode *topNode)
{
  /* Get PYTHIA Header */
  phpythiaheader = findNode::getClass<PHPythiaHeader>(topNode,"PHPythiaHeader");
  if(!phpythiaheader)
  {
    cout << PHWHERE << "Unable to get PHPythiaHeader, is Node missing?" << endl;
    return ABORTEVENT;
  }

  int proc_id = phpythiaheader->GetProcessid();
  double pt_event = phpythiaheader->GetPt();
  h2_proc_pt->Fill(pt_event, (double)proc_id);

  /* Uncomment when only calculate the weights */
  //return EVENT_OK;

  /* Get PYTHIA Particles */
  phpythia = findNode::getClass<PHPythiaContainer>(topNode,"PHPythia");
  if(!phpythia)
  {
    cout << PHWHERE << "Unable to get PHPythia, is Node missing?" << endl;
    return ABORTEVENT;
  }

  /* Loop over all particles */
  int npart = phpythia->size();
  for (int ipart=0; ipart<npart; ipart++)
  {
    TMCParticle *part = phpythia->getParticle(ipart);
    TMCParticle *parent = phpythia->getParent(part);

    /* Put particle's momentum into TVector3 */
    TVector3 v3_part(part->GetPx(), part->GetPy(), part->GetPz());
    double pt = v3_part.Pt();

    /* Only consider high-pT particles 
     * Test if particle passes energy threshold
     * and in Central Arm acceptance */
    if( pt < ptMin ||
        part->GetEnergy() < eMin ||
        !InAcceptance(v3_part) )
      continue;

    /* Test if particle is stable photon */
    if( part->GetKF() == PY_GAMMA &&
        part->GetKS() == 1 )
    {
      /* Test if particle is direct photon */
      int direct = 0;
      if( parent &&
          part->GetKF() == PY_GAMMA &&
          parent->GetKF() == PY_GAMMA &&
          parent->GetKS() == 21 &&
          !parent->GetParent() )
        direct = 1;

      /* Vary isolation cone angles and cut energy ratios */
      for(int icone=0; icone<11; icone++)
        for(int ie=0; ie<20; ie++)
        {
          double rcone = icone * 0.1;
          double econe = SumETruth(part, rcone);
          double re = (ie+1) * 0.02;
          int isolated = 0;
          if( econe < re * part->GetEnergy() )
            isolated = 1;

          double fill_hn_photon[] = {pt, rcone, re, (double)isolated, (double)direct};
          hn_photon->Fill(fill_hn_photon);
        }

      /* Fill two particle phi correlations */
      FillCorrelation(part, 1-direct);
    } // photon

    /* Test if particle is pi0 */
    if( part->GetKF() == PY_PIZERO )
      /* Fill two particle phi correlations */
      FillCorrelation(part, 2);
  }

  return EVENT_OK;
}

int AnaPHPythiaHistos::End(PHCompositeNode *topNode)
{
  /* Write histograms */
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void AnaPHPythiaHistos::BookHistograms()
{
  /* Initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  h2_proc_pt = new TH2F("h2_proc_pt", "Pt and process id", 500,0.,50., 120,0.5,120.5);
  hm->registerHisto(h2_proc_pt);

  /* Uncomment when only calculate the weights */
  return;

  /* pT bins */
  const int npT = 30;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  const int nbins_hn_photon[] = {npT, 11, 20, 2, 2};
  const double xmin_hn_photon[] = {0., -0.05, 0.01, -0.5, -0.5};
  const double xmax_hn_photon[] = {0., 1.05, 0.41, 1.5, 1.5};
  hn_photon = new THnSparseF("hn_photon", "Photon counts;p_{T} [GeV];rcone;renergy;isolated;direct;",
      5, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, pTbin);
  hm->registerHisto(hn_photon);

  const int nbins_hn_corr[] = {npT, npT, 30, 3};
  const double xmin_hn_corr[] = {0., 0., 0., -0.5};
  const double xmax_hn_corr[] = {0., 0., PI, 2.5};
  hn_corr = new THnSparseF("hn_corr", "Tow-particle correlation;p^{trig}_{T} [GeV];p^{part}_{T} [GeV];#Delta#phi [rad];type;",
      4, nbins_hn_corr, xmin_hn_corr, xmax_hn_corr);
  hn_corr->SetBinEdges(0, pTbin);
  hn_corr->SetBinEdges(1, pTbin);
  hm->registerHisto(hn_corr);

  return;
}

double AnaPHPythiaHistos::SumETruth(const TMCParticle *pref, double rcone)
{
  /* Sum up all energy in cone around particle */
  double econe = 0;

  /* Get reference vector */
  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());
  TVector2 v2_pref = v3_pref.EtaPhiVector();

  /* Loop over all particles */
  int npart = phpythia->size();
  for(int ipart2=0; ipart2<npart; ipart2++)
  {
    TMCParticle *part2 = phpythia->getParticle(ipart2);

    /* Only consider stable and interacting particles
     * skip if pointer identical to 'reference' particle */
    int id = abs(part2->GetKF());
    if( part2 == pref || part2->GetKS() != 1 ||
        id == PY_MU || id == PY_NU_E || id == PY_NU_MU )
      continue;

    /* Get particle vector */
    TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());
    TVector2 v2_part2 = v3_part2.EtaPhiVector();

    /* Test if particle is in Central Arm acceptance
     * and passes energy threshold */
    if( part2->GetEnergy() < eClusMin ||
        !InAcceptance(v3_part2) )
      continue;

    /* Check if particle within cone */
    TVector2 v2_diff(v2_part2 - v2_pref);
    if( v2_diff.Y() > PI ) v2_diff -= v2_2PI;
    else if( v2_diff.Y() < -PI ) v2_diff += v2_2PI;
    if( v2_diff.Mod() < rcone )
      econe += part2->GetEnergy();
  } // ipart2

  return econe;
}

void AnaPHPythiaHistos::FillCorrelation(const TMCParticle *pref, int type)
{
  /* Get reference vector */
  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());

  /* Loop over all particles */
  int npart = phpythia->size();
  for(int ipart2=0; ipart2<npart; ipart2++)
  {
    TMCParticle *part2 = phpythia->getParticle(ipart2);

    /* Only consider stable and interacting particles
     * skip if pointer identical to 'reference' particle */
    int id = abs(part2->GetKF());
    if( part2 == pref || part2->GetKS() != 1 ||
        id == PY_MU || id == PY_NU_E || id == PY_NU_MU )
      continue;

    /* Get particle vector */
    TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());

    /* Test if particle is in Central Arm acceptance
     * and passes energy threshold */
    if( part2->GetEnergy() < eMin ||
        !InAcceptance(v3_part2) )
      continue;

    double fill_hn_corr[] = {v3_pref.Pt(), v3_part2.Pt(), fabs(v3_part2.DeltaPhi(v3_pref)), (double)type};
    hn_corr->Fill(fill_hn_corr);
  } // ipart2

  return;
}

bool AnaPHPythiaHistos::InAcceptance(const TVector3 &v3_part)
{
  if(  fabs(v3_part.Eta()) < 0.35 &&
      (fabs(v3_part.Phi()) < PI/4. ||
       fabs(v3_part.Phi()) > PI*3./4.) )
    return true;
  return false;
}
