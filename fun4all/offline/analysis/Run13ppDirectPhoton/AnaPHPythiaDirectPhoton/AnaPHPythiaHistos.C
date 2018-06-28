#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <TPythia6.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8) 
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TVector3.h>
#include <TH1.h>
#include <THnSparse.h>

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHTimeStamp.h>
#include <Fun4AllReturnCodes.h>
#include <Fun4AllHistoManager.h>
#include <getClass.h>
#include <TOAD.h>

#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>
#include "AnaPHPythiaHistos.h"

using namespace std;

AnaPHPythiaHistos::AnaPHPythiaHistos(const string &name, const char *filename): SubsysReco(name)
{
  phpythia = 0;		// array of pythia particles
  phpythiaheader = 0;	// pythia header

  // Get PDGDatabase object
  pdg_db = new TDatabasePDG();

  // Construct output file names
  outFileName = "histos/AnaPHPythiaHistos-";
  outFileName.append(filename);

  hm = 0;
  hn_photon = 0;
  hn_corr = 0;
  for(int ih=0; ih<3; ih++)
    h_photon[ih] = 0;
}

AnaPHPythiaHistos::~AnaPHPythiaHistos()
{
}

// Done at the beginning of processing
int AnaPHPythiaHistos::Init(PHCompositeNode *topNode)
{
  // Create and register histograms
  BookHistograms();

  return EVENT_OK;
}

// Done at the end of processing
int AnaPHPythiaHistos::End(PHCompositeNode *topNode)
{
  // Write histogram output to ROOT file
  hm->dumpHistos();
  delete hm;
  delete pdg_db;

  return EVENT_OK;
}

// Process each event
int AnaPHPythiaHistos::process_event(PHCompositeNode *topNode)
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

  //static int event = 0;
  //bool none = true;
  //if(event > 1000) exit(0);

  // Loop over all particles & Fill output histograms
  int npart = phpythia->size();
  for (int ipart=0; ipart<npart; ipart++)
  {
    TMCParticle *part = phpythia->getParticle(ipart);
    TMCParticle *parent = phpythia->getParent(part);

    // Convert particle into TVector3
    double px = part->GetPx();
    double py = part->GetPy();
    double pz = part->GetPz();
    TVector3 v3_part(px,py,pz);

    //double eta = v3_part.Eta();
    //double phi = v3_part.Phi();
    //if(   abs(eta) < 0.35 &&
    //    ( abs(phi-PI/16.) < PI/4. ||
    //      phi > PI*11./16. ||
    //      phi < -PI*13./16. ) )
    //{
    //  if(none)
    //  {
    //    cout << "\n\nEvent " << event++ << endl;
    //    none = false;
    //  }
    //  cout << part->GetKF() << "(" << part->GetKS() << "), ";
    //}
    //continue;

    // Test if particle is in Central Arm acceptance
    // and passes energy threshold
    if( part->GetEnergy() < 0.3 ||
        v3_part.Pt() < 0.01 ||
        abs(v3_part.Eta()) > 0.35 ||
        (abs(v3_part.Phi()) > PI/4. &&
         abs(v3_part.Phi()) < PI*3./4.) )
      continue;

    // Test if particle is stable photon
    if( part->GetKF() == 22 &&
        part->GetKS() == 1 )
    {
      h_photon[0]->Fill(v3_part.Pt());

      // Test if particle is prompt photon
      if( parent &&
          part->GetKF() == 22 &&
          abs(parent->GetKF()) < 100 )
        h_photon[1]->Fill(v3_part.Pt());

      // Test if particle is direct photon
      int prompt = 0;
      if( parent &&
          part->GetKF() == 22 &&
          parent->GetKF() == 22 &&
          parent->GetKS() == 21 &&
          !(parent->GetParent()) )
      {
        h_photon[2]->Fill(v3_part.Pt());
        prompt = 1;
      }

      for(int icone=0; icone<11; icone++)
        for(int ie=0; ie<20; ie++)
        {
          double rcone = icone * 0.1;
          double econe = SumETruth(part, rcone);
          double re = (ie+1) * 0.02;
          int isolated = 0;
          if( econe < re * part->GetEnergy() )
            isolated = 1;

          double fill_hn_photon[] = {v3_part.Pt(), rcone, re, (double)isolated, (double)prompt};
          hn_photon->Fill(fill_hn_photon);
        }

      FillCorrelation(part, 1-prompt);
    } // photon

    // Test if particle is pi0
    if( part->GetKF() == 111 )
      FillCorrelation(part, 2);
  }

  return EVENT_OK;
}

void AnaPHPythiaHistos::BookHistograms()
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

  for(int ih=0; ih<3; ih++)
  {
    h_photon[ih] = new TH1F(Form("h_photon_%d",ih), "Photon counts;p_{T} [GeV];", npT,pTbin);
    hm->registerHisto(h_photon[ih]);
  }

  const int nbins_hn_photon[] = {npT, 11, 20, 2, 2};
  const double xmin_hn_photon[] = {0., -0.05, 0.01, -0.5, -0.5};
  const double xmax_hn_photon[] = {0., 1.05, 0.41, 1.5, 1.5};
  hn_photon = new THnSparseF("hn_photon", "Photon counts;p_{T} [GeV];rcone;renergy;isolated;prompt;",
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
  // Sum up all energy in cone around particle
  double econe = 0;

  // Get reference vector
  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());
  if( v3_pref.Pt() < 0.01 ) return econe;
  TVector2 v2_pref = v3_pref.EtaPhiVector();

  int npart = phpythia->size();

  for (int ipart2=0; ipart2<npart; ipart2++)
  {
    TMCParticle *part2 = phpythia->getParticle(ipart2);

    // Only consider stable particles, skip if pointer identical to 'reference' particle
    if( part2 != pref && part2->GetKS() == 1 )
    {
      // Get particle vector
      TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());
      if( v3_part2.Pt() < 0.01 ) continue;
      TVector2 v2_part2 = v3_part2.EtaPhiVector();

      // Test if particle is in Central Arm acceptance
      // and passes energy threshold
      if( part2->GetEnergy() < 0.15 ||
          abs(v2_part2.X()) > 0.35 ||
          (abs(v2_part2.Y()) > PI/4. &&
           abs(v2_part2.Y()) < PI*3./4.) )
        continue;

      // Check if particle within cone
      if( (v2_part2-v2_pref).Mod() < rcone )
        econe += part2->GetEnergy();
    }
  }

  return econe;
}

void AnaPHPythiaHistos::FillCorrelation(const TMCParticle *pref, int type)
{
  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());

  int npart = phpythia->size();
  for (int ipart2=0; ipart2<npart; ipart2++)
  {
    TMCParticle *part2 = phpythia->getParticle(ipart2);

    // Only consider stable particles, skip if pointer identical to 'reference' particle
    if( part2 != pref && part2->GetKS() == 1 )
    {
      TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());

      // Test if particle is in Central Arm acceptance
      // and passes energy threshold
      if( part2->GetEnergy() < 0.3 ||
          v3_part2.Pt() < 0.01 ||
          abs(v3_part2.Eta()) > 0.35 ||
          (abs(v3_part2.Phi()) > PI/4. &&
           abs(v3_part2.Phi()) < PI*3./4.) )
        continue;

      double fill_hn_corr[] = {v3_pref.Pt(), v3_part2.Pt(), abs(v3_pref.DeltaPhi(v3_part2)), (double)type};
      hn_corr->Fill(fill_hn_corr);
    }
  }

  return;
}
