#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <TPythia6.h>
#include <TVector3.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHTimeStamp.h>
#include <Fun4AllReturnCodes.h>
#include <getClass.h>

#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>
#include <AnaPHPythiaDirectPhoton.h>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

AnaPHPythiaDirectPhoton::AnaPHPythiaDirectPhoton(const std::string &name): SubsysReco(name),
									   phpythiaheader(0),
									   phpythia(0),
									   _tree_event_truth(nullptr),
									   _ievent(0),
									   _truth_pid(0),
									   _truth_parentid(0),
									   _truth_anclvl(0),
									   _truth_ptot(0),
									   _truth_pt(0),
									   _truth_eta(0),
									   _truth_phi(0),
									   _output_file_name("test.root"),
									   _fout(NULL)
{
}

AnaPHPythiaDirectPhoton::~AnaPHPythiaDirectPhoton()
{
}

// Done at the beginning of processing
int AnaPHPythiaDirectPhoton::Init(PHCompositeNode *topNode)
{
  _fout = new TFile(_output_file_name.c_str(),"RECREATE");

  /* Create tree for information about full event truth */
  _tree_event_truth = new TTree("event_truth", "a Tree with global event information and EM truth");

  /* Add event branches */
  _tree_event_truth->Branch( "event", &_ievent );
  _tree_event_truth->Branch( "pid", &_truth_pid );
  _tree_event_truth->Branch( "parentid", &_truth_parentid );
  _tree_event_truth->Branch( "anclvl", &_truth_anclvl );
  _tree_event_truth->Branch( "ptot", &_truth_ptot );
  _tree_event_truth->Branch( "pt", &_truth_pt );
  _tree_event_truth->Branch( "eta", &_truth_eta );
  _tree_event_truth->Branch( "phi", &_truth_phi );


//  int ndim_hn_photon = 5;
//  int nbins_hn_photon[] =   {100 , 100 ,  10000 ,   11  ,  2  };
//  double xmin_hn_photon[] = {  0.,   0.,      0.,  -0.05, -0.5};
//  double xmax_hn_photon[] = {100., 100.,    100.,   1.05,  1.5};
//
//  _hECone = new THnSparseF("hn_EConePhoton",
//                           "Energy in cone around Photon; E [GeV]; E_cone [GeV]; f_cone; r_cone [rad]; isPromptPhoton;",
//                           ndim_hn_photon,
//                           nbins_hn_photon,
//                           xmin_hn_photon,
//                           xmax_hn_photon );
//
  return EVENT_OK;
}

// Done at the end of processing
int AnaPHPythiaDirectPhoton::End(PHCompositeNode *topNode)
{
  _fout->cd();

  if ( _tree_event_truth )
    _tree_event_truth->Write();

  _fout->Close();

  return EVENT_OK;
}

// Process each event
int AnaPHPythiaDirectPhoton::process_event(PHCompositeNode *topNode)
{
  /* increment event counter */
  _ievent++;

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

  // Print some header information
  //  cout << "Event: " << phpythiaheader->GetEvt() << "\t" << phpythiaheader->GetNpart() << endl;

  //  cout << "KS\tKF\tpid\tname\thistory" << endl;

  /* Loop over all particles & Fill output tree */
  int npart = phpythia->size();
  for (int ipart=0; ipart<npart; ipart++)
    {
      TMCParticle *part = phpythia->getParticle(ipart);
      TMCParticle *parent = phpythia->getParent(part);

      /* test if particle is stable photon */
      if ( part->GetKF() != 22 ||
           part->GetKS() != 1 )
        continue;

      /* test if particle passed energy threshold */
      double photon_minEnergy = 1.0;
      if ( part->GetEnergy() < photon_minEnergy )
        continue;

      /* test if particle is prompt photon */
      bool isPromptPhoton = false;
      if ( parent )
        {
          if ( part->GetKF() == 22 &&
               parent->GetKF() == 22 &&
               parent->GetKS() == 21 &&
               !(parent->GetParent()) )
            {
              isPromptPhoton = true;
            }
        }

      /* Convert particle into TLorentzVector */
      Float_t px = part->GetPx();
      Float_t py = part->GetPy();
      Float_t pz = part->GetPz();
      Float_t energy = part->GetEnergy();
      TLorentzVector v_part(px,py,pz,energy);

      /* test if particle is in Central Arm acceptance */
      Double_t Eta = v_part.Eta();
      if ( fabs(Eta)<0.35 )
        continue;

      /* Set tree varables to store particle information */
      _truth_pid = part->GetKF();
      _truth_parentid = 0;
      if ( parent )
	_truth_parentid = parent->GetKF();
      _truth_anclvl = (!isPromptPhoton);
      _truth_ptot = v_part.P();
      _truth_pt = v_part.Perp();
      _truth_eta = v_part.Eta();
      _truth_phi = v_part.Phi();

      _tree_event_truth->Fill();

      //      cout << InCentralArmAcceptance(v_part) << endl;

      //        {
      //          cout << setw(4) << phpythia->getLineNumber(part)
      //               << "\t" << part->GetKS()
      //               << setw(8)  <<  part->GetKF() // << "\t" << part->GetEnergy()
      //               << setw(12) << part->GetName()
      //               << "\t" << part->GetEnergy()
      //            //<< "\t" << ipart
      //            //<< "\t" << part->GetParent()
      //            //<< "\t" << part->GetFirstChild()
      //            //<< "\t" << part->GetLastChild()
      //            //<< "\t" << phpythia->getChildNumber(part)
      //               << "\t" << phpythia->getHistoryString(part)
      //               << endl;


      /* loop over different cone radii */
//      for ( int round = 0; round < 10; round ++ )
//        {
//          double rcone = 0.1 + round * 0.1;
//
//          /* sum up all energy in cone around particle */
//          double econe = 0;
//          TVector3 v3_gamma(part->GetPx(), part->GetPy(), part->GetPz());
//
//          for (int ipart2=0; ipart2<npart; ipart2++)
//            {
//              TMCParticle *part2 = phpythia->getParticle(ipart2);
//
//              /* only consider stable particles, skip if pointer identical to 'reference' particle */
//              if ( part2->GetKS() == 1 && part2 != part )
//                {
//                  TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());
//
//                  /* test if particle is in Central Arm acceptance */
//                  Float_t px = part2->GetPx();
//                  Float_t py = part2->GetPy();
//                  Float_t pz = part2->GetPz();
//                  Float_t energy = part2->GetEnergy();
//                  TLorentzVector v_part2(px,py,pz,energy);
//
//                  //Double_t Eta = v_part2.Eta();
//                  //if ( fabs(Eta)<0.35 )
//                  //  continue;
//                  //
//                  //if ( part2->GetEnergy() < otherParticle_minEnergy )
//                  //  continue;
//
//                  if ( v3_gamma.Angle( v3_part2 ) < rcone )
//                    econe += part2->GetEnergy();
//                }
//            }
//          double econe_frac = econe / part->GetEnergy();
//          //cout << "Energy in cone, fraction: " << econe << ", " << econe_frac << endl;
//
//          double fill[] = {part->GetEnergy(), econe, econe_frac, rcone, (double)isPromptPhoton};
//
//          _hECone->Fill( fill );
//        }
//    }


  //  cout << endl << "Example 1 : Get ancestor Ks0(pid=310)'s energy" << endl;
  //  for (int ipart=0; ipart<npart; ipart++)  {
  //    TMCParticle *part = phpythia->getParticle(ipart);
  //    if (part->GetKS() != 1) continue;
  //
  //    TMCParticle *anc = phpythia->hasAncestor(part, 310);
  //    float energy = 0.0;
  //    if (anc) {
  //      energy = anc->GetEnergy();
  //      cout << setw(12) << phpythia->getLineNumber(part)
  //           << setw(12) << part->GetName()
  //           << " 's anc "
  //           << setw(12) << anc->GetName() << setw(12) << energy << endl;
  //    }
    }
  //
  //  cout << endl << "Example 2 : Get stable particle pi-, which has ancestor K0(311) and Ks0(pid=310),  get K0's energy" << endl;
  //  for (int ipart=0; ipart<npart; ipart++)  {
  //    TMCParticle *part = phpythia->getParticle(ipart);
  //    if (part->GetKF() != -211) continue;
  //
  //    if (phpythia->hasAncestor(part, 311, 310)) {
  //      TMCParticle *anc = phpythia->hasAncestor(part, 311);
  //      float energy = 0.0;
  //      if (anc) {
  //        energy = anc->GetEnergy();
  //        cout << setw(12) << phpythia->getLineNumber(part)
  //             << setw(12) << part->GetName()
  //             << " 's anc "
  //             << setw(12) << anc->GetName() << setw(12) << energy << endl;
  //      }
  //    }
  //  }
  //  cout << endl;

  // some useful functions
  // phpythia->getHistoryString(part)               // return a TString with part's all ancestor
  // phpythia->getParent(part)                      // return pointer to part's parent TMCParticle
  // phpythia->getChild(part, i)                    // return pointer to part's ith child TMCParticle
  // phpythia->getChildNumber(part)                 // return part's child number
  // phpythia->hasAncestor(part, pid1, pid2, pid3)  // whether part has all ancestor pid1, pid2, pid3
  // phpythia->hasAncestorAbs(part, pid)            // return the nearest ancestor pointer with abs(KF)==pid
  // phpythia->hasAncestorRange(part, pid1, pid2)   // return the nearest ancestor whose pid1<=KF<=pid2


  return EVENT_OK;
}
