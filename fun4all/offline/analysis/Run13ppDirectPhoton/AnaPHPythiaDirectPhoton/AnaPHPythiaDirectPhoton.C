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
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

using namespace std;

AnaPHPythiaDirectPhoton::AnaPHPythiaDirectPhoton(const std::string &name): SubsysReco(name),
									   _pdg_db(nullptr),
                                                                           phpythiaheader(0),
                                                                           phpythia(0),
                                                                           _tree_event_truth(nullptr),
                                                                           _ievent(0),
                                                                           _pythia_event(0),
                                                                           _pythia_processid(0),
                                                                           _truth_pid(0),
                                                                           _truth_parentid(0),
                                                                           _truth_ispromptphoton(0),
                                                                           _truth_ptot(0),
                                                                           _truth_pt(0),
                                                                           _truth_etot(0),
                                                                           _truth_eta(0),
                                                                           _truth_phi(0),
                                                                           _output_file_name("test.root"),
                                                                           _fout(NULL)
{
  /* define set of isolation cone sizes */
  _v_iso_conesize.push_back( 0.1 );
  _v_iso_conesize.push_back( 0.3 );
  _v_iso_conesize.push_back( 0.5 );
  _v_iso_conesize.push_back( 1.0 );

  /* set 0 as initial values for cone sums */
  for ( unsigned i = 0; i < _v_iso_conesize.size(); i++ )
    {
      _v_iso_eemcal.push_back( 0.0 );
      _v_iso_ptrack.push_back( 0.0 );
    }

  /* get PDGDatabase object */
  _pdg_db = new TDatabasePDG();
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
  _tree_event_truth->Branch( "pythia_event", &_pythia_event );
  _tree_event_truth->Branch( "pythia_processid", &_pythia_processid );

  /* Add particle branches */
  _tree_event_truth->Branch( "pid", &_truth_pid );
  _tree_event_truth->Branch( "parentid", &_truth_parentid );
  _tree_event_truth->Branch( "ispromptphoton", &_truth_ispromptphoton );
  _tree_event_truth->Branch( "ptot", &_truth_ptot );
  _tree_event_truth->Branch( "pt", &_truth_pt );
  _tree_event_truth->Branch( "etot", &_truth_etot );
  _tree_event_truth->Branch( "eta", &_truth_eta );
  _tree_event_truth->Branch( "phi", &_truth_phi );

  /* Add islolation cut branches */
  _tree_event_truth->Branch( "iso_conesize", &_v_iso_conesize );
  _tree_event_truth->Branch( "iso_eemcal", &_v_iso_eemcal );
  _tree_event_truth->Branch( "iso_ptrack", &_v_iso_ptrack );

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
  //  cout << "Process: " << phpythiaheader->GetProcessid() << endl;
  //  cout << "KS\tKF\tpid\tname\thistory" << endl;

  /* Update event parameter information */
  _pythia_event = phpythiaheader->GetEvt();
  _pythia_processid = phpythiaheader->GetProcessid();

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

      /* test if photon passes energy threshold */
      double photon_minEnergy = 1.0;
      if ( part->GetEnergy() < photon_minEnergy )
	continue;

      /* Set tree varables to store particle information */
      _truth_pid = part->GetKF();
      _truth_parentid = 0;
      if ( parent )
        _truth_parentid = parent->GetKF();
      _truth_ispromptphoton = isPromptPhoton;
      _truth_ptot = v_part.P();
      _truth_pt = v_part.Perp();
      _truth_etot = part->GetEnergy();
      _truth_eta = v_part.Eta();
      _truth_phi = v_part.Phi();

      /* Get eneries and momenta in different size isolation cones */
      for ( unsigned icone = 0; icone < _v_iso_conesize.size(); icone++ )
        {
          _v_iso_eemcal.at( icone ) = SumEEmcal( part , _v_iso_conesize.at( icone ) );
          _v_iso_ptrack.at( icone ) = SumPTrack( part , _v_iso_conesize.at( icone ) );
        }

      /* Fill tree */
      _tree_event_truth->Fill();

    }

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


float AnaPHPythiaDirectPhoton::SumEEmcal( TMCParticle* pref, float rcone )
{
  /* sum up all energy in cone around particle */
  double econe = 0;

  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());

  int npart = phpythia->size();
  for (int ipart2=0; ipart2<npart; ipart2++)
    {
      TMCParticle *part2 = phpythia->getParticle(ipart2);

      /* only consider stable particles, skip if pointer identical to 'reference' particle */
      if ( part2->GetKS() == 1 && part2 != pref )
        {
          /* skip particles that are not photons and not electrons because they'll give
	   * lower signal in EMCal */
          if ( ( part2->GetKF() != 22 ) &&
	       ( part2->GetKF() != 11 ) )
	    continue;

          TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());

          /* check if energy within cone */
          if ( v3_pref.Angle( v3_part2 ) < rcone )
            econe += part2->GetEnergy();
        }
    }

  return econe;
}


float AnaPHPythiaDirectPhoton::SumPTrack( TMCParticle* pref, float rcone )
{
  /* sum up all energy in cone around particle */
  double pcone = 0;

  TVector3 v3_pref(pref->GetPx(), pref->GetPy(), pref->GetPz());

  int npart = phpythia->size();
  for (int ipart2=0; ipart2<npart; ipart2++)
    {
      TMCParticle *part2 = phpythia->getParticle(ipart2);

      /* only consider stable particles, skip if pointer identical to 'reference' particle */
      if ( part2->GetKS() == 1 && part2 != pref )
        {
          /* skip particles that do not have an electric charge */
	  TParticlePDG* part2_pdg = _pdg_db->GetParticle( part2->GetKF() );
          if ( part2_pdg->Charge() == 0 )
	    continue;

          TVector3 v3_part2(part2->GetPx(), part2->GetPy(), part2->GetPz());

          /* check if particle within cone */
          if ( v3_pref.Angle( v3_part2 ) < rcone )
            pcone += part2->GetEnergy();
        }
    }

  return pcone;
}
