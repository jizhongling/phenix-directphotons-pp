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
									   _select_photons(false),
									   _pdg_db(nullptr),
                                                                           phpythiaheader(0),
                                                                           phpythia(0),
                                                                           _tree_event_truth(nullptr),
									   _ievent(0),
                                                                           _output_file_name("test.root"),
                                                                           _fout(nullptr)
{
  /* define set of isolation cone sizes in mrad */
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
  if(!_pdg_db)
  {
    cerr << "No pdg_db" << endl;
    exit(1);
  }
}

AnaPHPythiaDirectPhoton::~AnaPHPythiaDirectPhoton()
{
}

// Done at the beginning of processing
int AnaPHPythiaDirectPhoton::Init(PHCompositeNode *topNode)
{
  /* create output file */
  _fout = new TFile(_output_file_name.c_str(),"RECREATE");

  /* Add cluster properties to map that defines output tree */
  float dummy = 0;
  vector< float > vdummy;

  _branchmap_event.insert( make_pair( "eventcounter" , dummy ) );
  _branchmap_event.insert( make_pair( "pythia_event" , dummy ) );
  _branchmap_event.insert( make_pair( "pythia_processid" , dummy ) );

  _branchmap_mcparticles.insert( make_pair( "t_pid" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_parentpid" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_ispromptphoton" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_ptot" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_pt" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_etot" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_eta" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_phi" , vdummy ) );

  /* Create tree for information about full event truth */
  _tree_event_truth = new TTree("event_truth", "a Tree with global event information and EM truth");

  /* Add event branches */
  for ( map< string , float >::iterator iter = _branchmap_event.begin();
        iter != _branchmap_event.end();
        ++iter )
    {
      _tree_event_truth->Branch( (iter->first).c_str(),
				     &(iter->second) );
    }

  /* Add particle branches */
  for ( map< string , vector<float> >::iterator iter = _branchmap_mcparticles.begin();
        iter != _branchmap_mcparticles.end();
        ++iter )
    {
      _tree_event_truth->Branch( (iter->first).c_str(),
				 &(iter->second) );
    }

  /* Add islolation cut branches */
  //  _tree_event_truth->Branch( "iso_conesize", &_v_iso_conesize );
  //  _tree_event_truth->Branch( "iso_eemcal", &_v_iso_eemcal );
  //  _tree_event_truth->Branch( "iso_ptrack", &_v_iso_ptrack );

  return EVENT_OK;
}

// Done at the end of processing
int AnaPHPythiaDirectPhoton::End(PHCompositeNode *topNode)
{
  _fout->cd();

  if ( _tree_event_truth )
    _tree_event_truth->Write();

  if ( _pdg_db )
    delete _pdg_db;

  _fout->Close();

  return EVENT_OK;
}

// Process each event
int AnaPHPythiaDirectPhoton::process_event(PHCompositeNode *topNode)
{
  /* increment event counter */
  _ievent++;

  /* reset variables that store cluster properties */
  ResetBranchVariables();

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

  /* Set event parameters for putput trees */
  ( _branchmap_event.find("eventcounter") )->second  = _ievent;
  ( _branchmap_event.find("pythia_event") )->second  = phpythiaheader->GetEvt();
  ( _branchmap_event.find("pythia_processid") )->second  = phpythiaheader->GetProcessid();

  /* Loop over all particles & Fill output tree */
  int npart = phpythia->size();
  for (int ipart=0; ipart<npart; ipart++)
    {
      TMCParticle *part = phpythia->getParticle(ipart);
      TMCParticle *parent = phpythia->getParent(part);

      /* test if particle is stable photon */
      if ( _select_photons & ( part->GetKF() != 22 ||
			      part->GetKS() != 1 ) )
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

      /* test if particle is in Central Arm acceptance- continue if not */
      Double_t Eta = v_part.Eta();
      if ( _select_photons & ( ! ( fabs(Eta)<0.35 ) ) )
        continue;

      /* test if photon passes energy threshold- continue of not */
      double photon_minEnergy = 0.15;
      if ( _select_photons & ( ! ( part->GetEnergy() > photon_minEnergy ) ) )
      	continue;

      /* Set tree varables to store particle information */
      ( _branchmap_mcparticles.find("t_pid") )->second.push_back( part->GetKF() );

      int truth_parentpid = 0;
      if ( parent )
	truth_parentpid =  parent->GetKF();
      ( _branchmap_mcparticles.find("t_parentpid") )->second.push_back( truth_parentpid );
      ( _branchmap_mcparticles.find("t_ispromptphoton") )->second.push_back( isPromptPhoton );
      ( _branchmap_mcparticles.find("t_ptot") )->second.push_back(  v_part.P() );
      ( _branchmap_mcparticles.find("t_pt") )->second.push_back( v_part.Perp() );
      ( _branchmap_mcparticles.find("t_etot") )->second.push_back( part->GetEnergy() );
      ( _branchmap_mcparticles.find("t_eta") )->second.push_back( v_part.Eta() );
      ( _branchmap_mcparticles.find("t_phi") )->second.push_back( v_part.Phi() );

      /* Get eneries and momenta in different size isolation cones */
      for ( unsigned icone = 0; icone < _v_iso_conesize.size(); icone++ )
        {
          _v_iso_eemcal.at( icone ) = SumEEmcal( part , _v_iso_conesize.at( icone ) );
          _v_iso_ptrack.at( icone ) = SumPTrack( part , _v_iso_conesize.at( icone ) );
        }
    }
  /* Fill tree */
  _tree_event_truth->Fill();

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
          /* skip particles that are not photons - they'll be accounted for in momentum */
          if ( ( part2->GetKF() != 22 ) )
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

void AnaPHPythiaDirectPhoton::ResetBranchVariables( )
{
  /* Event branches */
  for ( map< string , float >::iterator iter = _branchmap_event.begin();
        iter != _branchmap_event.end();
        ++iter)
    {
      (iter->second) = NAN;
    }

  /* Particle branches */
  for ( map< string , vector<float> >::iterator iter = _branchmap_mcparticles.begin();
        iter != _branchmap_mcparticles.end();
        ++iter)
    {
      (iter->second).clear();
    }

  return;
}
