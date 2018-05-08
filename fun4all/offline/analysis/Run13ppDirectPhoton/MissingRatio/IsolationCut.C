#include "IsolationCut.h"

#include <PHGlobal.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <emcClusterContainer.h>

#include <PHCentralTrack.h>

#include <TOAD.h>
#include <phool.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

/* ROOT includes */
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TVector3.h>
//#include <TLorentzVector.h>
//#include <TDatabasePDG.h>

/* STL includes */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/foreach.hpp>

using namespace std;

IsolationCut::IsolationCut(const char *filename) : _ievent(0),
                                                   _event_nphotons(0),
                                                   _hn_energy_cone( NULL ),
                                                   _hn_energy_cone_reco( NULL ),
                                                   _tree_event_cluster(nullptr),
                                                   _tree_event_truth(nullptr),
						   _truth_pid(0),
						   _truth_parentid(0),
						   _truth_anclvl(0),
						   _truth_ptot(0),
						   _truth_pt(0),
						   _truth_eta(0),
						   _truth_phi(0),
                                                   _output_file_name("IsolationCut_output.root"),
                                                   _file_output( NULL )
{
  /* construct output file names */
  _output_file_name = "histos/IsolationCut-";
  _output_file_name.append(filename);

  /* initialize array for tower status */
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        _tower_status[isector][ibiny][ibinz] = 0;

}

IsolationCut::~IsolationCut()
{
}

int IsolationCut::Init(PHCompositeNode *topNode)
{
  /* create output file */
  _file_output = new TFile( _output_file_name.c_str(), "RECREATE" );

  /* Add cluster properties to map that defines output tree */
  float dummy = 0;

  _map_cluster_branches.insert( make_pair( "cluster_ecore" , dummy ) ); // cluster energy
  _map_cluster_branches.insert( make_pair( "cluster_pt" , dummy ) ); // cluster transverse momentum
  _map_cluster_branches.insert( make_pair( "cluster_rcone_r01" , dummy ) ); // cone radius
  _map_cluster_branches.insert( make_pair( "cluster_rcone_r02" , dummy ) ); // cone radius
  _map_cluster_branches.insert( make_pair( "cluster_rcone_r03" , dummy ) ); // cone radius
  _map_cluster_branches.insert( make_pair( "cluster_rcone_r04" , dummy ) ); // cone radius
  _map_cluster_branches.insert( make_pair( "cluster_rcone_r05" , dummy ) ); // cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_emcal_r01" , dummy ) ); // energy in EMCal clusters within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_emcal_r02" , dummy ) ); // energy in EMCal clusters within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_emcal_r03" , dummy ) ); // energy in EMCal clusters within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_emcal_r04" , dummy ) ); // energy in EMCal clusters within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_emcal_r05" , dummy ) ); // energy in EMCal clusters within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_tracks_r01" , dummy ) ); // charged tracks momenta within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_tracks_r02" , dummy ) ); // charged tracks momenta within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_tracks_r03" , dummy ) ); // charged tracks momenta within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_tracks_r04" , dummy ) ); // charged tracks momenta within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_econe_tracks_r05" , dummy ) ); // charged tracks momenta within cone radius
  _map_cluster_branches.insert( make_pair( "cluster_truth_pid" , dummy ) ); // pid of associated truth particle
  _map_cluster_branches.insert( make_pair( "cluster_truth_parentid" , dummy ) ); // parent id of associated truth particle
  _map_cluster_branches.insert( make_pair( "cluster_truth_anclvl" , dummy ) ); // ancestry level of associated truth particle

  /* Create tree for information about full event */
  _tree_event_cluster = new TTree("event_cluster", "a Tree with global event information and EM cluster");

  /* Add event branches */
  _tree_event_cluster->Branch( "event", &_ievent );
  _tree_event_cluster->Branch( "event_nphotons", &_event_nphotons );

  /* Add cluster branches */
  for ( map< string , float >::iterator iter = _map_cluster_branches.begin();
        iter != _map_cluster_branches.end();
        ++iter )
    {
      _tree_event_cluster->Branch( (iter->first).c_str(),
                                   &(iter->second) );
    }

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

  /* create output histogram */
  //int ndim_hn_photon = 5;
  //int nbins_hn_photon[] =   {100 , 100 ,  10000 ,   11  ,  2  };
  //double xmin_hn_photon[] = {  0.,   0.,      0.,  -0.05, -0.5};
  //double xmax_hn_photon[] = {100., 100.,    100.,   1.05,  1.5};
  //
  //_hn_energy_cone = new THnSparseF("hn_energy_cone",
  //                                 "Energy in cone around Photon; E [GeV]; E_cone [GeV]; f_cone; r_cone [rad]; isPromptPhoton;",
  //                                 ndim_hn_photon,
  //                                 nbins_hn_photon,
  //                                 xmin_hn_photon,
  //                                 xmax_hn_photon );
  //
  //_hn_energy_cone_reco = (THnSparseF*)_hn_energy_cone->Clone("hn_energy_cone_reco");

  /* create vector with neutral particle PID's */
  /* create vector of PID's of neutral particles */
  //_v_pid_neutral.push_back( PHOTON );
  //_v_pid_neutral.push_back( NEUTRINO );
  //_v_pid_neutral.push_back( PIZERO );
  //_v_pid_neutral.push_back( KAON_0_LONG );
  //_v_pid_neutral.push_back( NEUTRON );
  //_v_pid_neutral.push_back( KAON_0_SHORT );
  //_v_pid_neutral.push_back( ETA );
  //_v_pid_neutral.push_back( LAMBDA );
  //_v_pid_neutral.push_back( SIGMA_0 );
  //_v_pid_neutral.push_back( XI_0 );
  //_v_pid_neutral.push_back( ANTINEUTRON );
  //_v_pid_neutral.push_back( ANTILAMBDA );
  //_v_pid_neutral.push_back( ANTISIGMA_0 );
  //_v_pid_neutral.push_back( ANITXI_0 );

  return EVENT_OK;
}

int IsolationCut::process_event(PHCompositeNode *topNode)
{
  _ievent++;
  _event_nphotons = 0;

  /* global event info */
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
    {
      cout << "Cannot find PHGlobal" << endl;
      return DISCARDEVENT;
    }

  /* TRUTH track info (truth particles container) */
  emcGeaTrackContainer *truth_particles = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!truth_particles)
    {
      cout << "Cannot find emcGeaTrackContainer" << endl;
      return DISCARDEVENT;
    }

  /* TRUTH EMC cluster info */
  emcGeaClusterContainer *truth_emcclusters = truth_particles->GetClusters();
  if(!truth_emcclusters)
    {
      cout << "Cannot find emcGeaClusterContainer" << endl;
      return DISCARDEVENT;
    }

  /* Reco tracks from Drift Chamber */
  PHCentralTrack *reco_tracks = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!reco_tracks)
    {
      cout << "Cannot find PHCentralTrack" << endl;
      return DISCARDEVENT;
    }

  /* Reco cluster from Electromagnetic Calorimeter */
  emcClusterContainer* reco_emcclusters = findNode::getClass<emcClusterContainer> (topNode, "emcClusterContainer");
  if(!reco_emcclusters)
    {
      cout << "Cannot find emcClusterContainer" << endl;
      return DISCARDEVENT;
    }

  /* Loop over all turth tracks and store photon information */
  for( unsigned iparticle=0; iparticle < truth_particles->size(); iparticle++ )
    {
      /* get pointer to reco cluster */
      emcGeaTrackContent *truth_particle_i = truth_particles->get( iparticle );
      //      emcGeaTrackContent *truth_parent_i = truth_particles->get_common_parent( truth_particle_i, truth_particle_i );

      /* fill tree variables */
      _truth_pid = truth_particle_i->get_pid();
      _truth_parentid = 0;
      //      if ( truth_parent_i )
      //	_truth_parentid = truth_parent_i->get_pid();
      if ( truth_particle_i->get_parent_trkno() != 0 )
	_truth_parentid = truth_particles->find( truth_particle_i->get_parent_trkno() )->get_pid();
      _truth_anclvl = truth_particle_i->get_anclvl();
      _truth_ptot = truth_particle_i->get_ptot();
      _truth_pt = truth_particle_i->get_pt();

      TVector3 v( truth_particle_i->get_px(), truth_particle_i->get_py(), truth_particle_i->get_pz() );//, truth_particle_i->get_ekin() );

      _truth_eta = v.Eta(); // -2 * log( atan2( truth_particle_i->get_pt() , truth_particle_i->get_ptot() ) / 2.0 );
      _truth_phi = v.Phi(); // atan2( truth_particle_i->get_py() , truth_particle_i->get_px() );

      /* fill tree */
      _tree_event_truth->Fill();
    }

  /* store id's of photon candidate clusters */
  vector<unsigned> cluster_photons;

  /* cuts for photon candidate selection */
  float cluster_ecore_min = 0.3; //GeV
  float photon_ecore_min = 1; //GeV
  float photon_prob_min = 0.02;

  /* loop over all cluster in EMCAL - add to photon candidate collection if it meets photon criteria */
  for( unsigned icluster=0; icluster < reco_emcclusters->size(); icluster++ )
    {
      /* get pointer to reco cluster */
      emcClusterContent *reco_emc_cluster_i = reco_emcclusters->getCluster( icluster );

      /* is photon candidate? */
      if ( reco_emc_cluster_i->ecore() < cluster_ecore_min ||
           reco_emc_cluster_i->ecore() < photon_ecore_min ||
           reco_emc_cluster_i->prob_photon() < photon_prob_min )
        continue;

      /* charge veto? use truth information for charge veto? */
      // ...

      /* check tower- is in guard ring of 10(12) towers ob PbSc(PbGl)? */
      // ...

      /* append cluster id to photon candiadte collection */
      cluster_photons.push_back( icluster );
      _event_nphotons++;
    }

  /* loop over photon candidate clusters */
  for( unsigned i=0; i < cluster_photons.size(); i++ )
    {
      /* reset variables that store cluster properties */
      ResetBranchVariables();

      /* get cluster id from photons vector */
      unsigned icluster = cluster_photons.at( i );

      /* get pointer to reco cluster corresponding to this reco cluster */
      emcClusterContent *reco_emc_cluster_i = reco_emcclusters->getCluster( icluster );

      /* calculate cluster transverse momentum */
      float cluster_ecore = reco_emc_cluster_i->ecore();
      TVector3 cluster_pos( reco_emc_cluster_i->x(), reco_emc_cluster_i->y(), reco_emc_cluster_i->z() );
      float cluster_pt = cluster_ecore * ( cluster_pos.Perp() / cluster_pos.Mag() );

      /* get pointer to truth cluster corresponding to this reco cluster */
      emcGeaClusterContent *truth_emc_cluster_i = truth_emcclusters->getCluster( icluster );

      /* get truth particle associated with this cluster */
      emcGeaTrackContent *truth_particle_i = FindTruthParticle( truth_emc_cluster_i );

      /* get particle ID of truth particle */
      unsigned pid_i = 0;
      unsigned anclvl_i = 0;

      if ( truth_particle_i )
        {
          pid_i = truth_particle_i->get_pid();
          anclvl_i = truth_particle_i->get_anclvl();
        }
      else
        {
          cerr << "ERROR: Truth particle not found for this cluster" << endl;
        }

      /* get particle ID of parent of truth particle (if any) */
      unsigned parent_id = 0;
      emcGeaTrackContent *truth_parent_i = NULL;

      if ( truth_particle_i )
        {
          truth_parent_i = truth_particle_i->get_trackcontainer()->find( truth_particle_i->get_parent_trkno() );
        }

      if ( truth_parent_i )
        {
          parent_id = truth_parent_i->get_pid();
        }

      /* cut for calorimter cluster in cone */
      float cluster_emin = 0.3;

      /* cut for track selection in cone */
      float track_pmin = 0.2; //GeV
      float track_pmax = 15; //GeV

      /* cone radius definitions */
      float rcone_r01 = 0.1;
      float rcone_r02 = 0.2;
      float rcone_r03 = 0.3;
      float rcone_r04 = 0.4;
      float rcone_r05 = 0.5;

      /* add up calorimeter energy in cone around cluster */
      double econe_emcal_r01 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r01 );
      double econe_emcal_r02 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r02 );
      double econe_emcal_r03 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r03 );
      double econe_emcal_r04 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r04 );
      double econe_emcal_r05 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r05 );

      /* add up tracking energy in cone around cluster */
      double econe_track_r01 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r01 );
      double econe_track_r02 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r02 );
      double econe_track_r03 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r03 );
      double econe_track_r04 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r04 );
      double econe_track_r05 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r05 );

      /* Update tree branch variables */
      _map_cluster_branches.find( "cluster_ecore" )->second = cluster_ecore;
      _map_cluster_branches.find( "cluster_pt" )->second = cluster_pt;
      _map_cluster_branches.find( "cluster_rcone_r01" )->second = rcone_r01;
      _map_cluster_branches.find( "cluster_rcone_r02" )->second = rcone_r02;
      _map_cluster_branches.find( "cluster_rcone_r03" )->second = rcone_r03;
      _map_cluster_branches.find( "cluster_rcone_r04" )->second = rcone_r04;
      _map_cluster_branches.find( "cluster_rcone_r05" )->second = rcone_r05;
      _map_cluster_branches.find( "cluster_econe_emcal_r01" )->second = econe_emcal_r01;
      _map_cluster_branches.find( "cluster_econe_emcal_r02" )->second = econe_emcal_r02;
      _map_cluster_branches.find( "cluster_econe_emcal_r03" )->second = econe_emcal_r03;
      _map_cluster_branches.find( "cluster_econe_emcal_r04" )->second = econe_emcal_r04;
      _map_cluster_branches.find( "cluster_econe_emcal_r05" )->second = econe_emcal_r05;
      _map_cluster_branches.find( "cluster_econe_tracks_r01" )->second = econe_track_r01;
      _map_cluster_branches.find( "cluster_econe_tracks_r02" )->second = econe_track_r02;
      _map_cluster_branches.find( "cluster_econe_tracks_r03" )->second = econe_track_r03;
      _map_cluster_branches.find( "cluster_econe_tracks_r04" )->second = econe_track_r04;
      _map_cluster_branches.find( "cluster_econe_tracks_r05" )->second = econe_track_r05;
      _map_cluster_branches.find( "cluster_truth_pid" )->second = pid_i;
      _map_cluster_branches.find( "cluster_truth_parentid" )->second = parent_id;
      _map_cluster_branches.find( "cluster_truth_anclvl" )->second = anclvl_i;

      /* Fill tree */
      _tree_event_cluster->Fill();

    }

  return EVENT_OK;
}

int IsolationCut::End(PHCompositeNode *topNode)
{
  _file_output->cd();

  if ( _tree_event_cluster )
    _tree_event_cluster->Write();

  if ( _tree_event_truth )
    _tree_event_truth->Write();

  if ( _hn_energy_cone )
    _hn_energy_cone->Write();

  if ( _hn_energy_cone_reco )
    _hn_energy_cone_reco->Write();

  _file_output->Close();

  cout << "---------------------------------------------------" << endl;
  cout << "IsolationCut module:" << endl;
  cout << "Number of events processed:            " << _ievent << endl;
  cout << "Number of events with photon in EMCAL: " << _event_nphotons << endl;
  cout << "---------------------------------------------------" << endl;

  return EVENT_OK;
}


emcGeaTrackContent* IsolationCut::FindTruthParticle( emcGeaClusterContent* cluster )
{
  /* pointer to truth particle associated with cluster */
  emcGeaTrackContent* truthparticle = NULL;

  /* loop over all truth particle and pick the one with maximum deposited energy in this cluster */
  float edepMax = 0.;
  emc_tracklist_t particles_list = cluster->get_track_list();
  BOOST_FOREACH(const emc_trkno_t &ipart, particles_list)
    {
      float edep = cluster->get_edep_bytrack(ipart);
      if( edep > edepMax )
        {
          truthparticle = cluster->get_trackcontainer()->find(ipart);
          edepMax = edep;
        }
    }

  return truthparticle;
}


float IsolationCut::SumEmcalEnergyInCone( emcClusterContent* emccluster_ref,
                                          emcClusterContainer* emcclusters,
                                          double cluster_emin,
                                          double rcone )
{
  float econe = 0;

  /* loop over all cluster in EMCAL */
  for( unsigned icluster=0; icluster < emcclusters->size(); icluster++ )
    {
      /* Get cluster */
      emcClusterContent *emccluster_check = emcclusters->getCluster( icluster );

      /* avoid double counting cluster with same id as reference cluster */
      if ( emccluster_check->id() == emccluster_ref->id() )
        continue;

      /* skip cluster below energy threshold */
      if ( emccluster_check->ecore() < cluster_emin )
        continue;

      /* check if cluster is within cone around reference cluster */
      float dr = ( sqrt( pow( emccluster_ref->theta() - emccluster_check->theta() , 2 ) +
                         pow( emccluster_ref->phi()   - emccluster_check->phi()   , 2 ) ) );

      if ( dr < rcone )
        {
          /* Still here? Add cluster energy to energy in cone. */
          econe += emccluster_check->ecore();
        }
    }

  return econe;
}


float IsolationCut::SumTrackEnergyInCone( emcClusterContent* emccluster_ref,
                                          PHCentralTrack* tracks,
                                          double track_pmin,
                                          double track_pmax,
                                          double rcone )
{
  float econe = 0;

  /* loop over all charged tracks that could be associated with a cluster in EMCAL */
  for ( unsigned itrack = 0; itrack < tracks->get_npart(); itrack++ )
    {
      /* track momentum below threshold? */
      if ( ( tracks->get_mom( itrack ) < track_pmin ) ||
           ( tracks->get_mom( itrack ) > track_pmax ) )
        continue;

      /* determine track theta and phi directions */
      float track_theta = tracks->get_the0( itrack ); //atan2( emctrk->get_pt() , emctrk->get_ptot() );
      float track_phi = tracks->get_phi0( itrack ); //atan2( emctrk->get_py() , emctrk->get_px() );

      /* check if track within cone */
      float dr = ( sqrt( pow( track_theta - emccluster_ref->theta() , 2 ) +
                         pow( track_phi   - emccluster_ref->phi()   , 2 ) ) );

      if ( dr < rcone )
        {
          /* Still here? Add track momentum to energy in cone. */
          econe += tracks->get_mom( itrack );
        }
    }

  return econe;
}


void IsolationCut::ResetBranchVariables( )
{
  /* Cluster branches */
  for ( map< string , float >::iterator iter = _map_cluster_branches.begin();
        iter != _map_cluster_branches.end();
        ++iter)
    {
      (iter->second) = NAN;
    }

  return;
}
