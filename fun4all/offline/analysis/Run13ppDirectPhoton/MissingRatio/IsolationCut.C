#include "IsolationCut.h"

#include <AnaTrk.h>

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
#include <TF1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

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
                                                   _events_photon(0),
                                                   _hn_energy_cone( NULL ),
                                                   _hn_energy_cone_reco( NULL ),
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
  _v_pid_neutral.push_back( PHOTON );
  _v_pid_neutral.push_back( NEUTRINO );
  _v_pid_neutral.push_back( PIZERO );
  _v_pid_neutral.push_back( KAON_0_LONG );
  _v_pid_neutral.push_back( NEUTRON );
  _v_pid_neutral.push_back( KAON_0_SHORT );
  _v_pid_neutral.push_back( ETA );
  _v_pid_neutral.push_back( LAMBDA );
  _v_pid_neutral.push_back( SIGMA_0 );
  _v_pid_neutral.push_back( XI_0 );
  _v_pid_neutral.push_back( ANTINEUTRON );
  _v_pid_neutral.push_back( ANTILAMBDA );
  _v_pid_neutral.push_back( ANTISIGMA_0 );
  _v_pid_neutral.push_back( ANITXI_0 );

  return EVENT_OK;
}

int IsolationCut::process_event(PHCompositeNode *topNode)
{
  _ievent++;
  unsigned nphotons = 0;

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

  /* number of tracks and clusters */
  //unsigned n_truth_particles = truth_particles->size();
  //unsigned n_truth_emcclusters = truth_emcclusters->size();
  //unsigned n_reco_tracks = reco_tracks->get_npart();
  unsigned n_reco_emcclusters = reco_emcclusters->size();
  //
  //cout << "Truth particles  found in event:   " << n_truth_particles << endl;
  //cout << "Truth EMC  cluster found in event: " << n_truth_emcclusters << endl;
  //cout << "Reco tracks  found in event:       " << n_reco_tracks << endl;
  //cout << "Reco cluster found in event:       " << n_reco_emcclusters << endl;

  ///* loop over all reconstructed tracks */
  //for ( unsigned i = 0; i < n_reco_tracks; i++ )
  //  {
  //    cout << "Track: p = " << reco_tracks->get_mom(i) << endl;
  //  }

  /* Associate truth ECMal cluster with truth particles for easy lookup between the two.
   * Use separate maps for charged and neutral particles.
   * Map key is cluster id. */
  //typedef map<int,AnaTrk*> map_Ana_t;
  //map_Ana_t map_cluster_track_charged;
  //map_Ana_t map_cluster_track_neutral;
  //
  //for(unsigned itrk=0; itrk<n_truth_particles; itrk++)
  //  {
  //    emcGeaTrackContent *emctrk = truth_particles->get(itrk);
  //    AnaTrk *track = new AnaTrk(emctrk, truth_emcclusters, (int*)_tower_status);
  //
  //    if(track)
  //    {
  //      /* only keep track if it's associated with a cluster in the EMCAL */
  //      if( track->emcclus )
  //        {
  //          /* charged truth particle? */
  //          if ( find( _v_pid_neutral.begin(), _v_pid_neutral.end(), track->pid ) == _v_pid_neutral.end() )
  //            {
  //              //cout << "Found a charged track from a track with PID " << track->pid << endl;
  //              map_cluster_track_charged.insert( make_pair( track->cid, track ) );
  //            }
  //          /* neutral truth particle? */
  //          else
  //            {
  //              //cout << "Found a neutral track from a track with PID " << track->pid << endl;
  //              map_cluster_track_neutral.insert( make_pair( track->cid, track ) );
  //            }
  //        }
  //    }
  //  }

  /* store id's of photon candidate clusters */
  vector<unsigned> cluster_photons;

  /* cuts for photon candidate selection */
  float cluster_ecore_min = 0.3; //GeV
  float photon_ecore_min = 1; //GeV
  float photon_prob_min = 0.02;

  /* loop over all cluster in EMCAL - add to photon candidate collection if it meets photon criteria */
  for( unsigned icluster=0; icluster < n_reco_emcclusters; icluster++ )
    {
      /* get pointer to reco cluster */
      emcClusterContent *reco_emc_cluster_i = reco_emcclusters->getCluster( icluster );

      /* is photon candidate? */
      if ( reco_emc_cluster_i->ecore() < cluster_ecore_min ||
           reco_emc_cluster_i->ecore() < photon_ecore_min ||
           reco_emc_cluster_i->prob_photon() < photon_prob_min )
        continue;

      /* charge veto? use truth information for charge veto? */
      //map_Ana_t::iterator track_cluster = map_cluster_track_charged.find( emc_cluster->id() );
      //if( track_cluster != map_cluster_track_charged.end() )
      //        {
      //          continue;
      //        }

      /* check tower- is in guard ring of 10(12) towers ob PbSc(PbGl)? */
      // ...

      /* append cluster id to photon candiadte collection */
      cluster_photons.push_back( icluster );
      nphotons++;
    }

  /* loop over photon candidate clusters */
  for( unsigned i=0; i < cluster_photons.size(); i++ )
    {
      /* reset variables that store cluster properties */
      ResetClusterTreeVariables();

      /* get cluster id from photons vector */
      unsigned icluster = cluster_photons.at( i );

      /* get pointer to reco cluster corresponding to this reco cluster */
      emcClusterContent *reco_emc_cluster_i = reco_emcclusters->getCluster( icluster );

      /* get pointer to truth cluster corresponding to this reco cluster */
      emcGeaClusterContent *truth_emc_cluster_i = truth_emcclusters->getCluster( icluster );

      /* get truth particle associated with this cluster */
      emcGeaTrackContent *truth_particle_i = FindTruthParticle( truth_emc_cluster_i );

      /* get particle ID of truth particle */
      unsigned pid_i = 0;

      if ( truth_particle_i )
        {
          pid_i = truth_particle_i->get_pid();
        }
      else
        {
          cerr << "ERROR: Truth particle not found for this cluster" << endl;
        }

      /* get particle ID of parent of truth particle (if any) */
      unsigned parent_i = 0;
      emcGeaTrackContent *truth_parent_i = NULL;

      if ( truth_particle_i )
        {
          truth_parent_i = truth_particle_i->get_trackcontainer()->find( truth_particle_i->get_parent_trkno() );
        }

      if ( truth_parent_i )
        {
          parent_i = truth_parent_i->get_pid();
        }

      /* what's the origin of this cluster- direct photon, pi0 photon, something else? */
      //      bool isPromptPhoton = true;
      //      if ( lvl0_trk->anclvl == 0 )
      //        isPromptPhoton = true;

      /* loop over different isolation cone radii */
      for ( int round = 0; round < 10; round ++ )
        {
          double rcone = 0.1 + round * 0.1;

          /* add up calorimeter energy in cone around cluster */
          double econe_emcal = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, 0.3, rcone);

          /* cut for track selection */
          float track_pmin = 0.2; //GeV
          float track_pmax = 15; //GeV

          /* add up tracking energy in cone around cluster */
          double econe_track = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone);

          cout << "econe EMCal: " << econe_emcal << " , econe Track: " << econe_track << endl;
        }

      cout << "Cluster: E = " << truth_emc_cluster_i->ecore() << ", truth PID: " << pid_i << " , truth Parent: " << parent_i << " trkno = " << truth_particle_i->get_trkno() << endl;
    }

  //
  //
  //      /* determine energy fraction of cone */
  //          double econe_frac = econe / emc_cluster->ecore();
  //          //cout << "Energy in cone, fraction: " << econe << ", " << econe_frac << endl;
  //
  //          /* fill histogram */
  //          double fill[] = {(double)emc_cluster->ecore(), econe, econe_frac, rcone, (double)isPromptPhoton};
  //          _hn_energy_cone->Fill( fill );
  //
  //        } // END: loop over isolation cone radius
  //    } // END: loop over all cluster in EMCAL

  if ( cluster_photons.size() > 0 )
    _events_photon++;

  return EVENT_OK;
}

int IsolationCut::End(PHCompositeNode *topNode)
{
  _file_output->cd();

  if ( _hn_energy_cone )
    _hn_energy_cone->Write();

  if ( _hn_energy_cone_reco )
    _hn_energy_cone_reco->Write();

  _file_output->Close();

  cout << "---------------------------------------------------" << endl;
  cout << "IsolationCut module:" << endl;
  cout << "Number of events processed:            " << _ievent << endl;
  cout << "Number of events with photon in EMCAL: " << _events_photon << endl;
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


void IsolationCut::ResetClusterTreeVariables( )
{
  return;
}
