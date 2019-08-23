#include "IsolationCut.h"

#include <AnaToolsTowerID.h>
#include <EMCWarnmapChecker.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <PHGlobal.h>
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
  _emcwarnmap( NULL ),
  _hn_energy_cone( NULL ),
  _hn_energy_cone_reco( NULL ),
  _tree_recocluster(nullptr),
  _tree_mcparticles(nullptr),
  _truth_pid(-9999.),
  _truth_parentpid(-9999.),
  _truth_anclvl(-9999.),
  _truth_ptot(-9999.),
  _truth_pt(-9999.),
  _truth_eta(-9999.),
  _truth_phi(-9999.),
  _output_file_name("IsolationCut_output.root"),
  _file_output( NULL )
{
  /* construct output file names */
  _output_file_name = "histos/IsolationCut-";
  _output_file_name.append(filename);
}

IsolationCut::~IsolationCut()
{
}

int IsolationCut::Init(PHCompositeNode *topNode)
{
  /* create output file */
  _file_output = new TFile( _output_file_name.c_str(), "RECREATE" );

  /* Initialize EMC warnmap checker */
  _emcwarnmap = new EMCWarnmapChecker();
  if(!_emcwarnmap)
  {
    cerr << "No emcwarnmap" << endl;
    exit(1);
  }

  /* Add cluster properties to map that defines output tree */
  float dummy = 0;
  vector< float > vdummy;

  _branchmap_event.insert( make_pair( "eventcounter" , dummy ) );

  _branchmap_mcparticles.insert( make_pair( "t_pid" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_parentpid" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_anclvl" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_ptot" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_pt" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_eta" , vdummy ) );
  _branchmap_mcparticles.insert( make_pair( "t_phi" , vdummy ) );

  /* Cluster information */
  _branchmap_cluster.insert( make_pair( "cluster_ecore" , vdummy ) ); // cluster energy
  _branchmap_cluster.insert( make_pair( "cluster_pt" , vdummy ) ); // cluster transverse momentum
  _branchmap_cluster.insert( make_pair( "cluster_prob" , vdummy ) ); // cluster em-like probability
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r01" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r02" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r03" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r04" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r05" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r06" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r07" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r08" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r09" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_rcone_r10" , vdummy ) ); // cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r01" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r02" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r03" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r04" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r05" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r06" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r07" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r08" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r09" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_emcal_r10" , vdummy ) ); // energy in EMCal clusters within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r01" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r02" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r03" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r04" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r05" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r06" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r07" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r08" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r09" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_econe_tracks_r10" , vdummy ) ); // charged tracks momenta within cone radius
  _branchmap_cluster.insert( make_pair( "cluster_truth_pdgpid" , vdummy ) ); // pid of associated truth particle
  _branchmap_cluster.insert( make_pair( "cluster_truth_pdgparentpid" , vdummy ) ); // parent pid of associated truth particle
  _branchmap_cluster.insert( make_pair( "cluster_truth_pid" , vdummy ) ); // pid of associated truth particle
  _branchmap_cluster.insert( make_pair( "cluster_truth_parentpid" , vdummy ) ); // parent pid of associated truth particle
  _branchmap_cluster.insert( make_pair( "cluster_truth_anclvl" , vdummy ) ); // ancestry level of associated truth particle

  /* Create tree for information about full event */
  _tree_recocluster = new TTree("recocluster", "EMCal cluster and isolation cone information");

  /* Add event branches */
  for ( map< string , float >::iterator iter = _branchmap_event.begin();
      iter != _branchmap_event.end();
      ++iter )
  {
    _tree_recocluster->Branch( (iter->first).c_str(),
        &(iter->second) );
  }

  /* Add cluster branches */
  for ( map< string , vector<float> >::iterator iter = _branchmap_cluster.begin();
      iter != _branchmap_cluster.end();
      ++iter )
  {
    _tree_recocluster->Branch( (iter->first).c_str(),
        &(iter->second) );
  }

  /* Create tree for information about full event truth */
  _tree_mcparticles = new TTree("mcparticles", "MC truth particles and global event information");

  /* Add event branches */
  for ( map< string , float >::iterator iter = _branchmap_event.begin();
      iter != _branchmap_event.end();
      ++iter )
  {
    _tree_mcparticles->Branch( (iter->first).c_str(),
        &(iter->second) );
  }

  /* Add particle branches */
  for ( map< string , vector<float> >::iterator iter = _branchmap_mcparticles.begin();
      iter != _branchmap_mcparticles.end();
      ++iter )
  {
    _tree_mcparticles->Branch( (iter->first).c_str(),
        &(iter->second) );
  }

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

  /* reset variables that store cluster properties */
  ResetBranchVariables();

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

  /* Set event parameters for putput trees */
  ( _branchmap_event.find("eventcounter") )->second  = _ievent;

  /* Loop over all truth tracks and store photon information */
  for( unsigned iparticle=0; iparticle < truth_particles->size(); iparticle++ )
  {
    /* get pointer to reco cluster */
    emcGeaTrackContent *truth_particle_i = truth_particles->get( iparticle );
    //      emcGeaTrackContent *truth_parent_i = truth_particles->get_common_parent( truth_particle_i, truth_particle_i );

    /* fill tree variables */
    ( _branchmap_mcparticles.find("t_pid") )->second.push_back( truth_particle_i->get_pid() );

    int truth_parentpid = 0;
    if ( truth_particle_i->get_parent_trkno() != 0 )
      truth_parentpid = truth_particles->find( truth_particle_i->get_parent_trkno() )->get_pid();
    ( _branchmap_mcparticles.find("t_parentpid") )->second.push_back(truth_parentpid);
    ( _branchmap_mcparticles.find("t_anclvl") )->second.push_back(truth_particle_i->get_anclvl());
    ( _branchmap_mcparticles.find("t_ptot") )->second.push_back(truth_particle_i->get_ptot());
    ( _branchmap_mcparticles.find("t_pt") )->second.push_back(truth_particle_i->get_pt());

    TVector3 v( truth_particle_i->get_px(), truth_particle_i->get_py(), truth_particle_i->get_pz() );

    ( _branchmap_mcparticles.find("t_eta") )->second.push_back(v.Eta());
    ( _branchmap_mcparticles.find("t_phi") )->second.push_back(v.Phi());
  }
  /* fill tree */
  _tree_mcparticles->Fill();

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

    /* check warn map and guard ring - is tower a good tower? */
    if ( !_emcwarnmap->IsGoodTower(reco_emc_cluster_i) )
      continue;

    /* is photon candidate? */
    if ( reco_emc_cluster_i->ecore() < cluster_ecore_min ||
        reco_emc_cluster_i->ecore() < photon_ecore_min ||
        reco_emc_cluster_i->prob_photon() < photon_prob_min )
      continue;

    /* charge veto? use truth information for charge veto? */
    // ...

    /* append cluster id to photon candiadte collection */
    cluster_photons.push_back( icluster );
    _event_nphotons++;
  }

  /* loop over photon candidate clusters */
  for( unsigned i=0; i < cluster_photons.size(); i++ )
  {
    /* get cluster id from photons vector */
    unsigned icluster = cluster_photons.at( i );

    /* get pointer to reco cluster corresponding to this reco cluster */
    emcClusterContent *reco_emc_cluster_i = reco_emcclusters->getCluster( icluster );

    /* calculate cluster transverse momentum */
    float cluster_ecore = reco_emc_cluster_i->ecore();
    TVector3 cluster_pos( reco_emc_cluster_i->x(), reco_emc_cluster_i->y(), reco_emc_cluster_i->z() );
    float cluster_pt = cluster_ecore * ( cluster_pos.Perp() / cluster_pos.Mag() );

    /* get cluster photon probability */
    float cluster_prob = reco_emc_cluster_i->prob_photon();

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
    float rcone_r06 = 0.6;
    float rcone_r07 = 0.7;
    float rcone_r08 = 0.8;
    float rcone_r09 = 0.9;
    float rcone_r10 = 1.0;

    /* add up calorimeter energy in cone around cluster */
    double econe_emcal_r01 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r01 );
    double econe_emcal_r02 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r02 );
    double econe_emcal_r03 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r03 );
    double econe_emcal_r04 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r04 );
    double econe_emcal_r05 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r05 );
    double econe_emcal_r06 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r06 );
    double econe_emcal_r07 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r07 );
    double econe_emcal_r08 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r08 );
    double econe_emcal_r09 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r09 );
    double econe_emcal_r10 = SumEmcalEnergyInCone(reco_emc_cluster_i, reco_emcclusters, cluster_emin, rcone_r10 );

    /* add up tracking energy in cone around cluster */
    double econe_track_r01 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r01 );
    double econe_track_r02 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r02 );
    double econe_track_r03 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r03 );
    double econe_track_r04 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r04 );
    double econe_track_r05 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r05 );
    double econe_track_r06 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r06 );
    double econe_track_r07 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r07 );
    double econe_track_r08 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r08 );
    double econe_track_r09 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r09 );
    double econe_track_r10 = SumTrackEnergyInCone(reco_emc_cluster_i, reco_tracks, track_pmin, track_pmax, rcone_r10 );

    /* Update tree branch variables */
    _branchmap_cluster.find( "cluster_ecore" )->second.push_back( cluster_ecore );
    _branchmap_cluster.find( "cluster_pt" )->second.push_back( cluster_pt );
    _branchmap_cluster.find( "cluster_prob" )->second.push_back( cluster_prob );
    _branchmap_cluster.find( "cluster_rcone_r01" )->second.push_back( rcone_r01 );
    _branchmap_cluster.find( "cluster_rcone_r02" )->second.push_back( rcone_r02 );
    _branchmap_cluster.find( "cluster_rcone_r03" )->second.push_back( rcone_r03 );
    _branchmap_cluster.find( "cluster_rcone_r04" )->second.push_back( rcone_r04 );
    _branchmap_cluster.find( "cluster_rcone_r05" )->second.push_back( rcone_r05 );
    _branchmap_cluster.find( "cluster_rcone_r06" )->second.push_back( rcone_r06 );
    _branchmap_cluster.find( "cluster_rcone_r07" )->second.push_back( rcone_r07 );
    _branchmap_cluster.find( "cluster_rcone_r08" )->second.push_back( rcone_r08 );
    _branchmap_cluster.find( "cluster_rcone_r09" )->second.push_back( rcone_r09 );
    _branchmap_cluster.find( "cluster_rcone_r10" )->second.push_back( rcone_r10 );
    _branchmap_cluster.find( "cluster_econe_emcal_r01" )->second.push_back( econe_emcal_r01 );
    _branchmap_cluster.find( "cluster_econe_emcal_r02" )->second.push_back( econe_emcal_r02 );
    _branchmap_cluster.find( "cluster_econe_emcal_r03" )->second.push_back( econe_emcal_r03 );
    _branchmap_cluster.find( "cluster_econe_emcal_r04" )->second.push_back( econe_emcal_r04 );
    _branchmap_cluster.find( "cluster_econe_emcal_r05" )->second.push_back( econe_emcal_r05 );
    _branchmap_cluster.find( "cluster_econe_emcal_r06" )->second.push_back( econe_emcal_r06 );
    _branchmap_cluster.find( "cluster_econe_emcal_r07" )->second.push_back( econe_emcal_r07 );
    _branchmap_cluster.find( "cluster_econe_emcal_r08" )->second.push_back( econe_emcal_r08 );
    _branchmap_cluster.find( "cluster_econe_emcal_r09" )->second.push_back( econe_emcal_r09 );
    _branchmap_cluster.find( "cluster_econe_emcal_r10" )->second.push_back( econe_emcal_r10 );
    _branchmap_cluster.find( "cluster_econe_tracks_r01" )->second.push_back( econe_track_r01 );
    _branchmap_cluster.find( "cluster_econe_tracks_r02" )->second.push_back( econe_track_r02 );
    _branchmap_cluster.find( "cluster_econe_tracks_r03" )->second.push_back( econe_track_r03 );
    _branchmap_cluster.find( "cluster_econe_tracks_r04" )->second.push_back( econe_track_r04 );
    _branchmap_cluster.find( "cluster_econe_tracks_r05" )->second.push_back( econe_track_r05 );
    _branchmap_cluster.find( "cluster_econe_tracks_r06" )->second.push_back( econe_track_r06 );
    _branchmap_cluster.find( "cluster_econe_tracks_r07" )->second.push_back( econe_track_r07 );
    _branchmap_cluster.find( "cluster_econe_tracks_r08" )->second.push_back( econe_track_r08 );
    _branchmap_cluster.find( "cluster_econe_tracks_r09" )->second.push_back( econe_track_r09 );
    _branchmap_cluster.find( "cluster_econe_tracks_r10" )->second.push_back( econe_track_r10 );
    _branchmap_cluster.find( "cluster_truth_pdgpid" )->second.push_back( 0 );
    _branchmap_cluster.find( "cluster_truth_pdgparentpid" )->second.push_back( 0 );
    _branchmap_cluster.find( "cluster_truth_pid" )->second.push_back( pid_i );
    _branchmap_cluster.find( "cluster_truth_parentpid" )->second.push_back( parent_id );
    _branchmap_cluster.find( "cluster_truth_anclvl" )->second.push_back( anclvl_i );
  }

  /* Fill tree */
  _tree_recocluster->Fill();

  return EVENT_OK;
}

int IsolationCut::End(PHCompositeNode *topNode)
{
  _file_output->cd();

  if ( _emcwarnmap )
    delete _emcwarnmap;

  if ( _tree_recocluster )
    _tree_recocluster->Write();

  if ( _tree_mcparticles )
    _tree_mcparticles->Write();

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

    /* skip clusters in towers flagged bad by warn map,
       but allow for edge towers (status 20) */
    if ( _emcwarnmap->IsBadTower(emccluster_check) )
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
  for ( map< string , vector<float> >::iterator iter = _branchmap_cluster.begin();
      iter != _branchmap_cluster.end();
      ++iter)
  {
    (iter->second).clear();
  }

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
