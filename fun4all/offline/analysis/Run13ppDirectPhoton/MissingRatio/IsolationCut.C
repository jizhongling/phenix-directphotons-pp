#include "IsolationCut.h"

#include <AnaTrk.h>

#include <PHGlobal.h>

#include <emcNodeHelper.h>
#include <emcGeaTrackContainer.h>
#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

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

// global constants
const int PHOTON_PID = 1;
const int POSITRON_PID = 2;
const int ELECTRON_PID = 3;
const int NEUTRINO_PID = 4;
const int PIZERO_PID = 7;
const int KAON_0_LONG = 10;
const int NEUTRON = 13;
const int KAON_0_SHORT = 16;
const int ETA = 17;
const int LAMBDA = 18;
const int SIGMA_0 = 20;
const int XI_0 = 22;
const int ANTINEUTRON = 25;
const int ANTILAMBDA = 26;
const int ANTISIGMA_0 = 28;
const int ANITXI_0 = 30;

IsolationCut::IsolationCut(const char *filename) : _ievent(0),
                                                   _events_photon(0),
                                                   _hn_energy_cone( NULL ),
                                                   _file_output( NULL )
{
  // construct output file names
  _output_file_name = "histos/IsolationCut-";
  _output_file_name.append(filename);

  // initialize array for tower status
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
  _file_output = new TFile( _output_file_name.c_str(), "RECREATE" );

  int ndim_hn_photon = 5;
  int nbins_hn_photon[] =   {100 , 100 ,  10000 ,   11  ,  2  };
  double xmin_hn_photon[] = {  0.,   0.,      0.,  -0.05, -0.5};
  double xmax_hn_photon[] = {100., 100.,    100.,   1.05,  1.5};

  _hn_energy_cone = new THnSparseF("hn_energy_cone",
                                   "Energy in cone around Photon; E [GeV]; E_cone [GeV]; f_cone; r_cone [rad]; isPromptPhoton;",
                                   ndim_hn_photon,
                                   nbins_hn_photon,
                                   xmin_hn_photon,
                                   xmax_hn_photon );

  return EVENT_OK;
}

int IsolationCut::process_event(PHCompositeNode *topNode)
{
  _ievent++;

  /* create vector of PID's of neutral particles */
  vector< int > v_pid_neutral;
  v_pid_neutral.push_back( PHOTON_PID );
  v_pid_neutral.push_back( NEUTRINO_PID );
  v_pid_neutral.push_back( PIZERO_PID );
  v_pid_neutral.push_back( KAON_0_LONG );
  v_pid_neutral.push_back( NEUTRON );
  v_pid_neutral.push_back( KAON_0_SHORT );
  v_pid_neutral.push_back( ETA );
  v_pid_neutral.push_back( LAMBDA );
  v_pid_neutral.push_back( SIGMA_0 );
  v_pid_neutral.push_back( XI_0 );
  v_pid_neutral.push_back( ANTINEUTRON );
  v_pid_neutral.push_back( ANTILAMBDA );
  v_pid_neutral.push_back( ANTISIGMA_0 );
  v_pid_neutral.push_back( ANITXI_0 );

  // global info
  PHGlobal *data_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!data_global)
    {
      cout << "Cannot find PHGlobal" << endl;
      return DISCARDEVENT;
    }

  // track info
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
    {
      cout << "Cannot find emcGeaTrackContainer" << endl;
      return DISCARDEVENT;
    }

  // cluster info
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
    {
      cout << "Cannot find emcGeaClusterContainer" << endl;
      return DISCARDEVENT;
    }

  // associate cluster with track
  // map key is trkno
  typedef map<int,AnaTrk*> map_Ana_t;
  map_Ana_t track_list;

  // store photon clusters
  vector<emcGeaClusterContent*> cluster_photons;

  // store direct photons
  vector<AnaTrk*> photons_direct;

  // store photons from pi0 (and other hadrons) decay
  vector<AnaTrk*> photons_decay;

  // store other photons
  vector<AnaTrk*> photons_other;

  // number of tracks and clusters
  unsigned nemctrk = emctrkcont->size();
  unsigned nemcclus = emccluscont->size();

  // cuts for photon candidate selection
  float cluster_ecore_min = 0.3; //GeV
  float photon_ecore_min = 1; //GeV
  float photon_prob_min = 0.02;

  /* loop over all cluster in EMCAL - add to AnaTrk collection if it meets photon criteria */
  unsigned nphotons = 0;
  for( unsigned icluster=0; icluster < nemcclus; icluster++ )
    {
      emcGeaClusterContent *emc_cluster = emccluscont->get( icluster );

      /* is photon candidate? */
      if ( emc_cluster->ecore() < cluster_ecore_min ||
           emc_cluster->prob_photon() < photon_prob_min )
	continue;

      /* charged track veto? */
      //track matching parameters not set, see cout statement below to test
      //cout << "Track dz: " << emc_cluster->emctrkdz() << " , dphi: " << emc_cluster->emctrkdphi() << " , quality: " << emc_cluster->emctrkquality() << endl;
      /* use truth info here to identify photon */
      //cout << "Cluster pid: " << emc_cluster->pid() << endl;
      //if ( emc_cluster->pid() != PHOTON_PID )
      //continue;

      /* append to photon candiadte collection */
      cluster_photons.push_back( emc_cluster );
      nphotons++;
    }

  /* loop only over photon candidate clusters */
  for( unsigned icluster=0; icluster < cluster_photons.size(); icluster++ )
    {
      emcGeaClusterContent *emc_cluster = cluster_photons.at( icluster );

      /* is cluster energy below threshold for DIRECT photon candidate? */
      if ( emc_cluster->ecore() < photon_ecore_min )
	continue;

      /* check tower- is in guard ring of 10(12) towers ob PbSc(PbGl)? */
      // ...

      /* what's the origin of this cluster- direct photon, pi0 photon, something else? */
      bool isPromptPhoton = true;
      //if ( lvl0_trk->anclvl == 0 )
      //  isPromptPhoton = true;

      /* loop over different cone radii */
      for ( int round = 0; round < 10; round ++ )
        {
          double rcone = 0.1 + round * 0.1;

          /* add up energy in cone around cluster */
          double econe = 0;

          /* loop over all cluster in EMCAL */
	  for( unsigned icluster2=0; icluster < cluster_photons.size(); icluster++ )
	    {
              if ( icluster2 == icluster )
                continue;

	      emcGeaClusterContent *emc_cluster2 = cluster_photons.at( icluster2 );

              if ( emc_cluster2->ecore() < cluster_ecore_min )
                continue;

              float dr = ( sqrt( pow( emc_cluster->theta() - emc_cluster2->theta() , 2 ) +
                                 pow( emc_cluster->phi() - emc_cluster2->phi() , 2 ) ) );

              if ( dr < rcone )
                {
                  econe += emc_cluster2->ecore();
                }
            }

          /* loop over all tracks - ideally would use DC tracks but not available  */
	  float track_pmin = 0.2; //GeV
	  float track_pmax = 15; //GeV
	  for(unsigned itrk=0; itrk<nemctrk; itrk++)
	    {
	      emcGeaTrackContent *emctrk = emctrkcont->get(itrk);

	      /* Check if particle is charged */
	      if ( find( v_pid_neutral.begin(), v_pid_neutral.end(), emctrk->get_pid() ) != v_pid_neutral.end() )
		{
		  //cout << "impz(neutral): " << emctrk->get_impz() << " , pid = " << emctrk->get_pid() << endl;
		  continue;
		}
	      //cout << "impz(charged): " << emctrk->get_impz() << " , pid = " << emctrk->get_pid() << endl;

	      /* track NOT hitting EMCAL? (substitute to test if track is outside of drift chamber acceptance) */
	      if ( ( emctrk->get_impx() == 0 ) &&
		   ( emctrk->get_impy() == 0 ) &&
		   ( emctrk->get_impz() == 0 ) )
		continue;

	      /* track momentum below threshold? */
	      if ( ( emctrk->get_ptot() < track_pmin ) ||
		   ( emctrk->get_ptot() > track_pmax ) )
		continue;

	      /* determine track theta and phi directions */
	      float track_theta = atan2( emctrk->get_pt() , emctrk->get_ptot() );
	      float track_phi = atan2( emctrk->get_py() , emctrk->get_px() );

	      /* check if track within cone */
              float dr = ( sqrt( pow( track_theta - emc_cluster->theta() , 2 ) +
                                 pow( track_phi   - emc_cluster->phi()   , 2 ) ) );

              if ( dr < rcone )
                {
                  econe += emctrk->get_ptot();
		  //cout << "Add track with pid " << pdg_p->GetName() << " and momentum " << emctrk->get_ptot() << " to cone." << endl;
                }
	      else
		{
		  //cout << "Track outside of cone." << endl;
		}
	    }

	  /* determine energy fraction of cone */
          double econe_frac = econe / emc_cluster->ecore();
          //cout << "Energy in cone, fraction: " << econe << ", " << econe_frac << endl;

          /* fill histogram */
          double fill[] = {(double)emc_cluster->ecore(), econe, econe_frac, rcone, (double)isPromptPhoton};
          _hn_energy_cone->Fill( fill );

        } // END: loop over isolation cone radius
    } // END: loop over all cluster in EMCAL

  if ( nphotons > 0 )
    _events_photon++;

  /////// old code below
  //  for(unsigned itrk=0; itrk<nemctrk; itrk++)
  //    {
  //      emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
  //      AnaTrk *track = new AnaTrk(emctrk, emccluscont, (int*)_tower_status);
  //      if(track)
  //        track_list.insert( make_pair(track->trkno,track) );
  //    }
  //
  //  // analyze tracks
  //  map_Ana_t::iterator it_track_ref;
  //  map_Ana_t::iterator it_track_comp;
  //
  //  for ( it_track_ref = track_list.begin(); it_track_ref != track_list.end(); it_track_ref++ )
  //    {
  //      AnaTrk *lvl0_trk = it_track_ref->second;
  //
  //      /* track is photon */
  //      if( lvl0_trk->pid == PHOTON_PID )
  //        {
  //
  //          /* test if particle passes energy threshold */
  //          double photon_minEnergy = 1.0;
  //          if ( lvl0_trk->ecore < photon_minEnergy )
  //            continue;
  //
  //          /* test if particle is prompt photon */
  //          bool isPromptPhoton = false;
  //          if ( lvl0_trk->anclvl == 0 )
  //            isPromptPhoton = true;
  //          //      if ( parent )
  //          //        {
  //          //          if ( part->GetKF() == 22 &&
  //          //               parent->GetKF() == 22 &&
  //          //               parent->GetKS() == 21 &&
  //          //               !(parent->GetParent()) )
  //          //            {
  //          //              isPromptPhoton = true;
  //          //            }
  //          //        }
  //
  //          /* Convert particle into TLorentzVector */
  //          Float_t px = lvl0_trk->emctrk->get_px();
  //          Float_t py = lvl0_trk->emctrk->get_py();
  //          Float_t pz = lvl0_trk->emctrk->get_pz();
  //          Float_t energy = lvl0_trk->ecore;
  //          TLorentzVector v_part(px,py,pz,energy);
  //
  //          /* loop over different cone radii */
  //          for ( int round = 0; round < 10; round ++ )
  //            {
  //              double rcone = 0.1 + round * 0.1;
  //
  //              /* sum up all energy in cone around particle */
  //              double econe = 0;
  //              TVector3 v3_gamma(lvl0_trk->emctrk->get_px(), lvl0_trk->emctrk->get_py(), lvl0_trk->emctrk->get_pz());
  //
  //              for ( it_track_comp = track_list.begin(); it_track_comp != track_list.end(); it_track_comp++ )
  //                {
  //                  AnaTrk *lvl0_trk2 = it_track_comp->second;
  //
  //                  /* only consider stable particles, skip if pointer identical to 'reference' particle */
  //                  if ( lvl0_trk != lvl0_trk2 )
  //                    {
  //                      TVector3 v3_part2(lvl0_trk2->emctrk->get_px(), lvl0_trk2->emctrk->get_py(), lvl0_trk2->emctrk->get_pz());
  //
  //                      /* test if particle is in Central Arm acceptance */
  //                      //Float_t px = part2->get_px();
  //                      //Float_t py = part2->get_py();
  //                      //Float_t pz = part2->get_pz();
  //                      //Float_t energy = part2->GetEnergy();
  //                      //TLorentzVector v_part2(px,py,pz,energy);
  //
  //                      if ( v3_gamma.Angle( v3_part2 ) < rcone )
  //                      {
  //                        //cout << "ecore: " << lvl0_trk2->ecore << endl;
  //                        econe += lvl0_trk2->ecore;
  //                      }
  //                    }
  //                }
  //              double econe_frac = econe / lvl0_trk->ecore;
  //              //cout << "Energy in cone, fraction: " << econe << ", " << econe_frac << endl;
  //
  //              double fill[] = {lvl0_trk->ecore, econe, econe_frac, rcone, (double)isPromptPhoton};
  //
  //              _hn_energy_cone->Fill( fill );
  //            }
  //
  //          /* photon is direct photon */
  //          if ( lvl0_trk->anclvl == 0 )
  //            {
  //              photons_direct.push_back(lvl0_trk);
  //            }
  //          /* photon is from pi0 decay */
  //          if ( lvl0_trk->anclvl == 1 )
  //            {
  //              photons_decay.push_back(lvl0_trk);
  //            }
  //          /* photon is from other source */
  //          else
  //            {
  //              photons_other.push_back(lvl0_trk);
  //            }
  //        } // end: track is photon
  //    } // end: track lvl0 loop

  //cout << "Direct photons: " << photons_direct.size() << endl;
  //cout << "Decay  photons: " << photons_decay.size() << endl;
  //cout << "Other  photons: " << photons_other.size() << endl;

  return EVENT_OK;
}

int IsolationCut::End(PHCompositeNode *topNode)
{
  _file_output->cd();

  if ( _hn_energy_cone )
    _hn_energy_cone->Write();

  _file_output->Close();

  cout << "---------------------------------------------------" << endl;
  cout << "IsolationCut module:" << endl;
  cout << "Number of events processed:            " << _ievent << endl;
  cout << "Number of events with photon in EMCAL: " << _events_photon << endl;
  cout << "---------------------------------------------------" << endl;

  return EVENT_OK;
}
