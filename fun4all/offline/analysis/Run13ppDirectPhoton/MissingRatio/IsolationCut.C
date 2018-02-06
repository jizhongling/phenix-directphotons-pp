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

#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <THnSparse.h>

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
const int PIZERO_PID = 7;


IsolationCut::IsolationCut(const char *filename)
{
  // construct output file names
  outFileName = "histos/IsolationCut-";
  outFileName.append(filename);
}

IsolationCut::~IsolationCut()
{
}

int IsolationCut::Init(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int IsolationCut::process_event(PHCompositeNode *topNode)
{
  cout << "New event! " << endl;

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

  // number of tracks and clusters
  int nemctrk = emctrkcont->size();
  int nemcclus = emccluscont->size();

  cout << "Found " << nemctrk << " tracks and " << nemcclus << " cluster in event!" << endl;

  if( nemctrk<=0 || nemcclus<=0 ) return DISCARDEVENT;

  /* Loop over all cluster */
  




  // associate cluster with track
  // map key is trkno
//  typedef map<int,AnaTrk*> map_Ana_t;
//  map_Ana_t track_list;

  // BBC ZVertex
  //float bbc_z = data_global->getBbcZVertex();
  //if( fabs(bbc_z) > 10. ) return DISCARDEVENT;

//  for(int itrk=0; itrk<nemctrk; itrk++)
//  {
//    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
//    AnaTrk *track = new AnaTrk(emctrk, emccluscont, (int*)tower_status);
//    if(track)
//      track_list.insert( make_pair(track->trkno,track) );
//  }
//
//  // initialize weight of this event
//  // it will depend on pion pT
//  double pionpt = -9999.;
//  double weight = 0.;
//
//  // store photon
//  // multimap key is criteria for missing ratio
//  // criteria = 4*part+2*prob+warnmap
//  multimap<int,AnaTrk*> photon;
//  multimap<int,AnaTrk*> target;
//
//  // positron and electron from decayed photon
//  vector<AnaTrk*> positron;
//  vector<AnaTrk*> electron;
//
//  // analyze tracks
//  BOOST_FOREACH( map_Ana_t::value_type &lvl0, track_list)
//  {// lvl0 loop
//    AnaTrk *lvl0_trk = lvl0.second;
//    if( lvl0_trk->pid == PIZERO_PID && lvl0_trk->anclvl == 0 )
//    {// pion
//      pionpt = lvl0_trk->trkpt;
//      weight = Get_wpT( pionpt );
//      emc_tracklist_t pion_daughter = lvl0_trk->daughter_list;
//      BOOST_FOREACH( const emc_trkno_t &lvl1_trkno, pion_daughter )
//      {// lvl1 loop
//        AnaTrk *lvl1_trk = track_list.find(lvl1_trkno)->second;
//        lvl1_trk->parent_trkpt = lvl0_trk->trkpt;  // fill parent_trkpt
//        if( lvl1_trk->pid == PHOTON_PID && lvl1_trk->anclvl == 1 )
//        {// photon
//          //float edep = lvl1_trk->trkedep;
//          //lvl1_trk->FillCluster(edep);  // fill associated cluster info
//          // count photons for photon conversion
//          h2_photon->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm, weight);
//          h2_photon->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm+2, weight);
//          if( lvl1_trk->decayed == false )
//          {// not decayed photon
//            h2_noconv->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm, weight);
//            if( lvl1_trk->cid >= 0 )
//              // store photon info for missing ratio
//              StorePhoton(&photon, &target, lvl1_trk);
//          }// not decayed photon
//          else
//          {// decayed photon
//            emc_tracklist_t photon_daughter = lvl1_trk->daughter_list;
//            bool is_fillconv = false;
//            BOOST_FOREACH( const emc_trkno_t &lvl2_trkno, photon_daughter )
//            {// lvl2 loop
//              AnaTrk *lvl2_trk = track_list.find(lvl2_trkno)->second;
//              lvl2_trk->parent_trkpt = lvl1_trk->trkpt;  // fill parent_trkpt
//              // radius of the birth position
//              double radius = lvl2_trk->trkrbirth;
//              if( !is_fillconv && abs(radius) < merge_radius )
//              {// photon conversion inside magnetic field
//                h2_vtxconv->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm, weight);
//                is_fillconv = true;
//              }
//              // store positron and electron from decayed photon
//              if( lvl2_trk->pid == POSITRON_PID && lvl2_trk->anclvl == 2 )
//                positron.push_back(lvl2_trk);
//              else if( lvl2_trk->pid == ELECTRON_PID && lvl2_trk->anclvl == 2 )
//                electron.push_back(lvl2_trk);
//            } // end of lvl2 loop
//          } // decayed photon
//        } // photon
//      } // end of lvl1 loop
//      break; // only one pion
//    } // pion
//  } // end of lvl0 loop
//
//  // analyze positron and electron from decayed photon
//  vector<emc_trkno_t> used;
//  BOOST_FOREACH( AnaTrk *pos, positron )
//    BOOST_FOREACH( AnaTrk *ele, electron )
//      if( find(used.begin(), used.end(), ele->trkno) == used.end() &&
//          pos->parent_trkno == ele->parent_trkno )
//      {
//        //float edep = pos->trkedep + ele->trkedep;
//        //pos->FillCluster(edep);  // fill associated cluster info
//        //ele->FillCluster(edep);  // fill associated cluster info
//
//        // get parent photon pT
//        double photonpt = pos->parent_trkpt;
//
//        // skip if no associated cluster
//        //if( pos->cid < 0 || ele->cid < 0 ) continue;
//
//        // birth position of photon conversion
//        double fill_hn_conversion_position[] = { photonpt, pos->trkposbirth.X(), pos->trkposbirth.Y() };
//        hn_conversion_position->Fill( fill_hn_conversion_position, weight );
//
//        // radius of the birth position
//        double radius = pos->trkrbirth;
//
//        // angle(rad) between positron and electron track
//        TVector3 posline = pos->trkposemcal - pos->trkposbirth;
//        TVector3 eleline = ele->trkposemcal - ele->trkposbirth;
//        double angle = posline.Angle( eleline );
//
//        if( pos->arm == 0 || ele->arm == 0 )
//        {// plus sign for west arm
//          radius = radius;
//          angle = angle;
//        }
//        else if( pos->arm == 1 || ele->arm == 1 )
//        {// minus sign for the east arm
//          radius = -radius;
//          angle = -angle;
//        }
//        else
//        {// pairs outside detector geometry
//          radius = -9999.;
//          angle = -9999.;
//        }
//
//        // fill histograms for radius and angle
//        h2_radius->Fill(photonpt, radius, weight);
//        h2_angle->Fill(photonpt, angle, weight);
//
//        if( abs(radius) < merge_radius )
//        {// photon conversion inside magnetic field
//          h2_conversion->Fill(photonpt, (double)pos->arm, weight);
//        }
//        else
//        {// photon conversion outside magnetic field
//          // store photon info for missing ratio
//          StorePhoton(&photon, &target, pos);
//          h2_conversion->Fill(photonpt, (double)pos->arm+2, weight);
//        }
//
//        // already matched, do not consider this electron again
//        used.push_back(ele->trkno);
//        break;
//      }
//
//  // consider different criterias for missing ratio
//  // criteria = 12*merge+4*part+2*prob+warnmap
//  for(int icr=0; icr<24; icr++)
//  {
//    pair<multimap<int,AnaTrk*>::iterator, multimap<int,AnaTrk*>::iterator> cr_photon_it;
//    pair<multimap<int,AnaTrk*>::iterator, multimap<int,AnaTrk*>::iterator> cr_target_it;
//    cr_photon_it = photon.equal_range(icr%12);
//    cr_target_it = target.equal_range(icr%12);
//    int nphoton = distance(cr_photon_it.first, cr_photon_it.second);
//    int ntarget = distance(cr_target_it.first, cr_target_it.second);
//    int ntarget0 = ntarget;
//
//    // check whether two photons from pi0 are merged to one cluster
//    if( icr >= 12 )
//    {
//      multimap<int,AnaTrk*>::iterator ph_it = cr_photon_it.first;
//      if( nphoton == 2 && ph_it->second->cid == (++ph_it)->second->cid )
//      {
//        cr_photon_it.second = ph_it;
//        nphoton = 1;
//      }
//      multimap<int,AnaTrk*>::iterator tar_it = cr_target_it.first;
//      if( ntarget == 2 && tar_it->second->cid == (++tar_it)->second->cid )
//      {
//        cr_target_it.second = tar_it;
//        ntarget = 1;
//        h2_merge_photon->Fill(tar_it->second->cluspt, icr, weight);  // merged
//        h2_merge_pion->Fill(pionpt, icr, weight);  // merged
//        h2_merge_pion->Fill(pionpt, 24+icr, 0.5*weight);  // total
//      }
//    }
//
//    // fill histograms with target photons
//    for( multimap<int,AnaTrk*>::iterator tar_it = cr_target_it.first; tar_it != cr_target_it.second; tar_it++ )
//    {
//      h2_incident->Fill(tar_it->second->trkpt, (double)icr, weight);
//      h2_measured->Fill(tar_it->second->cluspt, (double)icr, weight);
//      if( nphoton == 1 )
//        h2_missing->Fill(tar_it->second->cluspt, (double)icr, weight);
//      else if( nphoton == 2 )
//        h2_nomissing->Fill(tar_it->second->cluspt, (double)icr, weight);
//      if( icr==1 || icr==5 || icr==9 )
//        FillCluster(tar_it->second->emcclus, pionpt, weight);
//      if( ntarget0 == 2 )
//      {
//        h2_merge_photon->Fill(tar_it->second->cluspt, 24+icr, weight);  // total
//        h2_merge_pion->Fill(pionpt, 24+icr, 0.5*weight);  // total
//      }
//    }
//
//    if( nphoton == 1 )
//    {
//      h2_missing_pion->Fill(pionpt, (double)icr, weight);
//    }
//    else if( nphoton == 2 )
//    {
//      multimap<int,AnaTrk*>::iterator ph_it = cr_photon_it.first;
//      AnaTrk *trk1 = ph_it->second;
//      AnaTrk *trk2 = (++ph_it)->second;
//      h2_incident_pion->Fill(pionpt, (double)icr, weight);
//      double totpt = anatools::GetTot_pT(trk1->emcclus, trk2->emcclus);
//      h2_measured_pion->Fill(totpt, (double)icr, weight);
//      h2_nomissing_pion->Fill(pionpt, (double)icr, weight);
//      if( icr==15 || icr==19 || icr==23 )
//        FillTwoClusters(trk1->emcclus, trk2->emcclus, weight);
//    }
//  }
//
//  // clear track list
//  BOOST_FOREACH( map_Ana_t::value_type &trk, track_list )
//    delete trk.second;

  return EVENT_OK;
}

int IsolationCut::End(PHCompositeNode *topNode)
{
  return EVENT_OK;
}
