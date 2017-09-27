#include "MissingRatio.h"

#include "AnaTrk.h"
#include "AnaToolsCluster.h"

#include <PHGlobal.h>
//#include <PHCentralTrack.h>
//#include <PHSnglCentralTrack.h>
//#include <McEvalSingleList.h>

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

MissingRatio::MissingRatio(const char *filename) :
  hm(NULL),
  h2_radius(NULL),
  h2_angle(NULL),
  h2_missing(NULL),
  h2_nomissing(NULL),
  h2_missing_pion(NULL),
  h2_nomissing_pion(NULL),
  h2_incident(NULL),
  h2_measured(NULL),
  h2_incident_pion(NULL),
  h2_measured_pion(NULL),
  h2_conversion(NULL),
  h2_photon(NULL),
  h2_merge_photon(NULL),
  h2_merge_pion(NULL),
  hn_conversion_position(NULL),
  hn_photon(NULL),
  hn_pion(NULL)
{
  // construct output file names
  outFileName = "histos/MissingRatio-";
  outFileName.append(filename);

  // initialize array for tower status
  for(int isector=0; isector<8; isector++)
    for(int ibiny=0; ibiny<48; ibiny++)
      for(int ibinz=0; ibinz<96; ibinz++)
        tower_status[isector][ibiny][ibinz] = 0;

  // initialize pT bins
  double apT[npT] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };
  copy(apT, apT+npT, vpT);
  sort(vpT, vpT+npT);

  // initialize pT weight, scale for bin width, count tracks
  double awpT[npT] = {0,
    0, 0.307319, 0.186529, 0.137168, 0.102891, 0.0784738, 0.0574793, 0.0402374, 0.0271419, 0.0184652,
    0.0127012, 0.00886497, 0.00609304, 0.00424588, 0.00307189, 0.00218866, 0.00163727, 0.00120101, 0.000907576, 0.000483097,
    0.000187945, 8.33357e-05, 4.30282e-05, 2.11378e-05, 1.13421e-05, 6.3131e-06, 4.1332e-06, 1.9904e-06, 1.10784e-06, 0};
  copy(awpT, awpT+npT, wpT);

  cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);
}

MissingRatio::~MissingRatio()
{
}

int MissingRatio::Init(PHCompositeNode *topNode)
{
  // initialize histogram manager
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  // create and register histograms
  BookHistograms();

  // read warnmap
  //ReadTowerStatus("Warnmap_Run13pp510.txt");
  ReadSashaWarnmap("dead_eff_run13pp500gev.dat");

  return EVENT_OK;
}

int MissingRatio::process_event(PHCompositeNode *topNode)
{
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
  if( nemctrk<=0 || nemcclus<=0 ) return DISCARDEVENT;

  // associate cluster with track
  // map key is trkno
  typedef map<int,AnaTrk*> map_Ana_t;
  map_Ana_t track_list;

  // BBC ZVertex
  //float bbc_z = data_global->getBbcZVertex();
  //if( fabs(bbc_z) > 10. ) return DISCARDEVENT;

  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *track = new AnaTrk(emctrk, emccluscont, (int*)tower_status);
    if(track)
      track_list.insert( make_pair(track->trkno,track) );
  }

  // initialize weight of this event
  // it will depend on pion pT
  double pionpt = -9999.;
  double weight = 0.;

  // store photon
  // multimap key is criteria for missing ratio
  // criteria = 4*part+2*prob+warnmap
  multimap<int,AnaTrk*> photon;
  multimap<int,AnaTrk*> target;

  // positron and electron from decayed photon
  vector<AnaTrk*> positron;
  vector<AnaTrk*> electron;

  // analyze tracks
  BOOST_FOREACH( map_Ana_t::value_type &lvl0, track_list)
  {// lvl0 loop
    AnaTrk *lvl0_trk = lvl0.second;
    if( lvl0_trk->pid == PIZERO_PID && lvl0_trk->anclvl == 0 )
    {// pion
      pionpt = lvl0_trk->trkpt;
      //weight = Get_wpT( pionpt );
      weight = 1.;
      emc_tracklist_t pion_daughter = lvl0_trk->daughter_list;
      BOOST_FOREACH( const emc_trkno_t &lvl1_trkno, pion_daughter )
      {// lvl1 loop
        AnaTrk *lvl1_trk = track_list.find(lvl1_trkno)->second;
        lvl1_trk->parent_trkpt = lvl0_trk->trkpt;  // fill parent_trkpt
        if( lvl1_trk->pid == PHOTON_PID && lvl1_trk->anclvl == 1 )
        {// photon
          //float edep = lvl1_trk->trkedep;
          //lvl1_trk->FillCluster(edep);  // fill associated cluster info
          // count photons for photon conversion
          h2_photon->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm, weight);
          h2_photon->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm+2, weight);
          if( lvl1_trk->decayed == false && lvl1_trk->cid >= 0 )
          {// not decayed photon on EMCal
            // store photon info for missing ratio
            StorePhoton(&photon, &target, lvl1_trk);
          } // not decayed photon on EMCal
          else
          {// decayed photon
            emc_tracklist_t photon_daughter = lvl1_trk->daughter_list;
            BOOST_FOREACH( const emc_trkno_t &lvl2_trkno, photon_daughter )
            {// lvl2 loop
              AnaTrk *lvl2_trk = track_list.find(lvl2_trkno)->second;
              lvl2_trk->parent_trkpt = lvl1_trk->trkpt;  // fill parent_trkpt
              // store positron and electron from decayed photon
              if( lvl2_trk->pid == POSITRON_PID && lvl2_trk->anclvl == 2 )
                positron.push_back(lvl2_trk);
              else if( lvl2_trk->pid == ELECTRON_PID && lvl2_trk->anclvl == 2 )
                electron.push_back(lvl2_trk);
            } // end of lvl2 loop
          } // decayed photon
        } // photon
      } // end of lvl1 loop
      break; // only one pion
    } // pion
  } // end of lvl0 loop

  // analyze positron and electron from decayed photon
  //const double merge_angle = 0.03;
  const double merge_radius = 50.;
  vector<emc_trkno_t> used;
  BOOST_FOREACH( AnaTrk *pos, positron )
    BOOST_FOREACH( AnaTrk *ele, electron )
      if( find(used.begin(), used.end(), ele->trkno) == used.end() &&
          pos->parent_trkno == ele->parent_trkno )
      {
        //float edep = pos->trkedep + ele->trkedep;
        //pos->FillCluster(edep);  // fill associated cluster info
        //ele->FillCluster(edep);  // fill associated cluster info

        // get parent photon pT
        double photonpt = pos->parent_trkpt;

        // skip if no associated cluster
        //if( pos->cid < 0 || ele->cid < 0 ) continue;

        // birth position of photon conversion
        double fill_hn_conversion_position[] = { photonpt, pos->trkposbirth.X(), pos->trkposbirth.Y() };
        hn_conversion_position->Fill( fill_hn_conversion_position, weight );

        // radius of the birth position
        double radius = pos->trkrbirth;

        // angle(rad) between positron and electron track
        TVector3 posline = pos->trkposemcal - pos->trkposbirth;
        TVector3 eleline = ele->trkposemcal - ele->trkposbirth;
        double angle = posline.Angle( eleline );

        if( pos->arm == 0 || ele->arm == 0 )
        {// plus sign for west arm
          radius = radius;
          angle = angle;
        }
        else if( pos->arm == 1 || ele->arm == 1 )
        {// minus sign for the east arm
          radius = -radius;
          angle = -angle;
        }
        else
        {// pairs outside detector geometry
          radius = -9999.;
          angle = -9999.;
        }

        // fill histograms for radius and angle
        h2_radius->Fill(photonpt, radius, weight);
        h2_angle->Fill(photonpt, angle, weight);

        if( abs(radius) < merge_radius )
        {// photon conversion inside magnetic field
          h2_conversion->Fill(photonpt, (double)pos->arm, weight);
        }
        else
        {// photon conversion outside magnetic field
          // store photon info for missing ratio
          StorePhoton(&photon, &target, pos);
          h2_conversion->Fill(photonpt, (double)pos->arm+2, weight);
        }

        // already matched, do not consider this electron again
        used.push_back(ele->trkno);
        break;
      }

  // consider different criterias for missing ratio
  // criteria = 12*merge+4*part+2*prob+warnmap
  for(int icr=0; icr<24; icr++)
  {
    pair<multimap<int,AnaTrk*>::iterator, multimap<int,AnaTrk*>::iterator> cr_photon_it;
    pair<multimap<int,AnaTrk*>::iterator, multimap<int,AnaTrk*>::iterator> cr_target_it;
    cr_photon_it = photon.equal_range(icr%12);
    cr_target_it = target.equal_range(icr%12);
    int nphoton = distance(cr_photon_it.first, cr_photon_it.second);
    int ntarget = distance(cr_target_it.first, cr_target_it.second);
    int ntarget0 = ntarget;

    // check whether two photons from pi0 are merged to one cluster
    if( icr >= 12 )
    {
      multimap<int,AnaTrk*>::iterator ph_it = cr_photon_it.first;
      if( nphoton == 2 && ph_it->second->cid == (++ph_it)->second->cid )
      {
        cr_photon_it.second = ph_it;
        nphoton = 1;
      }
      multimap<int,AnaTrk*>::iterator tar_it = cr_target_it.first;
      if( ntarget == 2 && tar_it->second->cid == (++tar_it)->second->cid )
      {
        cr_target_it.second = tar_it;
        ntarget = 1;
        h2_merge_photon->Fill(tar_it->second->cluspt, icr, weight);  // merged
        h2_merge_pion->Fill(pionpt, icr, weight);  // merged
        h2_merge_pion->Fill(pionpt, 24+icr, 0.5*weight);  // total
      }
    }

    // fill histograms with target photons
    for( multimap<int,AnaTrk*>::iterator tar_it = cr_target_it.first; tar_it != cr_target_it.second; tar_it++ )
    {
      h2_incident->Fill(tar_it->second->trkpt, (double)icr, weight);
      h2_measured->Fill(tar_it->second->cluspt, (double)icr, weight);
      if( nphoton == 1 )
        h2_missing->Fill(tar_it->second->cluspt, (double)icr, weight);
      else if( nphoton == 2 )
        h2_nomissing->Fill(tar_it->second->cluspt, (double)icr, weight);
      if( icr==1 || icr==5 || icr==9 )
        FillCluster(tar_it->second->emcclus, pionpt, weight);
      if( ntarget0 == 2 )
      {
        h2_merge_photon->Fill(tar_it->second->cluspt, 24+icr, weight);  // total
        h2_merge_pion->Fill(pionpt, 24+icr, 0.5*weight);  // total
      }
    }

    if( nphoton == 1 )
    {
      h2_missing_pion->Fill(pionpt, (double)icr, weight);
    }
    else if( nphoton == 2 )
    {
      multimap<int,AnaTrk*>::iterator ph_it = cr_photon_it.first;
      AnaTrk *trk1 = ph_it->second;
      AnaTrk *trk2 = (++ph_it)->second;
      h2_incident_pion->Fill(pionpt, (double)icr, weight);
      double totpt = anatools::GetTot_pT(trk1->emcclus, trk2->emcclus);
      h2_measured_pion->Fill(totpt, (double)icr, weight);
      h2_nomissing_pion->Fill(pionpt, (double)icr, weight);
      if( icr==15 || icr==19 || icr==23 )
        FillTwoClusters(trk1->emcclus, trk2->emcclus, weight);
    }
  }

  // clear track list
  BOOST_FOREACH( map_Ana_t::value_type &trk, track_list )
    delete trk.second;

  return EVENT_OK;
}

int MissingRatio::End(PHCompositeNode *topNode)
{
  // write histogram output to ROOT file
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void MissingRatio::BookHistograms()
{
  // all new histograms will automatically activate Sumw2
  TH1::SetDefaultSumw2();

  // radius for birth position of electron-positron pair
  // positive for the west arm
  // negative for th east arm
  h2_radius= new TH2D("h2_radius", "Radius for birth position of the e^{+}e^{-} pair; p_{T} [GeV/c]; Radius [cm];", npT-1,vpT, 600,-300.,300.);

  // angle between electron-positron pair
  h2_angle= new TH2D("h2_angle", "Angle distribution of e^{+}e^{-} pair; p_{T} [GeV/c]; rad;", npT-1,vpT, 100,-0.05,0.05);

  // 2D histograms for missing ratio
  // criteria = 12*merge+4*part+2*prob+warnmap
  h2_missing = new TH2D("h2_missing", "EMCal missing count for photon; p_{T} [GeV/c]; criteria;", npT-1, vpT, 24, -0.5, 23.5);
  h2_nomissing = (TH2*)h2_missing->Clone("h2_nomissing");
  h2_nomissing->SetTitle("EMCal nomissing count for photon");
  h2_missing_pion = (TH2*)h2_missing->Clone("h2_missing_pion");
  h2_missing_pion->SetTitle("EMCal missing count for pion");
  h2_nomissing_pion = (TH2*)h2_missing->Clone("h2_nomissing_pion");
  h2_nomissing_pion->SetTitle("EMCal nomissing count for pion");

  // 2D histograms for smearing
  // criteria = 12*merge+4*part+2*prob+warnmap
  h2_incident = (TH2*)h2_missing->Clone("h2_incident");
  h2_incident->SetTitle("Incident photon p_{T}");
  h2_measured = (TH2*)h2_missing->Clone("h2_measured");
  h2_measured->SetTitle("Measured photon p_{T}");
  h2_incident_pion = (TH2*)h2_missing->Clone("h2_incident_pion");
  h2_incident_pion->SetTitle("Incident pion p_{T}");
  h2_measured_pion = (TH2*)h2_missing->Clone("h2_measured_pion");
  h2_measured_pion->SetTitle("Measured pion p_{T}");

  // 2D histograms for photon conversion
  // 0 - west arm conversion inside magnetic field
  // 1 - east arm conversion inside magnetic field
  // 2 - west arm conversion outside magnetic field
  // 3 - east arm conversion outside magnetic field
  h2_conversion = new TH2D("h2_conversion", "EMCal converstion count; p_{T} [GeV/c]; criteria;", npT-1, vpT, 4, -0.5, 3.5);
  h2_photon = new TH2D("h2_photon", "EMCal photon count; p_{T} [GeV/c]; criteria;", npT-1, vpT, 4, -0.5, 3.5);

  // 2D histograms for photon merging
  // category = criteria * (0:merged or 1:total)
  h2_merge_photon = new TH2D("h2_merge_photon", "EMCal photons merging; p_{T} [GeV/c]; category;", npT-1, vpT, 48, -0.5, 47.5);
  h2_merge_pion = (TH2*)h2_merge_photon->Clone("h2_merge_pion");

  int nbins_hn_conversion_position[] = {npT-1, 300, 300};
  double xmin_hn_conversion_position[] = {0., -300., -300.};
  double xmax_hn_conversion_position[] = {0., 300., 300.};
  hn_conversion_position = new THnSparseD("hn_conversion_position", "Photon conversion position;p^{#pi^{0}}_{T} [GeV];X [cm];Y [cm];",
      3, nbins_hn_conversion_position, xmin_hn_conversion_position, xmax_hn_conversion_position);
  hn_conversion_position->SetBinEdges(0, vpT);

  int nbins_hn_photon[] = {npT-1, npT-1, 8, 70, 157};
  double xmin_hn_photon[] = {0., 0., -0.5, -0.35, -1.57};
  double xmax_hn_photon[] = {0., 0., 7.5, 0.35, 4.71};
  hn_photon = new THnSparseD("hn_photon", "EMCal photon count;p_{T} [GeV/c];p^{#pi^{0}}_{T} [GeV/c];sector;#eta;#phi [rad];",
      5, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, vpT);
  hn_photon->SetBinEdges(1, vpT);

  int nbins_hn_pion[] = {npT-1, 700, 8};
  double xmin_hn_pion[] = {0., 0., -0.5};
  double xmax_hn_pion[] = {0., 0.7, 7.5};
  hn_pion = new THnSparseD("hn_pion", "EMCal pion count;p^{#pi^{0}}_{T} [GeV];m_{inv} [GeV];sector;",
      3, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(0, vpT);

  // register histograms with histogram manager
  hm->registerHisto(h2_radius);
  hm->registerHisto(h2_angle);
  hm->registerHisto(h2_missing);
  hm->registerHisto(h2_nomissing);
  hm->registerHisto(h2_missing_pion);
  hm->registerHisto(h2_nomissing_pion);
  hm->registerHisto(h2_incident);
  hm->registerHisto(h2_measured);
  hm->registerHisto(h2_incident_pion);
  hm->registerHisto(h2_measured_pion);
  hm->registerHisto(h2_conversion);
  hm->registerHisto(h2_photon);
  hm->registerHisto(h2_merge_photon);
  hm->registerHisto(h2_merge_pion);
  hm->registerHisto(hn_conversion_position);  hn_conversion_position->Sumw2();
  hm->registerHisto(hn_photon);  hn_photon->Sumw2();
  hm->registerHisto(hn_pion);  hn_pion->Sumw2();
}

void MissingRatio::ReadTowerStatus(const string& filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  unsigned int sector = 0;
  unsigned int biny = 0;
  unsigned int binz = 0;
  unsigned int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> biny >> binz >> status )
  {
    // count tower with bad status for PbSc and PbGl
    if ( status > 10 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status[sector][biny][binz] = status;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void MissingRatio::ReadSashaWarnmap(const string &filename)
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  unsigned int sector = 0;
  unsigned int biny = 0;
  unsigned int binz = 0;

  TOAD *toad_loader = new TOAD("MissingRatio");
  string file_location = toad_loader->location(filename);
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> binz >> biny )
  {
    // count tower with bad status for PbSc and PbGl
    if( sector < 6 ) nBadSc++;
    else nBadGl++;

    tower_status[sector][biny][binz] = 1;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

double MissingRatio::Get_wpT(double pT)
{
  if(pT > 1.)
    return cross->Eval(pT);
  else
    return 0.;

  int ipT = static_cast<int>( upper_bound(vpT, vpT+npT, pT) - vpT ); ipT--;
  return wpT[ipT];
}

// store photon info for one criteria
// criteria = 4*part+2*prob+warnmap
void MissingRatio::StorePhoton(multimap<int,AnaTrk*> *photon, multimap<int,AnaTrk*> *target, AnaTrk *trk, int part, int prob, int warnmap)
{
  if( trk->part != part ) return;

  // criteria
  int icr = 4*part + 2*prob + warnmap;

  if( trk->prob_photon > 0.02*prob )
  {// prob_photon > 0 or 0.02
    if(!warnmap)
    {
      if( trk->ecore > 0.3 )
        photon->insert( make_pair(icr,trk) );
      if( trk->ecore > 0.5 )
        target->insert( make_pair(icr,trk) );
    }
    else
    {
      if( trk->status == 0 && trk->ecore > 0.3 )
        photon->insert( make_pair(icr,trk) );
      if( trk->status == 0 && trk->ecore > 0.3 )
        target->insert( make_pair(icr,trk) );
    }
  }

  return;
}

// store photon info for all criterias
// criteria = 4*part+2*prob+warnmap
void MissingRatio::StorePhoton(multimap<int,AnaTrk*> *photon, multimap<int,AnaTrk*> *target, AnaTrk *trk)
{
  for(int part=0; part<3; part++)
    for(int prob=0; prob<2; prob++)
      for(int warnmap=0; warnmap<2; warnmap++)
        StorePhoton(photon, target, trk, part, prob, warnmap);

  return;
}

void MissingRatio::FillCluster(emcGeaClusterContent *emccluster, double pionpT, double weight)
{
  int arm = emccluster->arm();
  int sector = arm==0 ? emccluster->sector() : 7-emccluster->sector();

  TLorentzVector pE = anatools::Get_pE(emccluster);
  double px = pE.Px();
  double py = pE.Py();
  double pz = pE.Pz();
  double pT = anatools::Get_pT(emccluster);
  double mom = pE.P();

  double eta = mom > 0. ? atan(pz/mom) : 9999.;
  double phi = px > 0. ? atan(py/px) : 3.1416+atan(py/px);

  double fill_hn_photon[] = {pT, pionpT, sector, eta, phi};
  hn_photon->Fill(fill_hn_photon, weight);

  return;
}

void MissingRatio::FillTwoClusters(emcGeaClusterContent *emccluster1, emcGeaClusterContent *emccluster2, double weight)
{
  int arm = emccluster1->arm();
  int sector = arm==0 ? emccluster1->sector() : 7-emccluster1->sector();

  double tot_pT = anatools::GetTot_pT(emccluster1, emccluster2);
  double minv = anatools::GetInvMass(emccluster1, emccluster2);

  double fill_hn_pion[] = {tot_pT, minv, sector};
  hn_pion->Fill(fill_hn_pion, weight);

  return;
}
