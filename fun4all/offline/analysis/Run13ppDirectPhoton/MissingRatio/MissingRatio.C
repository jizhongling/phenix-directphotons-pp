#include "MissingRatio.h"

#include <AnaToolsTowerID.h>
#include <AnaToolsCluster.h>
#include <EMCWarnmapChecker.h>
#include "AnaTrk.h"

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

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/foreach.hpp>

using namespace std;

double MissingRatio::vpT[] = { 0.0,
  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
  5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

/* Global constants */
const int PHOTON_PID = 1;
const int POSITRON_PID = 2;
const int ELECTRON_PID = 3;
const int PIZERO_PID = 7;

const double eMin = 0.3;
const double AsymCut = 0.8;
const double merge_angle = 0.03;
const double merge_radius = 50.;

MissingRatio::MissingRatio(const string &name):
  SubsysReco(name),
  emcwarnmap(nullptr),
  hm(nullptr),
  h_events(nullptr),
  hn_conversion_position(nullptr),
  h2_radius(nullptr),
  h2_angle(nullptr),
  h2_photon(nullptr),
  h2_noconv(nullptr),
  h2_vtxconv(nullptr),
  h2_eeinconv(nullptr),
  h2_eeoutconv(nullptr),
  hn_merge(nullptr),
  hn_photon(nullptr),
  hn_pion(nullptr)
{
  /* Construct output file names */
  cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);
}

MissingRatio::~MissingRatio()
{
}

int MissingRatio::Init(PHCompositeNode *topNode)
{
  /* Create and register histograms */
  BookHistograms();

  /* Initialize EMC warnmap checker */
  emcwarnmap = new EMCWarnmapChecker();
  if(!emcwarnmap)
  {
    cerr << "No emcwarnmap" << endl;
    exit(1);
  }

  return EVENT_OK;
}

int MissingRatio::process_event(PHCompositeNode *topNode)
{
  /* Count events */
  h_events->Fill(1.);

  /* EMC track info */
  emcGeaTrackContainer *emctrkcont = emcNodeHelper::getObject<emcGeaTrackContainer>("emcGeaTrackContainer", topNode);
  if(!emctrkcont)
  {
    cout << "Cannot find emcGeaTrackContainer" << endl;
    return DISCARDEVENT;
  }

  /* EMC cluster info */
  emcGeaClusterContainer *emccluscont = emctrkcont->GetClusters();
  if(!emccluscont)
  {
    cout << "Cannot find emcGeaClusterContainer" << endl;
    return DISCARDEVENT;
  }

  /* Initialize weight of this event
   * It will depend on pion pT */
  double pionpt = 0.;
  double weight = 1.;

  /* Associate cluster with track
   * Map key is trkno */
  typedef map<int,AnaTrk*> map_Ana_t;
  map_Ana_t track_list;

  /* Store photon */
  vector<AnaTrk*> photon;

  /* Positron and electron from decayed photon */
  vector<AnaTrk*> positron;
  vector<AnaTrk*> electron;

  /* Number of tracks and clusters */
  int nemctrk = emctrkcont->size();
  int nemcclus = emccluscont->size();
  if( nemctrk<=0 || nemcclus<=0 ) return DISCARDEVENT;

  for(int itrk=0; itrk<nemctrk; itrk++)
  {
    emcGeaTrackContent *emctrk = emctrkcont->get(itrk);
    AnaTrk *track = new AnaTrk(emctrk, emccluscont);
    if(track)
      track_list.insert( make_pair(track->trkno,track) );
  }

  /* Analyze tracks */
  BOOST_FOREACH( map_Ana_t::value_type &lvl0, track_list )
  {// lvl0 loop
    AnaTrk *lvl0_trk = lvl0.second;
    if( lvl0_trk->pid == PIZERO_PID && lvl0_trk->anclvl == 0 )
    {// pion
      pionpt = lvl0_trk->trkpt;
      if( pionpt > 1. )
        weight = cross->Eval(pionpt);
      else
        weight = cross->Eval(1.);

      emc_tracklist_t pion_daughter = lvl0_trk->daughter_list;
      BOOST_FOREACH( const emc_trkno_t &lvl1_trkno, pion_daughter )
      {// lvl1 loop
        AnaTrk *lvl1_trk = track_list.find(lvl1_trkno)->second;
        lvl1_trk->parent_trk = lvl0_trk;  // fill parent_trkpt
        if( lvl1_trk->pid == PHOTON_PID && lvl1_trk->anclvl == 1 )
        {// photon
          /* Count photons for photon conversion */
          h2_photon->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm);
          if( lvl1_trk->decayed == false )
          {// not decayed photon
            h2_noconv->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm);
            if( lvl1_trk->cid >= 0 && emcwarnmap->IsGoodTower(lvl1_trk->emcclus) )
              /* Store photon info for missing ratio */
              photon.push_back(lvl1_trk);
          } // not decayed photon

          else
          {// decayed photon
            emc_tracklist_t photon_daughter = lvl1_trk->daughter_list;
            bool is_fillconv = false;
            BOOST_FOREACH( const emc_trkno_t &lvl2_trkno, photon_daughter )
            {// lvl2 loop
              AnaTrk *lvl2_trk = track_list.find(lvl2_trkno)->second;
              lvl2_trk->parent_trk = lvl1_trk;  // fill parent_trk
              /* Radius of the birth position */
              double radius = lvl2_trk->trkrbirth;
              if( !is_fillconv && fabs(radius) < merge_radius )
              {// photon conversion inside magnetic field
                h2_vtxconv->Fill(lvl1_trk->trkpt, (double)lvl1_trk->arm);
                is_fillconv = true;
              }
              /* Store positron and electron from decayed photon */
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

  /* Analyze positron and electron from decayed photon */
  vector<emc_trkno_t> used;
  BOOST_FOREACH( AnaTrk *pos, positron )
    BOOST_FOREACH( AnaTrk *ele, electron )
    if( pos->parent_trkno == ele->parent_trkno &&
        find(used.begin(), used.end(), pos->trkno) == used.end() )
    {
      /* Consider one parent photon only once */
      used.push_back( pos->parent_trkno );

      /* Get parent photon pT */
      double photonpt = pos->parent_trk->trkpt;

      /* Skip if no associated cluster */
      //if(pos->cid < 0 || ele->cid < 0)
      //  continue;

      /* Birth position of photon conversion */
      double fill_hn_conversion_position[] = { photonpt, pos->trkposbirth.X(), pos->trkposbirth.Y() };
      hn_conversion_position->Fill( fill_hn_conversion_position );

      /* Radius of the birth position */
      double radius = pos->trkrbirth;

      /* Angle(rad) between positron and electron track */
      TVector3 posline = pos->trkposemcal - pos->trkposbirth;
      TVector3 eleline = ele->trkposemcal - ele->trkposbirth;
      double angle = posline.Angle( eleline );

      if( pos->arm == 0 && ele->arm == 0 )
      {// plus sign for west arm
      }
      else if( pos->arm == 1 && ele->arm == 1 )
      {// minus sign for the east arm
        radius = -radius;
        angle = -angle;
      }
      else
      {// pairs outside detector geometry
        radius = -9999.;
        angle = -9999.;
      }

      /* Fill histograms for radius and angle */
      h2_radius->Fill(photonpt, radius);
      h2_angle->Fill(photonpt, angle);

      if( fabs(radius) < merge_radius )
      {// photon conversion inside magnetic field
        h2_eeinconv->Fill(photonpt, (double)pos->arm);
      }
      else
      {// photon conversion outside magnetic field
        /* Store photon info for missing ratio */
        photon.push_back(pos->parent_trk);
        h2_eeoutconv->Fill(photonpt, (double)pos->arm);
      }

      /* Already matched, no need to consider other electrons */
      break;
    }

  /* BEGIN photon merging and missing ratio evaluation */

  /* Get number of photons and number of separated peaks 
   * nphoton=2: No missing partner for pi0
   * nphoton=1: Missing partner fro pi0
   * npeak=2: No merging for two photons
   * npeak=1: Merging for two photons */
  int nphoton = photon.size();
  int npeak = nphoton;

  /* Two photons merged */
  if( nphoton == 2 &&
      photon.at(0)->cid >= 0 &&
      photon.at(0)->cid == photon.at(1)->cid &&
      emcwarnmap->IsGoodTower(photon.at(0)->emcclus) )
  {
    npeak = 1;
    int sector = photon.at(0)->sector;
    int passed = photon.at(0)->prob_photon > 0.02 ? 1 : 0;
    double eta = fabs( photon.at(0)->trkvp.Eta() );
    double fill_hn_merge[] = {pionpt, (double)sector, (double)passed, eta};
    hn_merge->Fill(fill_hn_merge);
  }

  /* Loop over all clusters in calorimeter */
  BOOST_FOREACH( const AnaTrk *trk, photon )
  {
    /* Fill photon histogram */
    double pt_truth = trk->trkpt;
    double pt_reco = trk->cluspt;
    int sector = trk->sector;
    double fill_hn_photon[] = {pt_truth, pt_reco, (double)sector, (double)nphoton};
    hn_photon->Fill(fill_hn_photon, weight);
  }

  if( nphoton == 2 )
  {
    /* Parameters for two clusters from pi0 decay photons */
    int sec1 = photon.at(0)->sector;
    int sec2 = photon.at(1)->sector;
    double e1 = photon.at(0)->ecore;
    double e2 = photon.at(1)->ecore;

    /* Check in the same detector part
     * and pass the energy and asymmetry cuts */
    if( anatools::SectorCheck(sec1,sec2) &&
        e1 > eMin && e2 > eMin && fabs(e1-e2)/(e1+e2) < AsymCut )
    {
      /* Fill pi0 histogram */
      double ptsim = anatools::GetTot_pT( photon.at(0)->emcclus, photon.at(1)->emcclus );
      double minv = anatools::GetInvMass( photon.at(0)->emcclus, photon.at(1)->emcclus );
      double fill_hn_pion[] = {pionpt, ptsim, minv, (double)sec1, (double)npeak};
      hn_pion->Fill(fill_hn_pion, weight);
    }
  }

  /* Clear track list */
  BOOST_FOREACH( map_Ana_t::value_type &trk, track_list )
    delete trk.second;

  return EVENT_OK;
}

int MissingRatio::End(PHCompositeNode *topNode)
{
  /* Write histogram output to ROOT file */
  hm->dumpHistos();
  delete hm;

  return EVENT_OK;
}

void MissingRatio::BookHistograms()
{
  /* Initialize histogram manager */
  hm = new Fun4AllHistoManager("HistoManager");
  hm->setOutfileName(outFileName);

  /* Count events */
  h_events = new TH1F("h_events", "Events count", 1, 0.5, 1.5);
  hm->registerHisto(h_events);

  /* Photon conversion XY position */
  int nbins_hn_conversion_position[] = {npT, 300, 300};
  double xmin_hn_conversion_position[] = {0., -300., -300.};
  double xmax_hn_conversion_position[] = {0., 300., 300.};
  hn_conversion_position = new THnSparseF("hn_conversion_position", "Photon conversion position;p^{#pi^{0}}_{T} [GeV];X [cm];Y [cm];",
      3, nbins_hn_conversion_position, xmin_hn_conversion_position, xmax_hn_conversion_position);
  hn_conversion_position->SetBinEdges(0, vpT);
  hm->registerHisto(hn_conversion_position);

  /* Radius for birth position of electron-positron pair
   * Positive for the west arm
   * Negative for th east arm */
  h2_radius= new TH2F("h2_radius", "Radius for birth position of the e^{+}e^{-} pair;p_{T} [GeV/c];Radius [cm];", npT,vpT, 600,-300.,300.);
  hm->registerHisto(h2_radius);

  /* Angle between electron-positron pair */
  h2_angle= new TH2F("h2_angle", "Angle distribution of e^{+}e^{-} pair;p_{T} [GeV/c];rad;", npT,vpT, 100,-0.05,0.05);
  hm->registerHisto(h2_angle);

  /* 2D histograms for photon conversion */
  h2_photon = new TH2F("h2_photon", "EMCal photon count;p_{T} [GeV/c];arm;", npT,vpT, 2,-0.5,1.5);
  h2_noconv = new TH2F("h2_noconv", "No converted photon count;p_{T} [GeV/c];arm;", npT,vpT, 2,-0.5,1.5);
  h2_vtxconv = new TH2F("h2_vtxconv", "VTX Converted photon count;p_{T} [GeV/c];arm;", npT,vpT, 2,-0.5,1.5);
  h2_eeinconv = new TH2F("h2_eeinconv", "EMCal converstion inside magnetic field count;p_{T} [GeV/c];arm;", npT,vpT, 2,-0.5,1.5);
  h2_eeoutconv = new TH2F("h2_eeoutconv", "EMCal converstion outside magnetic field count;p_{T} [GeV/c];arm;", npT,vpT, 2,-0.5,1.5);
  hm->registerHisto(h2_photon);
  hm->registerHisto(h2_noconv);
  hm->registerHisto(h2_vtxconv);
  hm->registerHisto(h2_eeinconv);
  hm->registerHisto(h2_eeoutconv);

  int nbins_hn_merge[] = {60, 8, 2, 7};
  double xmin_hn_merge[] = {0., -0.5, -0.5, 0.};
  double xmax_hn_merge[] = {30., 7.5, 1.5, 0.35};
  hn_merge = new THnSparseF("hn_merge", "Merged photon count;p_{T} truth [GeV];sector;passed;eta;",
      4, nbins_hn_merge, xmin_hn_merge, xmax_hn_merge);
  hn_merge->SetBinEdges(0, vpT);
  hm->registerHisto(hn_merge);

  int nbins_hn_photon[] = {npT, npT, 8, 4};
  double xmin_hn_photon[] = {0., 0., -0.5, -0.5};
  double xmax_hn_photon[] = {0., 0., 7.5, 3.5};
  hn_photon = new THnSparseF("hn_photon", "EMCal photon count;p_{T} truth [GeV];p_{T} reco [GeV];sector;nphoton;",
      4, nbins_hn_photon, xmin_hn_photon, xmax_hn_photon);
  hn_photon->SetBinEdges(0, vpT);
  hn_photon->SetBinEdges(1, vpT);
  hn_photon->Sumw2();
  hm->registerHisto(hn_photon);

  int nbins_hn_pion[] = {npT, npT, 300, 8, 4};
  double xmin_hn_pion[] = {0., 0., 0., -0.5, -0.5};
  double xmax_hn_pion[] = {0., 0., 0.3, 7.5, 3.5};
  hn_pion = new THnSparseF("hn_pion", "EMCal pion count;Truth p_{T} truth [GeV];p_{T} reco [GeV];m_{inv} [GeV];sector;npeak;",
      5, nbins_hn_pion, xmin_hn_pion, xmax_hn_pion);
  hn_pion->SetBinEdges(0, vpT);
  hn_pion->SetBinEdges(1, vpT);
  hn_pion->Sumw2();
  hm->registerHisto(hn_pion);

  return;
}
