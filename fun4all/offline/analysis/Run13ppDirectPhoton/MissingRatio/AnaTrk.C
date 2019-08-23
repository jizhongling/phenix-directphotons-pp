#include "AnaTrk.h"

#include <AnaToolsTowerID.h>
#include "AnaToolsTrack.h"
#include <EMCWarnmapChecker.h>

#include <emcGeaTrackContent.h>
#include <emcGeaClusterContainer.h>
#include <emcGeaClusterContent.h>

#include <boost/foreach.hpp>

AnaTrk::AnaTrk(emcGeaTrackContent *trk, emcGeaClusterContainer *cluscont):
  trkno(-9999), pid(-9999), anclvl(-9999), parent_trkno(-9999), parent_trk(NULL),
  decayed(false), trkpt(-9999.), trkedep(-9999.), trkrbirth(-9999.),
  cid(-9999), arm(-9999), sector(-9999), ecore(-9999), cluspt(-9999.), prob_photon(-9999.),
  emctrk(trk), emccluscont(cluscont), emcclus(NULL), emcwarnmap(NULL)
{
  daughter_list.clear();
  trkvp.SetXYZ(-9999., -9999., -9999.);
  trkposbirth.SetXYZ(-9999., -9999., -9999.);
  trkposemcal.SetXYZ(-9999., -9999., -9999.);
  if(emctrk)
  {
    trkno = emctrk->get_trkno();
    pid = emctrk->get_pid();
    anclvl = emctrk->get_anclvl();
    parent_trkno = emctrk->get_parent_trkno();
    daughter_list = emctrk->get_daughter_list();
    decayed = daughter_list.empty() ? false : true;
    trkvp.SetXYZ( emctrk->get_px(), emctrk->get_py(), emctrk->get_pz() );
    trkpt = emctrk->get_pt();
    trkedep = emctrk->get_edep();
    trkposbirth.SetXYZ( emctrk->get_x(), emctrk->get_y(), emctrk->get_z() );
    trkposemcal.SetXYZ( emctrk->get_impx(), emctrk->get_impy(), emctrk->get_impz() );
    trkrbirth = emctrk->get_radius();
    arm = anatools::GetTrackArm(emctrk);
  }

  FillCluster();
}

AnaTrk::~AnaTrk()
{
  if(emcwarnmap)
    delete emcwarnmap;
}

void AnaTrk::FillCluster()
{
  /* find the cluster which has highest energy deposit */
  FindCluster();

  /* initialize EMC warnmap checker */
  emcwarnmap = new EMCWarnmapChecker();

  /* fill cluster info by this cluster */
  if( emcclus && emcwarnmap &&
      emcwarnmap->IsGoodTower(emcclus) )
  {
    cid = emcclus->id();
    sector = anatools::CorrectClusterSector(arm, emcclus->sector());
    ecore = emcclus->ecore();
    TVector3 vx( emcclus->x(), emcclus->y(), emcclus->z() );
    cluspt = ecore * ( vx.Perp() / vx.Mag() );
    prob_photon = emcclus->prob_photon();
  }

  return;
}

void AnaTrk::FindCluster()
{
  emcclus = NULL;

  float edepMax = 0.;
  emc_clusterlist_t clus_list = emctrk->get_cluster_list();
  BOOST_FOREACH(const emc_clusterid_t &iclus, clus_list)
  {
    float edep = emctrk->get_edep_bycluster(iclus);
    if( edep > edepMax )
    {
      emcclus = emccluscont->find(iclus);  // use find(clusterID)
      edepMax = edep;
    }
  }

  return;
}
