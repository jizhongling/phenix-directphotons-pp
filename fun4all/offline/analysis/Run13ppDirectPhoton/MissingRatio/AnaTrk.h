#ifndef __ANATRK_H__
#define __ANATRK_H__

/* associated track and cluster info */

#include <emctypes.h>
#include <TVector3.h>

class emcGeaTrackContent;
class emcGeaClusterContainer;
class emcGeaClusterContent;
class EMCWarnmapChecker;

class AnaTrk
{
  public:
    AnaTrk(emcGeaTrackContent *trk, emcGeaClusterContainer *cluscont); 
    virtual ~AnaTrk();

    int trkno;
    int pid;
    int anclvl;
    int parent_trkno;
    AnaTrk *parent_trk;
    emc_tracklist_t daughter_list;
    bool decayed;
    TVector3 trkvp;
    float trkpt;
    float trkedep;
    TVector3 trkposbirth;
    TVector3 trkposemcal;
    float trkrbirth;

    int cid;
    int arm;
    int sector;
    float ecore;
    float cluspt;
    float prob_photon;

    emcGeaTrackContent *emctrk;
    emcGeaClusterContainer *emccluscont;
    emcGeaClusterContent *emcclus;

  protected:
    /* associate a cluster which has highest energy deposit */
    void FillCluster();
    void FindCluster();

    EMCWarnmapChecker *emcwarnmap;
};

#endif /* __ANATRK__H__ */
