#ifndef __ANATRK_H__
#define __ANATRK_H__

/* associated track and cluster info */

#include <emctypes.h>
#include <TVector3.h>

class emcGeaTrackContent;
class emcGeaClusterContainer;
class emcGeaClusterContent;

class AnaTrk
{
  public:
    AnaTrk(emcGeaTrackContent *trk, emcGeaClusterContainer *cluscont, int *vstatus); 
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

    emcGeaTrackContent *emctrk;
    emcGeaClusterContainer *emccluscont;
    emcGeaClusterContent *emcclus;

  protected:
    // associate a cluster which has highest energy deposit
    void FillCluster();
    void FindCluster();

    int *vtower_status;
};

#endif /* __ANATRK__H__ */
