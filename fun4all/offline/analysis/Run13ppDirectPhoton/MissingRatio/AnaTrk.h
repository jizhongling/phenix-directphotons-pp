#ifndef __ANATRK_H__
#define __ANATRK_H__

/* associated track and cluster info */

#include <TVector3.h>
#include <emctypes.h>

class emcGeaTrackContent;
class emcGeaClusterContainer;
class emcGeaClusterContent;

class AnaTrk
{
  public:
    AnaTrk(emcGeaTrackContent *trk, emcGeaClusterContainer *cluscont, int *vstatus); 
    virtual ~AnaTrk();

    // associate a cluster which has closest energy to edep
    //void FillCluster(float edep);

    int trkno;
    int pid;
    int anclvl;
    int parent_trkno;
    emc_tracklist_t daughter_list;
    bool decayed;
    TVector3 trkvp;
    float trkpt;
    float parent_trkpt;
    float trkedep;
    TVector3 trkposbirth;
    TVector3 trkposemcal;
    float trkrbirth;

    int cid;
    int arm;
    int part;
    int sector;
    int status;
    float ecore;
    float cluspt;
    float prob_photon;

    emcGeaTrackContent *emctrk;
    emcGeaClusterContainer *emccluscont;
    emcGeaClusterContent *emcclus;

  protected:
    // associate a cluster which has highest energy deposit
    void FillCluster();
    void FindCluster();
    void FindCluster(float edep);

    int *vtower_status;
};

#endif /* __ANATRK__H__ */
