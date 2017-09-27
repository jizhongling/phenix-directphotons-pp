#ifndef __PHOTONNODE_H__
#define __PHOTONNODE_H__

#include <SubsysReco.h>

class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;
class PhotonContainer;
class SpinPattern;
class SpinDBContent;
class emcClusterContent;
class PHCentralTrack;
class PHCompositeNode;

class PhotonNode: public SubsysReco
{
  public:
    PhotonNode(const std::string &name = "PHOTONNODE");
    virtual ~PhotonNode();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    void EMCRecalibSetup();
    bool DispCut(const emcClusterContent *emccluster);
    bool TestPhoton(const emcClusterContent *emccluster, float bbc_t0);
    float GetTrackConeEnergy(const PHCentralTrack *tracks, const emcClusterContent *cluster, double cone_angle);
    void UpdateSpinPattern(SpinDBContent &spin_cont);

    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;
    PhotonContainer *photoncont;
    SpinPattern *spinpattern;

    int runnumber;
    int fillnumber;
};

#endif /* __PHOTONNODE_H__ */
