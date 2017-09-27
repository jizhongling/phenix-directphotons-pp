#ifndef __PHOTONNODEMB_H__
#define __PHOTONNODEMB_H__

#include <SubsysReco.h>

class EmcLocalRecalibrator;
class PhotonContainerMB;
class SpinPattern;
class SpinDBContent;
class emcClusterContent;
class PHCentralTrack;
class PHCompositeNode;

class PhotonNodeMB: public SubsysReco
{
  public:
    PhotonNodeMB(const std::string &name = "PHOTONNODEMB");
    virtual ~PhotonNodeMB();

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
    PhotonContainerMB *photoncont;
    SpinPattern *spinpattern;
};

#endif /* __PHOTONNODEMB_H__ */
