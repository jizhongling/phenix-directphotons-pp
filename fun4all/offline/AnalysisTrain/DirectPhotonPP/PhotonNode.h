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
    PhotonNode(const std::string &name = "PhotonNode");
    virtual ~PhotonNode();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SelectMB();
    void SelectERT();

  protected:
    void EMCRecalibSetup();
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    bool TestPhoton(const emcClusterContent *emccluster, float bbc_t0);
    bool DispCut(const emcClusterContent *emccluster);

    int GetStatus(const emcClusterContent *emccluster);
    int GetStatusSasha(const emcClusterContent *emccluster);
    float GetTrackConeEnergy(const PHCentralTrack *tracks, const emcClusterContent *cluster, double cone_angle);

    void UpdateSpinPattern(SpinDBContent &spin_cont);

    enum DataType {MB, ERT};
    DataType datatype;

    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;
    PhotonContainer *photoncont;
    SpinPattern *spinpattern;

    int runnumber;
    int fillnumber;

    // tower status for warnmap
    int tower_status[8][48][96];
    int tower_status_sasha[8][48][96];
};

#endif /* __PHOTONNODE_H__ */
