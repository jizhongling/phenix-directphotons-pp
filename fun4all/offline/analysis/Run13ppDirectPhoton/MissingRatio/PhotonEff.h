#ifndef __PHOTONEFF_H__
#define __PHOTONEFF_H__

#include <SubsysReco.h>
#include <string>

class emcClusterContent;
class PHCompositeNode;
class Fun4AllHistoManager;

class TF1;
class TH2;
class THnSparse;

class PhotonEff: public SubsysReco
{
  public:
    PhotonEff(const std::string &name = "PhotonEff", const char *filename = "histo.root");
    virtual ~PhotonEff();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    void BookHistograms();
    void ReadTowerStatus(const std::string& filename);
    void ReadSashaWarnmap(const std::string& filename);
    int GetStatus(const emcClusterContent *emccluster);

    // tower status for warnmap
    int tower_status[8][48][96];

    std::string outFileName;

    Fun4AllHistoManager *hm;
    TH2 *h2_photon_eta_phi[3];
    THnSparse *hn_1photon;

    TF1 *cross;
};

#endif /* __PHOTONEFF_H__ */
