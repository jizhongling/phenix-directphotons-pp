#ifndef __PHOTONEFF_H__
#define __PHOTONEFF_H__

#include <SubsysReco.h>
#include <string>

class EMCWarnmapChecker;

class Fun4AllHistoManager;
class PHCompositeNode;
class emcClusterContent;

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

    std::string outFileName;

    EMCWarnmapChecker *emcwarnmap;

    Fun4AllHistoManager *hm;
    TH2 *h2_photon_eta_phi[3];
    THnSparse *hn_1photon;

    TF1 *cross;
};

#endif /* __PHOTONEFF_H__ */
