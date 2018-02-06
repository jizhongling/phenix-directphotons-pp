#ifndef __MISSINGRATIO_H__
#define __MISSINGRATIO_H__

#include <SubsysReco.h>
#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;

class TF1;
class TFile;
class TH2;
class THnSparse;

class MissingRatio: public SubsysReco
{
  public:
    MissingRatio(const std::string &name = "MissingRatio", const char *filename = "histo.root");
    virtual ~MissingRatio();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    void BookHistograms();
    void ReadTowerStatus(const std::string& filename);
    void ReadSashaWarnmap(const std::string& filename);

    // number of pT bins
    static const int npT = 30;

    // pT bins
    static double vpT[npT+1];

    // tower status for warnmap
    int tower_status[8][48][96];

    std::string outFileName;
    Fun4AllHistoManager *hm;

    // 2D histograms for west and east arms with different criterias
    THnSparse *hn_conversion_position;
    TH2 *h2_radius;
    TH2 *h2_angle;
    TH2 *h2_photon;
    TH2 *h2_noconv;
    TH2 *h2_vtxconv;
    TH2 *h2_eeinconv;
    TH2 *h2_eeoutconv;
    THnSparse *hn_merge;
    THnSparse *hn_photon;
    THnSparse *hn_pion;

    TF1 *cross;
};

#endif /* __MISSINGRATIO_H__ */
