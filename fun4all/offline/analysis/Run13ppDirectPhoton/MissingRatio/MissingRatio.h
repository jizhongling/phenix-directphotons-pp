#ifndef __MISSINGRATIO_H__
#define __MISSINGRATIO_H__

#include <SubsysReco.h>
#include <string>

class EMCWarnmapChecker;

class PHCompositeNode;
class Fun4AllHistoManager;

class TF1;
class TFile;
class TH2;
class THnSparse;

class MissingRatio: public SubsysReco
{
  public:
    MissingRatio(const std::string &name = "MissingRatio");
    virtual ~MissingRatio();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void set_outfile(std::string filename) { outFileName = filename; }

  protected:
    /* Number of pT bins */
    static const int npT = 30;

    /* PT bins */
    static double vpT[npT+1];

    void BookHistograms();

    std::string outFileName;

    /* EMC warnmap checker */
    EMCWarnmapChecker *emcwarnmap;

    /* 2D histograms for west and east arms with different criterias */
    Fun4AllHistoManager *hm;
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
