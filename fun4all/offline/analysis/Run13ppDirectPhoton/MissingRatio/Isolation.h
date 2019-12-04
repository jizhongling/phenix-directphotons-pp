#ifndef __ISOLATION_H__
#define __ISOLATION_H__

#include <SubsysReco.h>
#include <string>

class AnaTrk;
class EMCWarnmapChecker;

class Fun4AllHistoManager;
class PHCompositeNode;
class PHCentralTrack;

class TH1;
class THnSparse;

class Isolation: public SubsysReco
{
  public:
    Isolation(const std::string &name = "Isolation");
    virtual ~Isolation();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void set_outfile(std::string filename) { outFileName = filename; }

  protected:
    void BookHistograms();

    /* Sum up energy within cone of fixed opening
     * Angle rcone around particle with vector anatrk */
    double SumEEmcal(const AnaTrk *anatrk, double rcone);
    double SumPTrack(const AnaTrk *anatrk, const PHCentralTrack *tracks, double rcone);

    std::string outFileName;

    /* EMC warnmap checker */
    EMCWarnmapChecker *emcwarnmap;

    Fun4AllHistoManager *hm;
    TH1 *h_events;
    THnSparse *hn_photon;
};

#endif /* __ISOLATION_H__ */
