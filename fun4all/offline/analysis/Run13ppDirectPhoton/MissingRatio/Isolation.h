#ifndef __ISOLATION_H__
#define __ISOLATION_H__

#include <SubsysReco.h>
#include <string>

class AnaTrk;

class Fun4AllHistoManager;
class PHCompositeNode;
class PHCentralTrack;

class THnSparse;

class Isolation: public SubsysReco
{
  public:
    Isolation(const std::string &name = "Isolation", const char *filename = "histo.root");
    virtual ~Isolation();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    void BookHistograms();

    /* Sum up energy within cone of fixed opening
     * Angle rcone around particle with vector anatrk */
    double SumEEmcal(const AnaTrk *anatrk, double rcone);
    double SumPTrack(const AnaTrk *anatrk, const PHCentralTrack *tracks, double rcone);

    std::string outFileName;

    Fun4AllHistoManager *hm;
    THnSparse *hn_photon;
};

#endif /* __ISOLATION_H__ */
