#ifndef __ISOLATION_H__
#define __ISOLATION_H__

#include <SubsysReco.h>
#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;
class THnSparse;

class AnaTrk;
class PHCentralTrack;

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
    void ReadSashaWarnmap(const std::string& filename);

    // Sum up energy within cone of fixed opening
    // angle rcone around particle with vector anatrk
    double SumEEmcal(const AnaTrk *anatrk, double rcone);
    double SumPTrack(const AnaTrk *anatrk, const PHCentralTrack *tracks, double rcone);

    // tower status for warnmap
    int tower_status[8][48][96];

    std::string outFileName;

    Fun4AllHistoManager *hm;
    THnSparse *hn_photon;
};

#endif /* __ISOLATION_H__ */
