#ifndef __ANAPHPYTHIAHISTOS_H__
#define __ANAPHPYTHIAHISTOS_H__

#include <string>
#include <SubsysReco.h>

class PHCompositeNode;
class PHPythiaContainer;
class Fun4AllHistoManager;
class TMCParticle;
class TVector3;

class TH1;
class THnSparse;

class AnaPHPythiaHistos: public SubsysReco
{
  public:
    AnaPHPythiaHistos(const std::string &name = "AnaPHPythiaHistos", const char *filename = "histo.root");
    virtual ~AnaPHPythiaHistos();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    /* Sum up truth energy within cone of fixed opening
     * angle rcone around particle with vector vref */
    double SumETruth(const TMCParticle* pref, double rcone);

    /* Fill two-particle correlation histogram */
    void FillCorrelation(const TMCParticle *pref, int type);

    /* Test if particle is in Central Arm acceptance */
    bool InAcceptance(const TVector3 &v3_part);

    void BookHistograms();

    std::string outFileName;

    PHPythiaContainer *phpythia;

    Fun4AllHistoManager *hm;
    THnSparse *hn_photon;
    THnSparse *hn_corr;
};

#endif	/* __ANAPHPYTHIAHISTOS_H__ */
