#ifndef __ANAPHPYTHIAHISTOS_H__
#define __ANAPHPYTHIAHISTOS_H__

#include <string>
#include <SubsysReco.h>

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;
class Fun4AllHistoManager;

class TDatabasePDG;
class TH1;
class THnSparse;

class AnaPHPythiaHistos: public SubsysReco
{
  public:
    AnaPHPythiaHistos(const std::string &name = "AnaPHPythiaHistos", const char *filename = "histo.root");
    virtual ~AnaPHPythiaHistos();

    // Methods Derived from SubsysReco
    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    // Sum up truth energy within cone of fixed opening
    // angle rcone around particle with vector vref
    double SumETruth(const TMCParticle* pref, double rcone);

    // Fill two-particle correlation histogram
    void FillCorrelation(const TMCParticle *pref, int type);

    void BookHistograms();

    // Some constants
    static const double PI = 3.1415927;

    // ROOT databsae with PDG properties
    TDatabasePDG* pdg_db;

    PHPythiaHeader *phpythiaheader;
    PHPythiaContainer *phpythia;

    std::string outFileName;
    Fun4AllHistoManager *hm;
    THnSparse *hn_photon;
    THnSparse *hn_corr;
};

#endif	/* __ANAPHPYTHIAHISTOS_H__ */
