#ifndef __ANAFASTMC_H__
#define __ANAFASTMC_H__

#include <SubsysReco.h>

#include <Rtypes.h>
#include <TLorentzVector.h>

#include <string>

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;

class Fun4AllHistoManager;

class TFile;
class TH1;
class TH2;
class THnSparse;
class TF1;

enum MCMethod {PHParticleGen, FastMC};
enum WarnMap {None, Nils, Sasha};

class AnaFastMC: public SubsysReco
{
  public:
    AnaFastMC(const std::string &name = "AnaFastMC");
    virtual ~AnaFastMC();

    // Methods Derived from SubsysReco
    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void set_outfile(std::string filename) { outFileName = filename; }
    void set_mcmethod(MCMethod method) { mcmethod = method; }
    void set_warnmap(WarnMap warn) { warnmap = warn; }

  protected:
    void BookHistograms();
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    bool pi0_sim(TLorentzVector *pG1, TLorentzVector *pG2, float& ptsim, float& mm, float& dist);
    bool photon_sim(TLorentzVector *pG1, float& ptsim);
    void ResetTowerEnergy();
    void FillTowerEnergy( int sec, int iy, int iz, float e );
    float GetETwr( int sec, int iy, int iz );
    int GetNpeak();
    bool GetImpactSectorTower(Double_t px, Double_t py, Double_t pz,  
        int& sec, int& iz, int& iy, float& zz, float& yy, 
        float& phi0, float& ximp, float& yimp, float& zimp );
    bool GetShower(Double_t px, Double_t py, Double_t pz, float& eout, int& itw );
    bool Gamma_En(Double_t px, Double_t py, Double_t pz, float& eout, int& itw,
        float& ximp, float& yimp, float& zimp);
    bool Gamma_Pos(Double_t& px, Double_t& py, Double_t& pz );

    std::string outFileName;
    MCMethod mcmethod;
    WarnMap warnmap;

    static const int MAXPEAK = 10;

    static const int NSEC = 8;
    static const int NY = 48;
    static const int NZ = 96;

    static const int nY_sc = 36;
    static const int nZ_sc = 72;
    static const int nY_gl = 48;
    static const int nZ_gl = 96;

    // tower status for warnmap
    int tower_status[NSEC][NY][NZ];

    int NPart;
    int NPeak;
    TLorentzVector Vpart[MAXPEAK];
    int itw_part[MAXPEAK];
    float eTwr[NSEC][NY][NZ];

    PHPythiaHeader *phpythiaheader;
    PHPythiaContainer *phpythia;

    Fun4AllHistoManager *hm;
    TH1 *h_pion;
    TH1 *h_photon;
    TH2 *h2_pion_eta_phi[3];
    TH2 *h2_photon_eta_phi[3];
    THnSparse* hn_missing;
    THnSparse* hn_pion;
    THnSparse *hn_photon;

    TF1 *cross_pi0;
    TF1 *cross_ph;
};

#endif	/* __ANAFASTMC_H__ */
