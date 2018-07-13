#ifndef __ANAFASTMC_H__
#define __ANAFASTMC_H__

#include <SubsysReco.h>

#include <Rtypes.h>
#include <TLorentzVector.h>

#include <string>

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;
class TMCParticle;

class Fun4AllHistoManager;
class TFile;
class TH1;
class TH2;
class TH3;
class THnSparse;
class TF1;

enum MCMethod {PHParticleGen, FastMC};

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

  protected:
    void PythiaInput(PHCompositeNode *topNode);
    void FastMCInput();

    void SumETruth(const TMCParticle *pref, double &econe_all, double &econe_acc);

    void BookHistograms();
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);
    void ReadSimWarnmap(const std::string &filename);

    void pi0_sim(const TLorentzVector beam_pi0);
    void photon_sim(const TLorentzVector beam_ph);
    void geom_sim(const TVector3 beam);
    void ResetClusters();
    void ResetTowerEnergy();
    void FillTowerEnergy( int sec, int iy, int iz, double e );
    double GetETwr( int sec, int iy, int iz );
    int GetNpeak();
    bool CheckWarnMap( int itower );
    bool InFiducial( int itower );
    bool IsGoodTower( int itower );
    bool IsHotTower( int itower );
    bool GetImpactSectorTower(Double_t px, Double_t py, Double_t pz,  
        int& sec, int& iz, int& iy, double& zz, double& yy, 
        double& phi0, double& ximp, double& yimp, double& zimp );
    bool GetShower(Double_t px, Double_t py, Double_t pz, double& eout, int& itw );
    bool Gamma_En(Double_t px, Double_t py, Double_t pz, double& eout, int& itw,
        double& ximp, double& yimp, double& zimp);
    bool Gamma_Pos(Double_t& px, Double_t& py, Double_t& pz );

    std::string outFileName;
    MCMethod mcmethod;

    static const int MAXPEAK = 2;

    static const int NSEC = 8;
    static const int NY = 48;
    static const int NZ = 96;

    static const int nY_sc = 36;
    static const int nZ_sc = 72;
    static const int nY_gl = 48;
    static const int nZ_gl = 96;

    static const int n_twrs = 24768;

    // Tower status for warnmap
    int tower_status_nils[NSEC][NY][NZ];
    int tower_status_sasha[NSEC][NY][NZ];
    int tower_status_sim[NSEC][NY][NZ];

    int NPart;
    int NPeak;
    TLorentzVector Vpart[MAXPEAK];
    int itw_part[MAXPEAK];
    int sec_part[MAXPEAK];
    double eTwr[NSEC][NY][NZ];

    PHPythiaContainer *phpythia;
    PHPythiaContainer *phpythia_bg;

    Fun4AllHistoManager *hm;
    TH1 *h_pion;
    TH1 *h_photon;
    TH2 *h2_pion_eta_phi[3];
    TH2 *h2_photon_eta_phi[3];
    TH3 *h3_isopi0;
    TH3 *h3_isoeta;
    THnSparse* hn_pion;
    THnSparse* hn_missing;
    THnSparse *hn_photon;
    THnSparse *hn_isoprompt;

    TF1 *cross_pi0;
    TF1 *cross_ph;
};

#endif	/* __ANAFASTMC_H__ */
