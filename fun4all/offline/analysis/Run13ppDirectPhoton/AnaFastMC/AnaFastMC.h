#ifndef __ANAFASTMC_H__
#define __ANAFASTMC_H__

#include <SubsysReco.h>

#include <Rtypes.h>
#include <TLorentzVector.h>

#include <string>

class EMCWarnmapChecker;
class DCDeadmapChecker;

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;
class TMCParticle;
class TDatabasePDG;

class Fun4AllHistoManager;
class TFile;
class TH1;
class TH2;
class TH3;
class THnSparse;
class TF1;

enum MCMethod {FastMC, PHParticleGen};

class AnaFastMC: public SubsysReco
{
  public:
    AnaFastMC(const std::string &name = "AnaFastMC");
    virtual ~AnaFastMC();

    // Methods Derived from SubsysReco
    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void InitBatch(int thread, int scale);
    void set_outfile(std::string filename) { outFileName = filename; }
    void set_mcmethod(MCMethod method) { mcmethod = method; }

  protected:
    static const int MAXPEAK = 2;

    static const int NSEC = 8;
    static const int NY = 48;
    static const int NZ = 96;

    static const int nY_sc = 36;
    static const int nZ_sc = 72;
    static const int nY_gl = 48;
    static const int nZ_gl = 96;

    static const int n_twrs = 24768;

    static const int nh_eta_phi = 3*3;

    void FastMCInput();
    void PythiaInput(PHCompositeNode *topNode);
    void SumETruth(const TMCParticle *pref, bool prefInAcc,
        double &econe_all, double &econe_emc, double econe_trk[]);

    void BookHistograms();
    void ReadSimWarnmap();

    void pi0_sim(const TLorentzVector &beam_pi0);
    void photon_sim(const TLorentzVector &beam_ph);
    void geom_sim(const TVector3 &beam);
    void ResetClusters();
    void ResetTowerEnergy();
    void FillTowerEnergy( int sec, int iy, int iz, double e );
    double GetETwr( int sec, int iy, int iz );
    int GetNpeak();
    bool CheckWarnMap( int itower );
    bool InDCAcceptance( const TVector3 &v3_part, int charge );
    double GetEMCResponse(double mom);
    bool GetImpactSectorTower(double px, double py, double pz,  
        int& sec, int& iz, int& iy, double& zz, double& yy, 
        double& phi0, double& ximp, double& yimp, double& zimp);
    bool GetShower(double px, double py, double pz, double& eout, int& itw);
    bool Gamma_En(double px, double py, double pz, double& eout, int& itw,
        double& ximp, double& yimp, double& zimp);
    bool Gamma_Pos(double& px, double& py, double& pz);

    std::string outFileName;
    MCMethod mcmethod;

    /* Tower status for sim warnmap */
    int tower_status_sim[NSEC][NY][NZ];

    int NPart;
    int NPeak;
    TLorentzVector Vpart[MAXPEAK];
    int itw_part[MAXPEAK];
    int sec_part[MAXPEAK];
    double eTwr[NSEC][NY][NZ];

    PHPythiaContainer *phpythia;

    TDatabasePDG *pdg_db;
    EMCWarnmapChecker *emcwarnmap;
    DCDeadmapChecker *dcdeadmap;

    Fun4AllHistoManager *hm;
    TH1 *h_events;
    TH1 *h_pion;
    TH1 *h_photon;
    TH1 *h_photon_eta050;
    TH1 *h_photon_eta025;
    TH1 *h_isophoton_eta050;
    TH1 *h_isophoton_eta025;
    TH2 *h2_pion_eta_phi[nh_eta_phi];
    TH2 *h2_photon_eta_phi[nh_eta_phi];
    TH3 *h3_isopi0;
    TH3 *h3_isoeta;
    THnSparse* hn_pion;
    THnSparse* hn_missing;
    THnSparse* hn_missing_eta;
    THnSparse *hn_photon;
    THnSparse *hn_geom;
    THnSparse *hn_isolated;

    TF1 *cross_pi0;
    TF1 *cross_ph;
    double weight_pythia;
};

#endif	/* __ANAFASTMC_H__ */
