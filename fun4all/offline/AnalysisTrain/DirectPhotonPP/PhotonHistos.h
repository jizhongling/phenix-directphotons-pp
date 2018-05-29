#ifndef __PHOTONHISTOS_H__
#define __PHOTONHISTOS_H__

#include <SubsysReco.h>

class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;

class PhotonContainer;
class SpinPattern;
class SpinDBContent;

class PHGlobal;
class PHCentralTrack;
class emcClusterContainer;
class emcClusterContent;
class emcTowerContainer;
class TrigLvl1;
class ErtOut;
class PHCompositeNode;

class Fun4AllHistoManager;
class TH1;
class TH2;
class TH3;
class THnSparse;

class PhotonHistos: public SubsysReco
{
  public:
    PhotonHistos(const std::string &name = "PHOTONHISTOS", const char *filename = "PhotonHistos.root");
    virtual ~PhotonHistos();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SelectMB();
    void SelectERT();

  protected:
    int FillClusterTofSpectrum(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const std::string &qualii = "");
    int FillPi0InvariantMass(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const std::string &quali = "");
    int FillBBCEfficiency(const emcClusterContainer *data_emccontainer, const TrigLvl1 *data_triggerlvl1);
    int FillERTEfficiency(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global,
        const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype);
    int FillTowerEnergy(const emcClusterContainer *data_emccontainer, const emcTowerContainer *data_emctwrcontainer,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1);
    int FillPi0Spectrum(const emcClusterContainer *data_emccontainer,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype);
    int FillPhotonSpectrum(const emcClusterContainer *data_emccontainer,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert, const int evtype);

    void BookHistograms();
    void EMCRecalibSetup();
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    bool IsEventType(const int evtype, const TrigLvl1 *data_triggerlvl1);
    bool TestPhoton(const emcClusterContent *cluster, double bbc_t0);
    bool DispCut(const emcClusterContent *cluster);

    int GetStatus(const emcClusterContent *cluster);
    double GetTrackConeEnergy(const PHCentralTrack *tracks, const emcClusterContent *cluster, double cone_angle);
    int GetPattern(int crossing);

    void UpdateSpinPattern(SpinDBContent &spin_cont);

    // some constants
    static const double PI = 3.14151927;
    static const int NSEC = 8;
    static const int NY = 48;
    static const int NZ = 96;

    // some cuts
    static const double eMin = 0.3;
    static const double probMin = 0.02;
    static const double tofMax = 10.;
    static const double AsymCut = 0.8;

    // triger bit
    static const unsigned bit_ppg = 0x70000000;
    static const unsigned bit_bbcnarrow = 0x00000010;
    static const unsigned bit_bbcnovtx = 0x00000002;
    static const unsigned bit_ert4x4[3];  // ert4x4a/b/c

    enum DataType {MB, ERT};
    DataType datatype;

    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;
    PhotonContainer *photoncont;
    SpinPattern *spinpattern;

    int runnumber;
    int fillnumber;

    std::string outFile;
    Fun4AllHistoManager *hm;
    TH1 *h_events;
    TH3 *h3_tof;
    TH3 *h3_tof_raw;
    TH3 *h3_minv;
    TH3 *h3_minv_raw;
    TH3 *h3_bbc;
    THnSparse *hn_bbc_pion;
    TH3 *h3_ert;
    THnSparse *hn_ert_pion;
    TH2 *h2_photon_eta_phi[3];
    TH2 *h2_cluster_eta_phi[3];
    THnSparse *hn_etwr;
    THnSparse *hn_pion;
    THnSparse *hn_1photon;
    THnSparse *hn_2photon;
    THnSparse *hn_photonbg;

    // tower status for warnmap
    int tower_status[8][48][96];
    int tower_status_sasha[8][48][96];
};

#endif /* __PHOTONHISTOS_H__ */
