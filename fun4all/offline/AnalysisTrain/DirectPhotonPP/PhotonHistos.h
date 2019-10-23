#ifndef __PHOTONHISTOS_H__
#define __PHOTONHISTOS_H__

#include <SubsysReco.h>

class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;
class EMCWarnmapChecker;
class DCDeadmapChecker;

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

class PhotonHistos: public SubsysReco
{
  public:
    PhotonHistos(const std::string &name = "PhotonHistos", const char *filename = "PhotonHistos.root");
    virtual ~PhotonHistos();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    /* Select to run on MinBias or ERT sample */
    void SelectMB();
    void SelectERT();

  protected:
    /* Number of histogram array */
    static const int nh_calib = 8*2;
    static const int nh_bbc = 3*2;
    static const int nh_ertsm = 8*2*3;
    static const int nh_ert = 3*2*3*2*2*3;
    static const int nh_dcpartqual = 2*2*3;
    static const int nh_dcgood = 2;
    static const int nh_pion = 3*3*2*2*2*2*3;
    static const int nh_pion_pol = 2*2;
    static const int nh_eta_phi = 3*2*2*3;
    static const int nh_1photon = 3*3*2*2*3;
    static const int nh_1photon_pol = 2*2*2;
    static const int nh_2photon = 3*3*2*2*2*3;
    static const int nh_2photon_pol = 2*2*2*2*2;
    static const int nh_mul_pion = 2;
    static const int nh_mul_photon = 2*2;

    /* pT bins for ALL */
    static const int npT_pol = 15;
    static double pTbin_pol[npT_pol+1];

    /* Event counts */
    int FillEventCounts(const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1);

    /* ToF and energy calibration */
    int FillClusterTofSpectrum(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const std::string &qualii = "");
    int FillPi0InvariantMass(const emcClusterContainer *data_emccontainer, const PHGlobal *data_global, const std::string &quali = "");

    /* BBC and ERT trigger efficiency */
    int FillBBCEfficiency(const emcClusterContainer *data_emccontainer, const TrigLvl1 *data_triggerlvl1);
    int FillERTEfficiency(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert);

    /* Study DC track quality */
    int FillTrackQuality(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert);

    /* Count pi0 yield */
    int FillPi0Spectrum(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert);

    /* Count direct photon yield */
    int FillPhotonSpectrum(const emcClusterContainer *data_emccontainer, const PHCentralTrack *data_tracks,
        const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1, const ErtOut *data_ert);

    /* Create histograms */
    void BookHistograms();

    /* Sum energy in cone around the reference particle
     * for isolated photon and isolated pair */
    void SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks, double bbc_t0, double econe[]);
    void SumEEmcal(const emcClusterContent *cluster, const emcClusterContent *cluster_part, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks, double bbc_t0, double econe[]);
    void SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks, double econe[]);
    void SumEPi0(const emcClusterContent *cluster1, const emcClusterContent *cluster2, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks, double bbc_t0, double econe[]);

    /* Check event type, BBC and photon cuts */
    bool IsEventType(int evtype, const TrigLvl1 *data_triggerlvl1);
    bool BBC10cm(const PHGlobal *data_global, const TrigLvl1 *data_triggerlvl1);
    bool TestPhoton(const emcClusterContent *cluster, double bbc_t0);
    bool PassChargeVeto(const emcClusterContent *cluster);

    /* Get spin pattern and EMCal associated track */
    int GetPattern(int crossing);

    /* Update spin pattern information and store in class */
    void UpdateSpinPattern(SpinDBContent &spin_cont);

    enum DataType {MB, ERT};
    DataType datatype;

    /* EMCal recalibrator, DC deadmap and spin information*/
    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;
    EMCWarnmapChecker *emcwarnmap;
    DCDeadmapChecker *dcdeadmap;
    SpinPattern *spinpattern;

    int runnumber;
    int fillnumber;

    /* Output histograms */
    std::string outFile;
    Fun4AllHistoManager *hm;
    TH1 *h_events;
    TH2 *h2_tof[nh_calib];
    TH2 *h2_minv[nh_calib];
    TH1 *h_bbc[nh_bbc];
    TH2 *h2_bbc_pion[nh_bbc];
    TH2 *h2_ertsm[nh_ertsm];
    TH1 *h_ert[nh_ert];
    TH2 *h2_ert_pion[nh_ert];
    TH3 *h3_dcdphiz[nh_dcpartqual];
    TH2 *h2_alphaboard[nh_dcpartqual];
    TH3 *h3_dclive[nh_dcgood];
    TH1 *h_prod;
    TH2 *h2_pion[nh_pion];
    TH2 *h2_pion_pol[nh_pion_pol];
    TH2 *h2_eta_phi[nh_eta_phi];
    TH1 *h_1photon[nh_1photon];
    TH1 *h_1photon_pol[nh_1photon_pol];
    TH2 *h2_2photon[nh_2photon];
    TH2 *h2_2photon2pt[nh_2photon];
    TH2 *h2_2photon_pol[nh_2photon_pol];
    TH2 *h2_2photon2pt_pol[nh_2photon_pol];
    TH2 *h2_mul_pion_sig[nh_mul_pion];
    TH2 *h2_mul_pion_bg[nh_mul_pion];
    TH2 *h2_mul_photon_sig[nh_mul_photon];
    TH2 *h2_mul_photon_bg[nh_mul_photon];
};

#endif /* __PHOTONHISTOS_H__ */
