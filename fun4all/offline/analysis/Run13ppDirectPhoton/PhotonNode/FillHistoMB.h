#ifndef __FILLHISTOMB_H__
#define __FILLHISTOMB_H__

#include <SubsysReco.h>

#include <string>

class PhotonMB;
class PhotonContainerMB;
class EmcLocalRecalibratorMB;

class PHCompositeNode;
class Fun4AllHistoManager;

class TH1;
class TH3;
class THnSparse;
class TGraphErrors;

class FillHistoMB: public SubsysReco
{
  public:
    FillHistoMB(const std::string &name = "FILLHISTOMB", const char *filename = "histo.root");
    virtual ~FillHistoMB();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int EndRun(const int runnumber); int End(PHCompositeNode *topNode);

  protected:
    int FillClusterTofSpectrum( const PhotonContainerMB *photoncont, const std::string &quali = "" );
    int FillPi0InvariantMass( const PhotonContainerMB *photoncont, const std::string &quali = "" );
    int FillTriggerEfficiency(const PhotonContainerMB *photoncont);
    int FillSinglePhotonSpectrum(const PhotonContainerMB *photoncont);
    int FillTwoPhotonSpectrum(const PhotonContainerMB *photoncont);
    int FillPi0Spectrum(const PhotonContainerMB *photoncont);
    int FillPileup(const PhotonContainerMB *photoncont);

    void BookHistograms();
    void EMCRecalibSetup();
    void ReadClockCounts(const std::string &filename);
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    unsigned long long GetClockLive(unsigned runnumber);
    unsigned long long GetBBCNovtxLive(unsigned runnumber);
    unsigned long long GetBBCNarrowLive(unsigned runnumber);
    unsigned long GetBBCNovtxScaledown(unsigned runnumber);
    unsigned long GetBBCNarrowScaledown(unsigned runnumber);
    int GetPattern(const PhotonContainerMB *photoncont);
    int GetStatus(const PhotonMB *photon);

    bool TestPhoton(const PhotonMB *photon, double bbc_t0);
    //bool TestTrackVeto(const PhotonMB *photon);

    std::string outFileName;
    Fun4AllHistoManager *hm;
    EmcLocalRecalibratorMB *emcrecalib;

    TH1 *h_events;
    TH3 *h3_tof;
    TH3 *h3_tof_raw;
    TH3 *h3_inv_mass_pi0calib;
    TH3 *h3_inv_mass_pi0calib_raw;
    TH3 *h3_trig;
    TH3 *h3_trig_pion;
    THnSparse *hn_1photon;
    THnSparse *hn_2photon;
    THnSparse *hn_pion;
    THnSparse *hn_minv;
    TGraphErrors *g_pileup_PbSc;
    TGraphErrors *g_pileup_PbGl;
    TGraphErrors *g_pileup_PbSc_notof;
    TGraphErrors *g_pileup_PbGl_notof;

    // tower status for warnmap
    int tower_status[8][48][96];

    int crossing_shift;
    int spinpattern_blue[120];
    int spinpattern_yellow[120];

    int irun;
    double runtime;
    // 0: Runnumber; 1: CLOCK live trigger count; 2: BBC novtx live count; 3: BBC narrow live count;
    // 4: BBC novtx scaledown; 5: BBC narrow scaledown
    unsigned long long n_clock_bbc[6][1020];
    unsigned long long nmb;
    unsigned long long npions_sig[2];
    unsigned long long npions_bg[2];
    unsigned long long npions_sig_notof[2];
    unsigned long long npions_bg_notof[2];
};

#endif /* __FILLHISTOMB_H__ */
