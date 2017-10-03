#ifndef __FILLHISTO_H__
#define __FILLHISTO_H__

#include <SubsysReco.h>

#include <string>

class Photon;
class PhotonContainer;
class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;

class PHCompositeNode;
class Fun4AllHistoManager;

class TH1;
class TH3;
class THnSparse;
class TGraphErrors;

class FillHisto: public SubsysReco
{
  public:
    FillHisto(const std::string &name = "FILLHISTO", const char *filename = "histo.root");
    virtual ~FillHisto();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int EndRun(const int runnumber);
    int End(PHCompositeNode *topNode);

  protected:
    int FillClusterTofSpectrum( const PhotonContainer *photoncont, const std::string &quali = "" );
    int FillPi0InvariantMass( const PhotonContainer *photoncont, const std::string &quali = "" );
    int FillTriggerEfficiency(const PhotonContainer *photoncont);
    int FillSinglePhotonSpectrum(const PhotonContainer *photoncont);
    int FillTwoPhotonSpectrum(const PhotonContainer *photoncont);
    int FillPi0Spectrum(const PhotonContainer *photoncont);

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
    int GetPattern(const PhotonContainer *photoncont);
    int GetStatus(const Photon *photon);

    bool TestPhoton(const Photon *photon, double bbc_t0);
    //bool TestTrackVeto(const Photon *photon);

    std::string outFileName;
    Fun4AllHistoManager *hm;
    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;

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
    THnSparse *hn_asym;
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
    int runnumber;
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

#endif /* __FILLHISTO_H__ */
