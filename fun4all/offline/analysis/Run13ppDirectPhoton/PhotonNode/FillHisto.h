#ifndef __FILLHISTO_H__
#define __FILLHISTO_H__

#include <SubsysReco.h>

#include <string>

class PhotonContainer;
class Photon;
class PhotonERT;
class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;

class PHCompositeNode;
class Fun4AllHistoManager;

class TH1;
class TH3;
class THnSparse;

class FillHisto: public SubsysReco
{
  public:
    FillHisto(const std::string &name = "FillHisto", const char *filename = "histo.root");
    virtual ~FillHisto();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int EndRun(const int runnumber);
    int End(PHCompositeNode *topNode);

    void SelectMB();
    void SelectERT();

  protected:
    int FillClusterTofSpectrum( const PhotonContainer *photoncont, const std::string &quali = "" );
    int FillPi0InvariantMass( const PhotonContainer *photoncont, const std::string &quali = "" );
    int FillBBCEfficiency(const PhotonContainer *photoncont);
    int FillERTEfficiency(const PhotonContainer *photoncont);
    int FillPhotonSpectrum(const PhotonContainer *photoncont);
    int FillPi0Spectrum(const PhotonContainer *photoncont);

    void BookHistograms();
    void EMCRecalibSetup();
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    bool TestPhoton(const Photon *photon, double bbc_t0);
    bool TestTrackVeto(const PhotonERT *photon);

    int GetStatus(const Photon *photon);
    int GetPattern(const PhotonContainer *photoncont);

    enum DataType {MB, ERT};
    DataType datatype;

    std::string outFileName;
    Fun4AllHistoManager *hm;
    EmcLocalRecalibrator *emcrecalib;
    EmcLocalRecalibratorSasha *emcrecalib_sasha;

    TH1 *h_events;
    TH3 *h3_tof;
    TH3 *h3_tof_raw;
    TH3 *h3_minv;
    TH3 *h3_minv_raw;
    TH3 *h3_bbc;
    THnSparse *hn_bbc_pion;
    TH3 *h3_ert;
    THnSparse *hn_ert_pion;
    THnSparse *hn_1photon;
    THnSparse *hn_2photon;
    THnSparse *hn_pion;

    // tower status for warnmap
    int tower_status[8][48][96];

    int crossing_shift;
    int spinpattern_blue[120];
    int spinpattern_yellow[120];

    int runnumber;
    double runtime;
};

#endif /* __FILLHISTO_H__ */
