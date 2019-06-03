#ifndef __HADRONRESPONSE_H__
#define __HADRONRESPONSE_H__

#include <SubsysReco.h>
#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;
class PHCentralTrack;
class emcClusterContainer;
class emcClusterContent;

class TFile;
class TH1;
class TH2;
class THnSparse;

class HadronResponse: public SubsysReco
{
  public:
    HadronResponse(const std::string &name = "HadronResponse",
        const char *filename = "histo.root",
        const double pt_init = 0.);
    virtual ~HadronResponse();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void InitBatch();

  protected:
    void BookHistograms();
    void ReadTowerStatus(const std::string& filename);
    void ReadSashaWarnmap(const std::string& filename);

    /* Sum energy in cone around the reference particle
     * for isolated photon and isolated pair */
    double SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks);
    double SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks);

    /* Check event type, photon cuts, charge veto and tower status */
    bool DCChargeVeto(const emcClusterContent *cluster, const PHCentralTrack *data_tracks);
    bool InFiducial(const emcClusterContent *cluster);
    bool IsGoodTower(const emcClusterContent *cluster);
    bool IsBadTower(const emcClusterContent *cluster);

    /* EMCal associated track */
    int GetEmcMatchTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks);

    /* Pythia weight normalization */
    TH1 *h_pt_weight;
    double pt_start;
    int pt_count;
    double weight_pythia;

    // number of pT bins
    static const int npT = 30;

    // pT bins
    static double vpT[npT+1];

    // tower status for warnmap
    int tower_status_nils[8][48][96];
    int tower_status_sasha[8][48][96];

    std::string outFileName;
    Fun4AllHistoManager *hm;

    THnSparse *hn_dc;
    THnSparse *hn_emcal;
    THnSparse *hn_cluster;
};

#endif /* __HADRONRESPONSE_H__ */
