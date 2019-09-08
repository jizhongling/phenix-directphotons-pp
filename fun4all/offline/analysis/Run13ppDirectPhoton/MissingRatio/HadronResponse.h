#ifndef __HADRONRESPONSE_H__
#define __HADRONRESPONSE_H__

#include <SubsysReco.h>
#include <string>

class ERTSimTrigger;
class EMCWarnmapChecker;
class DCDeadmapChecker;

class Fun4AllHistoManager;
class PHCompositeNode;
class emcClusterContainer;
class emcClusterContent;
class PHCentralTrack;

class TFile;
class TH1;
class TH2;
class THnSparse;

class HadronResponse: public SubsysReco
{
  public:
    HadronResponse(const std::string &name = "HadronResponse", const char *filename = "histo.root");
    virtual ~HadronResponse();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SetWeightPythia(double weight) { weight_pythia = weight; }
    void set_outfile(std::string filename) { outFileName = filename; }

  protected:
    /* Number of histogram array */
    static const int nh_eta_phi = 3*2*2*2;

    /* Create histograms */
    void BookHistograms();

    /* Sum energy in cone around the reference particle
     * for isolated photon and isolated pair */
    double SumEEmcal(const emcClusterContent *cluster, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks);
    void SumEEmcal(const emcClusterContent *cluster1, const emcClusterContent *cluster2, const emcClusterContainer *cluscont,
        const PHCentralTrack *data_tracks, double &econe);
    double SumPTrack(const emcClusterContent *cluster, const PHCentralTrack *data_tracks);

    /* Check charge veto and tower status */
    bool TestPhoton(const emcClusterContent *cluster);

    /* Pythia pT weights */
    double weight_pythia;

    /* ERT sim trigger */
    ERTSimTrigger *ertsim;

    /* EMC warnmap checker */
    EMCWarnmapChecker *emcwarnmap;

    /* DC deadmap checker */
    DCDeadmapChecker *dcdeadmap;

    std::string outFileName;
    Fun4AllHistoManager *hm;

    TH1 *h_events;
    TH2 *h2_eta_phi[nh_eta_phi];
    THnSparse *hn_1photon;
    THnSparse *hn_2photon;
    THnSparse *hn_prob_photon;
    THnSparse *hn_dcdphiz;
    THnSparse *hn_alphaboard;
    THnSparse *hn_dclive;
    TH1 *h_prod;
};

#endif /* __HADRONRESPONSE_H__ */
