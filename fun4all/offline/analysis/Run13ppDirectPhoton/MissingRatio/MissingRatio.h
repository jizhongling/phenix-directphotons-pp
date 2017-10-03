#ifndef __MISSINGRATIO_H__
#define __MISSINGRATIO_H__

#include <SubsysReco.h>
#include <string>
#include <map>

class AnaTrk;
class emcGeaClusterContent;

class PHCompositeNode;
class Fun4AllHistoManager;

class TFile;
class TF1;
class TH2;
class THnSparse;

class MissingRatio: public SubsysReco
{
  public:
    MissingRatio( const char *filename = "histo.root");
    virtual ~MissingRatio();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:
    void BookHistograms();
    void ReadTowerStatus(const std::string& filename);
    void ReadSashaWarnmap(const std::string& filename);

    // store photon info for different criterias
    // criteria = 4*part+2*prob+warnmap
    void StorePhoton(std::multimap<int,AnaTrk*> *photon, std::multimap<int,AnaTrk*> *target, AnaTrk *trk, int part, int prob, int warnmap);
    void StorePhoton(std::multimap<int,AnaTrk*> *photon, std::multimap<int,AnaTrk*> *target, AnaTrk *trk);
    void FillCluster(emcGeaClusterContent *emccluster, double pionpT, double weight);
    void FillTwoClusters(emcGeaClusterContent *emccluster1, emcGeaClusterContent *emccluster2, double weight);

    // get weight by pT
    double Get_wpT(double pT);

    // tower status for warnmap
    int tower_status[8][48][96];

    // number of pT bins
    static const int npT = 31;

    // pT bins
    double vpT[npT];

    // pT weight
    double wpT[npT];

    std::string outFileName;
    Fun4AllHistoManager *hm;

    // 2D histograms for west and east arms with different criterias
    TH2 *h2_radius;
    TH2 *h2_angle;
    TH2 *h2_missing;
    TH2 *h2_nomissing;
    TH2 *h2_missing_pion;
    TH2 *h2_nomissing_pion;
    TH2 *h2_incident;
    TH2 *h2_measured;
    TH2 *h2_incident_pion;
    TH2 *h2_measured_pion;
    TH2 *h2_conversion;
    TH2 *h2_photon;
    TH2 *h2_noconv;
    TH2 *h2_vtxconv;
    TH2 *h2_merge_photon;
    TH2 *h2_merge_pion;
    THnSparse *hn_conversion_position;
    THnSparse *hn_photon;
    THnSparse *hn_pion;

    TF1 *cross;
};

#endif /* __MISSINGRATIO_H__ */
