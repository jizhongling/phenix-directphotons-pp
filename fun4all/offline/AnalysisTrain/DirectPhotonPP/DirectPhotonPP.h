#ifndef __DIRECT_PHOTON_PP_H__
#define __DIRECT_PHOTON_PP_H__

#include <SubsysReco.h> // always put includes with "" before those with <>

#include <string>
#include <vector>
#include <fstream>

/* Fun4All Classes */
//class Fun4AllServer;
//class Fun4AllHistoManager;

/* Analysis Classes */
class PHCentralTrack;
//class McEvalSingleList_v1;
//class PHCompositeNode;
//class EventHeader;
class emcClusterContainer;
class emcClusterContent;
//class SvxClusterList;
class PHGlobal;
class ErtOut;
class TrigLvl1;
//class EventHeader;
//class RunHeader;
//class SvxCentralTrackList;
//class SvxCentralTrack;
class Fun4AllHistoManager;

class EmcLocalRecalibrator;

/* Root classes */
class TH1;
class TH2;
class TH3;
class THnSparse;
class TFile;

class DirectPhotonPP: public SubsysReco
{
  public:

    /**
     * Default constructor
     */
    DirectPhotonPP(const char* outputfilename);

    /**
     * Default destructor
     */
    virtual ~DirectPhotonPP();

    /**
     * Fun4All method
     */
    int Init(PHCompositeNode *topNode);

    /**
     * Fun4All method
     */
    int InitRun(PHCompositeNode *topNode);

    /**
     * Fun4All method
     */
    int process_event(PHCompositeNode *topNode);

    /**
     * Fun4All method
     */
    int End(PHCompositeNode *topNode);

    /**
     * Set local recalibrator
     */
    void SetEmcLocalRecalibrator( EmcLocalRecalibrator* emcrecalib )
    {
      _emcrecalib = emcrecalib;
    }

  protected:

    /**
     * Access tower status from tower status array (member variable)
     */
    int get_tower_status( int sector, int ybin, int zbin);

    /**
     * Check cluster is centered on good tower, i.e. NOT flagged bad on warnmap
     */
    bool testGoodTower( emcClusterContent *emccluster );

    /**
     * Check cluster is inside tight fiducial volume cut
     */
    bool testTightFiducial( emcClusterContent *emccluster );

    /**
     * Check if cluster matches criteria for photon candidate
     */
    bool testPhoton( emcClusterContent *emccluster , double bbct0 );

    /**
     * Check if cluster matches criteria for direct photon candidate
     */
    bool testDirectPhoton( emcClusterContent *emccluster , double bbct0 );

    /**
     * Check if cluster matches criteria for isolated photon candidate
     */
    bool testIsolatedPhoton( emcClusterContent *emccluster0 ,
        emcClusterContainer *emccontainer ,
        PHCentralTrack *tracks ,
        double coneangle ,
        double bbct0 );

    /**
     * Test cluster energy
     */
    bool testPhotonEnergy( emcClusterContent *emccluster );

    /**
     * Test cluster TOF
     */
    bool testPhotonTof( emcClusterContent *emccluster , double bbct0 );

    /**
     * Test cluster shower shape
     */
    bool testPhotonShape( emcClusterContent *emccluster );

    /**
     * Test cluster track veto
     */
    bool testPhotonTrackVeto( emcClusterContent *emccluster );


  private:

    /**
     * Read tower status from text file
     */
    void ReadTowerStatus(const std::string &filename);
    void ReadSashaWarnmap(const std::string &filename);

    /**
     * Fill histograms with cluster pT spectrum for trigger efficiency
     */
    int FillTriggerEfficiency( emcClusterContainer *data_emccontainer, PHGlobal *data_global, ErtOut *data_ert );

    /**
     * Fill histograms with cluster pT spectrum before and after applying bad tower map
     */
    int FillClusterPtSpectrum( emcClusterContainer *d_emcont, PHGlobal *d_gbl );

    /**
     * Fill histograms with cluster TOF spectrum before and after applying local TOF correction
     */
    int FillClusterTofSpectrum( emcClusterContainer *d_emcont, PHGlobal *d_gbl, std::string quali="" );

    /**
     * Fill histograms with invariant mass from two-photon pairs which are pi0 candidates
     * using tight cuts on pi0 identification before and after applying local tower energy correction
     */
    int FillPi0InvariantMass( emcClusterContainer *d_emcont, PHGlobal *d_gbl, TrigLvl1 *d_trig, ErtOut *data_ert, std::string quali="" );

    /**
     * Fill histograms with photon pT spectrum and invariant mass histogram
     * for combinations with all other photons
     */
    int FillPhotonPtSpectrum( emcClusterContainer *d_emccontainer, PHCentralTrack* d_tracks, PHGlobal *d_global );

    /**
     * Array providing status for each EMCal tower. Array indices are [sector][ytower][ztower]
     * 0=dead 1=good 10=iso fiducial 50=fiducial 100=hot
     */
    int _tower_status[8][48][96];

    /**
     * Event counter
     */
    int _ievent;

    /**
     * BBC z vertex range cut (in cm)
     */
    float _bbc_zvertex_cut;

    /**
     * minimum energy for cluster to be considered photon (in GeV)
     */
    float _photon_energy_min;

    /**
     * minimum EM shower shape probability for cluster to be considered as photon
     */
    float _photon_prob_min;

    /**
     * minimum TOF for cluster to be considered as photon (in ns)
     */
    float _photon_tof_min;

    /**
     * maximum TOF for cluster to be considered as photon (in ns)
     */
    float _photon_tof_max;

    /**
     * minimum energy for cluster to be considered direct photon (in GeV)
     */
    float _direct_photon_energy_min;

    /**
     * On-the-fly recalibration of EMCal towers
     */
    EmcLocalRecalibrator* _emcrecalib;

    /**
     * Name for output ROOT file for histograms
     */
    std::string _outfile_histos;

    /* histogram manager */
    Fun4AllHistoManager *_hm;

};
#endif
