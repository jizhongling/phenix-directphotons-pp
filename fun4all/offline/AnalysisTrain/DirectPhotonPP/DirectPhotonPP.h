#ifndef __DIRECT_PHOTON_PP_H__
#define __DIRECT_PHOTON_PP_H__

#include <SubsysReco.h> // always put includes with "" before those with <>

#include <string>
#include <vector>
#include <fstream>

/* Local analysis Classes */
class EmcLocalRecalibrator;
class EmcLocalRecalibratorSasha;

/* Fun4All classes */
class PHCentralTrack;
class emcClusterContainer;
class emcClusterContent;
class PHGlobal;
class ErtOut;
class TrigLvl1;
class Fun4AllHistoManager;

/* Root classes */
class TH1;
class TH2;
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

  /**
   * Set local recalibrator- Sasha style
   */
  void SetEmcLocalRecalibrator( EmcLocalRecalibratorSasha* emcrecalib )
  {
    _emcrecalib_sasha = emcrecalib;
  }


  /**
   * Choose type of DST input file- MinBias or ERT?
   */
  void SetDstDataType( std::string newtype )
  {
    _dsttype = newtype;
  }


  /**
   * Read tower status from text file with 4 columns (sector, y-idx, z-idx, status)
   */
  void ReadTowerStatus(const std::string &filename, unsigned ncols)
  {
    if ( ncols == 4 )
      ReadTowerStatus4Cols(filename);
    else if ( ncols == 2 )
      ReadTowerStatusSasha(filename);
  }


  /**
   * Set debug mode for detailed cluster information output
   */
  void SetClusterDebugMode( bool mode )
  {
    _debug_cluster = mode;
  }

  /**
   * Set debug mode for detailed event selection  information output
   */
  void SetTriggerDebugMode( bool mode )
  {
    _debug_trigger = mode;
  }

  /**
   * Set debug mode for detailed cluster pairing information output
   */
  void SetPi0DebugMode( bool mode )
  {
    _debug_pi0 = mode;
  }

protected:

  /**
   * Access tower status from tower status array (member variable)
   */
  int get_tower_status( int sector,
                        int ybin,
                        int zbin);

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
  bool testPhoton( emcClusterContent *emccluster,
                   double bbct0 );

  /**
   * Check if cluster matches criteria for direct photon candidate
   */
  bool testDirectPhoton( emcClusterContent *emccluster,
                         double bbct0 );

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
  bool testPhotonTof( emcClusterContent *emccluster,
                      double bbct0 );

  /**
   * Test cluster shower shape
   */
  bool testPhotonShape( emcClusterContent *emccluster );

  /**
   * Test cluster track veto
   */
  bool testPhotonTrackVeto( emcClusterContent *emccluster );

  /**
   * Select only clusters in good towers
   */
  int selectClusterGoodTower( emcClusterContainer *emccontainer );

  /**
   * Select only clusters which have EM like (photon or electron) shape
   */
  int selectClusterPhotonShape( emcClusterContainer *emccontainer );

  /**
   * Select only clusters above photon cutoff energy
   */
  int selectClusterPhotonEnergy( emcClusterContainer *emccontainer );

  /**
   * Select only clusters which have TOF within range for photons
   */
  int selectClusterPhotonTof( emcClusterContainer *emccontainer, double bbc_t0 );

  /**
   * Read tower status from text file with 4 columns (sector, y-idx, z-idx, status)
   */
  void ReadTowerStatus4Cols(const std::string &filename);

  /**
   * Read tower status from text file with 2 columns (tower id, status - e.g. Sasha)
   */
  void ReadTowerStatusSasha(const std::string &filename);

  /**
   * Fill histogram with trigger based event counts
   */
  int FillTriggerStats( std::string, TrigLvl1*, double );


  /**
   * Fill histograms with cluster pT spectrum for trigger efficiency
   */
  int FillTriggerEfficiency( emcClusterContainer *data_emccontainer,
                             PHGlobal *data_global,
                             ErtOut *data_ert );

  /**
   * Fill histograms with cluster pT spectrum before and after applying bad tower map
   */
  int FillClusterPtSpectrum( std::string,
			     emcClusterContainer* );

  /**
   * Fill histograms with cluster TOF spectrum before and after applying local TOF correction
   */
  int FillClusterTofSpectrum( std::string histname,
			      emcClusterContainer *data_emc,
			      PHGlobal *data_global,
			      double bbc_t0 );

  /**
   * Fill histograms with invariant mass from two-photon pairs which are pi0 candidates
   * using tight cuts on pi0 identification before and after applying local tower energy correction
   */
  int FillPi0InvariantMass( std::string histname,
			    emcClusterContainer *d_emcont,
			    TrigLvl1 *d_trig,
			    ErtOut *data_ert );

  /**
   * Fill histograms with photon pT spectrum and invariant mass histogram
   * for combinations with all other photons
   */
  int FillPhotonPtSpectrum( emcClusterContainer *d_emccontainer,
                            PHCentralTrack* d_tracks,
                            PHGlobal *d_global );

  /**
   * Print infomration of cluster container
   */
  void PrintClusterContainer( emcClusterContainer* , double );

  /**
   * Array providing status for each EMCal tower. Array indices are [sector][ytower][ztower]
   * 0=dead 1=good 10=iso fiducial 50=fiducial 100=hot
   */
  int _tower_status[8][48][96];

  /**
   * Type of DST file- MinBias or ERT?
   */
  std::string _dsttype;

  /**
   * Event counter
   */
  int _ievent;

  /**
   * Run number
   */
  int _runnumber;

  /**
   * BBC z vertex range cut (in cm)
   */
  double _bbc_zvertex_cut;

  /**
   * minimum energy for cluster to be considered photon (in GeV)
   */
  double _photon_energy_min;

  /**
   * minimum EM shower shape probability for cluster to be considered as photon
   */
  double _photon_prob_min;

  /**
   * minimum TOF for cluster to be considered as photon (in ns)
   */
  double _photon_tof_min;

  /**
   * maximum TOF for cluster to be considered as photon (in ns)
   */
  double _photon_tof_max;

  /**
   * minimum energy for cluster to be considered direct photon (in GeV)
   */
  double _direct_photon_energy_min;

  /**
   * On-the-fly recalibration of EMCal towers
   */
  EmcLocalRecalibrator* _emcrecalib;

  /**
   * On-the-fly recalibration of EMCal towers- based on Sasha's analysis
   */
  EmcLocalRecalibratorSasha* _emcrecalib_sasha;

  /**
   * Name for output ROOT file for histograms
   */
  std::string _outfile_histos;

  /**
   * histogram manager
   */
  Fun4AllHistoManager *_hm;

  /**
   * switch- set to TRUE to print detailed cluster information to log file
   */
  bool _debug_cluster;

  /**
   * switch- set to TRUE to print detailed trigger information to log file
   */
  bool _debug_trigger;

  /**
   * switch- set to TRUE to print detailed cluster pairing information to log file
   */
  bool _debug_pi0;

};
#endif
