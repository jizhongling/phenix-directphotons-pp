#ifndef __ISOLATIONCUT_H__
#define __ISOLATIONCUT_H__

#include <SubsysReco.h>
#include <string>
#include <map>
#include <vector>

class TFile;
class TTree;
class THnSparse;

class emcGeaTrackContent;
class emcGeaClusterContent;
class emcClusterContent;
class emcClusterContainer;
class PHCentralTrack;

class IsolationCut: public SubsysReco
{
public:
  IsolationCut( const char *filename = "isocut.root");
  virtual ~IsolationCut();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  /** Find truth particle with maximum deposited energy contribution to given cluster */
  emcGeaTrackContent* FindTruthParticle( emcGeaClusterContent* cluster );

  /** Sum and return all the energies of clusters found inside of cone of given radius around
      a photon candidate */
  float SumEmcalEnergyInCone( emcClusterContent*,
			      emcClusterContainer*,
			      double, double );

  /** Sum and return all the energies (momenta) of charged tracks found inside of cone of given radius around
      a photon candidate */
  float SumTrackEnergyInCone( emcClusterContent*,
			      PHCentralTrack*,
			      double, double, double );

  /** Reset global variables that store cluster information for filling output tree */
  void ResetBranchVariables();

  /** Read warnmap */
  void ReadTowerStatus(const std::string& filename);

  /** Read warnmap (Sasha's format) */
  void ReadSashaWarnmap(const std::string& filename);

  /** Get tower status for cluster from warnmap */
  int GetStatus(const emcClusterContent *emccluster);

  /** event counter */
  unsigned _ievent;

  /** count events with 1+ photon candidate cluster */
  unsigned _event_nphotons;

  /** tower status for warnmap */
  int _tower_status[8][48][96];

  /** vector with PID's of neutral particles */
  std::vector< int > _v_pid_neutral;

  /** histogram storing cone energy information */
  THnSparse*  _hn_energy_cone;

  /** histogram storing cone energy information */
  THnSparse*  _hn_energy_cone_reco;

  /** output tree with cluster information */
  TTree* _tree_event_cluster;

  /** output tree with truth information */
  TTree* _tree_event_truth;

  /** map with cluster variables */
  std::map< std::string , float > _map_cluster_branches;

  /** truth tree variables */
  float _truth_pid;
  float _truth_parentpid;
  float _truth_anclvl;
  float _truth_ptot;
  float _truth_pt;
  float _truth_eta;
  float _truth_phi;

  /** output file name */
  std::string _output_file_name;

  /** output file */
  TFile *_file_output;

  /** enum for PISA PID */
  enum PisaPid
    {
      PHOTON = 1,
      POSITRON = 2,
      ELECTRON = 3,
      NEUTRINO = 4,
      PIZERO = 7,
      KAON_0_LONG = 10,
      NEUTRON = 13,
      KAON_0_SHORT = 16,
      ETA = 17,
      LAMBDA = 18,
      SIGMA_0 = 20,
      XI_0 = 22,
      ANTINEUTRON = 25,
      ANTILAMBDA = 26,
      ANTISIGMA_0 = 28,
      ANITXI_0 = 30
    };

};

#endif /* __ISOLATIONCUT_H__ */
