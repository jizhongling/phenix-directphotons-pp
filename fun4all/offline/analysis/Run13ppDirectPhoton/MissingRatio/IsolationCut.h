#ifndef __ISOLATIONCUT_H__
#define __ISOLATIONCUT_H__

#include <SubsysReco.h>
#include <string>
#include <map>
#include <vector>

class TFile;
class THnSparse;

class IsolationCut: public SubsysReco
{
public:
  IsolationCut( const char *filename = "isocut.root");
  virtual ~IsolationCut();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  /** event counter */
  unsigned _ievent;

  /** count events with 1+ photon candidate cluster */
  unsigned _events_photon;

  /** tower status for warnmap */
  int _tower_status[8][48][96];

  /** vector with PID's of neutral particles */
  std::vector< int > _v_pid_neutral;

  /** histogram storing cone energy information */
  THnSparse*  _hn_energy_cone;

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
