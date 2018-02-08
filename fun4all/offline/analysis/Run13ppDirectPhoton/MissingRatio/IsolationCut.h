#ifndef __ISOLATIONCUT_H__
#define __ISOLATIONCUT_H__

#include <SubsysReco.h>
#include <string>
#include <map>

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

  /** histogram storing cone energy information */
  THnSparse*  _hn_energy_cone;

  /** output file name */
  std::string _output_file_name;

  /** output file */
  TFile *_file_output;

};

#endif /* __ISOLATIONCUT_H__ */
