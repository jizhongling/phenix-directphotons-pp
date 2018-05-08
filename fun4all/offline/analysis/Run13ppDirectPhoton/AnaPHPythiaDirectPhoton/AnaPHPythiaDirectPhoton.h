#ifndef __ANAPHPYTHIADIRECTPHOTON_H__
#define __ANAPHPYTHIADIRECTPHOTON_H__

#include "SubsysReco.h"

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;

class TFile;
class TTree;

class AnaPHPythiaDirectPhoton: public SubsysReco
{
public:
  AnaPHPythiaDirectPhoton(const std::string &name = "AnaPHPythiaDirectPhoton");
  virtual ~AnaPHPythiaDirectPhoton();

  // Methods Derived from SubsysReco
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  PHPythiaHeader *phpythiaheader;
  PHPythiaContainer *phpythia;

  /** output tree with truth information */
  TTree* _tree_event_truth;

  /** truth tree variables */
  unsigned int _ievent;
  float _truth_pid;
  float _truth_parentid;
  float _truth_anclvl;
  float _truth_ptot;
  float _truth_pt;
  float _truth_eta;
  float _truth_phi;

  /** output file name */
  std::string _output_file_name;

  /** output file */
  TFile *_fout;

};

#endif  /* __ANAPHPYTHIADIRECTPHOTON_H__ */

