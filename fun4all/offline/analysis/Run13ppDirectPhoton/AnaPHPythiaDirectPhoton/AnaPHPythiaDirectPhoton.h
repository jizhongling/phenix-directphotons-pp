#ifndef __ANAPHPYTHIADIRECTPHOTON_H__
#define __ANAPHPYTHIADIRECTPHOTON_H__

#include "SubsysReco.h"

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;

class TFile;
class TTree;
class TMCParticle;
class TDatabasePDG;

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

  /** sum up electromagnetic energy within cone of fixed opening
      angle rcone around particle with vector vref */
  float SumEEmcal( TMCParticle* pref, float rcone );

  /** sum up charged traack momenta within cone of fixed opening
      angle rcone around particle with vector vref */
  float SumPTrack( TMCParticle* pref, float rcone );

  /** ROOT databsae with PDG properties */
  TDatabasePDG* _pdg_db;

  PHPythiaHeader *phpythiaheader;
  PHPythiaContainer *phpythia;

  /** output tree with truth information */
  TTree* _tree_event_truth;

  /** truth tree variables */
  unsigned int _ievent;
  int _pythia_event;
  int _pythia_processid;
  float _truth_pid;
  float _truth_parentid;
  float _truth_ispromptphoton;
  float _truth_ptot;
  float _truth_pt;
  float _truth_etot;
  float _truth_eta;
  float _truth_phi;

  std::vector < float > _v_iso_conesize;
  std::vector < float > _v_iso_eemcal;
  std::vector < float > _v_iso_ptrack;

  /** output file name */
  std::string _output_file_name;

  /** output file */
  TFile *_fout;

};

#endif  /* __ANAPHPYTHIADIRECTPHOTON_H__ */

