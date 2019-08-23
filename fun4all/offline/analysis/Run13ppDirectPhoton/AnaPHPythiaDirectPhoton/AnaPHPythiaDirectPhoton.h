#ifndef __ANAPHPYTHIADIRECTPHOTON_H__
#define __ANAPHPYTHIADIRECTPHOTON_H__

#include "SubsysReco.h"

#include <map>

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

  /** Turn on/off photon and PHENIX acceptance cut */
  void SelectPhotons( bool yesno ) { _select_photons = yesno; }

protected:

  /** sum up electromagnetic energy within cone of fixed opening
      angle rcone around particle with vector vref */
  float SumEEmcal( TMCParticle* pref, float rcone );

  /** sum up charged traack momenta within cone of fixed opening
      angle rcone around particle with vector vref */
  float SumPTrack( TMCParticle* pref, float rcone );

  /** Reset global variables that store cluster information for filling output tree */
  void ResetBranchVariables();

  /** boolean to turn on/off photon and acceptance selection */
  bool _select_photons;

  /** ROOT databsae with PDG properties */
  TDatabasePDG* _pdg_db;

  PHPythiaHeader *phpythiaheader;
  PHPythiaContainer *phpythia;

  /** output tree with truth information */
  TTree* _tree_event_truth;

  /** Map of Event properties that will be written to
   * output ROOT Tree */
  std::map< std::string , float > _branchmap_event;

  /** Map of Particle (or cluster) properties that will be written to
   * output ROOT Tree */
  std::map< std::string , std::vector< float > > _branchmap_mcparticles;

  unsigned int _ievent;

  std::vector < float > _v_iso_conesize;
  std::vector < float > _v_iso_eemcal;
  std::vector < float > _v_iso_ptrack;

  /** output file name */
  std::string _output_file_name;

  /** output file */
  TFile *_fout;

};

#endif  /* __ANAPHPYTHIADIRECTPHOTON_H__ */
