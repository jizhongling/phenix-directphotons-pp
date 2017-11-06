#ifndef __ANAPHPYTHIADIRECTPHOTON_H__
#define __ANAPHPYTHIADIRECTPHOTON_H__

#include "SubsysReco.h"

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;

class AnaPHPythiaDirectPhoton: public SubsysReco
{
public:
  AnaPHPythiaDirectPhoton(const std::string &name = "AnaPHPythiaDirectPhoton");
  virtual ~AnaPHPythiaDirectPhoton();

  // Methods Derived from SubsysReco
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  PHPythiaHeader *phpythiaheader;
  PHPythiaContainer *phpythia;

};

#endif  /* __ANAPHPYTHIADIRECTPHOTON_H__ */

