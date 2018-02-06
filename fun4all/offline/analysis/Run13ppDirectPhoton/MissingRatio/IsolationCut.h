#ifndef __ISOLATIONCUT_H__
#define __ISOLATIONCUT_H__

#include <SubsysReco.h>
#include <string>
#include <map>

class IsolationCut: public SubsysReco
{
  public:
    IsolationCut( const char *filename = "isocut.root");
    virtual ~IsolationCut();

    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

  protected:

    std::string outFileName;
};

#endif /* __ISOLATIONCUT_H__ */
