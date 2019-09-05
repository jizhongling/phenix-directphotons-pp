#ifndef __ERTSIMTRIGGER_H__
#define __ERTSIMTRIGGER_H__

class emcClusterContent;

class TRandom3;
class TFile;
class TTree;

class ERTSimTrigger
{
  public:
    ERTSimTrigger();
    virtual ~ERTSimTrigger();

    bool Fired(const emcClusterContent *cluster);

  protected:
    void ReadTriggerInfo();

    TRandom3 *rnd;
    TFile *f_ertsm;
    TTree *t_ertsm;
};

#endif /* __ERTSIMTRIGGER_H__ */
