#ifndef __DCDEADMAPCHECKER_H__
#define __DCDEADMAPCHECKER_H__

#include <string>
#include <vector>
#include <map>

class PHCentralTrack;
class emcClusterContent;

struct KBB
{
  double k;
  double b1;
  double b2;
  double alpha1;
  double alpha2;

  KBB(double _k = 0., double _b1 = 0., double _b2 = 0.,
      double _alpha1= 0., double _alpha2 = 0.):
    k(_k), b1(_b1), b2(_b2), alpha1(_alpha1), alpha2(_alpha2)
  {}
};

class DCDeadmapChecker
{
  public:
    DCDeadmapChecker(int eventsmod = 10000);

    void SetMapByIndex(int mapindex) { imap = mapindex; }
    void SetMapByRunnumber(int runnumber); 
    void SetMapByEvent(); 

    bool IsDead(std::string nswe, double board, double alpha);
    bool IsDead(const PHCentralTrack *tracks, int itrk);

    bool ChargeVeto(const emcClusterContent *cluster, const PHCentralTrack *tracks);
    int GetEmcMatchTrack(const emcClusterContent *cluster, const PHCentralTrack *tracks);

  protected:
    static const int nmap = 15;
    int imap;
    int nevents;
    int nmod;
    int cum_nRuns[nmap];

    std::map<std::string,KBB> m_kbb;
    std::vector<std::string> v_deadmap[nmap];
};

#endif /* __DCDEADMAPCHECKER_H__ */
