#ifndef __EMCWARNMAPCHECKER_H__
#define __EMCWARNMAPCHECKER_H__

class emcClusterContent;

class EMCWarnmapChecker
{
  public:
    EMCWarnmapChecker();

    bool InFiducial(const emcClusterContent *cluster);
    bool IsGoodTower(const emcClusterContent *cluster);
    bool IsBadTower(const emcClusterContent *cluster);

    bool InFiducial(int itower);
    bool IsGoodTower(int itower);
    bool IsBadTower(int itower);

    bool PassCut(const emcClusterContent *cluster);

    int GetStatusNils(int sector, int iypos, int izpos);
    int GetStatusNils(const emcClusterContent *cluster);
    int GetStatusSasha(const emcClusterContent *cluster);

  protected:
    /* Number of warnmap array */
    static const int NSEC = 8;
    static const int NY = 48;
    static const int NZ = 96;
    static const int n_twrs = 24768;

    int tower_status_nils[NSEC][NY][NZ];
    int tower_status_sasha[NSEC][NY][NZ];
    int tower_status_inseok[NSEC][NY][NZ];

    void ReadWarnmapNils();
    void ReadWarnmapSasha();
    void ReadWarnmapInseok();
};

#endif /* __EMCWARNMAPCHECKER_H__ */
