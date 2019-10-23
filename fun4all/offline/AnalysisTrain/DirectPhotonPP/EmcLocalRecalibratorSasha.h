#ifndef __EMCLOCALRECALIBRATORSASHA_H__
#define __EMCLOCALRECALIBRATORSASHA_H__

class emcClusterContainer;
class emcClusterContent;

#define NMAXTWR 24768
#define NRUN 2000

/**
 * Class to access calibration correction for ECal from local files.
 *
 */
class EmcLocalRecalibratorSasha
{
  public:
    /**
     * Constructor
     */
    EmcLocalRecalibratorSasha();

    /**
     * Destructor
     */
    ~EmcLocalRecalibratorSasha() {}

    /**
     * Setup
     */
    void Setup();

    /**
     * Correct data for all cluster in cluster container
     */
    void ApplyClusterCorrection( const int runno, emcClusterContainer* data_emccontainer );

    void anaGetCorrCal(const char* fname);
    void anaGetCorrCal_run(const char* fname);
    void anaGetCorrTof(const char* fname);

  protected:
    int FindRun(int run, int* rlist, int nn);
    float GetTowerCorrCal(const emcClusterContent* emccl);
    float GetTowerTofCorr(int run, int id, float en);

    int fnRun_cal;
    int runlist_cal[NRUN];
    float fCorrCalRun[NRUN][8];

    float fCorrCal[NMAXTWR];
    float fCorrTof[NMAXTWR];
};

#endif /* __EMCLOCALRECALIBRATORSASHA_H__ */
