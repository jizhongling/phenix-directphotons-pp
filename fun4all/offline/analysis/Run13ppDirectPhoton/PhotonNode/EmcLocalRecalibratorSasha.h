#ifndef __EMC_LOCAL_RECALIBRATOR_SASHA_H__
#define __EMC_LOCAL_RECALIBRATOR_SASHA_H__

class PhotonContainer;
class Photon;

#define NMAXTWR 24768

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
     * Correct data for all cluster in cluster container
     */
    void ApplyClusterCorrection( const int runno, PhotonContainer* photoncont );
  
    void anaGetCorrTof(const char* fname);

  protected:
    float GetTowerTofCorr(int run, int id, float en);

    float fCorrTof[NMAXTWR];
};

#endif /* __EMC_LOCAL_RECALIBRATOR_SASHA_H__ */
