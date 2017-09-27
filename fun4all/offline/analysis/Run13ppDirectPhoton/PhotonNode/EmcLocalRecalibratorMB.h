#ifndef __EMC_LOCAL_RECALIBRATORMB_H__
#define __EMC_LOCAL_RECALIBRATORMB_H__

#include <string>

class PhotonContainerMB;
class PhotonMB;
class TF1;

/**
 * Class to access calibration correction for ECal from local files.
 *
 */
class EmcLocalRecalibratorMB
{
  public:

    /**
     * Constructor
     */
    EmcLocalRecalibratorMB();

    /**
     * Destructor
     */
    ~EmcLocalRecalibratorMB() {}

    /**
     * Correct data for all cluster in cluster container
     */
    void ApplyClusterCorrection( PhotonContainerMB* photoncont );

    /**
     * Read
     */
    void ReadEnergyCorrection(const int& a_runnumber);

    /**
     * Read
     */
    void ReadTofCorrection( const int& a_fillnumber );

    /**
     * Set source file
     */
    void SetEnergyCorrectionFile( const std::string a_energycalibrationfile )
    {
      _file_energycalibration = a_energycalibrationfile;
    }

    /**
     * Set source file
     */
    void SetTofCorrectionFile( const std::string a_tofmapfile )
    {
      _file_tofmap = a_tofmapfile;
    }

  protected:
    /**
     * Get energy corrected for non-linearity and run-by-run calibration for each sector
     */
    double GetCorrectedEcore( const PhotonMB *photon );

    /**
     * Get TOF corrected for tower-by-tower calibration
     */
    double GetCorrectedTof( const PhotonMB *photon );

    std::string _file_tofmap;
    std::string _file_energycalibration;

    double _tofmap[25000];
    double _energycalibration[8];

    TF1* _pbsc_cor_func;
    TF1* _pbgl_cor_func;
    TF1* _pbsc_recor_func;
};

#endif /* __EMC_LOCAL_RECALIBRATORMB_H__ */
