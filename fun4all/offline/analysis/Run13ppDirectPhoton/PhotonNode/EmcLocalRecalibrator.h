#ifndef __EMC_LOCAL_RECALIBRATOR_H__
#define __EMC_LOCAL_RECALIBRATOR_H__

#include <string>

class PhotonContainer;
class Photon;
class TF1;

/**
 * Class to access calibration correction for ECal from local files.
 *
 */
class EmcLocalRecalibrator
{
  public:

    /**
     * Constructor
     */
    EmcLocalRecalibrator();

    /**
     * Destructor
     */
    ~EmcLocalRecalibrator() {}

    /**
     * Correct data for all cluster in cluster container
     */
    void ApplyClusterCorrection( PhotonContainer* photoncont );

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

    void SelectMB();
    void SelectERT();

  protected:
    /**
     * Get energy corrected for non-linearity and run-by-run calibration for each sector
     */
    double GetCorrectedEcore( const Photon *photon );

    /**
     * Get TOF corrected for tower-by-tower calibration
     */
    double GetCorrectedTof( const Photon *photon );

    enum DataType {MB, ERT};
    DataType datatype;

    std::string _file_tofmap;
    std::string _file_energycalibration;

    double _tofmap[25000];
    double _energycalibration[8];

    TF1* _pbsc_cor_func;
    TF1* _pbgl_cor_func;
    TF1* _pbsc_recor_func;
};

#endif /* __EMC_LOCAL_RECALIBRATOR_H__ */
