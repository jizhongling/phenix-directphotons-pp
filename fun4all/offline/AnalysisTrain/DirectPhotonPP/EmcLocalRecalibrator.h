#ifndef __EMCLOCALRECALIBRATOR_H__
#define __EMCLOCALRECALIBRATOR_H__

#include "AnaToolsTowerID.h"

#include <TOAD.h>
#include <getClass.h>

#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

class emcClusterContainer;
class emcClusterContent;
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
     * Setup
     */
    void Setup();

    /**
     * Correct data for all cluster in cluster container
     */
    void ApplyClusterCorrection( emcClusterContainer* data_emccontainer );

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
    double GetCorrectedEcore( const emcClusterContent *cluster );

    /**
     * Get TOF corrected for tower-by-tower calibration
     */
    double GetCorrectedTof( const emcClusterContent *cluster );

    std::string _file_warnmap;
    std::string _file_tofmap;
    std::string _file_energycalibration;

    double _tofmap[25000];
    double _energycalibration[8];

    TF1* _pbsc_cor_func;
    TF1* _pbgl_cor_func;
};

#endif /* __EMCLOCALRECALIBRATOR_H__ */
