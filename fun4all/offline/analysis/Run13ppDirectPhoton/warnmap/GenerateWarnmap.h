#ifndef _GENERATE_WARNMAP_H_
#define _GENERATE_WARNMAP_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2I.h"

namespace direct_photon_pp
{

  class GenerateWarnmap
  {

  public:

    /**
     * Default constructor
     */
    GenerateWarnmap( unsigned nsigma , std::string output_file_plots );

    /**
     * Default destructor
     */
    virtual ~GenerateWarnmap(){}

    /**
     * Fill hits from input histogram
     */
    int FillHitsFromHistogram( std::string input_file , std::string input_histogram );

    /**
     * Fill hits from input 2D histogram
     */
    int FillHitsFrom2DHistogram( std::string input_file , std::string input_histogram , unsigned ybin_min , unsigned ybin_max );


    /**
     * Read list of uncalibrated towers
     */
    int ReadUncalibratedTowers( std::string input_file );


    /**
     * Find dead towers which are towers below minimum hit frequency threshold
     */
    int FindDeadTowers( float threshold_alive = 1 );

    /**
     * Determine hot and dead channels; nsigma = cut for hot/dead in multiples of standard deviations
     */
    int FindHotTowers( );

    /**
     * Repeat determine hot towers for niterations = number of iterations for excluding hot/dead towers
     */
    int FindHotTowers( int niterations );

    /**
     * Set towers around hot towers to appropriate status for fiducial volume cut
     */
    int FiducialCutHotTowers();

    /**
     * Set towers on the edge of sectors to appropriate status for fiducial volume cut
     */
    int FiducialCutSectorEdges();

    /**
     * Write Warnmap to textfile
     */
    int WriteWarnmap( std::string filename );

    void GeneratePlots();

    void Statistics();
    void Finish();

  protected:

    static const double ntowers_ = 24768;
    static const double ntowers_pbsc_ = 2592;
    static const double ntowers_pbgl_ = 4608;

    static const unsigned n_sector_ = 8;
    static const unsigned n_tower_z_ = 96;
    static const unsigned n_tower_y_ = 48;

    /**
     * Map of number of hits in each tower
     */
    long hitmap_[n_sector_][n_tower_z_][n_tower_y_];

    /**
     * Status codes for warnmap:
     * GOOD = good tower
     * FIDUCIAL_ARM = outside fiducial volume (edge of arm, still useful for finding partner photons)
     * FIDUCIAL_EDGE = outside fiducial volume (edge of sector)
     * FIDUCIAL_HOT = outside fiducial volume (neighbor of hot tower)
     * HOT = hot tower
     * DEAD = dead tower
     * UNCALIBRATED = uncalibrated tower
     */
    enum status_tower {
      GOOD = 0,
      FIDUCIAL_ARM = 10,
      FIDUCIAL_SECTOR = 20,
      FIDUCIAL_HOT = 40,
      HOT = 50,
      DEAD = 100,
      UNCALIBRATED = 150
    };

    /**
     * Warnmap array with status code for each tower
     */
    status_tower warnmap_[n_sector_][n_tower_z_][n_tower_y_];

    // TGraph "Threshold" vs "Iteration" for hot tower location
    std::vector< unsigned > v_iteration_;
    std::vector< std::vector< long long > > v_thresholds_;
    std::vector< std::vector< float > > v_means_;
    std::vector< std::vector< float > > v_sdevs_;
    std::vector< std::vector< unsigned > > v_newhot_;
    std::vector< std::vector< unsigned > > v_allhot_;

    unsigned nsigma_;
    bool generate_plots_;
    TFile *fout_;

  };

} // end namespace

#endif /* __GENERATEWARNMAP_H__ */
