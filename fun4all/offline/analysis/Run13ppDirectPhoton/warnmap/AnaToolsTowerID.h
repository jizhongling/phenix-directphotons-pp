#ifndef ANATOOLSTOWERID_H__
#define ANATOOLSTOWERID_H__

/*! Namespace with various functions for analysis.
 * This file provides functions to encode / decode EMCal tower ID's.
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
namespace anatools
{

  /*! Maximum tower id in EMCal
   */
  const unsigned int kMaxTowerID = 25000;

  /*! Return tower ID for given sector, yrow, and zrow location of tower
   *
   * EMCAL: 8 sectors (2 PbGl, 6 PbSc)
   *
   * PbSc: each sector has 72 * 36 towers
   *
   * PbGl: each sector has 96 * 48 towers
   *
   * Sector counting scheme: Sectors 0 - 5 are PbSc, sectors 6 - 7 are PbGl
   *
   */
  unsigned int TowerID( const int sector, const int ytower, const int ztower )
  {
    unsigned int id = 0;
    unsigned int itower = 0;

    if( sector < 6 ){ // PbSc

      if( ytower < 0 || ytower > 35 || ztower < 0 || ztower > 71 )
	{
	  std::cerr << "Exception: Bad PbSc tower location " << sector << " " << ztower << " " << ytower << std::endl;
	  throw 1;
	}

      itower = ytower * 72 + ztower;
      id = sector * 2592 + itower;
    }

    else{ // PbGl

      if( ytower < 0 || ytower > 47 || ztower < 0 || ztower > 95 )
	{
	  std::cerr << "Exception: Bad PbGl tower location " << sector << " " << ztower << " " << ytower << std::endl;
	  throw 1;
	}

      itower = ytower * 96 + ztower;
      id = 15552 + ( sector - 6 ) * 4608 + itower;
    }

    return id;

  }


  /*! calculate sector, yrow, and zrow location of tower for given tower ID
   */
  void TowerLocation( const unsigned int towerid, int &sector, int &ytower, int &ztower )
  {
    int itower=0;

    if( towerid < 15552 )
      { // PbSc
	sector = towerid / 2592;
	itower = towerid % 2592;
	ztower = itower % 72;
	ytower = itower / 72;
      }
    else
      { // PbGl
	sector = 6 + ( towerid - 15552 ) / 4608;
	itower = ( towerid - 15552 ) % 4608;
	ztower = itower % 96;
	ytower = itower / 96;
      }

    return;
  }



  /*! Correct counting of sectors for emcClusterContent, return corrected sector number
   */
  int CorrectClusterSector( const int arm , const int rawsector )
  {
    int sector = rawsector;

    if ( arm == 1 )
      sector = 7 - rawsector;

    return sector;
  }

};

#endif
