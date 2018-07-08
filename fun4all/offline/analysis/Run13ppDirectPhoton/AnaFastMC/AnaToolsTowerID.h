#ifndef __ANATOOLSTOWERID_H__
#define __ANATOOLSTOWERID_H__

#include <TMCParticle.h>
#include <iostream>
#include <cmath>

namespace anatools
{
  inline bool Edge_cg(int sec, int iy, int iz)
  {
    // Check if shower center of gravity (maximal energy tower) belongs 
    // to edge tower defined by parameter "ncut"
    const int ncut = 2; // Number of edge towers

    if( sec < 0 || sec > 7 ) return true;
    if( sec<6 ) {
      if( iy<ncut || iy>35-ncut ) return true;
      if( iz<ncut || iz>71-ncut ) return true;
    }
    else {
      if( iy<ncut || iy>47-ncut ) return true;
      if( iz<ncut || iz>95-ncut ) return true;
    }
    return false;
  }

  /*
   * Check whether the tower is in fiducial arm
   */
  inline bool ArmEdge_cg(int sec, int iy, int iz)
  {
    // Check if shower center of gravity (maximal energy tower) belongs 
    // to edge tower defined by parameter "ncut"
    const int ncut = 10; // Number of edge towers

    if( sec < 0 || sec > 7 ) return true;
    if( sec < 6 )
    {
      if( iz < ncut || iz > 71-ncut ) return true;
      if( sec == 0 && iy < ncut ) return true;
      if( ( sec == 3 || sec == 4 ) && iy > 35-ncut ) return true;
    }
    else
    {
      if( iz < ncut || iz > 95-ncut ) return true;
      if( sec == 7 && iy < ncut ) return true;
    }
    return false;
  }

  inline unsigned int TowerID( const int sector, const int ytower, const int ztower )
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
  inline void TowerLocation( const unsigned int towerid, int &sector, int &ytower, int &ztower )
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

  /*
   * Check whether the two sectors are in the same part
   */
  inline bool SectorCheck(int pi0_sec1, int pi0_sec2)
  {
    if( pi0_sec1<4 && pi0_sec2<4 ) // W0,1,2,3
      return true;
    else if( (pi0_sec1==4 || pi0_sec1==5) && (pi0_sec2==4 || pi0_sec2==5) ) // E2,3
      return true;
    else if( (pi0_sec1==6 || pi0_sec1==7) && (pi0_sec2==6 || pi0_sec2==7) ) // PbGl
      return true;
    else
      return false;
  }

  inline float GetPt(TMCParticle *part)
  {
    float px = part->GetPx();
    float py = part->GetPy();
    float pt = std::sqrt(px*px + py*py);

    return pt;
  }
}

# endif /* __ANATOOLSTOWERID_H__ */
