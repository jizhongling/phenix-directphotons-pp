#include "PhotonContainer.h"
#include <TClass.h>

ClassImp(PhotonContainer)

void PhotonContainer::Reset()
{
  Class()->IgnoreTObjectStreamer();
  photon_list.clear();

  bbc_z = -9999.;
  bbc_t0 = -9999.;
  //crossing = -9999;
  trig = 0;
}

void PhotonContainer::set_trigger(unsigned lvl1_live, unsigned lvl1_scaled)
{
  const unsigned bit_bbcwide = 0x00000001;
  const unsigned bit_bbcnovtx = 0x00000002;
  const unsigned bit_bbcnarrow = 0x00000010;

  const unsigned bit_4x4b = 0x00000040;
  const unsigned bit_4x4a = 0x00000080;
  const unsigned bit_4x4c = 0x00000100;

  trig = 0;

  if( lvl1_live & bit_4x4a )
    trig += 0x01;
  if( lvl1_live & bit_4x4b )
    trig += 0x02;
  if( lvl1_live & bit_4x4c )
    trig += 0x04;

  if( lvl1_scaled & bit_4x4a )
    trig += 0x10;
  if( lvl1_scaled & bit_4x4b )
    trig += 0x20;
  if( lvl1_scaled & bit_4x4c )
    trig += 0x40;

  if( lvl1_live & bit_bbcnovtx )
    trig += 0x0100;
  if( lvl1_live & bit_bbcwide )
    trig += 0x0200;
  if( lvl1_live & bit_bbcnarrow )
    trig += 0x0400;

  if( lvl1_scaled & bit_bbcnovtx )
    trig += 0x1000;
  if( lvl1_scaled & bit_bbcwide )
    trig += 0x2000;
  if( lvl1_scaled & bit_bbcnarrow )
    trig += 0x4000;
}
