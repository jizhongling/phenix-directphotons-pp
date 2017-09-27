#include "PhotonMB.h"
#include "AnaToolsTrigger.h"
#include <TClass.h>

ClassImp(PhotonMB)

PhotonMB::PhotonMB()
{
  Class()->IgnoreTObjectStreamer();

  towerid = -9999;

  x = -9999.;
  y = -9999.;
  z = -9999.;
  E = -9999.;

  tof = -9999.;
  trig = 0;
}

PhotonMB::PhotonMB(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_tof)
{
  towerid = a_towerid;

  x = a_x;
  y = a_y;
  z = a_z;
  E = a_E;

  tof = a_tof;
  trig = 0;
}

void PhotonMB::set_trig(ErtOut *ertout, emcClusterContent *cluster)
{
  trig = 0;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4a) )
    trig += 0x1;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4b) )
    trig += 0x2;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4c) )
    trig += 0x4;

  return;
}
