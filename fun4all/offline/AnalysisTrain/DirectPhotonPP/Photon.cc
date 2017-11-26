#include "Photon.h"
#include "AnaToolsTrigger.h"
#include <TClass.h>

ClassImp(Photon)

Photon::Photon()
{
  Class()->IgnoreTObjectStreamer();

  towerid = -9999;

  x = -9999.;
  y = -9999.;
  z = -9999.;
  E = -9999.;
  Ecorr = -9999.;

  tof = -9999.;
  trig = 0;
}

Photon::Photon(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_Ecorr, float a_tof)
{
  towerid = a_towerid;

  x = a_x;
  y = a_y;
  z = a_z;
  E = a_E;
  Ecorr = a_Ecorr;

  tof = a_tof;
  trig = 0;
}

void Photon::set_trig(ErtOut *ertout, emcClusterContent *cluster)
{
  trig &= 0x8888;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4a) )
    trig |= 0x0001;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4b) )
    trig |= 0x0002;
  if( anatools::PassERT(ertout, cluster, anatools::ERT_4x4c) )
    trig |= 0x0004;

  return;
}

void Photon::set_prob(bool is_prob)
{
  if( is_prob )
    trig |= 0x0008;
  else
    trig &= 0xFFF7;

  return;
}
