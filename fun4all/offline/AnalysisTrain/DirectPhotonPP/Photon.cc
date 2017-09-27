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

  tof = -9999.;
  Ecorr = -9999.;
  //theta_cv = -9999.;
  //cone_energy = -9999.;

  trig = 0;
}

//Photon::Photon(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_tof, float a_theta_cv, float a_cone_energy)
Photon::Photon(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_tof, float a_Ecorr)
{
  towerid = a_towerid;

  x = a_x;
  y = a_y;
  z = a_z;
  E = a_E;

  tof = a_tof;
  Ecorr = a_Ecorr;
  //theta_cv = a_theta_cv;
  //cone_energy = a_cone_energy;

  trig = 0;
}

void Photon::set_trig(ErtOut *ertout, emcClusterContent *cluster)
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
