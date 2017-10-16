#include "PhotonERT.h"
#include "AnaToolsTrigger.h"
#include <TClass.h>

ClassImp(Photon)

PhotonERT::PhotonERT()
{
  Class()->IgnoreTObjectStreamer();

  theta_cv = -9999.;
  cone_energy = -9999.;
}

PhotonERT::PhotonERT(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_theta_cv, float a_cone_energy)
{
  towerid = a_towerid;

  x = a_x;
  y = a_y;
  z = a_z;
  E = a_E;

  theta_cv = a_theta_cv;
  cone_energy = a_cone_energy;

  trig = 0;
}
