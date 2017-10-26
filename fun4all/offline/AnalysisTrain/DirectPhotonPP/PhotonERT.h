#ifndef __PHOTONERT_H__
#define __PHOTONERT_H__

#include "Photon.h"

class ErtOut;
class emcClusterContent;

class PhotonERT: public Photon
{
  public:
    PhotonERT();
    PhotonERT(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_theta_cv, float a_cone_energy);
    virtual ~PhotonERT() {}

    float get_theta_cv() const { return theta_cv; }
    float get_cone_energy() const { return cone_energy; }

    void set_theta_cv(float a_theta_cv) { theta_cv = a_theta_cv; }
    void set_cone_energy(float a_cone_energy) { cone_energy = a_cone_energy; }

  protected:
    float theta_cv;
    float cone_energy;

    ClassDef(PhotonERT, 1)
};

#endif /* __PHOTONERT_H__ */
