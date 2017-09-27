#ifndef __PHOTONMB_H__
#define __PHOTONMB_H__

#include <PHObject.h>

class ErtOut;
class emcClusterContent;

class PhotonMB: public PHObject
{
  public:
    PhotonMB();
    PhotonMB(short a_towerid, float a_x, float a_y, float a_z, float a_E, float a_tof);
    virtual ~PhotonMB() {}

    short get_towerid() const { return towerid; }
    float get_x() const { return x; }
    float get_y() const { return y; }
    float get_z() const { return z; }
    float get_E() const { return E; }
    float get_tof() const { return tof; }
    bool get_trg1() const { return ( trig & 0x1 ); }
    bool get_trg2() const { return ( trig & 0x2 ); }
    bool get_trg3() const { return ( trig & 0x4 ); }

    void set_towerid(short a_towerid) { towerid = a_towerid; }
    void set_x(float a_x) { x = a_x; }
    void set_y(float a_y) { y = a_y; }
    void set_z(float a_z) { z = a_z; }
    void set_E(float a_E) { E = a_E; }
    void set_tof(float a_tof) { tof = a_tof; }
    void set_trig(ErtOut *ertout, emcClusterContent *cluster);

  protected:
    short towerid;

    float x;
    float y;
    float z;
    float E;

    float tof;
    unsigned short trig;

    ClassDef(PhotonMB, 1)
};

#endif /* __PHOTONMB_H__ */
