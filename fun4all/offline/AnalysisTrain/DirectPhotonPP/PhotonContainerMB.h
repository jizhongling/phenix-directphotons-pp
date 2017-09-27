#ifndef __PHOTONCONTAINERMB_H__
#define __PHOTONCONTAINERMB_H__

#include "PhotonMB.h"
#include <PHObject.h>
#include <vector>

class PhotonContainerMB: public PHObject
{
  public:
    PhotonContainerMB() { Reset(); }
    virtual ~PhotonContainerMB() {}

    void AddPhoton(PhotonMB &photon) { photon_list.push_back(photon); }
    PhotonMB* GetPhoton(unsigned i) const { return const_cast<PhotonMB*>( &photon_list.at(i) ); }
    unsigned Size() const { return photon_list.size(); }
    void Reset(); 

    float get_bbc_t0() const { return bbc_t0; }
    short get_crossing() const { return crossing; }
    bool get_ert_a_live() const { return ( trig & 0x0001 ); }
    bool get_ert_b_live() const { return ( trig & 0x0002 ); }
    bool get_ert_c_live() const { return ( trig & 0x0004 ); }
    bool get_ert_a_scaled() const { return ( trig & 0x0010 ); }
    bool get_ert_b_scaled() const { return ( trig & 0x0020 ); }
    bool get_ert_c_scaled() const { return ( trig & 0x0040 ); }
    bool get_bbcnovtx_live() const { return ( trig & 0x0100 ); }
    bool get_bbcwide_live() const { return ( trig & 0x0200 ); }
    bool get_bbcnarrow_live() const { return ( trig & 0x0400 ); }
    bool get_bbcnovtx_scaled() const { return ( trig & 0x1000 ); }
    bool get_bbcwide_scaled() const { return ( trig & 0x2000 ); }
    bool get_bbcnarrow_scaled() const { return ( trig & 0x4000 ); }

    void set_bbc_t0(float a_bbc_t0) { bbc_t0 = a_bbc_t0; }
    void set_crossing(short a_crossing) { crossing = a_crossing; }
    void set_trigger(unsigned lvl1_live, unsigned lvl1_scaled);
    
  protected:
    std::vector<PhotonMB> photon_list;

    float bbc_t0;
    short crossing;
    unsigned short trig;

    ClassDef(PhotonContainerMB, 1)
};

#endif /* __PHOTONCONTAINERMB_H__ */
