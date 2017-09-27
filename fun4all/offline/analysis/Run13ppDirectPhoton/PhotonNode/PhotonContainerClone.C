#include "PhotonContainerClone.h"

PhotonContainerClone::PhotonContainerClone(PhotonContainer *photoncont)
{
  Reset();

  for(unsigned i=0; i<photoncont->Size(); i++)
  {
    Photon *photon = photoncont->GetPhoton(i);
    photon_list.push_back(*photon);
  }

  bbc_z = photoncont->get_bbc_z();
  bbc_t0 = photoncont->get_bbc_t0();
  crossing = photoncont->get_crossing();

  if( photoncont->get_ert_a_live() )
    trig += 0x01;
  if( photoncont->get_ert_b_live() )
    trig += 0x02;
  if( photoncont->get_ert_c_live() )
    trig += 0x04;

  if( photoncont->get_ert_a_scaled() )
    trig += 0x10;
  if( photoncont->get_ert_b_scaled() )
    trig += 0x20;
  if( photoncont->get_ert_c_scaled() )
    trig += 0x40;

  if( photoncont->get_bbcnovtx_live() )
    trig += 0x0100;
  if( photoncont->get_bbcwide_live() )
    trig += 0x0200;
  if( photoncont->get_bbcnarrow_live() )
    trig += 0x0400;

  if( photoncont->get_bbcnovtx_scaled() )
    trig += 0x1000;
  if( photoncont->get_bbcwide_scaled() )
    trig += 0x2000;
  if( photoncont->get_bbcnarrow_scaled() )
    trig += 0x4000;

}
