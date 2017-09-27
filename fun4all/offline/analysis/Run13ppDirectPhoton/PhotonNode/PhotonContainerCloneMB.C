#include "PhotonContainerCloneMB.h"

PhotonContainerCloneMB::PhotonContainerCloneMB(PhotonContainerMB *photoncont)
{
  Reset();

  for(unsigned i=0; i<photoncont->Size(); i++)
  {
    PhotonMB *photon = photoncont->GetPhoton(i);
    photon_list.push_back(*photon);
  }

  bbc_t0 = photoncont->get_bbc_t0();
  crossing = photoncont->get_crossing();

  if( photoncont->get_ert_a_live() )
    trig += 0x01;
  if( photoncont->get_ert_b_live() )
    trig += 0x02;
  if( photoncont->get_ert_a_live() )
    trig += 0x04;

  if( photoncont->get_ert_a_scaled() )
    trig += 0x10;
  if( photoncont->get_ert_b_scaled() )
    trig += 0x20;
  if( photoncont->get_ert_a_scaled() )
    trig += 0x40;
}
