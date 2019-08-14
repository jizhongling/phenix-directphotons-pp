#ifndef __ANATOOLSPHOTON_H__
#define __ANATOOLSPHOTON_H__

#include <Photon.h>
#include <TLorentzVector.h>
#include <cmath>

namespace anatools
{
  inline void TowerLocation(const Photon *photon, int &sector, int &ytower, int &ztower)
  {
    int towerid = photon->get_towerid();
    TowerLocation(towerid, sector, ytower, ztower);

    return;
  }

  inline int GetSector(const Photon *photon)
  {
    int towerid = photon->get_towerid();
    int sector, ytower, ztower;
    TowerLocation(towerid, sector, ytower, ztower);

    return sector;
  }

  inline TLorentzVector Get_pE(const Photon *photon)
  {
    double x = photon->get_x();
    double y = photon->get_y();
    double z = photon->get_z();
    double length = sqrt(x*x + y*y + z*z);

    double E = photon->get_E();
    double px = E*x/length;
    double py = E*y/length;
    double pz = E*z/length;

    TLorentzVector pE(px, py, pz, E);
    return pE;
  }

  inline double Get_pT(const Photon *photon)
  {
    TLorentzVector pE = Get_pE(photon);
    return pE.Pt();
  }

  inline double GetTot_pT(const Photon *photon1, const Photon *photon2)
  {
    TLorentzVector pE1 = Get_pE(photon1);
    TLorentzVector pE2 = Get_pE(photon2);
    return (pE1+pE2).Pt();
  }

  inline double GetInvMass(const Photon *photon1, const Photon *photon2)
  {
    TLorentzVector pE1 = Get_pE(photon1);
    TLorentzVector pE2 = Get_pE(photon2);
    return (pE1+pE2).M();
  }

  inline double GetAsymmetry_E(const Photon *photon1, const Photon *photon2)
  {
    double E1 = photon1->get_E();
    double E2 = photon2->get_E();
    double asymE = std::fabs(E1 - E2) / std::fabs(E1 + E2);
    return asymE;
  }
}

#endif /* __ANATOOLSPHOTON_H__ */
