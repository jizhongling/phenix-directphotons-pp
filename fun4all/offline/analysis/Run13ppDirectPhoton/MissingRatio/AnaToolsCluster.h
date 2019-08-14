#ifndef __ANATOOLSCLUSTER_H__
#define __ANATOOLSCLUSTER_H__

#include <emcClusterContent.h>
#include <emcGeaTrackContent.h>

#include <TMath.h>
#include <TLorentzVector.h>

namespace anatools
{
  const double PI = TMath::Pi();

  // get sector
  inline int GetSector(const emcClusterContent *cluster)
  {
    int arm = cluster->arm();
    int rawsector = cluster->sector();
    int sector = arm==0 ? rawsector : 7-rawsector;
    return sector;
  }

  // calculate difference vector between cluster and truth
  inline TVector3 Get_dR(const emcClusterContent* cluster, const emcGeaTrackContent* track)
  {
    double dx = cluster->x() - track->get_impx();
    double dy = cluster->y() - track->get_impy();
    double dz = cluster->z() - track->get_impz();
    return TVector3(dx, dy, dz);
  }

  // calculate momentum and energy of single cluster (=photon)
  inline TLorentzVector Get_pE(const emcClusterContent* cluster)
  {
    double x = cluster->x();
    double y = cluster->y();
    double z = cluster->z();
    double length = sqrt(x*x + y*y + z*z);

    double E = cluster->ecore();
    double px = E*x/length;
    double py = E*y/length;
    double pz = E*z/length;

    TLorentzVector pE(px, py, pz, E);
    return pE;
  }

  // calculate momentum and energy of single track
  inline TLorentzVector Get_pE(const emcGeaTrackContent* track)
  {
    //double x = track->get_impx();
    //double y = track->get_impy();
    //double z = track->get_impz();
    //double length = sqrt(x*x + y*y + z*z);

    //double E = track->get_edep();
    //double px = E*x/length;
    //double py = E*y/length;
    //double pz = E*z/length;

    double E = track->get_ekin();
    double px = track->get_px();
    double py = track->get_py();
    double pz = track->get_pz();

    TLorentzVector pE(px, py, pz, E);
    return pE;
  }

  // calculate transverse momentum pT of single cluster (=photon)
  inline double Get_pT(const emcClusterContent* cluster)
  {
    TLorentzVector pE = Get_pE(cluster);
    return pE.Pt();
  }

  // calculate total transverse momentum pT of cluster pair
  inline double GetTot_pT(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    TLorentzVector pE1 = Get_pE(cluster1);
    TLorentzVector pE2 = Get_pE(cluster2);
    return (pE1+pE2).Pt();
  }

  // calculate total transverse momentum pT of track pair
  inline double GetTot_pT(const emcGeaTrackContent* track1, const emcGeaTrackContent* track2)
  {
    TLorentzVector pE1 = Get_pE(track1);
    TLorentzVector pE2 = Get_pE(track2);
    return (pE1+pE2).Pt();
  }

  // calculate pseudorapidity of a single track
  inline double GetEta(const emcGeaTrackContent* track)
  {
    TLorentzVector pE = Get_pE(track);
    return pE.Eta();
  }

  // calculate phi angle of a single track
  inline double GetPhi(const emcGeaTrackContent* track)
  {
    double px = track->get_px();
    double py = track->get_py();
    if(px==0.) return -9999.;
    double phi = px > 0. ? atan(py/px) : PI + atan(py/px);
    return phi;
  }

  // calculate EMCal arm of a single track
  inline int GetTrackArm(const emcGeaTrackContent* track)
  {
    double phi = GetPhi(track);
    if( phi > -PI*3/16 && phi < PI*5/16 )
      return 0;
    else if( phi > PI*11/16 && phi < PI*19/16 )
      return 1;
    else
      return -9999;
  }

  // calculate invariant mass of cluster pair
  inline double GetInvMass(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    TLorentzVector pE1 = Get_pE(cluster1);
    TLorentzVector pE2 = Get_pE(cluster2);
    return (pE1+pE2).M();
  }

  // calculate invariant mass of track pair
  inline double GetInvMass(const emcGeaTrackContent* track1, const emcGeaTrackContent* track2)
  {
    TLorentzVector pE1 = Get_pE(track1);
    TLorentzVector pE2 = Get_pE(track2);
    return (pE1+pE2).M();
  }

  // calculate energy asymmetry of cluster pair
  inline double GetAsymmetry_E(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    double E1 = cluster1->ecore();
    double E2 = cluster2->ecore();
    double asymE = std::fabs(E1 - E2) / std::fabs(E1 + E2);
    return asymE;
  }
} /* namespace anatools */

#endif /* __ANATOOLSCLUSTER_H__ */
