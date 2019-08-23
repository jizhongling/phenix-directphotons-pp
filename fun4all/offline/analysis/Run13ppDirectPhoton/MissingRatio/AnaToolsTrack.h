#ifndef __ANATOOLSTRACK_H__
#define __ANATOOLSTRACK_H__

#include <emcGeaTrackContent.h>
#include <emcClusterContent.h>

#include <TMath.h>
#include <TLorentzVector.h>

namespace anatools
{
  const double PI = TMath::Pi();
  const double epsilon = TMath::Limits<float>::Epsilon();

  /* Calculate difference vector between cluster and truth */
  inline TVector3 Get_dR(const emcClusterContent* cluster, const emcGeaTrackContent* track)
  {
    double dx = cluster->x() - track->get_impx();
    double dy = cluster->y() - track->get_impy();
    double dz = cluster->z() - track->get_impz();
    return TVector3(dx, dy, dz);
  }

  /* Calculate momentum and energy of single track */
  inline TLorentzVector Get_pE(const emcGeaTrackContent* track)
  {
    double E = track->get_ekin();
    double px = track->get_px();
    double py = track->get_py();
    double pz = track->get_pz();

    TLorentzVector pE(px, py, pz, E);
    return pE;
  }

  /* Calculate total transverse momentum pT of track pair */
  inline double GetTot_pT(const emcGeaTrackContent* track1, const emcGeaTrackContent* track2)
  {
    TLorentzVector pE1 = Get_pE(track1);
    TLorentzVector pE2 = Get_pE(track2);
    return (pE1+pE2).Pt();
  }

  /* Calculate pseudorapidity of a single track */
  inline double GetEta(const emcGeaTrackContent* track)
  {
    TLorentzVector pE = Get_pE(track);
    return pE.Eta();
  }

  /* Calculate phi angle of a single track */
  inline double GetPhi(const emcGeaTrackContent* track)
  {
    double px = track->get_px();
    double py = track->get_py();
    if( fabs(px) < epsilon ) return -9999.;
    double phi = px > 0. ? atan(py/px) : PI + atan(py/px);
    return phi;
  }

  /* Calculate EMCal arm of a single track */
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

  /* Calculate invariant mass of track pair */
  inline double GetInvMass(const emcGeaTrackContent* track1, const emcGeaTrackContent* track2)
  {
    TLorentzVector pE1 = Get_pE(track1);
    TLorentzVector pE2 = Get_pE(track2);
    return (pE1+pE2).M();
  }
} /* namespace anatools */

#endif /* __ANATOOLSTRACK_H__ */
