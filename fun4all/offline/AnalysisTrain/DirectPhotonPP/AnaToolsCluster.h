/*
 * AnaToolsCluster.h
 *
 *  Created on: Jan 14, 2016
 *      Author: nfeege
 */

#ifndef OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSCLUSTER_H_
#define OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSCLUSTER_H_

#include <emcClusterContent.h>

#include <TVector3.h>
#include <TLorentzVector.h>

#include <cmath>

/*! Namespace with various functions for analysis.
 * This file provides functions to calculate values from cluster information
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
namespace anatools
{
  /*!
   * Get the cluster sector
   */
  inline int GetSector(const emcClusterContent* cluster)
  {
    int arm = cluster->arm();
    int sector = cluster->sector();
    if( arm == 1 )
      sector = 7 - sector;

    return sector;
  }

  /*!
   * Calculate momentum and energy of single cluster (=photon)
   */
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

  /*!
   * Calculate transverse momentum pT of single cluster (=photon)
   */
  inline double Get_pT(const emcClusterContent* cluster)
  {
    TLorentzVector pE = Get_pE(cluster);
    return pE.Pt();
  }

  /*!
   * Calculate total transverse momentum pT of cluster pair
   */
  inline double GetTot_pT(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    TLorentzVector pE1 = Get_pE(cluster1);
    TLorentzVector pE2 = Get_pE(cluster2);
    return (pE1+pE2).Pt();
  }

  /*!
   * Calculate energy asymmetry of cluster pair
   */
  inline double GetAsymmetry_E(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    double E1 = cluster1->ecore();
    double E2 = cluster2->ecore();

    double asym_E = std::abs( E1 - E2 ) / ( E1 + E2 );

    return asym_E;
  }

  /*!
   * Calculate invariant mass of cluster pair
   */
  inline double GetInvMass(const emcClusterContent* cluster1, const emcClusterContent* cluster2)
  {
    TLorentzVector pE1 = Get_pE(cluster1);
    TLorentzVector pE2 = Get_pE(cluster2);
    return (pE1+pE2).M();
  }

  /*!
   * Calculate angle between EMC cluster and PC3 track
   */
  inline double GetTheta_CV(const emcClusterContent *cluster)
  {
    //coordinate of emc hit
    double emc_phi = cluster->phi();
    double emc_z = cluster->z();

    int sector = anatools::CorrectClusterSector( cluster->arm() , cluster->sector() );
    double emc_r;
    if(sector>5) emc_r = 527.7;
    else emc_r = 507.7;

    //coodinate of pc3 hit                                           
    double pc3_dphi = cluster->emcpc3dphi();
    double pc3_dz = cluster->emcpc3dphi();

    //can't find nearest pc3 hit                                     
    if(pc3_dphi > 9990) return -9999.;

    double pc3_phi = emc_phi + pc3_dphi;
    double pc3_z = emc_z + pc3_dz;
    double pc3_r = 491.03;

    TVector3 v1(1,1,1);
    v1.SetPerp(emc_r);
    v1.SetPhi(emc_phi);
    v1.SetZ(emc_z);

    TVector3 v2(1,1,1);
    v2.SetPerp(pc3_r);
    v2.SetPhi(pc3_phi);
    v2.SetZ(pc3_z);

    // Angle between EMC cluster and PC3 track
    double theta_cv = v1.Angle(v2);

    return theta_cv;
  }
} /* namespace anatools */

#endif /* OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSCLUSTER_H_ */
