/*
 * AnaToolsTrigger.h
 *
 *  Created on: Feb 26, 2016
 *      Author: nfeege
 */

#ifndef OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSTRIGGER_H_
#define OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSTRIGGER_H_

#define NARMSECT 8

#include <TrigLvl1.h>
#include <TriggerHelper.h>
#include <ErtOut.h>
#include <emcClusterContent.h>

/*! Namespace with various functions for analysis.
 * This file provides functions to check trigger information.
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
namespace anatools
{
  /*!
   * Trigger mode
   *
   * (copied from /offline/AnalysisTrain/Run13_Pi0Ana_YIS)
   */
  enum TriggerMode {ERT_4x4a, ERT_4x4b, ERT_4x4c, ERT_2x2, RICH, ERT_Or, ERT_And};

  /*!
   * Trigger mask bits
   *
   * (copied from /offline/AnalysisTrain/Run13_Pi0Ana_YIS and extended)
   */
  enum TriggerMaskBits {
    Mask_ERT_4x4a=0x00000080,
    Mask_ERT_4x4b=0x00000040,
    Mask_ERT_4x4c=0x00000100,
    Mask_ERT_4x4or=0x000001C0,
    Mask_BBC_novtx=0x00000002,
    Mask_BBC_widevtx=0x00000001,
    Mask_BBC_narrowvtx=0x00000010,
    Mask_Clock=0x00000800};


  /*!
   * Check if cluster triggered ERT.
   *
   * (copied from /offline/AnalysisTrain/Run13_Pi0Ana_YIS)
   */
  inline Int_t PassERT( ErtOut* ertout, emcClusterContent* cluster, TriggerMode triggermode )
  {
    /*
      int ERThit_N Number of ERT hits in a event
      int ERTtrigmode[ERThit_N], Trigger mode 0:4x4a, 1:4x4b, 2:4x4c, 3:2x2, 4:RICH
      int ERTsm[ERThit_N], super module(sm) id, SM counts from left to right, from bottom to top with outside of the detector.

      0~17 for PbSc 0~31 for PbGl, 0~31 for RICH
    */

    /* Tis does not belong here... */
    Int_t nsm[NARMSECT];
    for(Int_t i=0;i<NARMSECT;++i)
      {
        if(i==4||i==5) nsm[i]=4;
        else nsm[i]=3;
      }

    /* Here's where the actual function starts */
    Int_t arm = cluster->arm();
    Int_t sector = cluster->sector();
    Int_t armsect = 4*arm + sector;
    Int_t y = cluster->iypos();
    Int_t z = cluster->izpos();
    Int_t sm = (Int_t)((y/12)*2*nsm[armsect]+z/12);
    Int_t trigger = ertout->get_ERTbit(triggermode, arm, sector, sm);

    return trigger;
  }


} /* namespace anatools */

#endif /* OFFLINE_ANALYSISTRAIN_DIRECTPHOTONPP_ANATOOLSTRIGGER_H_ */
