/*
 * AnaToolsTrigger.h
 *
 *  Created on: Feb 26, 2016
 *      Author: nfeege
 */

#ifndef __ANATOOLSTRIGGER_H__
#define __ANATOOLSTRIGGER_H__

#include "AnaToolsCluster.h"
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
  enum TriggerMode {ERT_4x4a, ERT_4x4b, ERT_4x4c, ERT_2x2, RICH, ERT_4x4or};

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
  inline int PassERT( const ErtOut* ertout, const emcClusterContent* cluster, const TriggerMode triggermode )
  {
    /*
       int ERThit_N Number of ERT hits in a event
       int ERTtrigmode[ERThit_N], Trigger mode 0:4x4a, 1:4x4b, 2:4x4c, 3:2x2, 4:RICH
       int ERTsm[ERThit_N], super module(sm) id, SM counts from left to right, from bottom to top with outside of the detector.

       0~17 for PbSc, 0~31 for PbGl, 0~31 for RICH
       */

    /* Here's where the actual function starts */
    int arm = cluster->arm();
    int sector = cluster->sector();
    int sm = GetSM(cluster);
    if(triggermode != ERT_4x4or)
    {
      int trigger = ertout->get_ERTbit(triggermode, arm, sector, sm);
      return trigger;
    }
    else
    {
      int trigger[3] = {};
      for(int i=0; i<3; i++)
        trigger[i] = ertout->get_ERTbit((TriggerMode)i, arm, sector, sm);
      return trigger[0] || trigger[1] || trigger[2];
    }
  }

} /* namespace anatools */

#endif /* __ANATOOLSTRIGGER_H__ */
