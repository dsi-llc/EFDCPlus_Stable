! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

SUBROUTINE WELCOME  
    
  ! DATE MODIFIED     BY           
  ! 06/25/2006        Paul M. Craig
  !                   Updated Code to Fortran 90
  USE GLOBAL,ONLY:EFDC_VER
  
  ! *** EFDC_VER is set in the main routine, aaefdc.f90
  WRITE(6,1) EFDC_VER      

1 FORMAT('***********************************************************************************'  &
      ,/,'*                                                                                 *'  &  
      ,/,'*                                                                                 *'  &
      ,/,'*              EEEEEEEEE    FFFFFFFFF    DDDDDDDD       CCCCCCCC                  *'  &
      ,/,'*             EEE          FFF          DDD     DD    CCC      CC   +             *'  &
      ,/,'*            EEE          FFF          DDD     DD    CCC            +             *'  &
      ,/,'*           EEEEEEEE     FFFFFFFF     DDD     DD    CCC         +++++++++         *'  &
      ,/,'*          EEE          FFF          DDD     DD    CCC              +             *'  &  
      ,/,'*         EEE          FFF          DDD     DD    CCC      CC       +             *'  &  
      ,/,'*        EEEEEEEEE    FFF          DDDDDDDDDD      CCCCCCCCC                      *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*                     ENVIRONMENTAL FLUID DYNAMICS CODE (PLUS)                    *'  &  
      ,/,'*                     ORIGINALLY DEVELOPED BY JOHN M. HAMRICK                     *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*        EFDC+   BY DSI, LLC,  EDMONDS WA, USA                                    *'  &  
      ,/,'*                   VERTICAL LAYERING WITH SIGMA-STRETCHED OR SIGMA-ZED           *'  &  
      ,/,'*                   MULTI-THREADED VERSION USING OpenMP                           *'  &  
      ,/,'*                   MULTI-PROCESSOR VERSION USING OpenMPI                         *'  &  
      ,/,'*                   PROPELLER WASH, LAGRANGIAN PARTICLE TRACKING,                 *'  &  
      ,/,'*                   NEW EUTROPHICATION WITH ZOOPLANKTON AND RPEM SUB-MODELS       *'  &  
      ,/,'*                   WITH SEDZLJ, HYDROKINETIC DEVICES (SNL)                       *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*        EFDC+ HAS BEEN CUSTOMIZED TO WORK WITH DSI''s EFDC_EXPLORER 11.2.0        *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*                         VERSION DATE: MPI ',A10,'                            *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'***********************************************************************************') 
  RETURN  
END  
