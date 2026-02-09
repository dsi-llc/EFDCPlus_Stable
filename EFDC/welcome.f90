! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

SUBROUTINE WELCOME  
    
  ! DATE MODIFIED     BY           
  ! 06/25/2006        Paul M. Craig
  !                   Updated Code to Fortran 90
  use GLOBAL,only:EFDC_VER
  
  ! *** EFDC_VER is set in the main routine, aaefdc.f90
  write(6,1) EFDC_VER      

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
      ,/,'*                INCLUDING:                                                       *'  &  
      ,/,'*                   GENERAL OCEAN TURBULENCE MODEL (GOTM)                         *'  &  
      ,/,'*                   VERTICAL LAYERING WITH SIGMA-STRETCHED OR SIGMA-ZED           *'  &  
      ,/,'*                   PROPELLER WASH, LAGRANGIAN PARTICLE TRACKING,                 *'  &  
      ,/,'*                   NEW EUTROPHICATION WITH ZOOPLANKTON AND RPEM SUB-MODELS       *'  &  
      ,/,'*                   WITH SEDZLJ, HYDROKINETIC DEVICES (SNL)                       *'  &  
      ,/,'*                   MULTI-THREADED VERSION USING OpenMP                           *'  &  
      ,/,'*                   MULTI-PROCESSOR VERSION USING OpenMPI                         *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*        EFDC+ HAS BEEN CUSTOMIZED TO WORK WITH DSI''s EFDC+ EXPLORER 12.4         *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'*                         VERSION DATE: MPI ',A10,'                            *'  &  
      ,/,'*                                                                                 *'  &  
      ,/,'***********************************************************************************') 
  return 
END  
