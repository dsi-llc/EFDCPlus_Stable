! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Maps flow type boundary condition values to local ones
! @author Zander Mausolff, adapted from O'Donncha's code
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_River

  USE GLOBAL
  USE HYDSTRUCMOD

  Use MPI
  Use Variables_MPI
  Use Variables_MPI_Mapping
  USE MPI_All_Reduce

  Implicit None

  ! *** Local
  Real(8) :: TTDS, TWAIT
  Integer :: II, LL, MS, K, III, IERR, JJJ, L, M
  Integer :: MMAX, MMIN
  Integer :: NQSIJ_GL
  Integer :: LLSave(NQSIJM)

  NQSIJ_GL = NQSIJ
  LLSave = 0

  ! *** Determine the local NQSIJ value
  NGRPID = 0
  NQSIJ = 0
  II = 0
  DO LL = 1, NQSIJ_GL
    LQS_GL(LL) = LIJ_Global(IQS_GL(LL),JQS_GL(LL))
    IF( LQS_GL(LL) < 2 )THEN
      WRITE(6,*) 'ERROR! FLOW CELL IS NOT IN GLOBAL DOMAIN FOR FLOW BOUNDARY: ',LL
      CALL STOPP('.')
    ENDIF

    III = IG2IL(IQS_GL(LL))
    JJJ = JG2JL(JQS_GL(LL))
    IF( III > 0 .AND. III <= IC )THEN     ! *** Allow ghost cells containing inflow
      IF( JJJ > 0 .AND. JJJ <= JC )THEN   ! *** Allow ghost cells containing inflow
        NQSIJ   = NQSIJ + 1

        II = II + 1
        LLSave(II) = LL
        IQS(II)    = III
        JQS(II)    = JJJ

        NQSMUL(II)   = NQSMUL_GL(LL)
        NQSMF(II)    = NQSMF_GL(LL)
        NQSERQ(II)   = NQSERQ_GL(LL)
        NCSERQ(II,:) = NCSERQ_GL(LL,:)
        QFACTOR(II)  = QFACTOR_GL(LL)
        GRPID(II)    = GRPID_GL(LL)
        
        QSSE(II)    = QSSE_GL(LL)
        DO K=1,KC
          QSS(K,II) = QSSE(II)*DZCK(K)
        ENDDO

        ! *** Constituent Constants
        MMAX = 3 + NDYM + NTOX
        DO MS = 1,MMAX
          DO K =1,KC
            CQS(K,II,MS) = CQSE_GL(LL,MS)
          END DO
        END DO

        MMIN = MMAX + 1
        MMAX = MMAX+NSED+NSND
        DO MS = MMIN,MMAX
          DO K =1,KC
            CQS(K,II,MS) = CQSE_GL(LL,MS)
          END DO
        END DO
        
        ! *** GET NUMBER OF GROUPS.  GROUP NUMBERS ARE NOT REQUIRED TO BE CONTINGUOUS.
        IF( GRPID(II) > NGRPID ) NGRPID = GRPID(II)
      END IF
    END IF
  END DO
  Call DSI_All_Reduce(NGRPID, M, MPI_Max, TTDS, 1, TWAIT)
  NGRPID = M
  
  !! *** Constituent Series
  !DO L = 1,NQSIJ
  !  DO N=1,NTOX
  !    M=MSVTOX(N)
  !    NCSERQ(L,M)=NTOXSRQ(L)
  !  ENDDO
  !  DO N=1,NSED
  !    M=MSVSED(N)
  !    NCSERQ(L,M)=NSEDSRQ(L)
  !  ENDDO
  !  DO N=1,NSND
  !    M=MSVSND(N)
  !    NCSERQ(L,M)=NSNDSRQ(L)
  !  ENDDO
  !END DO

  Call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' FLOW BOUNDARIES '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3','N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3'
  DO II = 1,NQSIJ
    LL = LLSave(II)
    write(mpi_mapping_unit,'(2(I5,2I8,F8.1,F8.4,F8.1,I8,1X))') LL,IQS_GL(LL),JQS_GL(LL),QSSE_GL(LL),QFACTOR_GL(LL),CQSE_GL(LL,3),NCSERQ_GL(LL,3),  &
                                                                II,IQS(II),   JQS(II),   QSSE(II),   QFACTOR(II),   CQS(KC,II,3), NCSERQ(II,3)
  ENDDO

  Call WriteBreak(mpi_mapping_unit)

  RETURN

End Subroutine Map_River

!---------------------------------------------------------------------------!
! @details Maps hydraulic structure boundary condition values to local ones
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_Hydraulic_Structures

  USE GLOBAL
  USE HYDSTRUCMOD
  Use Variables_MPI
  Use Variables_MPI_Mapping

  Implicit None

  ! *** Local
  Integer :: NC, LL, MS, K, III, JJJ, L, M
  Integer :: MMAX, MMIN
  Integer :: NQCTL_GL
  Integer :: LLSave(NQCTL)

  NQCTL_GL = NQCTL
  
  ALLOCATE(HYD_STR(NQCTL))
  ALLOCATE(HSCTL(NQCTL))
  DO NC = 1,NQCTL
    HYD_STR(NC).IQCTLU = 0
    HYD_STR(NC).JQCTLU = 0
    HYD_STR(NC).IQCTLD = 0
    HYD_STR(NC).JQCTLD = 0
    HYD_STR(NC).NQCTYP = 0
    
    HSCTL(NC).IREFUP = 0
    HSCTL(NC).JREFUP = 0
    HSCTL(NC).IREFDN = 0
    HSCTL(NC).JREFDN = 0
  ENDDO
  LLSave = 0

  ! *** Determine the local NQCTL value
  NQCTL = 0
  NC = 0
  DO LL = 1,NQCTL_GL
    III = HYD_STR_GL(LL).IQCTLU
    JJJ = HYD_STR_GL(LL).JQCTLU
    IF( LIJ_Global(III,JJJ) < 2 )THEN
      WRITE(6,*) 'ERROR! UPSTREAM CELL IS NOT VALID FOR HYDRAULIC STRUCTURE ',LL
      CALL STOPP('.')
    ENDIF

    III = IG2IL(HYD_STR_GL(LL).IQCTLU)
    JJJ = JG2JL(HYD_STR_GL(LL).JQCTLU)

    IF( III > 0 .AND. III <= IC )THEN     ! *** Allow ghost cells containing inflow
      IF( JJJ > 0 .AND. JJJ <= JC )THEN   ! *** Allow ghost cells containing inflow
        NQCTL = NQCTL + 1
        NC = NC + 1
        
        LLSave(NC) = LL
        HYD_STR(NC) = HYD_STR_GL(LL)
        HSCTL(NC)   = HSCTL_GL(LL)

        HYD_STR(NC).IQCTLU = III
        HYD_STR(NC).JQCTLU = JJJ
        HSCTL(NC).IREFUP   = IG2IL(HSCTL_GL(LL).IREFUP)   ! *** Assume reference cell is in the save subdomain
        HSCTL(NC).JREFUP   = JG2JL(HSCTL_GL(LL).JREFUP)
        
        ! *** CONVERT DOWNSTREAM
        IF( HYD_STR_GL(LL).IQCTLD > 0 .AND. HYD_STR_GL(LL).IQCTLD > 0 )THEN
          III = IG2IL(HYD_STR_GL(LL).IQCTLD)
          JJJ = JG2JL(HYD_STR_GL(LL).JQCTLD)
          IF( III == 0 .OR. JJJ == 0 )THEN
            WRITE(6,*) 'ERROR! DOWNSTREAM CELL IS NOT IN DOMAIN FOR HYDRAULIC STRUCTURE ',LL
            CALL STOPP('.')
          ENDIF
          IF( LIJ(III,JJJ) < 2 )THEN
            WRITE(6,*) 'ERROR! DOWNSTREAM CELL IS NOT VALID FOR HYDRAULIC STRUCTURE ',LL
            CALL STOPP('.')
          ENDIF
          
          HYD_STR(NC).IQCTLD = III
          HYD_STR(NC).JQCTLD = JJJ
          HSCTL(NC).IREFDN   = IG2IL(HSCTL_GL(LL).IREFDN)   ! *** Assume reference cell is in the save subdomain
          HSCTL(NC).JREFDN   = JG2JL(HSCTL_GL(LL).JREFDN)
        ENDIF
      END IF
    END IF
  END DO

  Call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' HYDR STRUCTURES '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IU','JU','ID','JD','NQCTYP','TABLEID','N','IU','JU','ID','JD','NQCTYP','TABLEID'
  DO NC = 1,NQCTL
    LL = LLSave(NC)
    write(mpi_mapping_unit,'(2(I5,6I8,1X))')                                              &
          LL, HYD_STR_GL(LL).IQCTLU, HYD_STR_GL(LL).JQCTLU, HYD_STR_GL(LL).IQCTLD, HYD_STR_GL(LL).JQCTLD, HYD_STR_GL(LL).NQCTYP, HYD_STR_GL(LL).NQCTLQ,  &
          NC, HYD_STR(NC).IQCTLU,    HYD_STR(NC).JQCTLU,    HYD_STR(NC).IQCTLD,    HYD_STR(NC).JQCTLD,    HYD_STR(NC).NQCTYP,    HYD_STR(NC).NQCTLQ
  ENDDO

  Call WriteBreak(mpi_mapping_unit)

  RETURN

End Subroutine Map_Hydraulic_Structures

!---------------------------------------------------------------------------!
! @details Maps withdrawal-return boundary condition values to local ones
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_Withdrawal_Return

  USE GLOBAL
  USE HYDSTRUCMOD

  Use Variables_MPI
  Use Variables_MPI_Mapping

  Implicit None

  ! *** Local
  Integer :: NC, LL, MS, K, III, JJJ, L, M
  Integer :: MMAX, MMIN
  Integer :: NQWR_GL
  Integer :: LLSave(NQWR)
  Real, Allocatable :: CQWR_GL(:,:)
  
  ALLOCATE(CQWR_GL(NQWRM,NSTVM))
  CQWR_GL = 0.0
  
  NQWR_GL = NQWR
  
  ALLOCATE(WITH_RET(NQWR))
  ALLOCATE(WRCTL(NQWR))
  DO NC = 1,NQWR
    WITH_RET(NC).IQWRU = 0
    WITH_RET(NC).JQWRU = 0
    WITH_RET(NC).KQWRU = 0
    WITH_RET(NC).IQWRD = 0
    WITH_RET(NC).JQWRD = 0
    WITH_RET(NC).KQWRD = 0
    
    WRCTL(NC).IREFUP = 0
    WRCTL(NC).JREFUP = 0
    WRCTL(NC).IREFDN = 0
    WRCTL(NC).JREFDN = 0
    
    DO MS = 1,NSTVM
      CQWR_GL(NC,MS) = CQWR(NC,MS)
      CQWR(NC,MS) = 0.
    ENDDO
  ENDDO

  LLSave = 0

  ! *** Determine the local NQWR value
  NQWR = 0
  NC = 0
  DO LL = 1,NQWR_GL
    III = WITH_RET_GL(LL).IQWRU
    JJJ = WITH_RET_GL(LL).JQWRU
    IF( LIJ_Global(III,JJJ) < 2 )THEN
      WRITE(6,*) 'ERROR! UPSTREAM CELL IS NOT VALID FOR WITHDRAWAL-RETURN BOUNDARY ',LL
      CALL STOPP('.')
    ENDIF
    III = IG2IL(WITH_RET_GL(LL).IQWRU)
    JJJ = JG2JL(WITH_RET_GL(LL).JQWRU)
    
    IF( III > 0 .AND. III <= IC )THEN     ! *** Allow ghost cells containing inflow
      IF( JJJ > 0 .AND. JJJ <= JC )THEN   ! *** Allow ghost cells containing inflow
        NQWR = NQWR + 1
        NC = NC + 1
        
        LLSave(NC) = LL
        WITH_RET(NC) = WITH_RET_GL(NC)
        WRCTL(NC)    = WRCTL_GL(NC)

        WITH_RET(NC).IQWRU = III
        WITH_RET(NC).JQWRU = JJJ
        WRCTL(NC).IREFUP   = IG2IL(WRCTL_GL(LL).IREFUP)   ! *** Assume reference cell is in the save subdomain
        WRCTL(NC).JREFUP   = JG2JL(WRCTL_GL(LL).JREFUP)

        ! *** CONVERT DOWNSTREAM
        IF( WITH_RET_GL(LL).IQWRD > 0 .AND. WITH_RET_GL(LL).IQWRD > 0 )THEN
          III = IG2IL(WITH_RET_GL(LL).IQWRD)
          JJJ = JG2JL(WITH_RET_GL(LL).JQWRD)
          IF( III == 0 .OR. JJJ == 0 )THEN
            WRITE(6,*) 'ERROR!  FLOW RETURN CELLS NOT IN SAME DOMAIN, CHECK CARD 33 and DECOMP.INP ',LL
            CALL STOPP('.')
          ENDIF
          WITH_RET(NC).IQWRD = III
          WITH_RET(NC).JQWRD = JJJ
          WRCTL(NC).IREFDN   = IG2IL(WRCTL_GL(LL).IREFDN)   ! *** Assume reference cell is in the save subdomain
          WRCTL(NC).JREFDN   = JG2JL(WRCTL_GL(LL).JREFDN)
        ENDIF
        
        ! *** CONVERT CONSTANT RISE/FALL CONCENTRATIONS
        DO MS = 1,NSTVM
          CQWR(NC,MS) = CQWR_GL(LL,MS)
        ENDDO
        
      END IF
    END IF
  END DO

  Call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' WR/RET BOUNDRY '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IU','JU','KU','ID','JD','TABLEID','N','IU','JU','KU','ID','JD','TABLEID'
  DO NC = 1,NQWR
    LL = LLSave(NC)
    write(mpi_mapping_unit,'(2(I5,6I8,1X))')                                              &
          LL, WITH_RET_GL(LL).IQWRU, WITH_RET_GL(LL).JQWRU, WITH_RET_GL(LL).KQWRU, WITH_RET_GL(LL).IQWRD, WITH_RET_GL(LL).JQWRD, WITH_RET_GL(LL).NQWRSERQ,  &
          NC, WITH_RET(NC).IQWRU,    WITH_RET(NC).JQWRU,    WITH_RET(NC).KQWRU,    WITH_RET(NC).IQWRD,    WITH_RET(NC).JQWRD,    WITH_RET(NC).NQWRSERQ
  ENDDO

  Call WriteBreak(mpi_mapping_unit)

  RETURN

End Subroutine Map_Withdrawal_Return

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Maps water quality inflow boundary condition values to local ones
! @date 2020-03-13
!---------------------------------------------------------------------------!

Subroutine Map_WQ_PointSource

  USE GLOBAL
  Use Variables_WQ   !, Only:IWQPSC, IWQPSV
  Use Variables_MPI
  Use Variables_MPI_Mapping

  Implicit None

  ! *** Local
  Integer :: II, LL, MS, K, III, JJJ, L, M, NT, NW
  Integer :: MMAX, MMIN
  Integer :: IWQPS_GL
  Integer :: KCPSL_GL(NWQPSM)
  Integer :: MVPSL_GL(NWQPSM)
  Real    :: XPSL_GL(0:NWQPSM,NWQVM)
  Integer :: LLSave(NWQPSM)

  IWQPS_GL = NWQPS
  KCPSL_GL = KCPSL
  MVPSL_GL = MVPSL
  XPSL_GL  = WQWPSLC
  LLSave = 0
  
  ! *** Determine the local NWQPS value
  NWQPS = 0
  KCPSL = 0
  MVPSL = 0
  WQWPSLC = 0.
  II = 0
  DO LL = 1, IWQPS_GL
  
    III = Map2Local(LQS_GL(LL)).IL           ! *** Assumes WQ point source cells use the same order and number
    JJJ = Map2Local(LQS_GL(LL)).JL           ! *** Assumes WQ point source cells use the same order and number
    IF( III > 0 .AND. III <= IC )THEN        ! *** Allow ghost cells containing inflow
      IF( JJJ > 0 .AND. JJJ <= JC )THEN      ! *** Allow ghost cells containing inflow
        NWQPS   = NWQPS + 1

        II = II + 1
        LLSave(II) = LL

        ICPSL(II) = III
        JCPSL(II) = JJJ
        KCPSL(II) = KCPSL_GL(LL)
        MVPSL(II) = MVPSL_GL(LL)
      
        ! *** ASSIGN GLOBAL CONCENTRATION TIME SERIES INDEX
        IF( IWQPSL == 2 )THEN
          ! *** CONSTAND CONCENTRATIONS
          NCSERQ(II,8) = NSERWQ(LL)          ! *** ALL WQ VARIABLES USE SAME TIME SERIES
          DO NW=1,NWQV
            IF( ISTRWQ(NW) > 0 )THEN
              NT = MSVWQV(NW)
              DO K =1,KC
                CQS(K,II,NT) = XPSL_GL(LL,NW)
              ENDDO
            ENDIF
          ENDDO
        ELSE
          ! *** CONSTANT MASS FLUXES
          ! *** NW 1 TO 19 ARE ALREADY IN MASS (G/DAY) FROM WQ3DCONTROL [CONC (mg/l) AND Q (m3/s)]
          ! *** NW = 20 ALREADY IN moles
          WQWPSLC(II,1:NWQV) = XPSL_GL(LL,1:NWQV)
          
          ! *** CONVERT FROM MPN/L TO MPN/DAY
          WQWPSLC(II,IFCB) = XPSL_GL(LL,IFCB) * 1000.
      
          ! *** Assign loading by layer
          L = LIJ(III,JJJ)
          IF( IWQPSL == 0 ) FORALL(K=1:KC) WQWPSL(L,K,1:NWQV) = WQWPSL(L,K,1:NWQV) + WQWPSLC(LL,1:NWQV)/DZI

          ! *** Uniform load in horizontal cell stack over all layers
          DO K=1,KC
            IWQPSC(L,K) = II
            IWQPSV(L,K) = MVPSL(II)
          ENDDO
        ENDIF

      END IF
    END IF
  END DO

  Call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' WQ POINT SRCS  '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3','N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3'
  DO II = 1,NWQPS
    LL = LLSave(II)
    write(mpi_mapping_unit,'(2(I5,2I8,F8.1,F8.4,F8.1,I8,1X))') LL,IQS_GL(LL),JQS_GL(LL),QSSE_GL(LL),QFACTOR_GL(LL),CQSE_GL(LL,3),NCSERQ_GL(LL,3),  &
                                                                II,IQS(II),   JQS(II),   QSSE(II),   QFACTOR(II),   CQS(KC,II,3), NCSERQ(II,3)
  ENDDO

  Call WriteBreak(mpi_mapping_unit)

  RETURN

End Subroutine Map_WQ_PointSource

  
