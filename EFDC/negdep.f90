! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE NEGDEP(NOPTIMAL1, LDMOPT1, QCHANUT, QCHANVT, SUB1, SVB1)

  ! *** SUBROUTINE NEGDEP CHECK EXTERNAL SOLUTION FOR NEGATIVE DEPTHS

  ! CHANGE RECORD

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !-------------------------------------------------------------------------------------------!
  !    2014-12       Paul M. Craig     Rewrote the NEGDEP function to be always called with
  !                                    a lot more information dumped to the EFDCLOG.OUT file,
  !                                    including the new ice sub-model information

  use GLOBAL
  use RESTART_MODULE
  use EFDCOUT
  use Mod_Map_Write_EE_Binary

  implicit none

  integer, intent(IN)    :: NOPTIMAL1,LDMOPT1

  real,intent(IN),dimension(:)    :: SUB1(LCM)
  real,intent(IN),dimension(:)    :: SVB1(LCM)

  integer :: INEGFLG, L, NMD, LHOST, IHOST, JHOST, LCHNU, LCHNV, ICHNU, JCHNU, ICHNV, JCHNV, K, ND, NNEG, LF, LL, LS, LN, LW, LE
  real    :: QCHANUT(NCHANM), QCHANVT(NCHANM), SRFCHAN, SRFHOST, SRFCHAN1, SRFHOST1, SURFTMP, DELTD2

  INEGFLG = 0
  
  ! ***  CHECK FOR NEGATIVE DEPTHS
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L)
  do ND = 1,NOPTIMAL1
    LF = 2+(ND-1)*LDMOPT1
    LL = MIN(LF+LDMOPT1-1,LA)
    do L = LF,LL
      if( HP(L) <=  0. )then
        INEGFLG = 1
        EXIT
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO
  
  if( INEGFLG == 0 )then
    return
  endif

  ! *** Count the negative depths
  NNEG = 0
  do L = 2,LA
    if( HP(L) <=  0. )then
      NNEG = NNEG + 1
    endif
  enddo

  PRINT '(A,I5,A,F15.5)', 'Found ', NNEG, ' CELLS WITH NEGATIVE DEPTHS @ ', TIMEDAY
  
  open(mpi_error_unit,FILE = OUTDIR//mpi_error_file,POSITION = 'APPEND')

  do L = 2,LA
    if( HP(L) <= 0.0 )then
      LN = LNC(L)
      LS = LSC(L)
      LW = LWC(L)
      LE = LEC(L)
      write(6,1111) TIMEDAY, NITER, ISTL, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, process_id
      WRITE (6,6060) Map2Global(L).IG, Map2Global(L).JG, HP(L), H1P(L), H2P(L)
      WRITE (6,6061) Map2Global(L).IG, Map2Global(L).JG, HU(L), H1U(L)
      WRITE (6,6062) Map2Global(L).IG, Map2Global(L).JG, HU(LE), H1U(LE)
      WRITE (6,6063) Map2Global(L).IG, Map2Global(L).JG, HV(L), H1V(L)
      WRITE (6,6064) Map2Global(L).IG, Map2Global(L).JG, HV(LN), H1V(LN)
      WRITE (6,6065) Map2Global(L).IG, Map2Global(L).JG, QSUME(L), QSUM1E(L)
      WRITE (6,6066) Map2Global(L).IG, Map2Global(L).JG, SUB(L), SUB(LE), SVB(L), SVB(LN)
      if( IEVAP > 0 )then
        WRITE (6,6067) Map2Global(L).IG, Map2Global(L).JG,RAINT(L),EVAPT(L)
      endif
      if( ISICE > 0 )then
        WRITE (6,6068) Map2Global(L).IG, Map2Global(L).JG, ICETHICK(L), ICETHICK1(L), TEM(L,KC), TATMT(L), SOLSWRT(L)
      endif

      if( ISDYNSTP == 0 )then
        DELT = DT
        DELTD2 = 0.5*DT
      else
        DELT = DTDYN
        DELTD2 = 0.5*DTDYN
      endif

      write(mpi_error_unit,1111) TIMEDAY, NITER, ISTL, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, process_id

      ! *** EE7.2 DIAGNOSTICS
      WRITE (mpi_error_unit,'(2X,A14,5I14)')'L  CWESN', Map2Global(L).LG, Map2Global(LW).LG, Map2Global(LE).LG, Map2Global(LS).LG, Map2Global(LN).LG
      WRITE (mpi_error_unit,'(A)')'DEPTHS'
      WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'HP CWESN',  HP(L), HP(LW), HP(LE), HP(LS), HP(LN)
      WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'H1P CWESN', H1P(L), H1P(LW), H1P(LE), H1P(LS), H1P(LN)
      if( IS2TIM == 0 ) WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'H2P CWESN', H2P(L), H2P(LW), H2P(LE), H2P(LS), H2P(LN)

      WRITE (mpi_error_unit,'(A)')'WATER SURFACE ELEVATIONS'
      WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'WS CWESN',  BELV(L)+HP(L), BELV(LW)+HP(LW), BELV(LE)+HP(LE), BELV(LS)+HP(LS), BELV(LN)+HP(LN)
      WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'WS1 CWESN', BELV(L)+H1P(L),BELV(LW)+H1P(LW),BELV(LE)+H1P(LE),BELV(LS)+H1P(LS),BELV(LN)+H1P(LN)
      if( IS2TIM == 0 ) WRITE (mpi_error_unit,'(2X,A14,5E14.6)')'WS2 CWESN', BELV(L)+H2P(L),BELV(LW)+H2P(LW),BELV(LE)+H2P(LE),BELV(LS)+H2P(LS),BELV(LN)+H2P(LN)

      WRITE (mpi_error_unit,'(A)')'FACE DEPTHS'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'HU/HV WESN',   HU(L),  HU(LE),  HV(L),  HV(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'H1U/H1V WESN', H1U(L), H1U(LE), H1V(L), H1V(LN)

      WRITE (mpi_error_unit,'(A)')'MASKS'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'SUB WESN',  SUB(L), SUB(LE), SVB(L), SVB(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'SUB1 WESN', SUB1(L),SUB1(LE),SVB1(L),SVB1(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'SUBO WESN', SUBO(L),SUBO(LE),SVBO(L),SVBO(LN)

      WRITE (mpi_error_unit,'(A)')'DEPTH AVERAGE VELOCITIES'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'U/V WESN',     UHDYE(L)/HU(L)/DYU(L),  UHDYE(LE)/HU(LE)/DYU(LE),  VHDXE(L)/HV(L)/DXV(L),  VHDXE(LN)/HV(LN)/DXV(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'U1/V1 WESN',   UHDY1E(L)/H1U(L)/DYU(L),UHDY1E(LE)/H1U(LE)/DYU(LE),VHDX1E(L)/H1V(L)/DXV(L),VHDX1E(LN)/H1V(LN)/DXV(LN)

      WRITE (mpi_error_unit,'(A)')'BOTTOM SHEAR VELOCITIES'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'UV/VU WESN',   UV(L), UV(LE), VU(L), VU(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'U1V/V1U WESN', U1V(L),U1V(LE),V1U(L),V1U(LN)

      WRITE (mpi_error_unit,'(A)')'FLUX TERMS (W/E AND S/N)'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'UHDYE/VHDXE',  UHDYE(L), UHDYE(LE), VHDXE(L), VHDXE(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'UHDY1E/VHDX1E',UHDY1E(L),UHDY1E(LE),VHDX1E(L),VHDX1E(LN)
      if( IS2TIM == 0 ) WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'UHDY2E/VHDX2E',UHDY2E(L),UHDY2E(LE),VHDX2E(L),VHDX2E(LN)
      WRITE (mpi_error_unit,'(A)')'MOMENTUM TERMS (W/E AND S/N)'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'FUHDYE/FVHDXE',FUHDYE(L),FUHDYE(LE),FVHDXE(L),FVHDXE(LN)
      WRITE (mpi_error_unit,'(A)')'SOURCE/SINK TERMS'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'QSUME/QSUM1E',QSUME(L),QSUM1E(L)

      ! ***                       FUHDYE(L) = UHDYE(L)-DELTD2*SUB(L)*HRUO(L)*HU(L)*(P(L)-P(LW))+SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L))+FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      ! ***                                       ***               WEST                            ******                 EAST
      WRITE (mpi_error_unit,'(16X,A)')'     WEST          EAST'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X P GRADIENT',  SUBO(L)*DELTD2*HRUO(L)*H1U(L)*(P1(LW)-P1(L)),     SUBO(LE)*DELTD2*HRUO(LE)*H1U(LE)*(P1(L)-P1(LE))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X SHEARS TOT',  SUBO(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L))+FCAXE(L)+FPGXE(L)-SNLT*FXE(L)), SUBO(LE)*DELT*DXIU(LE)*(DXYU(LE)*(TSX(LE)-RITB1*TBX(LE))+FCAXE(LE)+FPGXE(LE)-SNLT*FXE(LE))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X SHEARS T/B',  SUBO(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L))), SUBO(LE)*DELT*DXIU(LE)*(DXYU(LE)*(TSX(LE)-RITB1*TBX(LE)))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X SHEARS FCAXE',SUBO(L)*DELT*DXIU(L)*FCAXE(L),                        SUBO(LE)*DELT*DXIU(L)*FCAXE(LE)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X SHEARS FPGXE',SUBO(L)*DELT*DXIU(L)*FPGXE(L),                        SUBO(LE)*DELT*DXIU(L)*FPGXE(LE)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'X SHEARS FXE',  SUBO(L)*DELT*DXIU(L)*SNLT*FXE(L),                     SUBO(LE)*DELT*DXIU(L)*SNLT*FXE(LE)

      ! ***                       FVHDXE(L) = VHDXE(L)-DELTD2*SVB(L)*HRVO(L)*HV(L)*(P(L)-P(LS ))+SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L))-FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ! ***                                       ***               SOUTH                           ******                 NORTH
      WRITE (mpi_error_unit,'(16X,A)')'     SOUTH        NORTH'
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y P GRADIENT',  SVBO(L)*DELTD2*HRVO(L)*HV(L)*(P1(L)-P1(LS )),         SVBO(LN)*DELTD2*HRVO(LN)*HV(LN)*(P1(LN)-P1(L))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y SHEARS TOT',  SVBO(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L))-FCAYE(L)+FPGYE(L)-SNLT*FYE(L)), SVBO(LN)*DELT*DYIV(LN)*(DXYV(LN)*(TSY(LN)-RITB1*TBY(LN))-FCAYE(LN)+FPGYE(LN)-SNLT*FYE(LN))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y SHEARS T/B',  SVBO(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L))), SVBO(LN)*DELT*DYIV(LN)*(DXYV(LN)*(TSY(LN)-RITB1*TBY(LN)))
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y SHEARS FCAYE',SVBO(L)*DELT*DYIV(L)*FCAYE(L),                        SVBO(LN)*DELT*DYIV(LN)*FCAYE(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y SHEARS FPGYE',SVBO(L)*DELT*DYIV(L)*FPGYE(L),                        SVBO(LN)*DELT*DYIV(LN)*FPGYE(LN)
      WRITE (mpi_error_unit,'(2X,A14,4E14.6)')'Y SHEARS FYE',  SVBO(L)*DELT*DYIV(L)*SNLT*FYE(L),                     SVBO(LN)*DELT*DYIV(LN)*SNLT*FYE(LN)

      !RCX(L) = 1./( 1.+DELT*FXVEGE(L) )
      !RCY(L) = 1./( 1.+DELT*FYVEGE(L) )
      if( ISVEG > 0 )then
        WRITE (mpi_error_unit,'(/2X,A14,4E14.6)'),'FXVEGE X/Y',FXVEGE(L),FXVEGE(L),FYVEGE(L),FYVEGE(LN)
        WRITE (mpi_error_unit,'( 2X,A14,4E14.6)'),'RCX/RCY', 1./( 1.+DELT*FXVEGE(L) ),1./( 1.+DELT*FYVEGE(L) ),1./( 1.+DELT*FXVEGE(LE) ),1./( 1.+DELT*FYVEGE(LN) )
      endif

      if( ISICE > 0 )then
        WRITE (mpi_error_unit,'(A)')'ICE CONDITIONS'
        WRITE (mpi_error_unit,'(2X,A,6E14.6)') 'TEM,TATM,SOLAR,ICETEMP,FRAZILICE,ICEVOL',TEM(L,KC),TATMT(L),SOLSWRT(L),ICETEMP(L),FRAZILICE(L,KC),ICEVOL(L)
        WRITE (mpi_error_unit,'(A)')'ICE THICKNESS'
        WRITE (mpi_error_unit,'(2X,A14,5E14.6)') 'ICE  CWESN',  ICETHICK(L), ICETHICK(LE), ICETHICK(LW), ICETHICK(LS), ICETHICK(LN)
        WRITE (mpi_error_unit,'(2X,A14,5E14.6)') 'ICE1 CWESN',  ICETHICK1(L),ICETHICK1(LE),ICETHICK1(LW),ICETHICK1(LS),ICETHICK1(LN)
      endif
    endif
  enddo

  do L = 2,LA
    if( (HU(L) < 0. .and. SUBO(L) > 0.5) .or. (HV(L) < 0. .and. SVBO(L) > 0.5 ) )then
      LN = LNC(L)
      write(6,1112)
      WRITE (6,6060) Map2Global(L).IG, Map2Global(L).JG, HP(L), H1P(L), H2P(L)
      WRITE (6,6061) Map2Global(L).IG, Map2Global(L).JG, HU(L), H1U(L)
      WRITE (6,6062) Map2Global(L).IG, Map2Global(L).JG, HU(LEC(L)), H1U(LEC(L))
      WRITE (6,6063) Map2Global(L).IG, Map2Global(L).JG, HV(L), H1V(L)
      WRITE (6,6064) Map2Global(L).IG, Map2Global(L).JG, HV(LN), H1V(LN)
      WRITE (6,6065) Map2Global(L).IG, Map2Global(L).JG, QSUME(L), QSUM1E(L)
        
      write(mpi_error_unit,1112)
      WRITE (mpi_error_unit,6060) Map2Global(L).IG, Map2Global(L).JG, HP(L), H1P(L), H2P(L)
      WRITE (mpi_error_unit,6061) Map2Global(L).IG, Map2Global(L).JG, HU(L), H1U(L)
      WRITE (mpi_error_unit,6062) Map2Global(L).IG, Map2Global(L).JG, HU(LEC(L)), H1U(LEC(L))
      WRITE (mpi_error_unit,6063) Map2Global(L).IG, Map2Global(L).JG, HV(L), H1V(L)
      WRITE (mpi_error_unit,6064) Map2Global(L).IG, Map2Global(L).JG, HV(LN), H1V(LN)
      WRITE (mpi_error_unit,6065) Map2Global(L).IG, Map2Global(L).JG, QSUME(L), QSUM1E(L)
    endif
  enddo

  if( ISPPH == 1 .and. num_Processors == 1 )then
    WRITE (6,6069) mpi_log_file

    call Map_Write_EE_Binary
    if( process_id == master_id )then
      write(mpi_log_unit,*) 'NEGDEP: WRITING EE LINKAGE',TIMEDAY
      call EE_LINKAGE(-1)  
    endif
  endif
  
  if( ISRESTO == 0 .or. ISGREGOR > 0 )then
    call Restart_Out(TIMEDAY, 1)
  endif

  if( MDCHH > 0 )then
    write(mpi_error_unit,8000)
    do NMD = 1,MDCHH
      LHOST = LMDCHH(NMD)
      IHOST = IL(LHOST)
      JHOST = JL(LHOST)
      LCHNU = LMDCHU(NMD)
      LCHNV = LMDCHV(NMD)

      ! *** X-DIRECTION CHANNEL
      if( HP(LHOST) < 0.0 .or. HP(LCHNU) < 0.0 )then
        if( MDCHTYP(NMD) == 1 )then
          ICHNU = IL(LCHNU)
          JCHNU = JL(LCHNU)
          SRFCHAN = HP(LCHNU)+BELV(LCHNU)
          SRFHOST = HP(LHOST)+BELV(LHOST)
          SRFCHAN1 = H1P(LCHNU)+BELV(LCHNU)
          SRFHOST1 = H1P(LHOST)+BELV(LHOST)
          write(mpi_error_unit,8001) N,NMD,MDCHTYP(NMD),ICHNU,JCHNU,ISCDRY(LCHNU),SRFCHAN,HP(LCHNU),P1(LCHNU), H1P(LCHNU)
          write(mpi_error_unit,8002) IHOST,JHOST,ISCDRY(LHOST),SRFHOST,HP(LHOST),P1(LHOST), H1P(LHOST)
          write(mpi_error_unit,8003) QCHANU(NMD),QCHANUT(NMD),CCCCHU(NMD),CCCCHV(NMD)
        endif
      endif

      ! *** Y-DIRECTION CHANNEL
      if( HP(LHOST) < 0.0 .or. HP(LCHNV) < 0.0 )then
        if( MDCHTYP(NMD) == 2 )then
          ICHNV = IL(LCHNV)
          JCHNV = JL(LCHNV)
          SRFCHAN = HP(LCHNV)+BELV(LCHNV)
          SRFHOST = HP(LHOST)+BELV(LHOST)
          SRFCHAN1 = H1P(LCHNV)+BELV(LCHNV)
          SRFHOST1 = H1P(LHOST)+BELV(LHOST)
          write(mpi_error_unit,8001) NITER,NMD,MDCHTYP(NMD),ICHNV,JCHNV,ISCDRY(LCHNV),SRFCHAN,HP(LCHNV),SRFCHAN1,H1P(LCHNV)
          write(mpi_error_unit,8002) IHOST,JHOST,ISCDRY(LHOST),SRFHOST,HP(LHOST),SRFHOST1,H1P(LHOST)
          write(mpi_error_unit,8003) QCHANV(NMD),QCHANVT(NMD),CCCCHU(NMD),CCCCHV(NMD)
        endif
      endif
      write(mpi_error_unit,8004)
    enddo
  endif

  close(1)
  open(1,FILE = OUTDIR//'EQTERM.OUT',STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  open(1,FILE = OUTDIR//'EQTERM.OUT',POSITION = 'APPEND',STATUS = 'UNKNOWN')
  write(1,1001)NITER,ISTL
  do L = 2,LA
    write(1,1001) Map2Global(L).IG, Map2Global(L).JG, SUB(L), SVB(L), HRUO(L), HRVO(L), HU(L), HV(L)
  enddo
  close(1)
    
  if( ISINWV == 1 )then
    open(1,FILE = OUTDIR//'CFLMAX.OUT')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'CFLMAX.OUT')
    do L = 2,LA
      write(1,1991) Map2Global(L).IG, Map2Global(L).JG,(CFLUUU(L,K),K = 1,KC)
      write(1,1992) (CFLWWW(L,K),K = 1,KC)
      write(1,1992) (CFLCAC(L,K),K = 1,KC)
      write(1,1992) (CFLVVV(L,K),K = 1,KC)
    enddo
    close(1)
  endif

  close(mpi_error_unit)
#ifdef GNU
  STOP 'ABORTING RUN DUE TO NEGATIVE DEPTHS'
#else
  ERROR STOP('ABORTING RUN DUE TO NEGATIVE DEPTHS')
#endif  

1001 FORMAT(2I5,10(1X,E12.4))
1002 FORMAT(3I4,10(1X,E9.2))
1991 FORMAT(2I5,12F8.3)
1992 FORMAT(10X,12F8.3)
1111 FORMAT(' *************************************************************************************',/, &
            ' *************************************************************************************',/, &
            ' NEG DEPTH AT CELL CENTER: Timeday = ',F14.6,'  NITER = ',I12,'  ISTL, L, I, J = '4I10,'  Rank = ',I5)
1112 FORMAT(' NEG DEPTH AT WEST FACE')
1113 FORMAT(' NEG DEPTH AT SOUTH FACE')
6060 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  HP,H1P,H2P    = ',3(2X,E12.4))
6061 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  HUW,H1UW      = ',2(2X,E12.4))
6062 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  HUE,H1UE      = ',2(2X,E12.4))
6063 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  HVS,H1VS      = ',2(2X,E12.4))
6064 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  HVN,H1VN      = ',2(2X,E12.4))
6065 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  QSUME,QSUM1E  = ',2(2X,E12.4))
6067 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  RAIN,EVAP     = ',2(2X,E12.4))
6066 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  SUB,SVB       = ',4F6.1)
6068 FORMAT('  NEG DEPTH AT I,J  = ',2I4,'  ICETH,ICETH1,TEM,TATM,SOLAR  = ',2F8.3,2F8.2,F8.1)
6069 FORMAT('  ************************************************************************************',/, &
    '  ***      MORE DETAILS CAN BE FOUND IN THE FILE:  ',A22,'          ***',/, &
    '  ***      THE LATEST MODEL RESULTS HAVE ALSO BEEN SAVED TO THE EE LINKAGE         ***',/, &
    '  ************************************************************************************')
8001 FORMAT(I7,5I5,4E13.4)
8002 FORMAT(17X,3I5,4E13.4)
8003 FORMAT(32X,4E13.4)
8000 FORMAT('    NITER    NMD  MTYP   I    J  IDRY      P           H           P1           H1')
8004 FORMAT('                                     QCHANU       QCHANUT      CCCCHU       CCCCHV ')

  END

