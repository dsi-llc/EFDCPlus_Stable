! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! @details subroutine mapboundary surgically maps serial version boundary
  !! cell location to child partition
  ! @author Zander Mausolff, adapted from O'Donncha's code
  ! @date 9/4/2019
  !---------------------------------------------------------------------------!

  Subroutine Map_OpenBC_Pressure

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  ! *** Local variables
  integer :: ii, I, J, LL, M

  ! *** Set global values explicitly.  The variables without *_GL will become the values local to a process
  NPBW_GL = NPBW
  NPBE_GL = NPBE
  NPBN_GL = NPBN
  NPBS_GL = NPBS

  call WriteBreak(mpi_log_unit)
  write(mpi_log_unit,'(a)') 'NUMBER OF SURFACE ELEVATION OR PRESSURE BOUNDARY CONDITIONS CELLS ON OPEN BOUNDARIES'
  write(mpi_log_unit, '(a,I5)') 'Global NPBW = ', NPBW_GL
  write(mpi_log_unit, '(a,I5)') 'Global NPBE = ', NPBE_GL
  write(mpi_log_unit, '(a,I5)') 'Global NPBN = ', NPBN_GL
  write(mpi_log_unit, '(a,I5)') 'Global NPBS = ', NPBS_GL

  ! *** These are going to be redefined in the next several hundred lines of code....
  NPBW = 0
  NPBE = 0
  NPBN = 0
  NPBS = 0

  ii = 0
  do LL = 1,NPBW_GL
    I = IG2IL(IPBW_GL(LL))
    J = JG2JL(JPBW_GL(LL))

    if( I > 0 .and. I <= IC )then
      if( J > 0 .and. J <= JC )then
        if( IJCT(I,J)  >  0  .and. IJCT(I,J).LT.9 )then
          II = II + 1
          IPBW(II) = I
          JPBW(II) = J
          do M = 1,MTIDE
            PCBW(II,M) = PCBW_GL(LL,M)
            PSBW(II,M) = PSBW_GL(LL,M)
          enddo
          ISPRW(II)   = ISPRW_GL(LL)
          NPSERW(II)  = NPSERW_GL(LL)
          ISPBW(II)   = ISPBW_GL(LL)
          if( NPFORT >= 1 ) NPSERW1(II) = NPSERW1_GL(LL)
          TPCOORDW(II)   = TPCOORDW_GL(LL)
          
          NPBW = NPBW + 1
        endif
      endif
    endif
  enddo

  II = 0
  do LL = 1,NPBS_GL
    I = IG2IL(IPBS_GL(LL))
    J = JG2JL(JPBS_GL(LL))

    if( I > 0 .and. I <= IC )then
      if( J > 0 .and. J <= JC )then
        if( IJCT(I,J) >  0 .and. IJCT(I,J) .LT. 9 )then
          II = II + 1
          IPBS(II) = I
          JPBS(II) = J
          do M  = 1,MTIDE
            PCBS(II,M) = PCBS_GL(LL,M)
            PSBS(II,M) = PSBS_GL(LL,M)
          enddo
          ISPRS(II)   = ISPRS_GL(LL)
          NPSERS(II)  = NPSERS_GL(LL)
          ISPBS(II)   = ISPBS_GL(LL)
          if( NPFORT >= 1 ) NPSERS1(II) = NPSERS1_GL(LL)
          TPCOORDS(II)   = TPCOORDS_GL(LL)

          NPBS = NPBS  + 1
        endif
      endif
    endif
  enddo


  II = 0
  do LL = 1,NPBE_GL
    I = IG2IL(IPBE_GL(LL))
    J = JG2JL(JPBE_GL(LL))

    if( I > 0 .and. I <= IC )then
      if( J > 0 .and. J <= JC )then
        if( IJCT(I,J)  >  0 .and. IJCT(I,J).LT.9 )then
          II = II + 1
          IPBE(II) = I
          JPBE(II) = J
          do M = 1, MTIDE
            PCBE(II,M) = PCBE_GL(LL,M)
            PSBE(II,M) = PSBE_GL(LL,M)
          enddo
          ISPRE(II)  = ISPRE_GL(LL)
          NPSERE(II) = NPSERE_GL(LL)
          ISPBE(II)  = ISPBE_GL(LL)
          if( NPFORT >= 1 ) NPSERE1(II) = NPSERE1_GL(LL)
          TPCOORDE(II)   = TPCOORDE_GL(LL)
        
          NPBE = NPBE  + 1
        endif
      endif
    endif
  enddo

  II = 0
  do LL = 1,NPBN_GL
    I = IG2IL(IPBN_GL(LL))
    J = JG2JL(JPBN_GL(LL))
    if( I > 0 .and. I <= IC )then
      if( J > 0 .and. J <= JC )then
        if( IJCT(I,J)  >  0  .and. IJCT(I,J).LT.9 )then
          II = II + 1
          IPBN(II) = I
          JPBN(II) = J
          do M = 1,MTIDE
            PCBN(II,M) = PCBN_GL(LL,M)
            PSBN(II,M) = PSBN_GL(LL,M)
          enddo
          ISPRN(II)  = ISPRN_GL(LL)
          NPSERN(II) = NPSERN_GL(LL)
          ISPBN(II)  = ISPBN_GL(LL)
          if( NPFORT >= 1 ) NPSERN1(II) = NPSERN1_GL(LL)
          TPCOORDN(II)   = TPCOORDN_GL(LL)

          NPBN = NPBN  + 1
        endif
      endif
    endif
  enddo

  write(mpi_log_unit, '(a,I5)') ' '
  write(mpi_log_unit, '(a,I5)') 'Local NPBW  = ', NPBW
  write(mpi_log_unit, '(a,I5)') 'Local NPBE  = ', NPBE
  write(mpi_log_unit, '(a,I5)') 'Local NPBN  = ', NPBN
  write(mpi_log_unit, '(a,I5)') 'Local NPBS  = ', NPBS
  call WriteBreak(mpi_log_unit)

  return

  End Subroutine Map_OpenBC_Pressure
