! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE XYIJCONV
!Author: Dang Huu Chung

use GLOBAL  
USE INFOMOD
Use Variables_MPI
Use Variables_MPI_Mapping
    use Broadcast_Routines

implicit none

contains

SUBROUTINE XY2IJ(CEL, VALID) 
  type(CELL),intent(INOUT) :: CEL
  integer, intent(IN), OPTIONAL :: VALID
  integer :: N, NPMAX, ISVALID
  
  if( PRESENT(VALID) )then
    ISVALID = 1
  else
    ISVALID = 0
  endif
  
  NPMAX = SIZE(CEL.XCEL)
  
  do N = 1,NPMAX 
    call CONTAINERIJ_GL(N, CEL.XCEL(N), CEL.YCEL(N), CEL.ICEL(N), CEL.JCEL(N), ISVALID)
  enddo
  
END SUBROUTINE

SUBROUTINE CONTAINERIJ(NCEL, XCLL, YCLL, ICLL, JCLL, ISVALID)   
  integer,intent(IN)  :: NCEL
  real(8), intent(IN) :: XCLL, YCLL
  integer, intent(IN) :: ISVALID
  integer, intent(OUT) :: ICLL, JCLL
  integer :: LMILOC(1), L, I, J, ILN, JLN
  integer :: I1, I2, J1, J2
  real(8) :: RADLA(LA)
  
  ! *** FOR THE FIRST CALL                     
  RADLA(2:LA) = SQRT((XCLL-XCOR(2:LA,5))**2+(YCLL-YCOR(2:LA,5))**2)     ! *** Delme - This is inefficient.  Todo - Make more efficient.  use LEC, LWC, LSC, LNC
  LMILOC = MINLOC(RADLA(2:LA))
  ILN = IL(LMILOC(1)+1)    !I OF THE NEAREST CELL FOR DRIFTER
  JLN = JL(LMILOC(1)+1)    !J OF THE NEAREST CELL FOR DRIFTER     

  ! *** DETERMINE THE CELL CONTAINING THE DRIFTER WITHIN 9 CELLS: LLA(NCEL)
  I1 = max(1,ILN-1)
  I2 = min(ILN+1,ICM)
  J1 = max(1,JLN-1)
  J2 = min(JLN+1,JCM)
  do J = J1,J2
    do I = I1,I2
      L = LIJ(I,J)
      if( L < 2 ) CYCLE
      if( INSIDECELL(L, XCLL, YCLL) )then
        ICLL = I
        JCLL = J
        return
      endif
    enddo
  enddo

  ! *** Cooridnates not found.  Determin the action
  if( ISVALID > 0 )then
    ! *** Flag as invalid
    ICLL = -999
    JCLL = -999
  else
    PRINT *,'THIS CELL IS OUTSIDE THE DOMAIN:',NCEL
    call STOPP('Invalid XY in SUBSET.INP')
  endif
  
END SUBROUTINE

SUBROUTINE CONTAINERIJ_GL(NCEL, XCLL, YCLL, ICLL, JCLL, ISVALID)   
  integer,intent(IN)  :: NCEL
  real(8), intent(IN) :: XCLL, YCLL
  integer, intent(IN) :: ISVALID
  integer, intent(OUT) :: ICLL, JCLL
  integer :: LMILOC(1), L, I, J, ILN, JLN
  integer :: I1, I2, J1, J2
  real(8) :: RADLA(LA_GLOBAL)
  
  ! *** FOR THE FIRST CALL                     
  RADLA(2:LA_GLOBAL) = SQRT((XCLL-XCOR_Global(2:LA_GLOBAL,5))**2+(YCLL-YCOR_Global(2:LA_GLOBAL,5))**2)     ! *** Delme - This is inefficient.  Todo - Make more efficient.  use LEC, LWC, LSC, LNC
  LMILOC = MINLOC(RADLA(2:LA_GLOBAL))
  ILN = IL_GL(LMILOC(1)+1)    !I OF THE NEAREST CELL FOR DRIFTER
  JLN = JL_GL(LMILOC(1)+1)    !J OF THE NEAREST CELL FOR DRIFTER     

  ! *** DETERMINE THE CELL CONTAINING THE DRIFTER WITHIN 9 CELLS: LLA(NCEL)
  I1 = max(1,ILN-1)
  I2 = min(ILN+1,ICM_Global)
  J1 = max(1,JLN-1)
  J2 = min(JLN+1,JCM_Global)
  do J = J1,J2
    do I = I1,I2
      L = LIJ_Global(I,J)
      if( L < 2 ) CYCLE
      if( INSIDECELL_GL(L, XCLL, YCLL) )then
        ICLL = I
        JCLL = J
        return
      endif
    enddo
  enddo

  ! *** Cooridnates not found.  Determin the action
  if( ISVALID > 0 )then
    ! *** Flag as invalid
    ICLL = -999
    JCLL = -999
  else
    PRINT *,'THIS CELL IS OUTSIDE THE DOMAIN:',NCEL
    call STOPP('Invalid XY in SUBSET.INP')
  endif
  
END SUBROUTINE

!---------------------------------------------------------------------------!
!> @details AREA CALCULATION OF A POLYGON WITH GIVEN VEXTICES (XC,YC)
!
!> @param[in] XC(:) - x vertices
!> @param[in] YC(:) - y vertices
!> @param[out] AREA - area of the cell
!
!> @author  Zander Mausolff (not original author, modified function for clarity)
!> @date  2/17/2020
!---------------------------------------------------------------------------!
 SUBROUTINE AREACAL(XC,YC,AREA) 
 
  implicit none
  ! *** Passed in variables
  real(8),intent(IN)  :: XC(:)
  real(8),intent(IN)  :: YC(:)
  real(8),intent(OUT) :: AREA
  ! *** Local variables
  real(8) :: XVEC(2),YVEC(2)
  integer :: NPOL,K
  
  ! *** Get the number of vertices in polygon
  NPOL = SIZE(XC)
  AREA = 0.0
  
  XVEC(1) = XC(2) - XC(1)
  YVEC(1) = YC(2) - YC(1)
  do K = 3,NPOL
    XVEC(2) = XC(K) - XC(1)
    YVEC(2) = YC(K) - YC(1)
    AREA    = AREA + 0.5*ABS( XVEC(1)*YVEC(2) - XVEC(2)*YVEC(1) )
    XVEC(1) = XVEC(2)
    YVEC(1) = YVEC(2)
  enddo
END SUBROUTINE AREACAL

!---------------------------------------------------------------------------!
!> @details Function returns true if the x,y coordinates fall 
!! within the specified cell.
!
!> @param[in] L - Cell index
!> @param[in] XM - x cell centroid coordinate
!> @param[in] YM - y cell centroid coordinate
!> @param[out] INSIDE - True/False if coordinates fall inside/outside of cell
!
!> @author  Zander Mausolff (not original author, modified function for clarity)
!> @date  2/17/2020
!---------------------------------------------------------------------------!
FUNCTION INSIDECELL(L,XM,YM) RESULT(INSIDE)   

  implicit none

  ! *** Passed in an return variables
  logical(4) :: INSIDE
  integer, intent(IN) :: L
  real(8), intent(IN) :: XM
  real(8), intent(IN) :: YM
  ! *** Local variables
  real(8) :: XC(6),YC(6),AREA2
  
  XC(1)   = XM 
  YC(1)   = YM 
  XC(2:5) = XCOR(L, 1:4)
  YC(2:5) = YCOR(L, 1:4)
  XC(6)   = XC(2)
  YC(6)   = YC(2)

  call AREACAL(XC,YC,AREA2)
  
  if( ABS(AREA2-AREA(L)) <= 1D-6 )then
    INSIDE = .TRUE.
  else 
    INSIDE = .FALSE.
  endif
  
 END FUNCTION INSIDECELL
!---------------------------------------------------------------------------!
FUNCTION INSIDECELL_GL(L,XM,YM) RESULT(INSIDE)   

  implicit none

  ! *** Passed in an return variables
  logical(4) :: INSIDE
  integer, intent(IN) :: L
  real(8), intent(IN) :: XM
  real(8), intent(IN) :: YM
  ! *** Local variables
  real(8) :: XC(6),YC(6),AREA2
  
  XC(1)   = XM 
  YC(1)   = YM 
  XC(2:5) = XCOR_Global(L, 1:4)
  YC(2:5) = YCOR_Global(L, 1:4)
  XC(6)   = XC(2)
  YC(6)   = YC(2)

  call AREACAL(XC,YC,AREA2)
  
  if( ABS(AREA2-AREA_Global(L)) <= 1D-6 )then
    INSIDE = .TRUE.
  else 
    INSIDE = .FALSE.
  endif
  
END FUNCTION INSIDECELL_GL
!---------------------------------------------------------------------------!

SUBROUTINE AREA_CENTRD

  implicit none
  
  ! *** DETERMINING CELL CENTROID OF ALL CELLS
  ! *** AND CALCULATING THE AREA OF EACH CELL
  integer :: I,J,K, IIN, JIN, ii
  real(8) :: XC(4),YC(4),AREA2
  integer :: Q
  integer :: l_global, l_local, l
  
  ! *** ****************************
  if( process_id == master_id )then
      write(*,'(A)')'READING CORNERS.INP'
      open(UCOR,FILE = 'corners.inp',ACTION = 'READ')
      call SKIPCOM(UCOR, '*')
  endif
  
  ! *** At this point 'LA' is local to each process
  ! *** Want to read in using 'global' LA
  allocate(XCOR_Global(LCM_Global,5),YCOR_Global(LCM_Global,5),AREA_Global(LCM_Global))
  ! *** Setup local arrays 
  allocate(XCOR(LCM, 5), YCOR(LCM, 5), Area(LCM))
  
  XCOR_Global = 0
  YCOR_Global = 0
  AREA_GLobal = 0 !*** @todo not sure we even need at this time
  XCOR = 0
  YCOR = 0
  AREA = 0
  do Q = 1,LA_Global-1
     if( process_id == master_id )then 
         read(UCOR,*,err = 998) IIN,JIN,(XCOR_Global(LIJ_Global(IIN,JIN),K),YCOR_Global(LIJ_Global(IIN,JIN),K),K = 1,4)
     endif
      
     ! *** Setup xcor_global, ycor_global cell centroid and area_global
     if(process_id == master_id )then
        l_global = LIJ_Global(IIN,JIN)
        ! *** Select the 4 corners for x/y directions
        xc(1:4) = xcor_global(l_global, 1:4)
        yc(1:4) = ycor_global(l_global, 1:4)
        ! *** Calculatet area of the cell
        call AREACAL(xc,yc,area2)
        ! *** Set to global value
        Area_Global(l_global) = area2
        ! ***  STORE THE CELL CENTROID IN INDEX = 5
        XCOR_Global(l_global, 5) = 0.25*SUM(XC)        
        YCOR_Global(l_global, 5) = 0.25*SUM(YC)
     endif
     
     ! *** Broadcast read in I,J from the file for local remapping to determine local xcor/ycor
     call Broadcast_Scalar(IIN, master_id)
     call Broadcast_Scalar(JIN, master_id)
     ! *** Map to the local I,J
     I = IG2IL(IIN)
     J = JG2JL(JIN)
    
     call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),1), master_id)
     call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),2), master_id)
     call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),3), master_id)
     call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),4), master_id)
     call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),5), master_id)
     
     call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),1), master_id)
     call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),2), master_id)
     call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),3), master_id)
     call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),4), master_id)
     call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),5), master_id)
     
     ! *** Make sure we are only operating on data local to a process
     ! *** This works because IG2IL and JG2JL put zeros or ICM/JCM if the cell is not part of the domain
     if( I >0 .and. I <= IC )then
         if( J > 0 .and. J <= JC )then 
             l_global = LIJ_Global(IIN,JIN)
             l_local = LIJ(I,J)
             
             XCOR(l_local,1:4) = XCOR_Global(l_global, 1:4)
             YCOR(l_local,1:4) = YCOR_Global(l_global, 1:4)
             
             XC(1:4) = XCOR(l_local,1:4)  ! Access local only of the global xcor
             YC(1:4) = YCOR(l_local,1:4)
             call AREACAL(XC,YC,AREA2)
             AREA(l_local) = AREA2
             ! ***  STORE THE CELL CENTROID IN INDEX = 5
             XCOR(l_local,5) = 0.25*SUM(XC)        
             YCOR(l_local,5) = 0.25*SUM(YC)
         endif
     endif

  enddo
  
  !Call WriteBreak(mpi_log_unit)
  !!DO ii = 1, LA
  !!    write(mpi_log_unit,'(a,f7.2)') 'XCOR:   ', l, XCOR(ii,5),YCOR(ii,5)
  !!End do
  !Call WriteBreak(mpi_log_unit)

  !Call Broadcast_Array(XCOR_Global, master_id)
  !Call Broadcast_Array(YCOR_Global, master_id)
  
  if( process_id == master_id )then 
      100 close(UCOR)
  endif
    
  return
  998 STOP 'CORNERS.INP READING ERROR!'
END SUBROUTINE AREA_CENTRD

SUBROUTINE DIST2LINE(L,IP,X0,Y0,IPOINT,OFFSET,D,X3,Y3)
  ! *** COMPUTES THE DISTANCE FROM A POINT TO A LINE SEGMENT (X1,Y1),(X2,Y2) TO A POINT (X0,Y0) 
  ! *** RETURNS THE DISTANCE D AND THE INTERSECTION POINT ON THE LINE (X3,Y3)

  ! *** IP IS THE SIDE NUMBER, C1-C4 ARE THE POINT LOCATIONS AND ORDER
  ! *** 
  ! ***    C2    2   C3
  ! *** 
  ! ***    1          3
  ! *** 
  ! ***    C1   4    C4
  
  integer(IK4),intent(IN)  :: L,IP,IPOINT
  real(RKD)   ,intent(IN)  :: X0,Y0,OFFSET
  real(RKD)   ,intent(OUT) :: D,X3,Y3
    
  integer(IK4) :: I1,I2
  real(RKD)    :: H,XDEL,YDEL,ANG,EPSILON

  if( IP == 4 )then
    I1 = 4
    I2 = 1
  else
    I1 = IP
    I2 = IP+1
  endif
  XDEL = XCOR(L,I2) - XCOR(L,I1)
  YDEL = YCOR(L,I2) - YCOR(L,I1)
    
  D = YDEL*X0 - XDEL*Y0 + XCOR(L,I2)*YCOR(L,I1) - YCOR(L,I2)*XCOR(L,I1)
  H = SQRT(XDEL*XDEL + YDEL*YDEL)
  if( H < 1E-6 )then
    D = 0
    return
  endif
      
  ! *** SIGNED DISTANCE  <0 - LEFT OF LINE, >0 - RIGHT OF LINE
  D = D/H
  
  if( IPOINT > 0 )then
    !IF( ABS(XDEL) < 1E-12 )then
    !  M = 1E32
    !ELSEIF( ABS(YDEL) > 1E-12 )then
    !  M = YDEL/XDEL
    !ELSE
    !  M = 1E-32
    !ENDIF
    ANG = ATAN2(YDEL,XDEL)
    ANG = ANG + 0.5*PI
    
    !DS = SIGN(1.,XDEL)
    X3 = X0 + COS(ANG)*D  !*DS
    Y3 = Y0 + SIN(ANG)*D  !*DS

    ! *** CHECK RANGES BUT ADD A SMALL BUFFER FOR ROUNDOFF (EPSILON)
    EPSILON = 1E-12*X3
    if( (X3 + EPSILON) < min(XCOR(L,I1),XCOR(L,I2)) .or. (X3 - EPSILON) > max(XCOR(L,I1),XCOR(L,I2)) )then
      ! *** ON THE LINE BUT NOT IN THE SEGMENT
      D = 1E32
      return
    elseif( (Y3 + EPSILON) < min(YCOR(L,I1),YCOR(L,I2)) .or. (Y3 - EPSILON) > max(YCOR(L,I1),YCOR(L,I2)) )then
      ! *** ON THE LINE BUT NOT IN THE SEGMENT
      D = 1E32
      return
    endif
  
    ! *** RECOMPUTE X3,Y3 WITH OFFSET
    X3 = X0 + COS(ANG)*(D-OFFSET)  !*DS
    Y3 = Y0 + SIN(ANG)*(D-OFFSET)  !*DS

  endif

  return
    
END SUBROUTINE DIST2LINE

!---------------------------------------------------------------------------!
!> @details Function returns true if the the line (x1,y1) to (x2,y2)
!! is blocked by any cell mask 
!
!> @param[in] X1,Y1 - line begin coordinates
!> @param[in] X2,Y2 - line end coordinates
!> @param[out] XSECT - True/False if the line is blocked by the cell masks
!
!> @author  Nghiem Tien Lam
!> @date  04/27/2021
!> @modified Paul M. Craig
!---------------------------------------------------------------------------!
function BLOCKED(X1, Y1, X2, Y2, icell, jcell) result(intersect) 
  implicit none
  real(RKD), intent(in) :: X1, Y1, X2, Y2
  integer,   intent(in) :: icell, jcell
  integer  :: i, j, imin, imax, jmin, jmax, l
  
  logical :: intersect
    
  intersect = .false.

  ! *** Search in a 3 by 3 cell area
  imin = max(icell - 1,1)
  imax = min(icell + 1,IC)
  jmin = max(jcell - 1,1)
  jmax = min(jcell + 1,JC)
  do j = jmin, jmax
    do i = imin, imax
      L = LIJ(i,j)
      if( UMASK(L)  ==  1 )then
        if( L > 0 )then
          intersect = isintersect(X1, Y1, X2, Y2, XCOR(L,1), YCOR(L,1), XCOR(L,2), YCOR(L,2))
          if( intersect) return
        endif
      endif
      if( VMASK(L)  ==  1 )then
        if( L > 0 )then
          intersect = isintersect(X1, Y1, X2, Y2, XCOR(L,1), YCOR(L,1), XCOR(L,4), YCOR(L,4))
          if( intersect) return
        endif
      endif
    enddo
  enddo
  
end function

FUNCTION ISINTERSECT(X1,Y1,X2,Y2,X3,Y3,X4,Y4) RESULT(XSECT) 
    implicit none
    real(kind = RKD),intent(IN)  :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
    logical :: XSECT
    logical :: L134,L234,L123,L124

    L134 = ISCCW(X1,Y1,X3,Y3,X4,Y4)
    L234 = ISCCW(X2,Y2,X3,Y3,X4,Y4)
    L123 = ISCCW(X1,Y1,X2,Y2,X3,Y3)
    L124 = ISCCW(X1,Y1,X2,Y2,X4,Y4)
#ifdef GNU    
    XSECT = (L134 .NEQV. L234) .AND. (L123 .NEQV. L124)
#else
    XSECT = (L134 /= L234) .AND. (L123 /= L124)
#endif
END FUNCTION

FUNCTION ISCCW(X1,Y1,X2,Y2,X3,Y3) RESULT(CCW) 
    implicit none
    logical :: CCW
    real(kind = RKD),intent(IN)  :: X1,Y1,X2,Y2,X3,Y3
    CCW = (Y3-Y1)*(X2-X1) > (Y2-Y1)*(X3-X1)
END FUNCTION

END MODULE
