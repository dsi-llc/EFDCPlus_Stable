! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE XYIJCONV
!Author: Dang Huu Chung

USE GLOBAL  
USE INFOMOD
Use Variables_MPI
Use Variables_MPI_Mapping
#ifdef _MPI
    Use Broadcast_Routines
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE XY2IJ(CEL, VALID) 
  TYPE(CELL),INTENT(INOUT) :: CEL
  INTEGER, INTENT(IN), OPTIONAL :: VALID
  INTEGER :: N, NPMAX, ISVALID
  
  IF( PRESENT(VALID) )THEN
    ISVALID = 1
  ELSE
    ISVALID = 0
  ENDIF
  
  NPMAX = SIZE(CEL.XCEL)
  
  DO N=1,NPMAX 
    CALL CONTAINERIJ_GL(N, CEL.XCEL(N), CEL.YCEL(N), CEL.ICEL(N), CEL.JCEL(N), ISVALID)
  ENDDO
  
END SUBROUTINE

SUBROUTINE CONTAINERIJ(NCEL, XCLL, YCLL, ICLL, JCLL, ISVALID)   
  INTEGER,INTENT(IN)  :: NCEL
  REAL(8), INTENT(IN) :: XCLL, YCLL
  INTEGER, INTENT(IN) :: ISVALID
  INTEGER, INTENT(OUT) :: ICLL, JCLL
  INTEGER :: LMILOC(1), L, I, J, ILN, JLN
  INTEGER :: I1, I2, J1, J2
  REAL(8) :: RADLA(LA)
  
  ! *** FOR THE FIRST CALL                     
  RADLA(2:LA) = SQRT((XCLL-XCOR(2:LA,5))**2+(YCLL-YCOR(2:LA,5))**2)     ! *** Delme - This is inefficient.  Todo - Make more efficient.  Use LEC, LWC, LSC, LNC
  LMILOC = MINLOC(RADLA(2:LA))
  ILN = IL(LMILOC(1)+1)    !I OF THE NEAREST CELL FOR DRIFTER
  JLN = JL(LMILOC(1)+1)    !J OF THE NEAREST CELL FOR DRIFTER     

  ! *** DETERMINE THE CELL CONTAINING THE DRIFTER WITHIN 9 CELLS: LLA(NCEL)
  I1 = MAX(1,ILN-1)
  I2 = MIN(ILN+1,ICM)
  J1 = MAX(1,JLN-1)
  J2 = MIN(JLN+1,JCM)
  DO J=J1,J2
    DO I=I1,I2
      L = LIJ(I,J)
      IF( L < 2 ) CYCLE
      IF( INSIDECELL(L, XCLL, YCLL) )THEN
        ICLL = I
        JCLL = J
        RETURN
      ENDIF
    ENDDO
  ENDDO

  ! *** Cooridnates not found.  Determin the action
  IF( ISVALID > 0 )THEN
    ! *** Flag as invalid
    ICLL = -999
    JCLL = -999
  ELSE
    PRINT *,'THIS CELL IS OUTSIDE THE DOMAIN:',NCEL
    CALL STOPP('Invalid XY in SUBSET.INP')
  ENDIF
  
END SUBROUTINE

SUBROUTINE CONTAINERIJ_GL(NCEL, XCLL, YCLL, ICLL, JCLL, ISVALID)   
  INTEGER,INTENT(IN)  :: NCEL
  REAL(8), INTENT(IN) :: XCLL, YCLL
  INTEGER, INTENT(IN) :: ISVALID
  INTEGER, INTENT(OUT) :: ICLL, JCLL
  INTEGER :: LMILOC(1), L, I, J, ILN, JLN
  INTEGER :: I1, I2, J1, J2
  REAL(8) :: RADLA(LA_GLOBAL)
  
  ! *** FOR THE FIRST CALL                     
  RADLA(2:LA_GLOBAL) = SQRT((XCLL-XCOR_Global(2:LA_GLOBAL,5))**2+(YCLL-YCOR_Global(2:LA_GLOBAL,5))**2)     ! *** Delme - This is inefficient.  Todo - Make more efficient.  Use LEC, LWC, LSC, LNC
  LMILOC = MINLOC(RADLA(2:LA_GLOBAL))
  ILN = IL_GL(LMILOC(1)+1)    !I OF THE NEAREST CELL FOR DRIFTER
  JLN = JL_GL(LMILOC(1)+1)    !J OF THE NEAREST CELL FOR DRIFTER     

  ! *** DETERMINE THE CELL CONTAINING THE DRIFTER WITHIN 9 CELLS: LLA(NCEL)
  I1 = MAX(1,ILN-1)
  I2 = MIN(ILN+1,ICM_Global)
  J1 = MAX(1,JLN-1)
  J2 = MIN(JLN+1,JCM_Global)
  DO J=J1,J2
    DO I=I1,I2
      L = LIJ_Global(I,J)
      IF( L < 2 ) CYCLE
      IF( INSIDECELL_GL(L, XCLL, YCLL) )THEN
        ICLL = I
        JCLL = J
        RETURN
      ENDIF
    ENDDO
  ENDDO

  ! *** Cooridnates not found.  Determin the action
  IF( ISVALID > 0 )THEN
    ! *** Flag as invalid
    ICLL = -999
    JCLL = -999
  ELSE
    PRINT *,'THIS CELL IS OUTSIDE THE DOMAIN:',NCEL
    CALL STOPP('Invalid XY in SUBSET.INP')
  ENDIF
  
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
 
  Implicit none
  ! *** Passed in variables
  REAL(8),INTENT(IN)  :: XC(:)
  REAL(8),INTENT(IN)  :: YC(:)
  REAL(8),INTENT(OUT) :: AREA
  ! *** Local variables
  REAL(8) :: XVEC(2),YVEC(2)
  INTEGER :: NPOL,K
  
  ! *** Get the number of vertices in polygon
  NPOL = SIZE(XC)
  AREA = 0.0
  
  XVEC(1) = XC(2) - XC(1)
  YVEC(1) = YC(2) - YC(1)
  DO K=3,NPOL
    XVEC(2) = XC(K) - XC(1)
    YVEC(2) = YC(K) - YC(1)
    AREA    = AREA + 0.5*ABS( XVEC(1)*YVEC(2) - XVEC(2)*YVEC(1) )
    XVEC(1) = XVEC(2)
    YVEC(1) = YVEC(2)
  ENDDO
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

  Implicit none

  ! *** Passed in an return variables
  LOGICAL(4) :: INSIDE
  INTEGER, INTENT(IN) :: L
  REAL(8), INTENT(IN) :: XM
  REAL(8), INTENT(IN) :: YM
  ! *** Local variables
  REAL(8) :: XC(6),YC(6),AREA2
  
  XC(1)   = XM 
  YC(1)   = YM 
  XC(2:5) = XCOR(L, 1:4)
  YC(2:5) = YCOR(L, 1:4)
  XC(6)   = XC(2)
  YC(6)   = YC(2)

  CALL AREACAL(XC,YC,AREA2)
  
  IF( ABS(AREA2-AREA(L)) <= 1D-6 )THEN
    INSIDE=.TRUE.
  ELSE 
    INSIDE=.FALSE.
  ENDIF
  
 END FUNCTION INSIDECELL
!---------------------------------------------------------------------------!
FUNCTION INSIDECELL_GL(L,XM,YM) RESULT(INSIDE)   

  Implicit none

  ! *** Passed in an return variables
  LOGICAL(4) :: INSIDE
  INTEGER, INTENT(IN) :: L
  REAL(8), INTENT(IN) :: XM
  REAL(8), INTENT(IN) :: YM
  ! *** Local variables
  REAL(8) :: XC(6),YC(6),AREA2
  
  XC(1)   = XM 
  YC(1)   = YM 
  XC(2:5) = XCOR_Global(L, 1:4)
  YC(2:5) = YCOR_Global(L, 1:4)
  XC(6)   = XC(2)
  YC(6)   = YC(2)

  CALL AREACAL(XC,YC,AREA2)
  
  IF( ABS(AREA2-AREA_Global(L)) <= 1D-6 )THEN
    INSIDE=.TRUE.
  ELSE 
    INSIDE=.FALSE.
  ENDIF
  
END FUNCTION INSIDECELL_GL
!---------------------------------------------------------------------------!

SUBROUTINE AREA_CENTRD

  Implicit None
  
  ! *** DETERMINING CELL CENTROID OF ALL CELLS
  ! *** AND CALCULATING THE AREA OF EACH CELL
  INTEGER :: I,J,K, IIN, JIN, ii
  REAL(8) :: XC(4),YC(4),AREA2
  Integer :: Q
  Integer :: l_global, l_local, l
  
  ! *** ****************************
  if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING CORNERS.INP'
      OPEN(UCOR,FILE='corners.inp',ACTION='READ')
      CALL SKIPCOM(UCOR, '*')
  end if
  
  ! *** At this point 'LA' is local to each process
  ! *** Want to read in using 'global' LA
  ALLOCATE(XCOR_Global(LCM_Global,5),YCOR_Global(LCM_Global,5),AREA_Global(LCM_Global))
  ! *** Setup local arrays 
  Allocate(XCOR(LCM, 5), YCOR(LCM, 5), Area(LCM))
  
  XCOR_Global = 0
  YCOR_Global = 0
  AREA_GLobal = 0 !*** @todo not sure we even need at this time
  XCOR = 0
  YCOR = 0
  AREA = 0
  DO Q=1,LA_Global-1
     if( process_id == master_id )THEN 
         READ(UCOR,*,ERR=998) IIN,JIN,(XCOR_Global(LIJ_Global(IIN,JIN),K),YCOR_Global(LIJ_Global(IIN,JIN),K),K=1,4)
     end if
      
     ! *** Setup xcor_global, ycor_global cell centroid and area_global
     if(process_id == master_id )then
        l_global = LIJ_Global(IIN,JIN)
        ! *** Select the 4 corners for x/y directions
        xc(1:4) = xcor_global(l_global, 1:4)
        yc(1:4) = ycor_global(l_global, 1:4)
        ! *** Calculatet area of the cell
        Call AREACAL(xc,yc,area2)
        ! *** Set to global value
        Area_Global(l_global) = area2
        ! ***  STORE THE CELL CENTROID IN INDEX=5
        XCOR_Global(l_global, 5) = 0.25*SUM(XC)        
        YCOR_Global(l_global, 5) = 0.25*SUM(YC)
     end if
     
     ! *** Broadcast read in I,J from the file for local remapping to determine local xcor/ycor
     call Broadcast_Scalar(IIN, master_id)
     call Broadcast_Scalar(JIN, master_id)
     ! *** Map to the local I,J
     I = IG2IL(IIN)
     J = JG2JL(JIN)
    
     Call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),1), master_id)
     Call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),2), master_id)
     Call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),3), master_id)
     Call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),4), master_id)
     Call Broadcast_Scalar(XCOR_Global(LIJ_Global(IIN,JIN),5), master_id)
     
     Call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),1), master_id)
     Call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),2), master_id)
     Call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),3), master_id)
     Call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),4), master_id)
     Call Broadcast_Scalar(YCOR_Global(LIJ_Global(IIN,JIN),5), master_id)
     
     ! *** Make sure we are only operating on data local to a process
     ! *** This works because IG2IL and JG2JL put zeros or ICM/JCM if the cell is not part of the domain
     IF( I >0 .AND. I <= IC )THEN
         IF( J > 0 .AND. J <= JC )THEN 
             l_global = LIJ_Global(IIN,JIN)
             l_local = LIJ(I,J)
             
             XCOR(l_local,1:4) = XCOR_Global(l_global, 1:4)
             YCOR(l_local,1:4) = YCOR_Global(l_global, 1:4)
             
             XC(1:4) = XCOR(l_local,1:4)  ! Access local only of the global xcor
             YC(1:4) = YCOR(l_local,1:4)
             CALL AREACAL(XC,YC,AREA2)
             AREA(l_local) = AREA2
             ! ***  STORE THE CELL CENTROID IN INDEX=5
             XCOR(l_local,5) = 0.25*SUM(XC)        
             YCOR(l_local,5) = 0.25*SUM(YC)
         ENDIF
     ENDIF

  ENDDO
  
  Call WriteBreak(mpi_log_unit)
  !DO ii = 1, LA
  !    write(mpi_log_unit,'(a,f7.2)') 'XCOR:   ', l, XCOR(ii,5),YCOR(ii,5)
  !End do
  Call WriteBreak(mpi_log_unit)

  !Call Broadcast_Array(XCOR_Global, master_id)
  !Call Broadcast_Array(YCOR_Global, master_id)
  
  if( process_id == master_id )THEN 
      100 CLOSE(UCOR)
  endif
    
  RETURN
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
  
  INTEGER(IK4),INTENT(IN)  :: L,IP,IPOINT
  REAL(RKD)   ,INTENT(IN)  :: X0,Y0,OFFSET
  REAL(RKD)   ,INTENT(OUT) :: D,X3,Y3
    
  INTEGER(IK4) :: I1,I2
  REAL(RKD)    :: H,XDEL,YDEL,ANG,EPSILON

  IF( IP == 4 )THEN
    I1 = 4
    I2 = 1
  ELSE
    I1 = IP
    I2 = IP+1
  ENDIF
  XDEL = XCOR(L,I2) - XCOR(L,I1)
  YDEL = YCOR(L,I2) - YCOR(L,I1)
    
  D = YDEL*X0 - XDEL*Y0 + XCOR(L,I2)*YCOR(L,I1) - YCOR(L,I2)*XCOR(L,I1)
  H = SQRT(XDEL*XDEL + YDEL*YDEL)
  IF( H < 1E-6 )THEN
    D = 0
    RETURN
  ENDIF
      
  ! *** SIGNED DISTANCE  <0 - LEFT OF LINE, >0 - RIGHT OF LINE
  D = D/H
  
  IF( IPOINT > 0 )THEN
    !IF( ABS(XDEL) < 1E-12 )THEN
    !  M = 1E32
    !ELSEIF( ABS(YDEL) > 1E-12 )THEN
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
    IF( (X3 + EPSILON) < MIN(XCOR(L,I1),XCOR(L,I2)) .OR. (X3 - EPSILON) > MAX(XCOR(L,I1),XCOR(L,I2)) )THEN
      ! *** ON THE LINE BUT NOT IN THE SEGMENT
      D = 1E32
      RETURN
    ELSEIF( (Y3 + EPSILON) < MIN(YCOR(L,I1),YCOR(L,I2)) .OR. (Y3 - EPSILON) > MAX(YCOR(L,I1),YCOR(L,I2)) )THEN
      ! *** ON THE LINE BUT NOT IN THE SEGMENT
      D = 1E32
      RETURN
    ENDIF
  
    ! *** RECOMPUTE X3,Y3 WITH OFFSET
    X3 = X0 + COS(ANG)*(D-OFFSET)  !*DS
    Y3 = Y0 + SIN(ANG)*(D-OFFSET)  !*DS

  ENDIF

  RETURN
    
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
      if ( UMASK(L)  ==  1 ) then
        if( L > 0 )then
          intersect = isintersect(X1, Y1, X2, Y2, XCOR(L,1), YCOR(L,1), XCOR(L,2), YCOR(L,2))
          if (intersect) return
        endif
      endif
      if ( VMASK(L)  ==  1 ) then
        if( L > 0 )then
          intersect = isintersect(X1, Y1, X2, Y2, XCOR(L,1), YCOR(L,1), XCOR(L,4), YCOR(L,4))
          if (intersect) return
        endif
      endif
    enddo
  enddo
  
end function

FUNCTION ISINTERSECT(X1,Y1,X2,Y2,X3,Y3,X4,Y4) RESULT(XSECT) 
    IMPLICIT NONE
    LOGICAL :: XSECT
    REAL(kind = RKD),INTENT(IN)  :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
    XSECT = (ISCCW(X1,Y1,X3,Y3,X4,Y4) /= ISCCW(X2,Y2,X3,Y3,X4,Y4)) &
      .AND. (ISCCW(X1,Y1,X2,Y2,X3,Y3) /= ISCCW(X1,Y1,X2,Y2,X4,Y4))
END FUNCTION

FUNCTION ISCCW(X1,Y1,X2,Y2,X3,Y3) RESULT(CCW) 
    IMPLICIT NONE
    LOGICAL :: CCW
    REAL(kind = RKD),INTENT(IN)  :: X1,Y1,X2,Y2,X3,Y3
    CCW = (Y3-Y1)*(X2-X1) > (Y2-Y1)*(X3-X1)
END FUNCTION

END MODULE
