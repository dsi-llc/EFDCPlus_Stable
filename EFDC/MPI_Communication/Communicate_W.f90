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
! @details Communicates part of the W solution --> UHDY, VHDX
! @author Zander Mausolff
! @date 9/4/2019
!---------------------------------------------------------------------------!
   
SUBROUTINE COMMUNICATE_W

    Use GLOBAL
    Use Variables_MPI
    Use MPI
    
    Implicit None

    !***Local variables
    Integer :: I, J, II, L, K
    Integer :: istatus
    Integer :: length_arg 
    INTEGER(4) :: IERR
    
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WSENDW
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WSENDE
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WRECVE
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WRECVW
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WSENDN
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WSENDS
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WRECVN
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WRECVS

    IF(.NOT.ALLOCATED(WSENDW))THEN
        ALLOCATE(WSENDW(max_width_y*4*KCM))
        ALLOCATE(WSENDE(max_width_y*4*KCM))
        ALLOCATE(WRECVE(max_width_y*4*KCM))
        ALLOCATE(WRECVW(max_width_y*4*KCM))
        ALLOCATE(WSENDN(max_width_x*4*KCM))
        ALLOCATE(WSENDS(max_width_x*4*KCM))
        ALLOCATE(WRECVN(max_width_x*4*KCM))
        ALLOCATE(WRECVS(max_width_x*4*KCM))
        WSENDW = 0.0
        WSENDE = 0.0
        WRECVE = 0.0
        WRECVW = 0.0
        WSENDN = 0.0
        WSENDS = 0.0
        WRECVN = 0.0
        WRECVS = 0.0
    END IF

    !***Send West
    IF (nbr_west.NE.-1)THEN
        II = 0
        DO K = 1, KC
            DO I =3,4
                DO J = 3,JC-2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        WSENDW(II) = UHDY(L,K)
                    endif
                    II = II + 1
                    if( l > 0) then
                        WSENDW(II) = VHDX(L,K)
                    endif
                END DO
            END DO
        END DO
        IERR = 0
        length_arg = (JC-4)*4*KC
        CALL MPI_SEND(WSENDW,length_arg,mpi_real,nbr_west,process_id,  &
                      comm_2d, IERR)
    END IF

    !***East recv
    IF (nbr_east.NE.-1)THEN
    
        length_arg = (JC-4)*4*KC
        
        CALL MPI_RECV(WRECVE,length_arg, mpi_real, nbr_east, nbr_east,  &
                      comm_2d, ISTATUS, IERR)
        
        II = 0
        DO K = 1,KC
            DO I =IC-1,IC
                DO J = 3,JC-2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        UHDY(L,K) = WRECVE(II)
                    endif
                    II = II + 1
                    if( l > 0) then
                        VHDX(L,K) = WRECVE(II)
                    endif
                END DO
            END DO
        END DO
    END IF

    !***East Send
    IF (nbr_east.NE.-1)THEN
        ii = 0
        DO K = 1,KC
            DO I =IC-3,IC-2
                DO J = 3,JC-2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        WSENDE(II) = UHDY(L,K)
                    endif
                    II = II + 1
                    if( l > 0) then
                        WSENDE(II) = VHDX(L,K)
                    endif
                END DO
            END DO
        END DO
        IERR = 0
        length_arg = (JC-4)*4*KC
        CALL MPI_SEND(WSENDE,length_arg,mpi_real,nbr_east,process_id,  &
            comm_2d, IERR)
    END IF

    !***West recv
    IF (nbr_west.NE.-1)THEN
    
        length_arg = (JC-4)*4*KC
        
        CALL MPI_RECV(WRECVW,length_arg,mpi_real,nbr_west,nbr_west, &
            comm_2d,ISTATUS,IERR)
        II = 0
        DO K = 1,KC
            DO I =1,2
                DO J = 3,JC-2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        UHDY(L,K) = WRECVW(II)
                    endif
                    II = II + 1
                    if( l > 0) then
                        VHDX(L,K) = WRECVW(II)
                    endif
                END DO
            END DO
        END DO
    END IF

    !***North send
    IF (nbr_north.NE.-1)THEN
        II = 0
        DO K = 1,KC
            DO I = 1, IC
                DO J = JC-3,JC-2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        WSENDN(II) = UHDY(L,K)
                    endif
                    II = II + 1
                    if( l > 0) then
                        WSENDN(II) = VHDX(L,K)
                    endif
                END DO
            END DO
        END DO
        IERR = 0
        length_arg = IC*4*KC
        CALL MPI_SEND(WSENDN,length_arg,mpi_real,nbr_north,process_id, &
            comm_2d, IERR)
    END IF

    !***South recv
    IF (nbr_south.NE.-1)THEN
        length_arg = IC*4*KC
        CALL MPI_RECV(WRECVS,length_arg,mpi_real,nbr_south,nbr_south, &
            comm_2d,ISTATUS,IERR)
        II = 0
        DO K = 1,KC
            DO I = 1,IC
                DO J = 1,2
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        UHDY(L,K) = WRECVS(II)
                    endif
                    II = II + 1
                    if( l > 0) then
                        VHDX(L,K) = WRECVS(II)
                    endif
                END DO
            END DO
        END DO
    END IF

    !***South Send
    IF (nbr_south.NE.-1)THEN
        II = 0
        DO K = 1,KC
            DO I = 1, IC
                DO J = 3,4
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        WSENDS(II) = UHDY(L,K)
                    endif
                    II = II + 1
                    if( l > 0) then
                        WSENDS(II) = VHDX(L,K)
                    endif
                    
                END DO
            END DO
        END DO
        IERR = 0
        length_arg = IC*4*KC
        CALL MPI_SEND(WSENDS,length_arg,mpi_real,nbr_south,process_id, &
            comm_2d, IERR)
    END IF

    !***North recv
    IF (nbr_north.NE.-1)THEN
        
        length_arg = IC*4*KC
        
        CALL MPI_RECV(WRECVN,length_arg,mpi_real,nbr_north,nbr_north, &
            comm_2d,ISTATUS,IERR)
        II = 0
        DO K = 1,KC
            DO I = 1,IC
                DO J = JC-1,JC
                    L = LIJ(I,J)
                    II =II + 1
                    if( l > 0) then
                        UHDY(L,K) = WRECVN(II)
                    endif
                    II =II + 1
                    if( l > 0) then
                        VHDX(L,K) = WRECVN(II)
                    endif
                END DO
            END DO
        END DO
    END IF


END SUBROUTINE COMMUNICATE_W