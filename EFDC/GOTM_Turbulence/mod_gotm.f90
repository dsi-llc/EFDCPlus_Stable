  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  !> @details Initiates and calls GOTM subroutines
  !> @author
  !---------------------------------------------------------------------------!
  Module Mod_GOTM

  use GLOBAL
  use Variables_MPI
  use OMP_LIB
  
  implicit none

  Contains

  Subroutine Init_GOTM

  use turbulence, only: init_turbulence_DSI
  use mtridiagonal, only: init_tridiagonal
  integer :: ks

  ks = kc - 1
  call init_turbulence_DSI(kc)
  call init_tridiagonal(kc)

  End Subroutine Init_GOTM

  Subroutine Advance_GOTM(ISTL_)

  use turbulence, only: do_turbulence
  use turbulence, only: eps_min,k_min
  use turbulence, only: tke1d => tke, eps1d => eps, GL1d => Lgs
  use turbulence, only: num1d => num, nuh1d => nuh
  use turbulence, only: z0s_min
  
  implicit none
  
  ! *** Local variables
  integer,intent(in) :: ISTL_
  integer :: K, L, LP, ND, LN, LE, LW, LS, LL
  integer :: i, j, itr, NLAYER, IT
  integer :: itz0b = 10 !< number of iterations for hydrodynamic bottom roughness
  real(rkd) :: rr, HPK2, rr_s
  real(rkd) :: charnock_val = 1400.
  real(rkd) :: GPO = 9.81
  
  ! *** Variable interfacing with GOTM
  real(rkd) :: DEPTH, DELT, NNC, NNW, NNE, NNN, NNS, SSU, SSV, DZA
  real(rkd) :: GTAUS, GTAUB, GUC, GVC, z0s_gotm, z0b_gotm, cnpar = 1.
  real(rkd), Allocatable, Dimension(:,:) :: NN
  real(rkd), Allocatable, Dimension(:,:) :: SS
  real(rkd)  :: NN1d(0:KC)
  real(rkd)  :: SS1d(0:KC)
  real(rkd)  :: HPK1d(0:KC)

  If( .Not. Allocated(NN) )then
    allocate(NN(LCM,0:KC))
    allocate(SS(LCM,0:KC))
  Endif
  
  GTAUS = 0.0
  GTAUB = 0.0
  NN = 0.0
  SS = 0.0
  IT = 1

  ! *** Set up Time Step [s]
  If( ISDYNSTP == 0 )then
    DELT = DT
  Else
    DELT = DTDYN
  endif
  IT = 1   ! *** Initialize thread for non-OMP compilations
  
  !$OMP PARALLEL DEFAULT(SHARED) 

  ! *** Calculate Buoyancy Frequency Squared (NN)
  If( (ISTRAN(1) + ISTRAN(2)) >= 1 )then 
    If( ICALNN == 0 )then
      !$OMP DO PRIVATE(ND,K,L,LP,DZA)
      Do ND = 1,NDM
        Do K = 1,KS ! *** Check the initial and final value of the loop
          Do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            DZA = 0.5*(HPK(L,K+1) + HPK(L,K))
            NN(L,K) = -GPO*(B(L,K+1)-B(L,K))/DZA
          enddo
        enddo
      enddo
      !$OMP END DO
    else
      !$OMP DO PRIVATE(ND,K,L,LP,LS,LE,LW,LN,NNC,NNS,NNN,NNE,NNW,DZA)
      Do ND = 1,NDM
        Do K = 1,KS ! *** Check the initial and final value of the loop
          Do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            ! *** Computing NN based on new approach from Hans Burchard, 2002 to improve the stability
            LS = LSC(L)
            LE = LEC(L)
            LW = LWC(L)
            LN = LNC(L)
            
            DZA = 0.5*(HPK(L,K+1) + HPK(L,K))
            NNC = -GPO*(B(L,K+1)-B(L,K))/DZA
            
            DZA = 0.5*(HPK(LS,K+1) + HPK(LS,K))
            if( DZA > 0. )then
              NNS = -GPO*(B(LS,K+1)-B(LS,K))/DZA
            else
              NNS = NNC 
            endif
            DZA = 0.5*(HPK(LN,K+1) + HPK(LN,K))
            if( DZA > 0. )then
              NNN = -GPO*(B(LN,K+1)-B(LN,K))/DZA
            else
              NNN = NNC
            endif
            DZA = 0.5*(HPK(LE,K+1) + HPK(LE,K))
            if( DZA > 0. )then
              NNE = -GPO*(B(LE,K+1)-B(LE,K))/DZA
            else
              NNE = NNC
            endif
            DZA = 0.5*(HPK(LW,K+1) + HPK(LW,K))
            if( DZA > 0. )then
              NNW = -GPO*(B(LW,K+1)-B(LW,K))/DZA
            else
              NNW = NNC
            endif
    
            NN(L,K) = 1./3.*NNC + 1./6.*(NNS + NNN + NNE + NNW)
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
  endif

  ! *** Set BC's as GOTM Does
  !$OMP DO PRIVATE(L)
  Do L = 2,LA
    NN(L,KSZ(L)-1) = 0.
    NN(L,KC)       = 0.
  Enddo
  !$OMP END DO
  
  If( ICALSS == 0 )then
    !$OMP DO PRIVATE(ND,K,L,LP,HPK2)
    Do ND = 1,NDM
      Do K = 1,KS ! *** Check the initial and final value of the loop
        Do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)
          HPK2 = HPK(L,K)*HPK(L,K)
          SS(L,K) = (UCTR(L,K+1)-UCTR(L,K))**2/HPK2 + (VCTR(L,K+1)-VCTR(L,K))**2/HPK2
        enddo
      enddo
    enddo
    !$OMP END DO
  else
    if( ISTL_ /= 2 )then
      !$OMP DO PRIVATE(ND,K,L,LP,SSU,SSV)
      Do ND = 1,NDM
        Do K = 1,KS ! *** Check the initial and final value of the loop
          Do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)
            ! *** Computing SS based on new approach from Hans Burchard, 2002 to improve the stability
            SSU = 0.5*(                                                     &
              (cnpar*    (UCTR(L,K+1) -UCTR(L,K))*(UCTR(L,K+1) -UCTR2(L,K))+           &
              (1.-cnpar)*(UCTR2(L,K+1)-UCTR2(L,K))*(UCTR2(L,K+1)-UCTR(L,K)))            &
                   /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K)                   &
              +(cnpar*   (UCTR(L,K+1)-UCTR(L,K))*(UCTR2(L,K+1)-UCTR(L,K))+            &
              (1.-cnpar)*(UCTR2(L,K+1)-UCTR2(L,K))*(UCTR(L,K+1)-UCTR2(L,K)))           &
                   /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K+1)                 &
              )
      
            SSV = 0.5*(                                                     &
              (cnpar*    (VCTR(L,K+1) -VCTR(L,K)) *(VCTR(L,K+1) -VCTR2(L,K))+           &
              (1.-cnpar)*(VCTR2(L,K+1)-VCTR2(L,K))*(VCTR2(L,K+1)-VCTR(L,K)))            &
                    /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K)                  &
              +(cnpar*   (VCTR(L,K+1) -VCTR(L,K)) *(VCTR2(L,K+1)-VCTR(L,K))+            &
              (1.-cnpar)*(VCTR2(L,K+1)-VCTR2(L,K))*(VCTR(L,K+1) -VCTR2(L,K)))           &
                    /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K+1)                &
              )
      
            SS(L,K) = SSU + SSV
          enddo
        enddo
      enddo
    !$OMP END DO
    else
      !$OMP DO PRIVATE(ND,K,L,LP,SSU,SSV)
      Do ND = 1,NDM
        Do K = 1,KS ! *** Check the initial and final value of the loop
          Do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)
            ! *** Computing SS based on new approach from Hans Burchard, 2002 to improve the stability
            SSU = 0.5*(                                                     &
              (cnpar*    (UCTR(L,K+1) -UCTR(L,K))*(UCTR(L,K+1) -UCTR1(L,K))+           &
              (1.-cnpar)*(UCTR1(L,K+1)-UCTR1(L,K))*(UCTR1(L,K+1)-UCTR(L,K)))            &
                   /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K)                   &
              +(cnpar*   (UCTR(L,K+1)-UCTR(L,K))*(UCTR1(L,K+1)-UCTR(L,K))+            &
              (1.-cnpar)*(UCTR1(L,K+1)-UCTR1(L,K))*(UCTR(L,K+1)-UCTR1(L,K)))           &
                   /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K+1)                 &
              )
      
            SSV = 0.5*(                                                     &
              (cnpar*    (VCTR(L,K+1) -VCTR(L,K)) *(VCTR(L,K+1) -VCTR1(L,K))+           &
              (1.-cnpar)*(VCTR1(L,K+1)-VCTR1(L,K))*(VCTR1(L,K+1)-VCTR(L,K)))            &
                    /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K)                  &
              +(cnpar*   (VCTR(L,K+1) -VCTR(L,K)) *(VCTR1(L,K+1)-VCTR(L,K))+            &
              (1.-cnpar)*(VCTR1(L,K+1)-VCTR1(L,K))*(VCTR(L,K+1) -VCTR1(L,K)))           &
                    /(0.5*(HPK(L,K+1)+ HPK(L,K)))/HPK(L,K+1)                &
              )
      
            SS(L,K) = SSU + SSV
          enddo
        enddo
      enddo
      !$OMP END DO
    Endif
  Endif

  ! *** Set BC's for SS as GOTM Does
  !$OMP DO PRIVATE(L)  
  Do L = 2,LA
    SS(L,KSZ(L)-1) = SS(L,KSZ(L))
    SS(L,KC)       = SS(L,KS)
  Enddo
  !$OMP END DO

  ! *** Main Loop Over Elements
  !$OMP DO PRIVATE(IT, ITR, L, HPK1D, NN1D, SS1D, Z0S_GOTM, z0b_gotm, RR, GUC, GVC, NLAYER, DEPTH)  &
  !$OMP    FIRSTPRIVATE(GTAUB, GTAUS)
  Do L = 2,LA
    !$  IT = OMP_GET_THREAD_NUM() + 1
  
    ! *** Initialize 1D array
    HPK1d = 0.
    SS1d  = 0.
    NN1d  = 0.
    
    ! *** Surface Roughness Height
    GUC = 0.5*(TSX(L) + TSX(LEC(L)))
    GVC = 0.5*(TSY(L) + TSY(LNC(L)))
    GTAUS = (GUC*GUC + GVC*GVC)**(1./4.)
    
    if( charnock == 1 )then
      z0s_gotm = charnock_val*GTAUS**2/GPO
      if( z0s_gotm.lt.z0s_min ) z0s_gotm = z0s_min ! *** GOTM formula
    elseif( charnock == 2 )then
      z0s_gotm = max(ZSRE(L),z0s_min)   ! *** Using surface roughness computed from COARE 3.6
    else
      z0s_gotm = z0s_min
    endif
  
    ! *** Bottom Roughness Height and Bottom Friction Velocity
    If( IFRICTION == 0 )then
      GUC = 0.5*(TBX(L) + TBX(LEC(L)))
      GVC = 0.5*(TBY(L) + TBY(LNC(L)))
      GTAUB = (GUC*GUC + GVC*GVC)**(1./4.)
      z0b_gotm = ZBR(L)
    else
      Do itr = 1, itz0b
        z0b_gotm = 0.1*AVO/max(AVO,GTAUB) + 0.03*ZBR(L)
        rr = VKC/(log((z0b_gotm + HPK(L,KSZ(L)))/z0b_gotm))
        GTAUB = rr*sqrt(UCTR(L,KSZ(L))*UCTR(L,KSZ(L)) + VCTR(L,KSZ(L))*VCTR(L,KSZ(L)))
      enddo
    endif

    ! *** Set up Depth [m] and Layer Thicknesses [m]
    NLAYER = KC - KSZ(L) + 1
    DEPTH = HP(L)
    HPK1d(1:NLAYER) = HPK(L,KSZ(L):KC)

    ! *** Set Up 1-D Arrays
    if( ISRESTI == 0 .and. NITER == 1 )then 
      num1d(1:NLAYER,IT) = AVOXY(L)                      !< Vertical turbulent eddy viscosity     [m^2/s]
      nuh1d(1:NLAYER,IT) = AVBXY(L)                      !< Vertical mass diffusion coefficient   [m^2/s]
    else
      num1d(1:NLAYER,IT) = Av(L,KSZ(L):KC)*DEPTH         !< Vertical turbulent eddy viscosity     [m^2/s]
      nuh1d(1:NLAYER,IT) = Ab(L,KSZ(L):KC)*DEPTH         !< Vertical mass diffusion coefficient   [m^2/s]
    endif
    SS1d (0:NLAYER)    = SS   (L,KSZ(L)-1:KC)            !< Shear Frequency Squared               [1/s^2]
    NN1d (0:NLAYER)    = NN   (L,KSZ(L)-1:KC)            !< Buoyancy Frequency Squared            [1/s^2]
    tke1d(0:NLAYER,IT) = tke3d(L,KSZ(L)-1:KC)            !< Turbulent Kinetic Energy              [m^2/s^2]
    eps1d(0:NLAYER,IT) = eps3d(L,KSZ(L)-1:KC)            !< Turbulence Dissipation Rate           [m^2/s^3]
    GL1d (0:NLAYER,IT) = GL3d (L,KSZ(L)-1:KC)            !< Turbulence Length Scale               [m]
    
    ! *** Update Turbulence Model
    If( NLAYER > 1 )then
      call do_turbulence(IT, NLAYER, DELT, DEPTH, GTAUS, GTAUB, z0s_gotm, z0b_gotm, HPK1d, NN1d, SS1d)

      ! *** Update 3D Fields of TKE, EPS, Av, Ab
      tke3d(L,KSZ(L)-1:KC) = tke1d(0:NLAYER,IT)
      eps3d(L,KSZ(L)-1:KC) = eps1d(0:NLAYER,IT)
      Av   (L,KSZ(L)  :KC) = num1d(1:NLAYER,IT)/DEPTH + AVOXY(L)/DEPTH 
      ! *** In GOTM, the diffusion coefficient takes the same value for temperature, salinity, etc
      Ab  (L,KSZ(L):KC)    = nuh1d(1:NLAYER,IT)/DEPTH + AVBXY(L)/DEPTH
      GL3d(L,KSZ(L):KC)    = GL1d(1:NLAYER,IT)
      
      ! *** Zeroing turbulence quantities at layer KC
      tke3d(L,KC) = 0.
      eps3d(L,KC) = 0.
      GL3d (L,KC) = 0.
      Av(L,KC) = 0.
      Ab(L,KC) = 0.
    Endif
  enddo
  ! $OMP END SINGLE
  !$OMP END DO  
  !$OMP END PARALLEL
    
  ! *** Zeroing Ab at open boundary cells
  if( NBCSOP > 0 )then
    do LL = 1,NBCSOP
      L = LOBCS(LL)
      do K = 1,KS  
        Ab(L,K) = 0.0
      enddo  
    enddo 
   
    do K = 1,KS
      do LL = 1,NPBS
        L  = LPBS(LL)
        LN = LNC(L)
        tke3d(L,K) = tke3d(LN,K)
        eps3d(L,K) = eps3d(LN,K)
        GL3d (L,K) = GL3d (LN,K)
      enddo
    enddo
    do K = 1,KS
      do LL = 1,NPBW
        L  = LPBW(LL)
        LE = LEC(L)
        tke3d(L,K) = tke3d(LE,K)
        eps3d(L,K) = eps3d(LE,K)
        GL3d (L,K) = GL3d (LE,K)
      enddo
    enddo
    do K = 1,KS
      do LL = 1,NPBE
        L = LPBE(LL)
        LW = LWC(L)
        tke3d(L,K) = tke3d(LW,K)
        eps3d(L,K) = eps3d(LW,K)
        GL3d (L,K) = GL3d (LW,K)
      enddo
    enddo
    
    do K = 1,KS
      do LL = 1,NPBN
        L  = LPBN(LL)
        LS = LSC(L)
        tke3d(L,K) = tke3d(LS,K)
        eps3d(L,K) = eps3d(LS,K)
        GL3d (L,K) = GL3d (LS,K)
      enddo
    enddo
  endif

  End Subroutine Advance_GOTM

  End Module Mod_GOTM
