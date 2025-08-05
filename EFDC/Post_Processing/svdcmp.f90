! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)

  ! ***  FROM NUMERICAL RECIPES
  ! CHANGE RECORD

  use GLOBAL,only: LWC
  implicit none
  
  integer :: M,N
  integer :: MP,NP,I,L,K,J,ITS,NM
  
  real :: A(MP,NP),W(NP),V(NP,NP)
  real,allocatable :: RV1(:)
  real :: G,S,SCALE,ANORM,F,SQRTSSS,H,C,Y,Z,X

  if( .not. allocated(RV1) )then
    allocate(RV1(N))
    RV1 = 0.
  endif
  G = 0.0
  SCALE = 0.0
  ANORM = 0.0
  do 25 I = 1,N
    L = I+1
    RV1(I) = SCALE*G
    G = 0.0
    S = 0.0
    SCALE = 0.0
    if( I <= M )then
      do 11 K = I,M
        SCALE = SCALE+ABS(A(K,I))
     11     continue
      if( SCALE /= 0.0 )then
        do 12 K = I,M
          A(K,I) = A(K,I)/SCALE
          S = S+A(K,I)*A(K,I)
     12       continue
        F = A(I,I)
        SQRTSSS = SQRT(S)
        G = -SIGN(SQRTSSS,F)
        H = F*G-S
        A(I,I) = F-G
        if( I /= N )then
          do 15 J = L,N
            S = 0.0
            do 13 K = I,M
              S = S+A(K,I)*A(K,J)
     13           continue
            F = S/H
            do 14 K = I,M
              A(K,J) = A(K,J)+F*A(K,I)
     14           continue
     15         continue
        endif
        do 16 K= I,M
          A(K,I) = SCALE*A(K,I)
     16       continue
      endif
    endif
    W(I) = SCALE *G
    G = 0.0
    S = 0.0
    SCALE = 0.0
    if( (I <= M) .and. (I /= N) )then
      do 17 K = L,N
        SCALE = SCALE+ABS(A(I,K))
     17     continue
      if( SCALE /= 0.0 )then
        do 18 K = L,N
          A(I,K) = A(I,K)/SCALE
          S = S+A(I,K)*A(I,K)
     18       continue
        F = A(I,L)
        G = -SIGN(SQRT(S),F)
        H = F*G-S
        A(I,L) = F-G
        do 19 K = L,N
          RV1(K) = A(I,K)/H
     19       continue
        if( I /= M )then
          do 23 J = L,M
            S = 0.0
            do 21 K = L,N
              S = S+A(J,K)*A(I,K)
     21           continue
            do 22 K = L,N
              A(J,K) = A(J,K)+S*RV1(K)
     22           continue
     23         continue
        endif
        do 24 K = L,N
          A(I,K) = SCALE*A(I,K)
     24       continue
      endif
    endif
    ANORM = MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
     25 continue
  do 32 I = N,1,-1
    if( I < N )then
      if( G /= 0.0 )then
        do 26 J = L,N
          V(J,I) = (A(I,J)/A(I,L))/G
     26       continue
        do 29 J = L,N
          S = 0.0
          do 27 K = L,N
            S = S+A(I,K)*V(K,J)
     27         continue
          do 28 K = L,N
            V(K,J) = V(K,J)+S*V(K,I)
     28         continue
     29       continue
      endif
      do 31 J = L,N
        V(I,J) = 0.0
        V(J,I) = 0.0
     31     continue
    endif
    V(I,I) = 1.0
    G = RV1(I)
    L = I
     32 continue
  do 39 I = N,1,-1
    L = I+1
    G = W(I)
    if( I < N )then
      do 33 J = L,N
        A(I,J) = 0.0
     33     continue
    endif
    if( G /= 0.0 )then
      G = 1.0/G
      if( I /= N )then
        do 36 J = L,N
          S = 0.0
          do 34 K = L,M
            S = S+A(K,I)*A(K,J)
     34         continue
          F = (S/A(I,I))*G
          do 35 K = I,M
            A(K,J) = A(K,J)+F*A(K,I)
     35         continue
     36       continue
      endif
      do 37 J = I,M
        A(J,I) = A(J,I)*G
     37     continue
    else
      do 38 J= I,M
        A(J,I) = 0.0
     38     continue
    endif
    A(I,I) = A(I,I)+1.0
     39 continue
  do 49 K = N,1,-1
    do 48 ITS = 1,30
      do 41 L = K,1,-1
        NM = LWC(L)
        if( (ABS(RV1(L))+ANORM) == ANORM)  GOTO 2
        if( (ABS(W(NM))+ANORM) == ANORM)  GOTO 1
     41     continue
      1     C = 0.0
      S = 1.0
      do 43 I = L,K
        F = S*RV1(I)
        if( (ABS(F)+ANORM) /= ANORM )then
          G = W(I)
          H = SQRT(F*F+G*G)
          W(I) = H
          H = 1.0/H
          C = (G*H)
          S = -(F*H)
          do 42 J = 1,M
            Y = A(J,NM)
            Z = A(J,I)
            A(J,NM) = (Y*C)+(Z*S)
            A(J,I) = -(Y*S)+(Z*C)
     42         continue
        endif
     43     continue
      2     Z = W(K)
      if( L == K )then
        if( Z < 0.0 )then
          W(K) = -Z
          do 44 J = 1,N
            V(J,K) = -V(J,K)
     44         continue
        endif
        GOTO 3
      endif
      if( ITS == 60 )then
        write(6,603)
        deallocate(RV1)
        return
      endif
    603 FORMAT('  NO CONVERGENCE IN 60 ITERATIONS')
      X = W(L)
      NM = K-1
      Y = W(NM)
      G = RV1(NM)
      H = RV1(K)
      F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
      G = SQRT(F*F+1.0)
      F = ((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
      C = 1.0
      S = 1.0
      do 47 J = L,NM
        I = J+1
        G = RV1(I)
        Y = W(I)
        H = S*G
        G = C*G
        Z = SQRT(F*F+H*H)
        RV1(J) = Z
        C = F/Z
        S = H/Z
        F= (X*C)+(G*S)
        G = -(X*S)+(G*C)
        H = Y*S
        Y = Y*C
        do 45 NM = 1,N
          X = V(NM,J)
          Z = V(NM,I)
          V(NM,J)= (X*C)+(Z*S)
          V(NM,I) = -(X*S)+(Z*C)
     45       continue
        Z = SQRT(F*F+H*H)
        W(J) = Z
        if( Z /= 0.0 )then
          Z = 1.0/Z
          C = F*Z
          S = H*Z
        endif
        F= (C*G)+(S*Y)
        X = -(S*G)+(C*Y)
        do 46 NM = 1,M
          Y = A(NM,J)
          Z = A(NM,I)
          A(NM,J)= (Y*C)+(Z*S)
          A(NM,I) = -(Y*S)+(Z*C)
     46       continue
     47     continue
      RV1(L) = 0.0
      RV1(K) = F
      W(K) = X
     48   continue
      3   continue
     49 continue
  deallocate(RV1)
  
  return
END

