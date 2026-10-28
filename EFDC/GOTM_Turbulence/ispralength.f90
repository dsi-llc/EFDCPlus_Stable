#include"cppdefs.h"
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Algebraic length-scale from ISPRAMIX \label{sec:ispramix}
!
! !INTERFACE:
   subroutine ispralength(nlev,NN,h,depth,ith)
!
! !DESCRIPTION:
!  This subroutine calculates the
!  lengthscale used in the ISPRAMIX model,
!  see \cite{EiflerSchrimpf92} and \cite{Demirovetal98}.
!  In both mixing regions (close to the surface and the bottom),
!  $l$ is obtained from the formula
!  \begin{equation}
!    \label {Lmixed}
!    l = \frac {\kappa \tilde z} {1+\frac {\kappa \tilde z} {c_2 \cdot h_m}}
!        (1-R_f)^e
!  \end{equation}
!  where $\tilde z$
!  is the distance from the interface (surface or bottom). The
!  fraction in (\ref{Lmixed})
!  predicts an approximation to a linear behavior of $l$ near boundaries
!  and a value proportional to the thickness of the mixed
!  layer far from the interface, $l = c_2 h_m$, where $c_2 = 0.065$
!  is estimated from experimental data as discussed in
!  \cite{EiflerSchrimpf92}.
!  The factor $(1-R_f)$, with the flux Richardson
!  number $R_f = -G/P$, accounts for the effect
!  of stratification on the length-scale.
!  The parameter $e$ is here a tuning parameter
!  (pers.\ comm.\ Walter Eifler, JRC, Ispra, Italy)
!  which is usually set to $e = 1$.
!
! !USES:
   use turbulence, only: Lgs,tke,k_min,eps_min,xRF,kappa,cde

   implicit none
!
! !INPUT parameterS:

!  number of vertical layers
   integer,  intent(in)                :: nlev
   
!  thread index
   integer,  intent(in)                :: ith

!  buoyancy frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  local depth (m)
   REALTYPE, intent(in)                :: depth

   ! !REVISION HISTORY:
!  Original author(s):  Manuel Ruiz Villarreal, Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer                    :: i,SLind,BLind,Index,Index2
  REALTYPE                   :: hms,hmb,db,ds
  REALTYPE                   :: kml,c2_i,c3_i
  REALTYPE                   :: l_min
!
!-------------------------------------------------------------------------
!BOC

   l_min = cde*k_min**1.5/eps_min

   kml   = 1.e-5
   c2_i  = 0.065

!  Calculation of surface mixed layer depth
   hms = 0.
   SLind = 1
   do i = nlev,1,-1
      hms = hms+h(i)
      if( tke(i,ith).le.kml )then
         SLind = i
         goto 500
      endif
   enddo
500  continue
!  Calculation of bottom mixed layer depth
   hmb = 0.
   BLind = nlev
   do i = 1,nlev
      hmb = hmb+h(i)
      if( tke(i,ith).le.kml )then
         BLind = i
         goto 501
      endif
   enddo
501  continue

! If there is no point where k < kml, the water column is assumed to be mixed.
   if( BLind.gt.SLind )then
      hms = 0.5*depth
      hmb = 0.5*depth
      BLind = int(nlev/2)
      SLind = int(nlev/2)+1
   endif

! Calculation of mixing length in bottom layer
   db = 0.
   do i = 1,BLind
      db = db+h(i)
      Lgs(i,ith) = kappa*db/(1.+kappa*db/(c2_i*hmb+L_min))*xRf(i,ith)**3
      if( Lgs(i,ith).lt.L_min) Lgs(i,ith) = L_min
   enddo

! Calculation of mixing length in surface layer
   ds = h(nlev)
   do i = nlev-1,SLind,-1
      ds = ds+h(i)
      Lgs(i,ith) = kappa*ds/(1.+kappa*ds/(c2_i*hms+L_min))*xRf(i,ith)**3
      if( Lgs(i,ith).lt.L_min) Lgs(i,ith) = L_min
   enddo

! Calculation of mixing length in the intermediate region

   c3_i = Lgs(SLind,ith)*sqrt(NN(SLind)/tke(SLind,ith))
   if( c3_i.lt.1e-10) c3_i = 0.
   Index = Slind-1
   do i = SLind-1,BLind+1,-1
      if( NN(i).le.0. )then
         Lgs(i,ith) = L_min
      else
         Lgs(i,ith) = max(c3_i*sqrt(tke(i,ith)/NN(i)),L_min)
         if( Lgs(i,ith).gt.Lgs(SLind,ith)) Lgs(i,ith) = Lgs(SLind,ith)
      endif
      if( Lgs(i,ith).eq.L_min )then
         Index = i
         goto 503
      endif
   enddo
503  continue
   c3_i = Lgs(BLind,ith)*sqrt(NN(BLind)/tke(BLind,ith))
   if( c3_i.lt.1e-10) c3_i = 0.
   Index2 = BLind+1
   do i = BLind+1,Index
      if( NN(i).le.0. )then
         Lgs(i,ith) = L_min
      else
         Lgs(i,ith) = max(c3_i*sqrt(tke(i,ith)/NN(i)),L_min)
         if(Lgs(i,ith).gt.Lgs(BLind,ith)) Lgs(i,ith) = Lgs(BLind,ith)
      endif
      if( Lgs(i,ith).eq.L_min )then
         Index2 = i
         goto 504
      endif
   enddo
504  continue
   do i = Index2+1,Index-1
      Lgs(i,ith) = L_min
   enddo

   return
   end subroutine ispralength
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
