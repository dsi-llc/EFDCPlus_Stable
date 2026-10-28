#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic velocity variances \label{sec:variances}
!
! !INTERFACE:
   subroutine variances(nlev,SSU,SSV,ith)
!
! !DESCRIPTION:

!  Using \eq{bijVertical} and the solution shown in \eq{b13} and
!  the variances of the turbulent velocity fluctations can be
!  evaluated according to
!  \begin{equation}
!    \label{variances}
!    \begin{array}{rcl}
!      \dfrac{\langle u'^2 \rangle}{k}
!      &=& \dfrac23 + \dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}+a_3\right) \nu_t \left( \partder{U}{z} \right)^2
!                          -\dfrac23 a_2 \nu_t \left( \partder{V}{z} \right)^2
!                          -\dfrac43 a_5 G \right)
!      \comma
!      \\[7mm]
!      \dfrac{\langle v'^2 \rangle}{k}
!      &=& \dfrac23
!      +\dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}+a_3\right) \nu_t \left(\partder{V}{z}\right)^2
!                          -\dfrac23 a_2 \nu_t \left(\partder{U}{z}\right)^2
!                          -\dfrac43 a_5 G \right)
!      \comma
!      \\[7mm]
!      \dfrac{\langle w'^2 \rangle}{k}
!      &=& \dfrac23
!      +\dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}-a_3\right) P
!        +\dfrac83 a_5 G
!      \right)
!      \comma
!    \end{array}
!  \end{equation}
!  where the diffusivities are computed according to \eq{nu}
!  (also see \sect{sec:cmueC} and \sect{sec:cmueD}),
!  and the buoyancy production, $G$, follows from \eq{computeG}.
!
! !USES:
  use turbulence,  only:     uu,vv,ww
  use turbulence,  only:     tke,eps,P,B,num
  use turbulence,  only:     cc1,ct1,a2,a3,a5
  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
  integer,  intent(in)                 :: nlev
  
! thread index
  integer,  intent(in)                 :: ith

! square of shear frequency (1/s^2)
! (from u- and v-component)
  REALTYPE, intent(in)                :: SSU(0:nlev),SSV(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

   integer                             :: i
   REALTYPE                            :: N,Nt
   REALTYPE                            :: fac1,fac2,fac3,fac4,fac5

!
!-----------------------------------------------------------------------
!BOC

   N    =   0.5*cc1
   Nt   =   ct1

   do i=0,nlev

      fac1 = 2./3.
      fac2 = 1.0/( N*eps(i,ith) )
      fac3 = a2/3.0 + a3
      fac4 = a2/3.0 - a3
      fac5 = 2./3.*a2

      uu(i) = tke(i,ith)*( fac1 + fac2*( fac3*num(i,ith)*SSU(i)                &
                          - fac5*num(i,ith)*SSV(i) - 4./3.*a5*B(i,ith) ) )

      vv(i) = tke(i,ith)*( fac1 + fac2*( fac3*num(i,ith)*SSV(i)                &
                          - fac5*num(i,ith)*SSU(i) - 4./3.*a5*B(i,ith) ) )

      ww(i) = tke(i,ith)*( fac1 + fac2*( fac4*P(i,ith) + 8./3.*a5*B(i,ith) ) )

   enddo

   return
   end subroutine variances

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
