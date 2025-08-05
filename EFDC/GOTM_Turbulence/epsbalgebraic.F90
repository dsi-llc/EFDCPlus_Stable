#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic epsilonb-equation\label{sec:epsbalgebraic}
!
! !INTERFACE:
   subroutine epsbalgebraic(nlev,ith)
!
! !DESCRIPTION:
! The algebraic equation for $\epsilon_b$, the molecular rate of
! destruction of buoyancy variance, see \eq{kbeq}, simply assumes a
! constant time scale ratio $r = c_b$, see \eq{DefR}. From
! this assumption, it follows immediately that
! \begin{equation}
!   \label{epsbAgebraic}
!     \epsilon_b = \dfrac{1}{c_b} \dfrac{\epsilon}{k} k_b
!   \point
! \end{equation}
!
! !USES:
  use turbulence,  only:     tke,eps,kb,epsb
  use turbulence,  only:     ctt,epsb_min

  implicit none
!
! !INPUT parameterS:

! number of vertical layers
  integer,  intent(in)                 :: nlev
  
! thread index
  integer,  intent(in)                 :: ith

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  REALTYPE                             :: one_over_ctt
  integer                              :: i
!
!-----------------------------------------------------------------------
!BOC

  one_over_ctt = 1.0D0/ctt

  do i = 0,nlev
     epsb(i,ith) = one_over_ctt*eps(i,ith)/tke(i,ith)*kb(i,ith)

!     clip at epsb_min
     epsb(i,ith) = max(epsb(i,ith),epsb_min)
  enddo

  return
  end subroutine epsbalgebraic
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
