#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic kb-equation\label{sec:kbalgebraic}
!
! !INTERFACE:
   subroutine kbalgebraic(nlev,ith)
!
! !DESCRIPTION:
! The algebraic equation for $k_b$ simply assumes equilibrium in \eq{kbeq},
! \begin{equation}
!   \label{kbEquilibrium}
!   P_b = \epsilon_b
!   \point
! \end{equation}
! This equation can be re-written as
! \begin{equation}
!   \label{kbAgebraic}
!   k_b = \dfrac{k_b \epsilon}{k \epsilon_b} \dfrac{k}{\epsilon} P_b
!       = r \dfrac{k}{\epsilon} P_b = c_b \dfrac{k}{\epsilon} P_b
!   \comma
! \end{equation}
! where we used the definition of the time scale ratio $r$ in
! \eq{DefR}, and assumed that $r = c_b$ is a constant.
!

!
! !USES:
   use turbulence,  only:     tke,eps,kb,Pb
   use turbulence,  only:     ctt,kb_min

  implicit none
!
! !INPUT parameterS:

! number of vertical layers
   integer,  intent(in)                 :: nlev

! thread index
   integer,  intent(in)                 :: ith
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

   integer                             :: i

!-----------------------------------------------------------------------
!BOC

   do i = 0,nlev
      kb(i,ith) = ctt*tke(i,ith)/eps(i,ith)*Pb(i,ith)

      !  clip at kb_min
      kb(i,ith) = max(kb(i,ith),kb_min)
   enddo

   return
   end subroutine kbalgebraic

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
