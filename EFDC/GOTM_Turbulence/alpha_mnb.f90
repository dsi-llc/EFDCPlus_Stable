#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update dimensionless alpha's\label{sec:alpha}
!
! !INTERFACE:
   subroutine alpha_mnb(nlev,NN,SS,ith)
!
! !DESCRIPTION:
! This subroutine updates the dimensionless numbers $\alpha_M$, $\alpha_N$,
! and $\alpha_b$ according to \eq{alphaMN}. Note that according to \eq{Nbar}
! and \eq{NbarVertical} the following identities are valid
! \begin{equation}
!  \label{alphaIdentities}
!    \alpha_M = \overline{S}^2 \comma
!    \alpha_N = \overline{N}^2 \comma
!    \alpha_b = \overline{T}   \point
! \end{equation}
!
!
! !USES:
  use turbulence,  only:     tke,eps,kb
  use turbulence,  only:     as,an,at
  implicit none
!
! !INPUT parameterS:
  integer,  intent(in)      :: nlev,ith
  REALTYPE, intent(in)      :: NN(0:nlev),SS(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer              :: i
  REALTYPE             :: tau2

!-----------------------------------------------------------------------
!BOC

  do i = 0,nlev
     tau2   = tke(i,ith)*tke(i,ith) / ( eps(i,ith)*eps(i,ith) )
     as(i,ith)  = tau2 * SS(i)
     an(i,ith)  = tau2 * NN(i)
     at(i,ith)  = tke(i,ith)/eps(i,ith) * kb(i,ith)/eps(i,ith)

!    clip negative values
     as(i,ith) = max(as(i,ith),1.e-10*_ONE_)
     at(i,ith) = max(at(i,ith),1.e-10*_ONE_)
  enddo

  return
end subroutine alpha_mnb

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
