#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mtridiagonal --- solving the system\label{sec:tridiagonal}
!
! !INTERFACE:
   MODULE mtridiagonal
!
! !DESCRIPTION:
!
!  Solves a linear system of equations with a tridiagonal matrix
!  using Gaussian elimination.
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_tridiagonal, tridiagonal, clean_tridiagonal
!
! !PUBLIC DATA MEMBERS:
   REALTYPE, dimension(:,:), allocatable     :: au,bu,cu,du
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
!  private data members
   REALTYPE, private, dimension(:,:),allocatable  ::  ru,qu
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory
!
! !INTERFACE:
   subroutine init_tridiagonal(N)
!
! !DESCRIPTION:
!  This routines allocates memory necessary to perform the Gaussian
!  elimination.
!
! !USES:
   use GLOBAL, only:NTHREADS
   implicit none
!
! !INPUT parameterS:
   integer, intent(in)                 :: N
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_tridiagonal'
   allocate(au(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating au)'
   au = _ZERO_

   allocate(bu(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating bu)'
   bu = _ZERO_

   allocate(cu(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating cu)'
   cu = _ZERO_

   allocate(du(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating du)'
   du = _ZERO_

   allocate(ru(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating ru)'
   ru = _ZERO_

   allocate(qu(0:N, NTHREADS),stat = rc)
   if( rc /= 0) stop 'init_tridiagonal: Error allocating qu)'
   qu = _ZERO_

   return
   end subroutine init_tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Simplified Gaussian elimination
!
! !INTERFACE:
   subroutine tridiagonal(ith,N,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix structure is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
! The method used here is the simplified Gauss elimination, also called
! \emph{Thomas algorithm}.
!
! !USES:
   use GLOBAL, only:NTHREADS
   implicit none
!
! !INPUT parameterS:
   integer, intent(in)                 :: ith, N, fi, lt
!
! !OUTPUT parameterS:
   REALTYPE                            :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ru(lt,ith) = au(lt,ith)/bu(lt,ith)
   qu(lt,ith) = du(lt,ith)/bu(lt,ith)

   do i = lt-1,fi+1,-1
      ru(i,ith) = au(i,ith)/(bu(i,ith)-cu(i,ith)*ru(i+1,ith))
      qu(i,ith) = (du(i,ith)-cu(i,ith)*qu(i+1,ith))/(bu(i,ith)-cu(i,ith)*ru(i+1,ith))
   enddo

   qu(fi,ith) = (du(fi,ith)-cu(fi,ith)*qu(fi+1,ith))/(bu(fi,ith)-cu(fi,ith)*ru(fi+1,ith))

   value(fi) = qu(fi,ith)
   do i = fi+1,lt
      value(i) = qu(i,ith)-ru(i,ith)*value(i-1)
   enddo


   return
   end subroutine tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: De-allocate memory
!
! !INTERFACE:
   subroutine clean_tridiagonal()
!
! !DESCRIPTION:
!  De-allocates memory allocated in init\_tridiagonal.
!
! !USES:
   implicit none
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_tridiagonal'

   if( allocated(au)) deallocate(au)

   if( allocated(bu)) deallocate(bu)

   if( allocated(cu)) deallocate(cu)

   if( allocated(du)) deallocate(du)

   if( allocated(ru)) deallocate(ru)

   if( allocated(qu)) deallocate(qu)

   return
   end subroutine clean_tridiagonal
!EOC
!-----------------------------------------------------------------------

   end module mtridiagonal

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
