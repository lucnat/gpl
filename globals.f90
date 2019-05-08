
MODULE globals
  implicit none

  integer, parameter :: prec = selected_real_kind(15,32)  
  real :: zero = 1e-15          ! values smaller than this count as zero
  integer :: GPLInfinity = 30   ! the default outermost expansion order for MPLs

END MODULE globals
