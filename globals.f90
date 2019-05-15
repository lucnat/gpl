
MODULE globals
  implicit none

  integer, parameter :: prec = selected_real_kind(15,32)  
  real, parameter :: zero = 1e-15          ! values smaller than this count as zero
  integer, parameter :: GPLInfinity = 30   ! the default outermost expansion order for MPLs
  real, parameter :: epsilon = 1e-15       ! used for the small imaginary part

END MODULE globals
