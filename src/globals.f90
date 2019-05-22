
MODULE globals
  implicit none

  integer, parameter :: prec = selected_real_kind(15,32)  
  integer, parameter :: GPLInfinity = 30   ! the default outermost expansion order for MPLs
  real, parameter :: epsilon = 1e-15       ! used for the small imaginary part
  real, parameter :: zero = 1e-15          ! values smaller than this count as zero
  real, parameter :: pi = 3.14159265358979323846

  integer :: verbosity = 100

END MODULE globals
