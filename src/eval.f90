PROGRAM eval
  use globals
  use utils
  use shuffle
  use maths_functions
  use mpl_module
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 

  res = G_flat(cmplx((/2,1/)),cmplx(3))
  print*, res

END PROGRAM eval