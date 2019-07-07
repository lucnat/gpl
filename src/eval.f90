 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  use shuffle

  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()

  res = GPL([1.1,3.0,1.0])
  print*, res
  
  ! res = pending_integral(cmplx([1,3]),2,cmplx([2,4]))
  ! print*, res

END PROGRAM eval

