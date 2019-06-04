 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()
  
  ! res = GPL([0,2,3,4])
  ! print*, res
  
  res = pending_integral(cmplx([2,3]),2,cmplx([0,4]))
  print*, res

  
END PROGRAM eval

