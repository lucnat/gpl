 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()
  
  ! res = GPL([0,1,3,2])
  ! print*, res
  
  res = pending_integral(cmplx([1,0]) + epsilon ,3, cmplx([0,0,2]) + epsilon)
  print*, res
  
END PROGRAM eval

