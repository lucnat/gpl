 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()
  
  res = GPL([2,1,3,4]) ! here's an example where we need to shuffle away from last place
  print*, res
  
  ! res = pending_integral(cmplx([1,0]) + epsilon ,1, cmplx([3,2]) + epsilon)
  ! print*, res
  
END PROGRAM eval

