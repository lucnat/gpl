 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()
  
  res = GPL([0,1,3,2])   ! should be 0.09593041677639341 - 0.8829351795197851*I
  print*, res
  
  ! res = pending_integral(cmplx([1,0]) + epsilon ,1, cmplx([3,2]) + epsilon)
  ! print*, res
  
END PROGRAM eval

