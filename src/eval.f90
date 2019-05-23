 
PROGRAM eval
  use globals
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()

  res = GPL([1,2,3])
  ! res = pending_integral(cmplx([1,0]),1,cmplx([2,3]))
  
  print*, res
  
END PROGRAM eval
