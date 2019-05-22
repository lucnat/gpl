 
PROGRAM eval
  use globals
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()

  res = GPL(cmplx([1.0,-0.618034,0.1]))
  print*, res
  
END PROGRAM eval
