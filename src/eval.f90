
PROGRAM eval
  use globals
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 

  call parse_cmd_args()
  print*, verb
  res = GPL(cmplx([1,2,5]))
  print*, res
  
END PROGRAM eval
