 
PROGRAM eval
  use GPL
  implicit none
  
  complex(kind=prec) :: res 

#ifdef DEBUG
  call parse_cmd_args()
#endif

  res = G([-1.,0.,0.,0.,1.])
  print*, res

  ! res = pending_integral(cmplx([1,3]),2,cmplx([2,4]))
  ! print*, res

END PROGRAM eval

