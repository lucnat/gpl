 
PROGRAM eval
  use globals
  use gpl_module
  use utils
  use shuffle

  implicit none
  
  integer :: i
  complex(kind=prec) :: res 
  ! complex(kind=prec) :: a(3), s(2)
  ! complex(kind=prec) :: alpha(product((/(i,i=1,size(a)+size(s))/))/  & 
  !     (product((/(i,i=1,size(a))/))*product((/(i,i=1,size(s))/))), & 
  !     size(a) + size(s))

  call parse_cmd_args()

  ! a = cmplx((/1,2,1/))
  ! s = cmplx((/4.0,42e50/))

  ! alpha = shuffle_product(a,s)
  ! call print_logical_matrix(alpha == 42e50)
  ! print*, find_first_true(alpha(6,:) == 42e50)

  ! res = GPL([2,1,3,4]) ! here's an example where we need to shuffle away from last place
  ! print*, res
  
  res = pending_integral(cmplx([1,3]),2,cmplx([2,4]))
  print*, res

END PROGRAM eval

