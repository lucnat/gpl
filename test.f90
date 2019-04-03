
! im terminal kann man den exit code bekommen via echo $?

PROGRAM TEST
  use functions
  implicit none
  real, parameter :: tol = 1.0e-15
  integer :: m(2)
  complex(kind=prec) :: x(2)
  complex(kind=prec) :: res, ref

  print*, 'testing multiple polylog...'
  
  m = (/ 1,1 /)
  x = cmplx((/ 0.3156498673740053, 0.3431255827785649/))
  ref = cmplx(0.022696600480693277651633)
  call check_multiple_polylog(m,x,ref)
  
  m = (/ 1,1 /)
  x = cmplx((/ 0.3156498673740053, 0.3431255827785649/))
  ref = cmplx(0.022696600480693277651633)
  call check_multiple_polylog(m,x,ref)

CONTAINS

  subroutine check_multiple_polylog(m,x,ref)
    ! checks multiple_polylog(m,x) against reference value ref
    integer :: m(2)
    complex(kind=prec) :: x(2)
    complex(kind=prec) :: res, ref
    res = multiple_polylog(m,x)
    if(abs((res-ref)/ref) < tol) then
      print*, 'passed'
    else 
      print*, 'm=',m,'x=',x,'failed'
      print*, 'res=',res,'ref=',ref
      stop 1
    end if

  end subroutine check_multiple_polylog
  
  
END PROGRAM TEST
 