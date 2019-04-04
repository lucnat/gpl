
! In terminal kann man den exit code bekommen via echo $?
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use functions
  implicit none
  real, parameter :: tol = 1.0e-14

  call do_multiple_polylog_tests() 

CONTAINS
  
  subroutine do_multiple_polylog_tests()
    integer :: m2(2), m3(3)
    integer :: bla
    complex(kind=prec) :: x2(2), x3(3)
    complex(kind=prec) :: res, ref
    
    print*, 'doing multiple polylog tests...'
    
    m2 = (/ 1,1 /)
    x2 = cmplx((/ 0.3156498673740053, 0.3431255827785649 /))
    ref = cmplx(0.022696600480693277651633)
    call check_multiple_polylog(m2,x2,ref)

    m2 = (/ 1,1 /)
    x2 = cmplx((/ 0.03, 0.5012562893380046 /))
    ref = cmplx(0.00023134615630308335448329926098409)
    call check_multiple_polylog(m2,x2,ref)
    
    m3 = (/ 2,1,2 /)
    x3 = cmplx((/ 0.03, 0.5012562893380046, 55.3832 /))
    ref = cmplx(0.000023446106415452030937059124671151)
    call check_multiple_polylog(m3,x3,ref)

  end subroutine do_multiple_polylog_tests

  subroutine check_multiple_polylog(m,x,ref)
    ! checks multiple_polylog(m,x) against reference value ref
    integer :: m(:)
    complex(kind=prec) :: x(:)
    complex(kind=prec) :: res, ref
    real(kind=prec) :: delta
    res = multiple_polylog(m,x)
    delta = abs((res-ref)/ref)
    if(delta < tol) then
      print*, 'passed with delta = ', delta
    else 
      print*, 'not passed with delta = ', delta, 'for arguments: '
      print*, 'm=',m,'x=',x,'failed'
      print*, 'res=',res,'ref=',ref
      stop 1
    end if

  end subroutine check_multiple_polylog
  
  
END PROGRAM TEST
 