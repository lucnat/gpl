
! In terminal kann man den exit code bekommen via echo $?
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use mpl_module
  use gpl_module
  implicit none

  complex(kind=prec) :: res 
  real, parameter :: tol = 1.0e-14
  integer :: testnr

  call do_MPL_tests() 
  call do_GPL_tests()

CONTAINS

  subroutine check(res, ref)
    complex(kind=prec) :: res, ref
    real(kind=prec) :: delta

    delta = abs((res-ref)/ref)
    if(delta < tol) then
      print*, testnr,'passed with delta = ', delta
    else 
      print*, testnr,'not passed with delta = ', delta
      stop 1
    end if
  end subroutine check
  
  subroutine do_MPL_tests()
    integer :: m2(2), m3(3)
    integer :: bla
    complex(kind=prec) :: x2(2), x3(3)
    complex(kind=prec) :: res, ref
    print*, 'doing MPL tests...'
    
    testnr = 1
    m2 = (/ 1,1 /)
    x2 = dcmplx((/ 0.3156498673740053, 0.3431255827785649 /))
    res = MPL(m2,x2)
    ref = dcmplx(0.022696600480693277651633)
    call check(res,ref)

    testnr = 2
    m2 = (/ 1,1 /)
    x2 = dcmplx((/ 0.03, 0.5012562893380046 /))
    res = MPL(m2,x2)
    ref = dcmplx(0.00023134615630308335448329926098409)
    call check(res,ref)
    
    testnr = 3
    m3 = (/ 2,1,2 /)
    x3 = dcmplx((/ 0.03, 0.5012562893380046, 55.3832 /))
    res = MPL(m3,x3)
    ref = dcmplx(0.000023446106415452030937059124671151)
    call check(res,ref)

  end subroutine do_MPL_tests

  subroutine check_one_MPL(m,x,ref)
    ! checks MPL(m,x) against reference value ref
    integer :: m(:)
    complex(kind=prec) :: x(:)
    complex(kind=prec) :: res, ref

    res = MPL(m,x)
    call check(res,ref)
  end subroutine check_one_MPL

  subroutine do_GPL_tests()
    integer :: m(2), k
    complex(kind=prec) :: z(2), y, res, ref
    print*, 'doing GPL tests...'

    testnr = 11
    m = (/ 1,1 /)
    z = dcmplx((/ 1.3d0, 1.1d0 /))
    y = 0.4
    k = 2
    res = GPL(m,z,y,k)
    ref = dcmplx(0.0819393734128676)
    call check(res,ref)

    testnr = 12
    m = (/ 3,2 /)
    z = dcmplx((/ 1.3d0, 1.1d0 /))
    y = 0.4
    k = 2
    res = GPL(m,z,y,k)
    ref = dcmplx(0.01592795952537145)
    call check(res,ref)

  end subroutine do_GPL_tests

  subroutine check_one_GPL(m,z,y,k,ref)
    ! checks one GPL(m,z,y,k) against reference value ref
    integer :: m(:), k
    complex(kind=prec) :: z(:), y, res, ref

    res = GPL(m,z,y,k)
    call check(res,ref)
  end subroutine check_one_GPL
END PROGRAM TEST
 