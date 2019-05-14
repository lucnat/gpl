
! In terminal kann man den exit code bekommen via echo $?
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use globals
  use utils
  use shuffle
  use mpl_module
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 
  real, parameter :: tol = 1.0e-14
  logical :: tests_successful = .true.

  ! call do_MPL_tests() 
  ! call do_GPL_tests()
  call do_shuffle_tests() ! put this somewhere else
  
  if(tests_successful) then
    print*, 'All tests passed. '
  else 
    print*, 'Some tests failed. '
    stop 1
  end if

CONTAINS
   
  subroutine check(res, ref)
    complex(kind=prec) :: res, ref
    real(kind=prec) :: delta

    delta = abs((res-ref)/ref)
    if(delta < tol) then
      print*, '  ',' passed with delta = ', delta
    else 
      print*, '  ',' failed with delta = ', delta
      tests_successful = .false.
    end if
  end subroutine check

  subroutine test_one_MPL(m,x,ref, test_id)
    integer :: m(:)
    complex(kind=prec) :: x(:), ref
    character(len=*) :: test_id
    
    print*, '  ', 'testing MPL ', test_id, ' ...'
    res = MPL(m,x)
    call check(res,ref)
  end subroutine test_one_MPL

  subroutine do_MPL_tests()
    complex(kind=prec) :: ref
    print*, 'doing MPL tests...'
    
    ref = dcmplx(0.022696600480693277651633)
    call test_one_MPL((/ 1,1 /),cmplx((/ 0.3156498673740053, 0.3431255827785649 /)),ref, '1.1')
    
    ref = dcmplx(0.00023134615630308335448329926098409)
    call test_one_MPL((/ 1,1 /),cmplx((/ 0.03, 0.5012562893380046 /)),ref, '1.2')
    
    ref = dcmplx(0.000023446106415452030937059124671151)
    call test_one_MPL((/ 2,1,2 /),cmplx((/ 0.03, 0.5012562893380046, 55.3832 /)),ref, '1.3')  
  end subroutine do_MPL_tests

  subroutine test_one_condensed(m,z,y,k,ref,test_id)
    integer :: m(:), k
    complex(kind=prec) :: z(:), y, res, ref
    character(len=*) :: test_id

    print*, '  ', 'testing GPL ', test_id, ' ...'
    res = G_condensed(m,z,y,k)
    call check(res,ref)
  end subroutine test_one_condensed

  subroutine test_one_flat(z,y,ref,test_id)
    complex(kind=prec) :: z(:), y, res, ref
    character(len=*) :: test_id

    print*, '  ', 'testing GPL ', test_id, ' ...'
    res = G_flat(z,y)
    call check(res,ref)
  end subroutine test_one_flat

  subroutine do_GPL_tests()
    complex(kind=prec) :: ref
    complex(kind=prec), parameter :: epsilon = 1E-14
    print*, 'doing GPL tests...'
    
    ref = dcmplx(0.0819393734128676)
    call test_one_condensed((/ 1,1 /),cmplx((/ 1.3d0, 1.1d0 /)),cmplx(0.4),2,ref,'2.1')
    
    ref = dcmplx(0.01592795952537145)
    call test_one_condensed((/ 3,2 /),cmplx((/ 1.3d0, 1.1d0 /)),cmplx(0.4),2,ref,'2.2')
    
    ref = dcmplx(0.0020332632172573974)
    call test_one_condensed((/ 4 /),cmplx((/ 0 /)),cmplx(1.6),1,ref,'2.3')

    ref = dcmplx(0.0020332632172573974)
    call test_one_flat(cmplx((/0.0,1.7,0.5/)),cmplx(1.1),ref,'2.5')
  end subroutine do_GPL_tests

  subroutine do_shuffle_tests() 
    complex(kind=prec) :: v(3) = cmplx((/1,2/))
    complex(kind=prec) :: w(2) = cmplx((/3,4/))

    call print_matrix(shuffle_product(v,w))
  end subroutine do_shuffle_tests

END PROGRAM TEST
 