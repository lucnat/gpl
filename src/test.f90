
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use globals
  use utils
  use shuffle
  use maths_functions
  use mpl_module
  use gpl_module
  implicit none
  
  complex(kind=prec) :: res 
  real, parameter :: tol = 1.0e-12
  logical :: tests_successful = .true. 

  call parse_cmd_args()
  call do_MPL_tests() 
  call do_GPL_tests()
  ! call do_shuffle_tests() ! put this somewhere else


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

    delta = abs(res-ref)
    if(delta < tol) then
      print*, '       ',' passed with delta = ', delta
    else 
      print*, '  ',' FAILED with delta = ', delta
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
    
    ref = cmplx(0.022696600480693277651633)
    call test_one_MPL((/ 1,1 /),cmplx((/ 0.3156498673740053, 0.3431255827785649 /)),ref, '1.1')
    
    ref = cmplx(0.00023134615630308335448329926098409)
    call test_one_MPL((/ 1,1 /),cmplx((/ 0.03, 0.5012562893380046 /)),ref, '1.2')
    
    ref = cmplx(0.000023446106415452030937059124671151)
    call test_one_MPL((/ 2,1,2 /),cmplx((/ 0.03, 0.5012562893380046, 55.3832 /)),ref, '1.3')  
  end subroutine do_MPL_tests

  subroutine test_one_condensed(m,z,y,k,ref,test_id)
    integer :: m(:), k
    complex(kind=prec) :: z(:), y, res, ref
    character(len=*) :: test_id

    print*, '  ', 'testing GPL ', test_id, ' ...'
    res = G_condensed(m,toinum(z),inum(y,di0),k)
    call check(res,ref)
  end subroutine test_one_condensed

  subroutine test_one_flat(z,ref,test_id)
    complex(kind=prec) :: z(:), res, ref
    character(len=*) :: test_id

    print*, '  ', 'testing GPL ', test_id, ' ...'
    res = GPL(z)
    call check(res,ref)
  end subroutine test_one_flat

  subroutine do_GPL_tests()
    complex(kind=prec) :: ref
    real(kind=prec) :: z, xchen
    print*, 'doing GPL tests...'
    
    ref = cmplx(0.0819393734128676)
    call test_one_condensed((/ 1,1 /),cmplx((/ 1.3, 1.1 /)),cmplx(0.4),2,ref,'2.1')
    
    ref = cmplx(0.01592795952537145)
    call test_one_condensed((/ 3,2 /),cmplx((/ 1.3, 1.1 /)),cmplx(0.4),2,ref,'2.2')
    
    ref = cmplx(0.0020332632172573974)
    call test_one_condensed((/ 4 /),cmplx((/ 0 /)),cmplx(1.6),1,ref,'2.3')

    ! requires making convergent
    ref = (0.09593041677639341, -0.8829351795197851)
    call test_one_flat(cmplx([0,1,3,2]),ref,'2.4')

    ref = (0.009947947789928621,0.0)
    call test_one_flat(cmplx([0.0, 0.0, -3.3333333333333335, -3.3333333333333335, 1.0]),ref,'2.5')

    ! requires hoelder convolution
    ref = (-0.012709942828250949,0.0)
    call test_one_flat(cmplx([0.0, 3.3333333333333335, 1.0, 3.3333333333333335, 1.0]),ref,'2.6')

    ! here the tests from the mathematica nb start
    ! --------------------------------------------------------------------------

    z = 1./200; xchen = 0.3;

    ref = (-0.0050125418235441935,0.0)
    call test_one_flat(cmplx([1./z,1.0]),ref,'3.1')

    ref = (-0.0015011261262671913,0.0)
    call test_one_flat(cmplx([1./(xchen*z),1.0]),ref,'3.2')

    ref = (-0.0007502860817810596,0.0)
    call test_one_flat(cmplx([(1+sqrt(1-z**2))/(xchen*z),1.0]),ref,'3.3')

    ref = (0.0074335969894765335,0.0)
    call test_one_flat(cmplx([-1./xchen,-1./xchen,1.,1.,1.0]),ref,'3.4')
    
    ref = (-8.403785974849544e-6,0.0)
    call test_one_flat(cmplx([-1./xchen,0.,-1./xchen,1./(xchen*z),1.0]),ref,'3.5')
    
    ref = cmplx((0.4925755847450199,2.6389214054743295))
    call test_one_flat(cmplx([-1.,-1.,z,z,1.]),ref,'3.6')
    
    ! ref = cmplx((-0.0015317713178859967,-0.00045003911367000565))
    ! call test_one_flat(cmplx([0.,0.,(1-sqrt(1-z**2))/(xchen*z), 1./(xchen*z),1.]),ref,'3.7')


    ! here the chen integral tests start
    ! ----------------------------------------------------------------------

    ref = (-1.2039728043259361,0.0)
    call test_one_flat(cmplx([0., 0.3]),ref,'4.1')

    ref = (10.603796220956799,0.0)
    call test_one_flat(cmplx([0., 0., 0.01]),ref,'4.2')

    ref = (0.0005042990065180765,0.0)
    call test_one_flat(cmplx([100., 1., 0.3]),ref,'4.3')

    ref = (0.05630861877185141,0.0)
    call test_one_flat(cmplx([1., 0., 0.01]),ref,'4.4')

    ref = (0.04032895150872735,0.2922709647568842)
    call test_one_flat([(0.01,0.9999499987499375), (0.3,0)],ref,'4.5')

    ref = cmplx(0.000025210645098340354)
    call test_one_flat(cmplx([100., 199.99499987499374, 1.]),ref,'4.6')

    ref = (0.0556470547632135,0.)
    call test_one_flat(cmplx([-1., 0.01, 0., 0.01]),ref,'4.7')

    ref = (0.00003794895846844715,0.)
    call test_one_flat(cmplx([100., 100., 1., 0., 1.]),ref,'4.8')

    ref = (0.00007574284433216497,0.)
    call test_one_flat(cmplx([100., 1., 100., 0., 1.]),ref,'4.9')

    ref = (1.8801911586012443,2.509434138598062)
    call test_one_flat(cmplx([0.01, -1.0, 0.01, 1.]),ref,'4.10')

    ref = (-0.012539108315054982, -0.015414250168437678)
    call test_one_flat(cmplx([0.01, 199.99499987499374, 0.01, 1.]),ref,'4.10')



  end subroutine do_GPL_tests
  
  ! subroutine do_shuffle_tests() 
  !   complex(kind=prec) :: v(2) = cmplx((/1,2/))
  !   complex(kind=prec) :: w(2) = cmplx((/3,4/))

  !   call print_matrix(shuffle_product(v,w))
  ! end subroutine do_shuffle_tests

END PROGRAM TEST

! In terminal kann man den exit code bekommen via echo $? 
