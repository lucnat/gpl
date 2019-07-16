
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use globals
  use utils
  use shuffle
  use maths_functions
  use mpl_module
  use gpl_module
  use chenreftest, only: do_chen_test
  implicit none
  
  complex(kind=prec) :: res 
  real :: tol = 8.0e-7
  logical :: tests_successful = .true. 
#ifdef HAVE_GINAC
  character(len=6) :: ginacwhat
#endif

#ifdef DEBUG
  call parse_cmd_args()
#endif

  tol = 8e-10
  call do_MPL_tests() 
  call do_GPL_tests()
  ! call do_shuffle_tests() ! put this somewhere else
  tol = 8.0e-7
  tests_successful = tests_successful .and. do_chen_test(cmplx(0.3),cmplx(0.1))

#ifdef HAVE_GINAC
  tol = 2.0e-5
  ginacwhat = 'values'
  call do_muone_tests(cmplx(0.4),cmplx(.7),"")
  ginacwhat = 'speed1'
  call do_muone_tests(cmplx(0.4),cmplx(.7),"using GiNaC")
  ginacwhat = 'speed2'
  call do_muone_tests(cmplx(0.4),cmplx(.7),"using GPL")
#endif

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

    ref = cmplx(-0.06565799418838372)
    call test_one_MPL((/1, 1/), cmplx((/-0.25,-2. /)), ref, '1.4')
    ref = cmplx(-0.03199896396564833)
    call test_one_MPL((/2, 1/), cmplx((/-0.25,-2. /)), ref, '1.4')
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
    res = G(z)
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
    
    ref = cmplx((-0.0015317713178859967,-0.00045003911367000565))
    call test_one_flat(cmplx([0.,0.,(1-sqrt(1-z**2))/(xchen*z), 1./(xchen*z),1.]),ref,'3.7')


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
    call test_one_flat(cmplx([0.01, 199.99499987499374, 0.01, 1.]),ref,'4.11')



  end subroutine do_GPL_tests
  



#ifdef HAVE_GINAC
  subroutine test_one_ginac(z,test_id)
    complex(kind=prec) :: z(:), res, ref, geval
    character(len=*) :: test_id
    if (ginacwhat=="values") then
      print*, '  ', 'testing GPL ', test_id, ' ...'
      ref = geval(z,size(z))
      res = G(z)
      call check(res,ref)
    elseif (ginacwhat=="speed1") then
      ref = geval(z,size(z))
    elseif (ginacwhat=="speed2") then
      res = G(z)
    endif
  end subroutine
  subroutine do_muone_tests(x,y,msg)
    use maths_functions, only:clearcache
    complex(kind=prec) x, y
    real(kind=prec) tstart, tend
    character(len=*) :: msg
    call clearcache
    call cpu_time(tstart)
    call test_one_ginac([(-1.,0.),x],'6.1')
    call test_one_ginac([(-1.,0.),(-1.,0.),x],'6.2')
    call test_one_ginac([(0.,0.),(-1.,0.),x],'6.3')
    call test_one_ginac([(-1.,0.),(-1.,0.),(-1.,0.),x],'6.4')
    call test_one_ginac([(-1.,0.),(0.,0.),(-1.,0.),x],'6.5')
    call test_one_ginac([(0.,0.),(-1.,0.),(-1.,0.),x],'6.6')
    call test_one_ginac([(0.,0.),(0.,0.),(-1.,0.),x],'6.7')
    call test_one_ginac([(-1.,0.),(-1.,0.),(-1.,0.),(-1.,0.),x],'6.8')
    call test_one_ginac([(-1.,0.),(-1.,0.),(0.,0.),(-1.,0.),x],'6.9')
    call test_one_ginac([(-1.,0.),(0.,0.),(-1.,0.),(-1.,0.),x],'6.10')
    call test_one_ginac([(-1.,0.),(0.,0.),(0.,0.),(-1.,0.),x],'6.11')
    call test_one_ginac([(0.,0.),(-1.,0.),(-1.,0.),(-1.,0.),x],'6.12')
    call test_one_ginac([(0.,0.),(-1.,0.),(0.,0.),(-1.,0.),x],'6.13')
    call test_one_ginac([(0.,0.),(0.,0.),(-1.,0.),(-1.,0.),x],'6.14')
    call test_one_ginac([(0.,0.),(0.,0.),(0.,0.),(-1.,0.),x],'6.15')
    call test_one_ginac([(0.,0.),y],'6.16')
    call test_one_ginac([(1.,0.),y],'6.17')
    call test_one_ginac([(0.,0.),(0.,0.),y],'6.18')
    call test_one_ginac([(0.,0.),(1.,0.),y],'6.19')
    call test_one_ginac([(1.,0.),(0.,0.),y],'6.20')
    call test_one_ginac([(1.,0.),(1.,0.),y],'6.21')
    call test_one_ginac([(0.,0.),(0.,0.),(0.,0.),y],'6.22')
    call test_one_ginac([(0.,0.),(0.,0.),(1.,0.),y],'6.23')
    call test_one_ginac([(0.,0.),(1.,0.),(0.,0.),y],'6.24')
    call test_one_ginac([(0.,0.),(1.,0.),(1.,0.),y],'6.25')
    call test_one_ginac([(1.,0.),(0.,0.),(0.,0.),y],'6.26')
    call test_one_ginac([(1.,0.),(0.,0.),(1.,0.),y],'6.27')
    call test_one_ginac([(1.,0.),(1.,0.),(0.,0.),y],'6.28')
    call test_one_ginac([(1.,0.),(1.,0.),(1.,0.),y],'6.29')
    call test_one_ginac([(0.,0.),(0.,0.),(0.,0.),(0.,0.),y],'6.30')
    call test_one_ginac([(0.,0.),(0.,0.),(0.,0.),(1.,0.),y],'6.31')
    call test_one_ginac([(0.,0.),(0.,0.),(1.,0.),(0.,0.),y],'6.32')
    call test_one_ginac([(0.,0.),(0.,0.),(1.,0.),(1.,0.),y],'6.33')
    call test_one_ginac([(0.,0.),(1.,0.),(0.,0.),(0.,0.),y],'6.34')
    call test_one_ginac([(0.,0.),(1.,0.),(0.,0.),(1.,0.),y],'6.35')
    call test_one_ginac([(0.,0.),(1.,0.),(1.,0.),(0.,0.),y],'6.36')
    call test_one_ginac([(0.,0.),(1.,0.),(1.,0.),(1.,0.),y],'6.37')
    call test_one_ginac([(1.,0.),(0.,0.),(0.,0.),(0.,0.),y],'6.38')
    call test_one_ginac([(1.,0.),(0.,0.),(0.,0.),(1.,0.),y],'6.39')
    call test_one_ginac([(1.,0.),(0.,0.),(1.,0.),(0.,0.),y],'6.40')
    call test_one_ginac([(1.,0.),(0.,0.),(1.,0.),(1.,0.),y],'6.41')
    call test_one_ginac([(1.,0.),(1.,0.),(0.,0.),(0.,0.),y],'6.42')
    call test_one_ginac([(1.,0.),(1.,0.),(0.,0.),(1.,0.),y],'6.43')
    call test_one_ginac([(1.,0.),(1.,0.),(1.,0.),(0.,0.),y],'6.44')
    call test_one_ginac([(1.,0.),(1.,0.),(1.,0.),(1.,0.),y],'6.45')
    call test_one_ginac([(-1.,0.),y],'6.46')
    call test_one_ginac([(-1.,0.),(0.,0.),(0.,0.),y],'6.47')
    call test_one_ginac([(-1.,0.),(0.,0.),(1.,0.),y],'6.48')
    call test_one_ginac([(-1.,0.),(-1.,0.),y],'6.49')
    call test_one_ginac([(-1.,0.),(0.,0.),y],'6.50')
    call test_one_ginac([(-1.,0.),(1.,0.),y],'6.51')
    call test_one_ginac([(0.,0.),(-1.,0.),y],'6.52')
    call test_one_ginac([(1.,0.),(-1.,0.),y],'6.53')
    call test_one_ginac([(-1.,0.),(-1.,0.),(0.,0.),(0.,0.),y],'6.54')
    call test_one_ginac([(-1.,0.),(-1.,0.),(0.,0.),(1.,0.),y],'6.55')
    call test_one_ginac([(-1.,0.),(0.,0.),(0.,0.),(0.,0.),y],'6.56')
    call test_one_ginac([(-1.,0.),(0.,0.),(0.,0.),(1.,0.),y],'6.57')
    call test_one_ginac([(-1.,0.),(0.,0.),(1.,0.),(0.,0.),y],'6.58')
    call test_one_ginac([(-1.,0.),(0.,0.),(1.,0.),(1.,0.),y],'6.59')
    call test_one_ginac([(-1.,0.),(1.,0.),(0.,0.),(0.,0.),y],'6.60')
    call test_one_ginac([(-1.,0.),(1.,0.),(0.,0.),(1.,0.),y],'6.61')
    call test_one_ginac([(0.,0.),(-1.,0.),(0.,0.),(0.,0.),y],'6.62')
    call test_one_ginac([(0.,0.),(-1.,0.),(0.,0.),(1.,0.),y],'6.63')
    call test_one_ginac([(1.,0.),(-1.,0.),(0.,0.),(0.,0.),y],'6.64')
    call test_one_ginac([(1.,0.),(-1.,0.),(0.,0.),(1.,0.),y],'6.65')
    call test_one_ginac([-(1/y),x],'6.66')
    call test_one_ginac([-y,x],'6.67')
    call test_one_ginac([-(1/y),(-1.,0.),x],'6.68')
    call test_one_ginac([-y,(-1.,0.),x],'6.69')
    call test_one_ginac([-(1/y),(-1.,0.),(-1.,0.),x],'6.70')
    call test_one_ginac([-(1/y),(0.,0.),(-1.,0.),x],'6.71')
    call test_one_ginac([-y,(-1.,0.),(-1.,0.),x],'6.72')
    call test_one_ginac([-y,(0.,0.),(-1.,0.),x],'6.73')
    call test_one_ginac([1 - 1/y - y,-(1/y),x],'6.74')
    call test_one_ginac([-(1/y),-(1/y),x],'6.75')
    call test_one_ginac([-y,-(1/y),x],'6.76')
    call test_one_ginac([1 - 1/y - y,x],'6.77')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(0.,0.),y],'6.78')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(1.,0.),y],'6.79')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),y],'6.80')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),y],'6.81')
    call test_one_ginac([1 - 1/y - y,(-1.,0.),x],'6.82')
    call test_one_ginac([1 - 1/y - y,-y,x],'6.83')
    call test_one_ginac([-(1/y),-y,x],'6.84')
    call test_one_ginac([-y,-y,x],'6.85')
    call test_one_ginac([1 - 1/y - y,-(1/y),(-1.,0.),x],'6.86')
    call test_one_ginac([1 - 1/y - y,-y,(-1.,0.),x],'6.87')
    call test_one_ginac([-(1/y),-(1/y),(-1.,0.),x],'6.88')
    call test_one_ginac([-(1/y),-y,(-1.,0.),x],'6.89')
    call test_one_ginac([-y,-(1/y),(-1.,0.),x],'6.90')
    call test_one_ginac([-y,-y,(-1.,0.),x],'6.91')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(0.,0.),(0.,0.),(0.,0.),y],'6.92')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(0.,0.),(0.,0.),(1.,0.),y],'6.93')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(0.,0.),(1.,0.),(0.,0.),y],'6.94')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(0.,0.),(1.,0.),(1.,0.),y],'6.95')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(1.,0.),(0.,0.),(0.,0.),y],'6.96')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(1.,0.),(0.,0.),(1.,0.),y],'6.97')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(1.,0.),(1.,0.),(0.,0.),y],'6.98')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),(1.,0.),(1.,0.),(1.,0.),y],'6.99')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(0.,0.),(0.,0.),y],'6.100')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(0.,0.),(1.,0.),y],'6.101')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(1.,0.),(0.,0.),y],'6.102')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(1.,0.),(1.,0.),y],'6.103')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),(0.,0.),(0.,0.),y],'6.104')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),(0.,0.),(1.,0.),y],'6.105')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),(1.,0.),(0.,0.),y],'6.106')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),(1.,0.),(1.,0.),y],'6.107')
    call test_one_ginac([1 - 1/y - y,(-1.,0.),(0.,0.),(-1.,0.),x],'6.108')
    call test_one_ginac([1 - 1/y - y,-(1/y),(-1.,0.),(-1.,0.),x],'6.109')
    call test_one_ginac([1 - 1/y - y,-(1/y),(0.,0.),(-1.,0.),x],'6.110')
    call test_one_ginac([1 - 1/y - y,-y,(-1.,0.),(-1.,0.),x],'6.111')
    call test_one_ginac([1 - 1/y - y,-y,(0.,0.),(-1.,0.),x],'6.112')
    call test_one_ginac([-(1/y),(-1.,0.),(-1.,0.),(-1.,0.),x],'6.113')
    call test_one_ginac([-(1/y),(-1.,0.),(0.,0.),(-1.,0.),x],'6.114')
    call test_one_ginac([-(1/y),(0.,0.),(-1.,0.),(-1.,0.),x],'6.115')
    call test_one_ginac([-(1/y),(0.,0.),(0.,0.),(-1.,0.),x],'6.116')
    call test_one_ginac([-(1/y),-(1/y),(-1.,0.),(-1.,0.),x],'6.117')
    call test_one_ginac([-(1/y),-(1/y),(0.,0.),(-1.,0.),x],'6.118')
    call test_one_ginac([-(1/y),-y,(-1.,0.),(-1.,0.),x],'6.119')
    call test_one_ginac([-(1/y),-y,(0.,0.),(-1.,0.),x],'6.120')
    call test_one_ginac([-y,(-1.,0.),(-1.,0.),(-1.,0.),x],'6.121')
    call test_one_ginac([-y,(-1.,0.),(0.,0.),(-1.,0.),x],'6.122')
    call test_one_ginac([-y,(0.,0.),(-1.,0.),(-1.,0.),x],'6.123')
    call test_one_ginac([-y,(0.,0.),(0.,0.),(-1.,0.),x],'6.124')
    call test_one_ginac([-y,-(1/y),(-1.,0.),(-1.,0.),x],'6.125')
    call test_one_ginac([-y,-(1/y),(0.,0.),(-1.,0.),x],'6.126')
    call test_one_ginac([-y,-y,(-1.,0.),(-1.,0.),x],'6.127')
    call test_one_ginac([-y,-y,(0.,0.),(-1.,0.),x],'6.128')
    call test_one_ginac([cmplx(0.5,-sqrt(3.)/2.),y],'6.129')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),y],'6.130')
    call test_one_ginac([(-1.,0.),-y,x],'6.131')
    call test_one_ginac([(-1.,0.),-(1/y),x],'6.132')
    call test_one_ginac([(-1.,0.),-(1/y),(-1.,0.),x],'6.133')
    call test_one_ginac([(-1.,0.),-y,(-1.,0.),x],'6.134')
    call test_one_ginac([(-1.,0.),-(1/y),(-1.,0.),(-1.,0.),x],'6.135')
    call test_one_ginac([(-1.,0.),-(1/y),(0.,0.),(-1.,0.),x],'6.136')
    call test_one_ginac([(-1.,0.),-y,(-1.,0.),(-1.,0.),x],'6.137')
    call test_one_ginac([(-1.,0.),-y,(0.,0.),(-1.,0.),x],'6.138')
    call test_one_ginac([(0.,0.),-(1/y),x],'6.139')
    call test_one_ginac([(0.,0.),-y,x],'6.140')
    call test_one_ginac([(0.,0.),-(1/y),(-1.,0.),x],'6.141')
    call test_one_ginac([(0.,0.),-y,(-1.,0.),x],'6.142')
    call test_one_ginac([(0.,0.),-(1/y),(-1.,0.),(-1.,0.),x],'6.143')
    call test_one_ginac([(0.,0.),-(1/y),(0.,0.),(-1.,0.),x],'6.144')
    call test_one_ginac([(0.,0.),-y,(-1.,0.),(-1.,0.),x],'6.145')
    call test_one_ginac([(0.,0.),-y,(0.,0.),(-1.,0.),x],'6.146')
    call test_one_ginac([(-1.,0.),(-1.,0.),(0.,0.),y],'6.147')
    call test_one_ginac([(0.,0.),(-1.,0.),(0.,0.),y],'6.148')
    call test_one_ginac([(-1.,0.),(-1.,0.),(-1.,0.),(0.,0.),y],'6.149')
    call test_one_ginac([(-1.,0.),(0.,0.),(-1.,0.),(0.,0.),y],'6.150')
    call test_one_ginac([(0.,0.),(-1.,0.),(-1.,0.),(0.,0.),y],'6.151')
    call test_one_ginac([(0.,0.),(0.,0.),(-1.,0.),(0.,0.),y],'6.152')
    call test_one_ginac([(0.,0.),(-1.,0.),(1.,0.),(0.,0.),y],'6.153')
    call test_one_ginac([(0.,0.),(1.,0.),(-1.,0.),(0.,0.),y],'6.154')
    call test_one_ginac([(1.,0.),(0.,0.),(-1.,0.),(0.,0.),y],'6.155')
    call test_one_ginac([(-1.,0.),(1.,0.),(0.,0.),y],'6.156')
    call test_one_ginac([(1.,0.),(-1.,0.),(0.,0.),y],'6.157')
    call test_one_ginac([(-1.,0.),(-1.,0.),(1.,0.),(0.,0.),y],'6.158')
    call test_one_ginac([(-1.,0.),(1.,0.),(-1.,0.),(0.,0.),y],'6.159')
    call test_one_ginac([(-1.,0.),(1.,0.),(1.,0.),(0.,0.),y],'6.160')
    call test_one_ginac([(1.,0.),(-1.,0.),(-1.,0.),(0.,0.),y],'6.161')
    call test_one_ginac([(1.,0.),(-1.,0.),(1.,0.),(0.,0.),y],'6.162')
    call test_one_ginac([(1.,0.),(1.,0.),(-1.,0.),(0.,0.),y],'6.163')
    call test_one_ginac([(0.,0.),(1.,0.),x],'6.164')
    call test_one_ginac([(0.,0.),(1.,0.),(0.,0.),(-1.,0.),x],'6.165')
    call test_one_ginac([(1.,0.),(0.,0.),(-1.,0.),x],'6.166')
    call test_one_ginac([(1.,0.),x],'6.167')
    call test_one_ginac([(-1.,0.),(1.,0.),(0.,0.),(-1.,0.),x],'6.168')
    call test_one_ginac([(1.,0.),(-1.,0.),(0.,0.),(-1.,0.),x],'6.169')
    call test_one_ginac([(1.,0.),(0.,0.),(-1.,0.),(-1.,0.),x],'6.170')
    call test_one_ginac([(1.,0.),(0.,0.),(0.,0.),(-1.,0.),x],'6.171')
    call test_one_ginac([(1.,0.),(1.,0.),(0.,0.),(-1.,0.),x],'6.172')
    call test_one_ginac([(-1.,0.),(1.,0.),x],'6.173')
    call test_one_ginac([(1.,0.),(-1.,0.),x],'6.174')
    call test_one_ginac([(1.,0.),(1.,0.),x],'6.175')
    call test_one_ginac([-(1/y),(1.,0.),x],'6.176')
    call test_one_ginac([-y,(1.,0.),x],'6.177')
    call test_one_ginac([-(1/y),(1.,0.),(0.,0.),(-1.,0.),x],'6.178')
    call test_one_ginac([-y,(1.,0.),(0.,0.),(-1.,0.),x],'6.179')
    call test_one_ginac([(-1 + y - y**2)/y,x],'6.180')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),y],'6.181')
    call test_one_ginac([-cmplx(-0.5,sqrt(3.)/2.),(0.,0.),y],'6.182')
    call test_one_ginac([(-1 + y - y**2)/y,(-1.,0.),x],'6.183')
    call test_one_ginac([(-1 + y - y**2)/y,-(1/y),x],'6.184')
    call test_one_ginac([(-1 + y - y**2)/y,-y,x],'6.185')
    call test_one_ginac([(-1 + y - y**2)/y,-(1/y),(-1.,0.),x],'6.186')
    call test_one_ginac([(-1 + y - y**2)/y,-y,(-1.,0.),x],'6.187')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(0.,0.),(0.,0.),y],'6.188')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(0.,0.),(1.,0.),(0.,0.),y],'6.189')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),(1.,0.),(0.,0.),(0.,0.),y],'6.190')
    call test_one_ginac([-cmplx(-0.5,sqrt(3.)/2.),(0.,0.),(0.,0.),(0.,0.),y],'6.191')
    call test_one_ginac([-cmplx(-0.5,sqrt(3.)/2.),(0.,0.),(1.,0.),(0.,0.),y],'6.192')
    call test_one_ginac([-cmplx(-0.5,sqrt(3.)/2.),(1.,0.),(0.,0.),(0.,0.),y],'6.193')
    call test_one_ginac([(-1 + y - y**2)/y,(-1.,0.),(0.,0.),(-1.,0.),x],'6.194')
    call test_one_ginac([(-1 + y - y**2)/y,-(1/y),(0.,0.),(-1.,0.),x],'6.195')
    call test_one_ginac([(-1 + y - y**2)/y,-y,(0.,0.),(-1.,0.),x],'6.196')
    call test_one_ginac([cmplx(0.5,sqrt(3.)/2.),y],'6.197')
    call test_one_ginac([-cmplx(-0.5,sqrt(3.)/2.),y],'6.198')
    call cpu_time(tend)
    write(*,900) msg,198./(tend-tstart)

900 format("Evaluating ",A," at ",F9.2,"G/s")
  end subroutine
#endif
  ! subroutine do_shuffle_tests() 
  !   complex(kind=prec) :: v(2) = cmplx((/1,2/))
  !   complex(kind=prec) :: w(2) = cmplx((/3,4/))

  !   call print_matrix(shuffle_product(v,w))
  ! end subroutine do_shuffle_tests

END PROGRAM TEST

! In terminal kann man den exit code bekommen via echo $? 
