
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
#ifdef HAVE_MM
  call do_ginac_tests
  call do_timing_tests(5)
#endif
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

  function evalt(arr, what)
    implicit none
    complex(kind=prec) :: arr(:), evalt, geval
    integer what, i, l

    evalt =0.

    do i=1,size(arr)
      if (arr(i) .ne. arr(i)) then ! isnan?
        l = i-1
        goto 123
      endif
    enddo
    l = size(arr)

123 continue
    if (l==0) return

    if (what .eq. 0) then
      evalt = G(arr(1:l))
    elseif (what.eq.1) then
      evalt = geval(arr(1:l),l)
    endif
  end function

  subroutine perform_ginacv(n, args)
    use maths_functions, only:clearcache
    complex(kind=prec) :: args(:,:)
    integer i,n
    call clearcache

    do i=1,size(args,1)
      write(*,900) n,i
      call check(evalt(args(i,:),0),evalt(args(i,:),1))
    enddo

900 format("   testing GPL ",I1,".",I4.4," ...")

  end subroutine

  subroutine do_ginac_tests
    use maths_functions, only:clearcache
    use gtestchen   , only: inichen   =>args
    use gtestmuone  , only: inimuone  =>args
    use gtestmuonenp, only: inimuonenp=>args
    implicit none
    tol = 6.0e-6
    call perform_ginacv( 6, inichen   (cmplx(0.3), cmplx(0.1)) )
    call perform_ginacv( 7, inimuone  (cmplx(0.5), cmplx(0.6)) )
    call perform_ginacv( 8, inimuonenp(cmplx(0.3), cmplx(0.6)) )

  end subroutine

  subroutine do_one_speed_test(args, u, msg)
    use maths_functions, only:clearcache
    implicit none
    complex(kind=prec) :: args(:,:,:),res
    real(kind=prec) :: tstart, tend, time(2), ttime(2)
    integer i,j, u
    character,parameter :: cr = achar(13)
    character(len=*) msg
    do j=1,size(args,1)
      ! try function a bunch of times
      call cpu_time(tstart)
      do i=1,size(args,3)
        res=evalt(args(j,:,i),0)
      enddo
      call cpu_time(tend)
      time(1) = (tend-tstart)/size(args,3)
      if (time(1).lt.zero) print*,j
      ttime(1) = ttime(1) + time(1)

      call cpu_time(tstart)
      do i=1,size(args,3)
        res=evalt(args(j,:,i),1)
      enddo
      call cpu_time(tend)
      time(2) = (tend-tstart)/size(args,3)
      if (time(2).lt.zero) print*,j
      ttime(2) = ttime(2) + time(2)

      write(u,*) time

      write( * , 900, advance='no' ) cr, j, size(args,1), msg
    enddo
    print*,
    write(*,901) msg, size(args,1)/ttime(2)/1000., size(args,1)/ttime(1)/1000., int(ttime(2)/ttime(1))

900 FORMAT(a , 'Function ',i4,'/',i4,' for ',a)
901 format('Evaluating ',A,' using GiNaC at ',F9.2,'kG/s and GPL at ',F9.2,'kG/s (',I3,'x)')
  end subroutine

  subroutine do_timing_tests(n)
    use gtestchen   , only: inichen   =>args
    use gtestmuone  , only: inimuone  =>args
    use gtestmuonenp, only: inimuonenp=>args
    implicit none
    integer, intent(in) :: n
    integer i
    complex(kind=prec) :: cargs( 1399,5,n)
    complex(kind=prec) :: pargs(  198,5,n)
    complex(kind=prec) :: nargs( 1733,5,n)
    real(kind=prec) :: z, x, y, w
    integer ranseed
    ranseed = 233123

    do i=1,n
      z = ran2(ranseed) / 2.
      x = ran2(ranseed)*(1-z) + z
      cargs(:,:,i) = inichen(cmplx(x), cmplx(z))

      w = ran2(ranseed) ! 0<w<1
      z = ran2(ranseed) * (sqrt(1-w+w**2)-sqrt(w)) + sqrt(w)
      nargs(:,:,i) = inimuonenp(cmplx(w), cmplx(z))

      x = ran2(ranseed)
      y = ran2(ranseed)
      pargs(:,:,i) = inimuone(cmplx(x), cmplx(y))
    enddo

    cargs(1181,:,:)=cargs(1181,:,:)/0


    open(unit=9, file="stats.txt")
    write(9,*) "Chen"
    call do_one_speed_test(cargs,9,"Chen")
    write(9,*) "MUonE-planar"
    call do_one_speed_test(pargs,9,"Muone planar")
    write(9,*) "MUonE-non-planar"
    call do_one_speed_test(nargs,9,"Muone non planar")
    close(unit=9)

  end subroutine

      FUNCTION RAN2(randy)

            ! This is the usual "random"

      implicit none
      real(kind=prec) :: MINV,RAN2
      integer m,a,Qran,r,hi,lo,randy
      PARAMETER(M=2147483647,A=16807,Qran=127773,R=2836)
      PARAMETER(MINV=0.46566128752458e-09)
      HI = RANDY/Qran
      LO = MOD(RANDY,Qran)
      RANDY = A*LO - R*HI
      IF(RANDY.LE.0) RANDY = RANDY + M
      RAN2 = RANDY*MINV
      END FUNCTION RAN2




#endif
  ! subroutine do_shuffle_tests() 
  !   complex(kind=prec) :: v(2) = cmplx((/1,2/))
  !   complex(kind=prec) :: w(2) = cmplx((/3,4/))

  !   call print_matrix(shuffle_product(v,w))
  ! end subroutine do_shuffle_tests

END PROGRAM TEST

! In terminal kann man den exit code bekommen via echo $? 
