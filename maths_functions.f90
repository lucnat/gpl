
MODULE maths_functions
  use globals
  use utils
  implicit none

  ! integer, parameter :: prec = selected_real_kind(15,32)  
  ! integer, parameter :: GPLInfinity = 30              ! the default outermost expansion order for MPLs
  ! real(kind=prec), parameter :: epsilon = 1e-15       ! used for the small imaginary part
  ! real(kind=prec), parameter :: zero = 1e-15          ! values smaller than this count as zero
  ! real(kind=prec), parameter :: pi = 3.14159265358979323846_prec

CONTAINS 

  FUNCTION naive_polylog(m,x,n_passed) result(res)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    integer, optional :: n_passed
    complex(kind=prec) :: x, res
    integer :: i,n
    integer, allocatable :: j(:)
    n = merge(n_passed,GPLInfinity,present(n_passed))  
    j = (/(i, i=1,n,1)/) 
    res = sum(x**j / j**m)
  END FUNCTION naive_polylog

  FUNCTION Li2(x)

   !! Dilogarithm for arguments x < = 1.0

   real (kind=prec):: X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
   real (kind=prec):: C(0:18),H,ALFA,B0,B1,B2,LI2_OLD
   real (kind=prec):: Li2
   integer :: i

   DATA ZERO /0.0_prec/, ONE /1.0_prec/
   DATA HALF /0.5_prec/, MALF /-0.5_prec/ 
   DATA MONE /-1.0_prec/, MTWO /-2.0_prec/
   DATA PI3 /3.289868133696453_prec/, PI6 /1.644934066848226_prec/

   DATA C( 0) / 0.4299669356081370_prec/
   DATA C( 1) / 0.4097598753307711_prec/
   DATA C( 2) /-0.0185884366501460_prec/
   DATA C( 3) / 0.0014575108406227_prec/
   DATA C( 4) /-0.0001430418444234_prec/
   DATA C( 5) / 0.0000158841554188_prec/
   DATA C( 6) /-0.0000019078495939_prec/
   DATA C( 7) / 0.0000002419518085_prec/
   DATA C( 8) /-0.0000000319334127_prec/
   DATA C( 9) / 0.0000000043454506_prec/
   DATA C(10) /-0.0000000006057848_prec/
   DATA C(11) / 0.0000000000861210_prec/
   DATA C(12) /-0.0000000000124433_prec/
   DATA C(13) / 0.0000000000018226_prec/
   DATA C(14) /-0.0000000000002701_prec/
   DATA C(15) / 0.0000000000000404_prec/
   DATA C(16) /-0.0000000000000061_prec/
   DATA C(17) / 0.0000000000000009_prec/
   DATA C(18) /-0.0000000000000001_prec/

   if(X > 1.00000000001_prec) then
     print*, 'crashes because Li called with bad arguments'
   elseif(X > 1.0_prec) then
     X = 1._prec
   endif    

   IF(X > 0.999999_prec) THEN
    LI2_OLD=PI6
    Li2 = Real(LI2_OLD,prec)
    RETURN
   ELSE IF(X .EQ. MONE) THEN
    LI2_OLD=MALF*PI6
    RETURN
   END IF
   T=-X
   IF(T .LE. MTWO) THEN
    Y=MONE/(ONE+T)
    S=ONE
    A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
   ELSE IF(T .LT. MONE) THEN
    Y=MONE-T
    S=MONE
    A=LOG(-T)
    A=-PI6+A*(A+LOG(ONE+ONE/T))
   ELSE IF(T .LE. MALF) THEN
    Y=(MONE-T)/T
    S=ONE
    A=LOG(-T)
    A=-PI6+A*(MALF*A+LOG(ONE+T))
   ELSE IF(T .LT. ZERO) THEN
    Y=-T/(ONE+T)
    S=MONE
    A=HALF*LOG(ONE+T)**2
   ELSE IF(T .LE. ONE) THEN
    Y=T
    S=ONE
    A=ZERO
   ELSE
    Y=ONE/T
    S=MONE
    A=PI6+HALF*LOG(T)**2
   END IF

   H=Y+Y-ONE
   ALFA=H+H
   B1=ZERO
   B2=ZERO
   DO  I = 18,0,-1
     B0=C(I)+ALFA*B1-B2
     B2=B1
     B1=B0
   ENDDO
   LI2_OLD=-(S*(B0-H*B2)+A)
         ! Artificial conversion           
   Li2 = Real(LI2_OLD,prec)

  END FUNCTION Li2

  FUNCTION dilog_in_unit_circle(x) result(res)
    ! evaluates for any argument x in unit circle
    complex(kind=prec) :: x, res
    res = naive_polylog(2,x)
  END FUNCTION dilog_in_unit_circle

  FUNCTION dilog(x) result(res)
    ! evaluates dilog for any argument
    complex(kind=prec) :: res
    complex(kind=prec) :: x
    if(abs(x) <= 1.0) then
     res = naive_polylog(2,x)
    else
     res = -dilog_in_unit_circle(1/x) - pi**2/6 - log(dcmplx(   add_ieps(-x)  ))**2 / 2
   end if
  END FUNCTION dilog

  FUNCTION polylog(m,x) result(res)
    ! computes the polylog for now naively (except for dilog half-naively)
    integer :: m
    complex(kind=prec) :: x,res
    print*, 'called polylog(',m,',',x,')'
    if(m == 2) then
      res = dilog(x)
    else 
      res = naive_polylog(m,x)
    end if
  END FUNCTION polylog

END MODULE maths_functions

! PROGRAM test
!   use maths_functions
!   implicit none
!   complex(kind=prec) :: res
!   res = dilog(dcmplx(0.4d0))
!   print*, res
! END PROGRAM test

