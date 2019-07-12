
MODULE maths_functions
  use globals
  use utils
  implicit none
  interface polylog
    module procedure polylog1, polylog2
  end interface polylog

  real(kind=prec), parameter :: zeta(2:10) = (/1.6449340668482262, 1.2020569031595942, 1.0823232337111381, &
               1.03692775514337, 1.0173430619844488, 1.008349277381923, &
               1.0040773561979441, 1.0020083928260821, 1.000994575127818/)

  type el
    type(inum) :: c
    complex(kind=prec) ans
  end type el

  type(el) :: cache(PolyLogCacheSize(1),PolyLogCacheSize(2))
  integer :: plcachesize(PolyLogCacheSize(1)) = 0
CONTAINS

  FUNCTION naive_polylog(m,x) result(res)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    complex(kind=prec) :: x, res
    integer :: i
    res=0.
    do i=1,PolylogInfinity
      if(i**m.lt.0) return ! roll over
      if(abs(x**i).lt.1.e-250) return
      res = res+x**i/i**m
    enddo
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
   ELSE IF(abs(x-MONE) < zero) THEN
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

  RECURSIVE FUNCTION dilog(x) result(res)
    ! evaluates dilog for any argument |x|<1
    complex(kind=prec) :: res
    complex(kind=prec) :: x

    if(abs(aimag(x)) < zero ) then
      res = Li2(real(x))
    else
      res = naive_polylog(2,x)
    endif
  END FUNCTION dilog

  FUNCTION Li3(x)
    ! Trilogarithm for arguments x < = 1.0
    ! This was hacked from LI2 to also follow C332
    ! In theory this could also produce Re[Li [x]] for x>1

    real (kind=prec):: X,S,A
    real (kind=prec):: CA(0:18),HA,ALFAA,BA0,BA1,BA2, YA
    real (kind=prec):: CB(0:18),HB,ALFAB,BB0,BB1,BB2, YB
    DATA CA(0) / 0.4617293928601208/
    DATA CA(1) / 0.4501739958855029/
    DATA CA(2) / -0.010912841952292843/
    DATA CA(3) / 0.0005932454712725702/
    DATA CA(4) / -0.00004479593219266303/
    DATA CA(5) / 4.051545785869334e-6/
    DATA CA(6) / -4.1095398602619446e-7/
    DATA CA(7) / 4.513178777974119e-8/
    DATA CA(8) / -5.254661564861129e-9/
    DATA CA(9) / 6.398255691618666e-10/
    DATA CA(10) / -8.071938105510391e-11/
    DATA CA(11) / 1.0480864927082917e-11/
    DATA CA(12) / -1.3936328400075057e-12/
    DATA CA(13) / 1.8919788723690422e-13/
    DATA CA(14) / -2.6097139622039465e-14/
    DATA CA(15) / 3.774985548158685e-15/
    DATA CA(16) / -5.671361978114946e-16/
    DATA CA(17) / 1.1023848202712794e-16/
    DATA CA(18) / -5.0940525990875006e-17/
    DATA CB(0) / -0.016016180449195803/
    DATA CB(1) / -0.5036424400753012/
    DATA CB(2) / -0.016150992430500253/
    DATA CB(3) / -0.0012440242104245127/
    DATA CB(4) / -0.00013757218124463538/
    DATA CB(5) / -0.000018563818526041144/
    DATA CB(6) / -2.841735345177361e-6/
    DATA CB(7) / -4.7459967908588557e-7/
    DATA CB(8) / -8.448038544563037e-8/
    DATA CB(9) / -1.5787671270014e-8/
    DATA CB(10) / -3.0657620579122164e-9/
    DATA CB(11) / -6.140791949281482e-10/
    DATA CB(12) / -1.2618831590198e-10/
    DATA CB(13) / -2.64931268635803e-11/
    DATA CB(14) / -5.664711482422879e-12/
    DATA CB(15) / -1.2303909436235178e-12/
    DATA CB(16) / -2.7089360852246495e-13/
    DATA CB(17) / -6.024075373994343e-14/
    DATA CB(18) / -1.2894320641440237e-14/
    real (kind=prec):: Li3
    real (kind=prec), parameter :: zeta2 = 1.6449340668482264365
    real (kind=prec), parameter :: zeta3 = 1.2020569031595942854
    integer :: i


    if(x > 1.00000000001_prec) then
      print*, 'need to crash Li3, since not convergent'
    elseif(x > 1.0_prec) then
      x = 1._prec
    endif

    IF(X > 0.999999_prec) THEN
      LI3=zeta3
    RETURN
    ELSE IF( abs(x+1) < zero) THEN
      LI3=-0.75_prec*zeta3
    RETURN
    END IF
    IF(X .LE. -1._prec) THEN
      YA=1._prec/x ; YB=0._prec
      S=-1._prec
      A=-LOG(-X)*(zeta2+LOG(-x)**2/6._prec)
    ELSE IF(X .LE. 0._prec) THEN
      YA=x ; YB=0._prec
      S=-1._prec
      A=0._prec
    ELSE IF(X .LE. 0.5_prec) THEN
      YA=0._prec ; YB=x
      S=-1._prec
      A=0._prec
    ELSE IF(X .LE. 1._prec) THEN
      YA=(x-1._prec)/x ; YB=1._prec-x
      S=1._prec
      A=zeta3 + zeta2*Log(x) - (Log(1._prec - X)*Log(X)**2)/2._prec + Log(X)**3/6._prec
    ELSE IF(X .LE. 2._prec) THEN
      YA=1._prec - X ; YB=(X-1._prec)/X
      S=1._prec
      A=zeta3 + zeta2*Log(x) - (Log(X - 1._prec)*Log(X)**2)/2._prec + Log(X)**3/6._prec
    ELSE
      YA=0._prec ; YB=1._prec/X
      S=-1._prec
      A=2*zeta2*Log(x)-Log(x)**3/6._prec
    END IF


    HA=-2._prec*YA-1._prec ; HB= 2._prec*YB
    ALFAA=HA+HA ; ALFAB = HB+HB

    BA0 = 0. ; BA1=0. ; BA2=0.
    BB0 = 0. ; BB1=0. ; BB2=0.
    DO  I = 18,0,-1
       BA0=CA(I)+ALFAA*BA1-BA2 ; BA2=BA1 ; BA1=BA0
       BB0=CB(I)+ALFAB*BB1-BB2 ; BB2=BB1 ; BB1=BB0
    ENDDO
    Li3 = A + S * (  (BA0 - HA*BA2) + (BB0 - HB*BB2) )
  END FUNCTION Li3

  FUNCTION trilog(x) result(res)
    ! evaluates trilog for any argument |x|<1
    complex(kind=prec) :: res
    complex(kind=prec) :: x
    if(abs(aimag(x)) < zero ) then
      res = Li3(real(x))
    else
      res = naive_polylog(3,x)
    endif
  END FUNCTION trilog

  FUNCTION BERNOULLI_POLYNOMIAL(n, x) result(res)
    integer n
    complex(kind=prec) :: x, res
    select case(n)
      case(1)
        res = -1/2. + x
      case(2)
        res = 1/6. - x + x**2
      case(3)
        res = x/2. - 3*x**2/2. + x**3
      case(4)
        res = -1/30. + x**2 - 2*x**3 + x**4
      case(5)
        res = -x/6. + 5*x**3/3 - 5*x**4/2 + x**5
      case(6)
        res = 1/42. - x**2/2 + 5*x**4/2 - 3*x**5 + x**6
      case(7)
        res = x/6. - 7*x**3/6 + 7*x**5/2 - 7*x**6/2 + x**7
      case(8)
        res = -1/30. + 2*x**2/3 - 7*x**4/3 + 14*x**6/3 - 4*x**7 + x**8
      case(9)
        res = -3*x/10 + 2*x**3 - 21*x**5/5 + 6*x**7 - 9*x**8/2 + x**9
      case(10)
        res = 5/66. - 3*x**2/2 + 5*x**4 - 7*x**6 + 15*x**8/2 - 5*x**9 + x**10
      case(11)
        res = 5*x/6 - 11*x**3/2 + 11*x**5 - 11*x**7 + 55*x**9/6 - 11*x**10/2 + x**11
      case(12)
        res = -691/2730. + 5*x**2 - 33*x**4/2 + 22*x**6 - 33*x**8/2 + 11*x**10 - 6*x**11 + x**12
      case(13)
        res = -691*x/210 + 65*x**3/3 - 429*x**5/10 + 286*x**7/7 - 143*x**9/6 + 13*x**11 - 13*x**12/2 + x**13
      case(14)
        res = 7/6. - 691*x**2/30 + 455*x**4/6 - 1001*x**6/10 + 143*x**8/2 - 1001*x**10/30 + 91*x**12/6 - 7*x**13 + x**14
      case(15)
        res = 35*x/2 - 691*x**3/6 + 455*x**5/2 - 429*x**7/2 + 715*x**9/6 - 91*x**11/2 + 35*x**13/2 - 15*x**14/2 + x**15
      case default
        print*,"Bernoulli beyond 15 is not implemented"
        stop
    end select


  END FUNCTION

  RECURSIVE FUNCTION polylog1(m,x) result(res)
    ! computes the polylog
    
    integer :: m
    type(inum) :: x, inv
    complex(kind=prec) :: res
    integer i

    
#ifdef DEBUG
    if(verb >= 70) print*, 'called polylog(',m,',',x%c,x%i0,')'
#endif
#ifndef NOCACHE
    if (m.le.5) then
      do i=1,plcachesize(m)
        if( abs(cache(m,i)%c%c-x%c).lt.zero ) then
          res = cache(m,i)%ans
          return
        endif
      enddo
    endif
#endif
    if ((m.le.9).and.(abs(x%c-1.).lt.zero)) then
      res = zeta(m)
    else if ((m.le.9).and.(abs(x%c+1.).lt.zero)) then
      res = -(1. - 2.**(1-m))*zeta(m)
    else if (abs(x) .gt. 1) then
      inv = inum(1./x%c, x%i0)
      res = (-1)**(m-1)*polylog(m,inv) &
          - cmplx(0,2*pi)**m * bernoulli_polynomial(m, 0.5-cmplx(0.,1.)*conjg(log(-x%c))/2/pi) / factorial(m)
    else if(m == 2) then
      res = dilog(x%c)
    else if(m == 3) then
      res = trilog(x%c)
    else
      res = naive_polylog(m,x%c)
    end if

#ifndef NOCACHE
    if (m.le.PolyLogCacheSize(1)) then
      if (plcachesize(m).lt.PolyLogCacheSize(2)) then
        plcachesize(m) = plcachesize(m) + 1
        cache(m,plcachesize(m)) = el(x,res)
      endif
    endif
#endif
  END FUNCTION polylog1




  RECURSIVE FUNCTION polylog2(m,x,y) result(res)
    type(inum) :: x, y
    integer m
    complex(kind=prec) :: res
    !TODO!!
    res=polylog1(m,inum(x%c/y%c,di0))
  END FUNCTION POLYLOG2


  FUNCTION PLOG1(a,b)
  ! calculates log(1-a/b)
  implicit none
  type(inum) :: a,b
  complex(kind=prec) plog1
  !TODO!!
  plog1 = log(1.-a%c/b%c)
  END FUNCTION

#ifndef NOCACHE
  SUBROUTINE CLEARCACHE
  plcachesize=0
  END SUBROUTINE
#endif

END MODULE maths_functions

! PROGRAM test
!   use maths_functions
!   implicit none
!   complex(kind=prec) :: res
!   res = Li3(0.4d0)
!   print*, res
! END PROGRAM test

