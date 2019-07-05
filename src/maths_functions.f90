
MODULE maths_functions
  use globals
  use utils
  implicit none

CONTAINS 

  FUNCTION naive_polylog(m,x,n_passed) result(res)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    integer, optional :: n_passed
    complex(kind=prec) :: x, res
    integer :: i,n
    integer, allocatable :: j(:)
    n = 1000
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

  RECURSIVE FUNCTION dilog(x) result(res)
    ! evaluates dilog for any argument
    complex(kind=prec) :: res
    complex(kind=prec) :: x

    if(abs(x) <= 1.0) then
      if(abs(aimag(x)) < zero ) then
        res = Li2(real(x))
      else
        res = naive_polylog(2,x)
      endif
    else
     res = -dilog(1/x) - pi**2/6 - log(add_ieps(-x))**2 / 2
   end if
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
    ELSE IF(X .EQ. -1._prec) THEN
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
    ! evaluates trilog for any argument
    complex(kind=prec) :: res
    complex(kind=prec) :: x
    if(abs(x) <= 1.0) then
      if(abs(aimag(x)) < zero ) then
        res = Li3(real(x))
      else
        res = naive_polylog(3,x)
      endif
    else
     res = naive_polylog(3,sub_ieps(x)**(-1)) - (log(-sub_ieps(x)))**3/6 - pi**2/6 * log(-sub_ieps(x))
   end if
  END FUNCTION trilog

  FUNCTION polylog(m,x) result(res)
    ! computes the polylog for now naively (except for dilog half-naively)
    integer :: m
    complex(kind=prec) :: x,res
      
    print*, 'called polylog with m = ', m

    if(verb >= 70) print*, 'called polylog(',m,',',x,')'
    if(m == 2) then
      res = dilog(x)
    else if(m == 3) then
      res = trilog(x)
    else
      res = naive_polylog(m,x)
    end if
  END FUNCTION polylog

END MODULE maths_functions

! PROGRAM test
!   use maths_functions
!   implicit none
!   complex(kind=prec) :: res
!   res = Li3(0.4d0)
!   print*, res
! END PROGRAM test

