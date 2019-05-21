
                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE FUNCTIONS
                          !!!!!!!!!!!!!!!!!!!!!!

  USE GLOBAL_DEF

  IMPLICIT NONE

  TYPE MLM   ! MassLess Momentum
    real(kind=prec), dimension(4) :: momentum
    complex(kind=prec), dimension(2) :: sp    ! <+|
    complex(kind=prec), dimension(2) :: sm    ! <-|
    complex(kind=prec), dimension(2) :: ep    ! |+>
    complex(kind=prec), dimension(2) :: em    ! |->
  END TYPE MLM

  TYPE PARTON
    character(len=1) :: flavour  
          ! choices for flavour:
          ! "j": light quarks or gluons
          ! "b": bottom 
          ! "B": antibottom
          ! "F": fat, bottom + antibottom
    real (kind=prec), dimension(4) :: momentum
  END TYPE PARTON

  INTERFACE OPERATOR(+)
    module procedure add_parton
  END INTERFACE

  INTERFACE OPERATOR(-)
    module procedure minus_mlm
  END INTERFACE

  INTERFACE A
    module procedure Aversion2, Aversion3, Aversion4, Aversion5
  END INTERFACE

  INTERFACE B
    module procedure Bversion2, Bversion3, Bversion4, Bversion5
  END INTERFACE

  INTERFACE S
    module procedure Sversion1, Sversion2, Sversion3, Sversion4
  END INTERFACE

  INTERFACE SQ
    module procedure SQversion1
  END INTERFACE

  INTERFACE EIK
    module procedure eik00, eik0M, eikM0, eikMM
  END INTERFACE

  INTERFACE IREG
    module procedure ireg_00, ireg_0M, ireg_M0, ireg_MM, ireg_SELF
  END INTERFACE
  INTERFACE ISIN
    module procedure isin_00, isin_0m, isin_m0, isin_MM, isin_SELF
  END INTERFACE
  INTERFACE ILIN
    module procedure ilin_MM, ilin_SELF
  END INTERFACE

  INTERFACE DISCB
    module procedure discbSMN, discbSMM
  END INTERFACE
  INTERFACE SCALARC0IR6
    module procedure scalarc0ir6SMM, scalarc0ir6SMN
  END INTERFACE


  CONTAINS


  FUNCTION MAKE_MLM(qq)
  real(kind=prec), intent(in) :: qq(4)
  type(mlm) :: make_mlm
  complex(kind=prec) :: sq

  if(qq(4) > 0._prec) then
    sq = sqrt(qq(4) + qq(1) + Tiny(qq(4)))
  else
    sq = imag*sqrt(-qq(4) - qq(1) + Tiny(qq(4)))
  endif

  make_mlm%momentum = qq
  make_mlm%sp = (/ (-imag*qq(2)-qq(3))/sq, sq /)
  make_mlm%sm = (/ (-imag*qq(2)+qq(3))/sq, sq /)
  make_mlm%ep = (/ sq, (imag*qq(2)-qq(3))/sq /)
  make_mlm%em = (/ sq, (imag*qq(2)+qq(3))/sq /)

  END FUNCTION MAKE_MLM


  FUNCTION MINUS_MLM(ki)
  type(mlm) :: minus_mlm
  type(mlm), intent(in) :: ki
  complex(kind=prec) :: factor
 
  minus_mlm%momentum = -ki%momentum

  if(ki%momentum(4) > 0._prec) then
    factor = imag
  else
    factor = -imag
  endif

  minus_mlm%sp = factor*ki%sp
  minus_mlm%sm = factor*ki%sm
  minus_mlm%ep = factor*ki%ep
  minus_mlm%em = factor*ki%em

  END FUNCTION MINUS_MLM



  FUNCTION ADD_PARTON(par1,par2)
  type(parton) :: add_parton
  type(parton), intent(in) :: par1, par2
  character(len=1) :: fl1,fl2
 
  fl1 = par1%flavour
  fl2 = par2%flavour

  select case(fl1)
  case("j")
    add_parton%flavour = fl2
  case("b")
    select case(fl2)
    case("j")
      add_parton%flavour = fl1
    case("B")
      add_parton%flavour = "F"
    case default
      call crash("ADD_PARTON_1")
    end select
  case("B")
    select case(fl2)
    case("j")
      add_parton%flavour = fl1
    case("b")
      add_parton%flavour = "F"
    case default
      call crash("ADD_PARTON_2")
    end select
  case("F")
    select case(fl2)
    case("j")
      add_parton%flavour = fl1
    case default
      call crash("ADD_PARTON_3")
    end select
  case default
    call crash("ADD_PARTON_0")
  end select

  add_parton%momentum = par1%momentum + par2%momentum

  END FUNCTION ADD_PARTON



  FUNCTION Sversion1(k1,k2)
     ! S(k1,k2) = 2 k1.k2

  type(mlm), intent(in) :: k1,k2
  real (kind=prec) :: Sversion1,dot_dot,q1(4),q2(4)
  q1 = k1%momentum
  q2 = k2%momentum

  dot_dot = q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3)
  Sversion1 =  2*dot_dot 

  END FUNCTION Sversion1  


  FUNCTION Sversion2(q1,q2)
     ! S(q1,q2) =  2 q1.q2

  real (kind=prec), intent(in) :: q1(4),q2(4)
  real (kind=prec) :: Sversion2,dot_dot
  dot_dot = q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3)
  Sversion2 =  2*dot_dot 

  END FUNCTION Sversion2  



  FUNCTION Sversion3(k1,q2)
     ! S(q1,q2) =  2 q1.q2

  type(mlm), intent(in) :: k1
  real (kind=prec), intent(in) :: q2(4)
  real (kind=prec) :: q1(4)
  real (kind=prec) :: Sversion3,dot_dot
  q1 = k1%momentum
  dot_dot = q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) 
  Sversion3 =  2*dot_dot 

  END FUNCTION Sversion3  



  FUNCTION Sversion4(q1,k2)
     ! S(q1,q2) = 2 q1.q2

  type(mlm), intent(in) :: k2
  real (kind=prec), intent(in) :: q1(4)
  real (kind=prec) :: q2(4)
  real (kind=prec) :: Sversion4,dot_dot
  q2 = k2%momentum
  dot_dot = q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) 
  Sversion4 =  2*dot_dot 

  END FUNCTION Sversion4  



  FUNCTION SQversion1(q1)
     ! SQ(q1) = (q1)^2

  real (kind=prec), intent(in) :: q1(4)
  real (kind=prec) :: SQversion1

  SQversion1  = q1(4)**2 - q1(1)**2 - q1(2)**2 - q1(3)**2

  END FUNCTION SQversion1  


  FUNCTION GRAM3(s1,s2,s3)

  real (kind=prec), intent(in) :: s1,s2,s3
  real (kind=prec) :: gram3
  gram3 = s1**2 + s2**2 + s3**2 - 2*(s1*s2 + s2*s3+ s3*s1)

  END FUNCTION GRAM3   


!==================  definition of angles and brackets ==============!

#define MakeMatrixP(q1_mat,q1) \
  q1_mat(1,:)  =  (/ -imag*q1(2) - q1(3), cmplx(q1(4) + q1(1)) /) ;\
  q1_mat(2,:)  =  (/  cmplx(q1(4) - q1(1)), imag*q1(2) - q1(3) /)

#define MakeMatrixM(q1_mat,q1) \
  q1_mat(1,:)  =  (/ -imag*q1(2) + q1(3), cmplx(q1(4) + q1(1)) /) ;\
  q1_mat(2,:)  =  (/  cmplx(q1(4) - q1(1)), imag*q1(2) + q1(3) /)


  FUNCTION Aversion2(k1,k2)

    complex(kind=prec) :: Aversion2
    type(mlm) :: k1, k2

    Aversion2 = dot_product(conjg(k1%sm),k2%ep)

  END FUNCTION Aversion2


  FUNCTION Aversion3(k1,q1,k2)

    complex(kind=prec) :: Aversion3, q1_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4)
    
    MakeMatrixP(q1_mat,q1)

    Aversion3 = dot_product(conjg(k1%sm),matmul(q1_mat,k2%em))

  END FUNCTION Aversion3


  FUNCTION Aversion4(k1,q1,q2,k2)

    complex(kind=prec) :: Aversion4, q1_mat(2,2), q2_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4), q2(4)
    
    MakeMatrixP(q1_mat,q1)
    MakeMatrixM(q2_mat,q2)

    Aversion4 = dot_product(conjg(k1%sm),matmul(q1_mat,matmul(q2_mat,k2%ep)))

  END FUNCTION Aversion4


  FUNCTION Aversion5(k1,q1,q2,q3,k2)

    complex(kind=prec) :: Aversion5, q1_mat(2,2), q2_mat(2,2), q3_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4), q2(4), q3(4)
    
    MakeMatrixP(q1_mat,q1)
    MakeMatrixM(q2_mat,q2)
    MakeMatrixP(q3_mat,q3)

    Aversion5 = dot_product(conjg(k1%sm), & 
                     matmul(q1_mat,matmul(q2_mat,matmul(q3_mat,k2%em))))

  END FUNCTION Aversion5


  FUNCTION Bversion2(k1,k2)

    complex(kind=prec) :: Bversion2
    type(mlm) :: k1, k2

    Bversion2 = dot_product(conjg(k1%sp),k2%em)

  END FUNCTION Bversion2


  FUNCTION Bversion3(k1,q1,k2)

    complex(kind=prec) :: Bversion3, q1_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4)
    
    MakeMatrixM(q1_mat,q1)

    Bversion3 = dot_product(conjg(k1%sp),matmul(q1_mat,k2%ep))

  END FUNCTION Bversion3


  FUNCTION Bversion4(k1,q1,q2,k2)

    complex(kind=prec) :: Bversion4, q1_mat(2,2), q2_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4), q2(4)
    
    MakeMatrixM(q1_mat,q1)
    MakeMatrixP(q2_mat,q2)

    Bversion4 = dot_product(conjg(k1%sp),matmul(q1_mat,matmul(q2_mat,k2%em)))

  END FUNCTION Bversion4


  FUNCTION Bversion5(k1,q1,q2,q3,k2)

    complex(kind=prec) :: Bversion5, q1_mat(2,2), q2_mat(2,2), q3_mat(2,2)
    type(mlm) :: k1, k2
    real(kind=prec) :: q1(4), q2(4), q3(4)
    
    MakeMatrixM(q1_mat,q1)
    MakeMatrixP(q2_mat,q2)
    MakeMatrixM(q3_mat,q3)

    Bversion5 = dot_product(conjg(k1%sp),   &
                     matmul(q1_mat,matmul(q2_mat,matmul(q3_mat,k2%ep))))

  END FUNCTION Bversion5


!==============  end definition of angles and brackets ==============!


!==============  eikonal factors  ==============!


  FUNCTION EIK00(k1,kg,k2)

    type(mlm), intent(in) :: k1,k2,kg
    real (kind=prec) :: eik00

    eik00 = 2*s(k1,k2)/(s(k1,kg)*s(kg,k2))

  END FUNCTION EIK00


  FUNCTION EIK0M(k1,kg,q2)

    type(mlm), intent(in) :: k1,kg
    real (kind=prec), intent(in) :: q2(4)
    real (kind=prec) :: eik0M

    eik0m = 2*s(k1,q2)/(s(k1,kg)*s(kg,q2)) - 2*sq(q2)/(s(kg,q2)**2)

  END FUNCTION EIK0M


  FUNCTION EIKM0(q1,kg,k2)

    type(mlm), intent(in) :: k2,kg
    real (kind=prec), intent(in) :: q1(4)
    real (kind=prec) :: eikM0

    eikm0 = 2*s(k2,q1)/(s(k2,kg)*s(kg,q1)) - 2*sq(q1)/(s(kg,q1)**2)

  END FUNCTION EIKM0


  FUNCTION EIKMM(q1,kg,q2)

    type(mlm), intent(in) :: kg
    real (kind=prec), intent(in) :: q1(4), q2(4)
    real (kind=prec) :: eikMM

    eikmm = 2*s(q1,q2)/(s(q1,kg)*s(kg,q2))     &
         - 2*sq(q1)/(s(kg,q1)**2) - 2*sq(q2)/(s(kg,q2)**2)

  END FUNCTION EIKMM


   
!===========  end of eikonal factors and AP splitting kernels ===========!



  FUNCTION LI2(x)

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

   if(x > 1.00000000001_prec) then
     call crash("LI2")
   elseif(x > 1.0_prec) then
     x = 1._prec
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

  END FUNCTION LI2

  FUNCTION LI3(x)
   ! This was hacked from LI2 to also follow C332
   !! Trilogarithm for arguments x < = 1.0
   ! In theory this could also produce Re[Li [x]] for x>1
   !                                        3

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
     call crash("LI3")
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

  END FUNCTION LI3




  FUNCTION LN(arg1, arg2)

    !! log(sij/mm) for positive and negative sij 

    complex (kind=prec) :: Ln
    real(kind=prec), intent(in) :: arg1, arg2
    real(kind=prec) :: abs_part

    abs_part = 1._prec
    
    if(arg1 > 0._prec) then
      ln = log(arg1/arg2)
    elseif(arg1 < 0._prec ) then
      ln = log(-arg1/arg2) - imag*abs_part*pi
    else
      call crash("LNc1")
    endif

    if(arg2 < 0._prec) call crash("LN")

  END FUNCTION LN



  FUNCTION DILOG(arg1,arg2)

    !! dilogarithm of 1 - arg1/arg2 for pos and neg arg2 

    complex (kind=prec) :: Dilog
    real(kind=prec), intent(in) :: arg1,arg2
    real(kind=prec) :: abs_part

    abs_part = 1._prec
   
    if(arg1/arg2 > 0._prec) then
      dilog = li2(1._prec - arg1/arg2)
    elseif(arg1 > 0._prec .and. arg2 < 0._prec) then
      dilog = pi**2/6 - li2(arg1/arg2)                             &
          + (log(-arg2/arg1) - imag*abs_part*pi)*log(1._prec - arg1/arg2) 
    else
     print*, "arg1, arg2: ",arg1,arg2
     call crash("DILOG")
    endif

  END FUNCTION DILOG
   

!======================    soft integrals    ======================!


  FUNCTION ISIN_00(xicut,epcmf,k1,k2) ! both massless
  real (kind=prec), intent(in) :: epcmf,xicut
  type(mlm), intent(in) :: k1,k2
  real (kind=prec) isin_00, q1(4), q2(4)
  q1 = k1%momentum
  q2 = k2%momentum

  isin_00 = -(log(s(k1,k2)/musq)-log(4*q1(4)*q2(4)/xicut**2/epcmf**2))
  END FUNCTION ISIN_00

  FUNCTION ISIN_M0(xicut,epcmf,q1,k2)  ! first massive second massless
  real (kind=prec), intent(in) :: epcmf,xicut,q1(4)
  type(mlm), intent(in) :: k2
  real (kind=prec) isin_m0,m12,q2(4)
  m12 = sq(q1)
  q2 = k2%momentum
  isin_m0 = -(log(s(q1,k2)/musq)-0.5*log(4.*m12**2*q2(4)**2/xicut**2/epcmf**2/musq))
  END FUNCTION ISIN_M0
  FUNCTION ISIN_0M(xicut,epcmf,k1,q2)
           ! first massless second massive, "same" as M0
  real (kind=prec), intent(in) :: epcmf,xicut,q2(4)
  type(mlm), intent(in) :: k1
  real (kind=prec) isin_0m,m12,q1(4)
  m12 = sq(q1)
  q1 = k1%momentum
  isin_0m = -(log(s(q2,k1)/musq)-0.5*log(4.*m12**2*q1(4)**2/xicut**2/epcmf**2/musq))
  END FUNCTION ISIN_0M

  FUNCTION ISIN_SELF(xicut,epcmf,q1)  ! massive self-eikonal
  real (kind=prec), intent(in) :: q1(4),epcmf,xicut
  real (kind=prec) :: isin_self

  isin_self = -1

  END FUNCTION ISIN_SELF

  FUNCTION ISIN_MM(xicut,epcmf,qk,ql)  ! both massive eq(A.11) of FFMS

  real (kind=prec), intent(in) :: qk(4),ql(4),epcmf,xicut
  real (kind=prec) :: isin_mm, vkl
  real (kind=prec) :: KdotL, mk2, ml2


  mk2 = sq(qk)
  ml2 = sq(ql)
  KdotL = 0.5_prec*s(ql,qk)
  vkl = sqrt(1. - mk2*ml2/KdotL**2)


  if(vkl > 1.0E-4) then
     isin_mm = - 0.5_prec/vkl*log((1.+vkl)/(1.-vkl))
  else
     isin_mm = - 0.5* (2 + (2*vkl**2)/3. + (2*vkl**4)/5. + (2*vkl**6)/7. )
  endif


  END FUNCTION ISIN_MM


  FUNCTION IREG_00(xicut,epcmf,k1,k2) ! both massless

  real (kind=prec), intent(in) :: epcmf,xicut
  type(mlm), intent(in) :: k1,k2
  real (kind=prec) :: ireg_00, e1, e2, arg1, arg2, arg24
  real (kind=prec) :: log1, log2, log24, twologtwo, q1(4), q2(4)

  q1 = k1%momentum     
  q2 = k2%momentum     

  e1 = q1(4); e2 = q2(4)
  arg1 = xicut**2*epcmf**2/musq
  arg2 = s(k1,k2)/e1/e2
  arg24 = arg2/4._prec
  log1 = log(arg1)
  log2 = log(arg2)
  log24 = log(arg24)
  twologtwo = 0.960906027836403_prec    ! 2 (log 2)^2

  if(arg2 > 4.000001_prec) call crash("IREG")
  if(arg2 > 3.999999_prec) arg2 = 3.999999_prec 

  ireg_00 = 0.5_prec*log1**2 + log1*log24 - li2(arg24)          &  
          + 0.5_prec*log2**2 - log(4._prec - arg2)*log24        &
          - twologtwo

  select case(cgamma)
  case("exp")
    ireg_00 = ireg_00 - pi**2/12.
  case("gam")
    ireg_00 = ireg_00 - 0. ! do nothing
  case default
    call crash("IREG_00")
  end select
 
    ! NOTE  -pi^2/12 term to change convention in prefactor
    ! Note that the overall factor 1/8/pi**2 is NOT included
    ! since the psi_mn have also a different normalization
  
  END FUNCTION IREG_00       


  FUNCTION IREG_M0(xicut,epcmf,q1,k2)  ! first massive second massless

  real (kind=prec), intent(in) :: q1(4),epcmf,xicut
  type(mlm), intent(in) :: k2
  real (kind=prec) :: ireg_m0, ml, ek, el, bl, q2(4)
  real (kind=prec) :: logX, logSQ, logKK, KdotL
     
  logX = log(xicut)
  logSQ = log(epcmf**2/musq)

  q2 = k2%momentum
  ek = q2(4); el = q1(4)
  KdotL = 0.5_prec*s(q1,k2)
  ml = sqrt(0.5_prec*s(q1,q1))
  bl = sqrt(q1(1)**2+q1(2)**2+q1(3)**2)/el
!  bl = sqrt(1._prec - (ml/el)**2)

  logKK = log(KdotL/ml/ek)
  ireg_m0 = logX*(logX+logSQ+2.*logKK) + 0.25*logSQ**2              &
          - 0.25*(log((1.+bl)/(1.-bl)))**2                          &
          + 0.5*(log(KdotL/(1.-bl)/el/ek))**2  + logSQ*logKK        &
          - li2(1. - (1.+bl)*ek*el/KdotL)                           &
          + li2(1. - KdotL/(1.-bl)/ek/el)

  select case(cgamma)
  case("exp")
    ireg_m0 = ireg_m0 - pi**2/8.
  case("gam")
    ireg_m0 = ireg_m0 - pi**2/12.
  case default
    call crash("IREG_M0")
  end select

    ! NOTE -pi^2/12 -> -pi^2/12 -pi^2/24 = -pi^2/8
    ! to change convention in prefactor
  
  END FUNCTION IREG_M0       



  FUNCTION IREG_0M(xicut,epcmf,k1,q2)  
           ! first massless second massive, "same" as M0

  real (kind=prec), intent(in) :: q2(4),epcmf,xicut
  type(mlm), intent(in) :: k1
  real (kind=prec) :: ireg_0m, ml, ek, el, bl
  real (kind=prec) :: logX, logSQ, logKK, KdotL
     
  logX = log(xicut)
  logSQ = log(epcmf**2/musq)

  ek = k1%momentum(4); el = q2(4)
  KdotL = 0.5_prec*s(k1,q2)
  ml = sqrt(0.5_prec*s(q2,q2))
  bl = sqrt(q2(1)**2+q2(2)**2+q2(3)**2)/el
!  bl = sqrt(1._prec - (ml/el)**2)

  logKK = log(KdotL/ml/ek)
  ireg_0m = logX*(logX+logSQ+2.*logKK)  + 0.25*logSQ**2             &
          - 0.25*(log((1.+bl)/(1.-bl)))**2                          &
          + 0.5*(log(KdotL/(1.-bl)/el/ek))**2  + logSQ*logKK        &
          - li2(1. - (1.+bl)*ek*el/KdotL)                           &
          + li2(1. - KdotL/(1.-bl)/ek/el)

  select case(cgamma)
  case("exp")
    ireg_0m = ireg_0m - pi**2/8.
  case("gam")
    ireg_0m = ireg_0m - pi**2/12.
  case default
    call crash("IREG_0M")
  end select

    ! NOTE -pi^2/12 -> -pi^2/12 -pi^2/24 = -pi^2/8
    ! to change convention in prefactor
  
  END FUNCTION IREG_0M       



  FUNCTION IREG_SELF(xicut,epcmf,q1)  ! massive self-eikonal

  real (kind=prec), intent(in) :: q1(4),epcmf,xicut
  real (kind=prec) :: ireg_self, ml, ek, el, bl
  real (kind=prec) :: logX, logSQ, logKK, KdotL
     
  logX = log(xicut)
  logSQ = log(epcmf**2/musq)

  el = q1(4)  
  bl = sqrt(q1(1)**2+q1(2)**2+q1(3)**2)/el

  if(bl > 1.0E-4) then
     ireg_self = 2.*logX + logSQ - 1._prec/bl*log((1.+bl)/(1.-bl))
  else
     ireg_self = 2.*logX + logSQ - (2 + (2*bl**2)/3. + (2*bl**4)/5. + (2*bl**6)/7. )
  endif

  END FUNCTION IREG_SELF       




  FUNCTION IREG_MM(xicut,epcmf,qk,ql)  ! both massive eq(A.12) of FFMS

  real (kind=prec), intent(in) :: qk(4),ql(4),epcmf,xicut
  real (kind=prec) :: ireg_mm, vkl, akl, bk, bl, nu, lam_nu_2
  real (kind=prec) :: KdotL, mk2, ml2, ek, el, LogV, LogX

  ek = qk(4)
  el = ql(4)
  bk = sqrt(qk(1)**2+qk(2)**2+qk(3)**2)/ek
  bl = sqrt(ql(1)**2+ql(2)**2+ql(3)**2)/el

  mk2 = sq(qk)
  ml2 = sq(ql)
  KdotL = 0.5_prec*s(ql,qk)
  vkl = sqrt(1. - mk2*ml2/KdotL**2)
  akl = (1.+vkl)*KdotL/mk2
  lam_nu_2 = akl**2*mk2 - ml2  ! 2 \lambda \nu
  nu = 0.5_prec*lam_nu_2/(akl*ek - el)
!  lam = akl*ek - el
!  nu = 0.5_prec*(akl**2*mk2 - ml2)/lam

  logX = log((xicut*epcmf)**2/musq)
  logV = log((1.+vkl)/(1.-vkl))

  ireg_mm = jaux(akl*ek, akl*ek*bk, nu, lam_nu_2)   &
          - jaux(el, el*bl, nu, lam_nu_2)
  ireg_mm = 0.5*(1.+vkl)*KdotL**2*ireg_mm/mk2
  ireg_mm = ireg_mm + 0.5*logV*logX/vkl

  END FUNCTION IREG_MM

  
  FUNCTION JAUX(xx,yy,nu,lam_nu_2)

  real(kind=prec) :: jaux, xx, yy, nu, lam_nu_2

  jaux = (log((xx-yy)/(xx+yy)))**2        &
       + 4*Li2(1._prec - (xx+yy)/nu)      &
       + 4*Li2(1._prec - (xx-yy)/nu)  

  jaux = jaux/lam_nu_2

  END FUNCTION JAUX  



  FUNCTION ILIN_SELF(xicut,epcmf,q1)  ! massive self-eikonal
  real (kind=prec), intent(in) :: q1(4),epcmf,xicut
  real (kind=prec) :: ilin_self, ml, ek, el, bl, x
  real (kind=prec) :: logX, logSQ, logC, log0, log10

  logX = log(xicut)
  logSQ = log(epcmf**2/musq)

  el = q1(4)
  bl = sqrt(q1(1)**2+q1(2)**2+q1(3)**2)/el
  logC = logSQ+2*logX

  if(bl > 1.0E-4) then
    x = (1-bl) / (1+bl)
    log0 = log(x)
    log10 = -log0*log(1.-x) - Li2(x)

    ilin_self = -(3*log0**2+12*log10+6*log0*logC+2*pi**2+3*logC**2*bl)/6./bl

  else
    ilin_self = - 4. - 0.5_prec*logC**2 + 2*logC*(1._prec + bl**2/3. + bl**4/5. + bl**6/7.)
  endif

  END FUNCTION ILIN_SELF

  FUNCTION ILIN_MM(xicut,epcmf,qk,ql)  ! both massive eq(A.12) of FFMS
  real (kind=prec), intent(in) :: qk(4),ql(4),epcmf,xicut
  real (kind=prec) :: ilin_mm, bl, mk2, ml2,x, log0, log10, logC,logStuff


  if (abs(qk(1))+abs(qk(2))+abs(qk(3)).gt.zero) call crash("ILIN_MM only for CMS of qk yet!")

  mk2 = sq(qk)
  ml2 = sq(ql)
  bl = sqrt(1. - mk2*ml2/(0.5*s(qk,ql))**2)

  logC = log(epcmf**2/musq) + 2*log(xicut)

  x = (1-bl) / (1+bl)
  log0 = log(x)
  log10 = -log0*log(1.-x) - Li2(x)

  logStuff = -0.5 * log0**2*log(1.-x) - Li3(x) - 2*Li3(1-x) + z3

  ilin_mm = (log0**3 + 3*log0**2*logC + 12*log10*logC + 3*log0*logC**2 +  12*logStuff + log0*Pi**2 + 2*logC*Pi**2)/(12.*bl)


  END FUNCTION ILIN_MM




!=====================  end soft integrals ====================!



  FUNCTION QTIL(cc,xicut,epcmf,q1)

  real (kind=prec), intent(in) :: q1(4),epcmf,xicut
  character (len=1), intent(in) :: cc
  real (kind=prec) :: qtil, e1, gamma, gammap, ccol

  e1 = q1(4)

  select case(cc)
  case("l")
    gamma = 1.5_prec
    gammap = (6.5_prec - 2*pi**2/3.)
    ccol = 1._prec
  case("q")
    gamma = 1.5_prec*Cf
    gammap = (6.5_prec - 2*pi**2/3.)*Cf
    ccol = Cf
  case("g")
    gamma = (11.*Nc - 4*Tf*5.)/6.  !!Nf -> 5   
    gammap = 67.*Nc/9. - 2*pi**2*Nc/3. - 23.*Tf*5./9.   !!Nf -> 5              
    ccol = Nc
  case default
    call crash("QTIL")
  end select

  qtil = gammap - gamma*log(2*delcut*e1**2/musq) +                 &
         2*ccol*log(2*e1/epcmf/xicut)*                             &
                 log(epcmf*delcut*e1*xicut/musq)

  END FUNCTION QTIL




  FUNCTION BOOST_BACK(rec,mo)   !!boosts to cms system  

  real (kind=prec), intent(in):: rec(4),mo(4)    
  real (kind=prec)  :: cosh_a, energy,  dot_dot,   &
                       n_vec(3), boost_back(4)

  energy = rec(4)**2 - rec(1)**2 - rec(2)**2 - rec(3)**2
  
  if(energy < 1.0E-16_prec) then
    throw_away = 1
    energy = 1.0E-12_prec   
  else
    energy = sqrt(energy)
  end if

  cosh_a = rec(4)/energy
  n_vec = - rec(1:3)/energy  ! 1/sinh_a omitted

  dot_dot = sum(n_vec*mo(1:3))  ! \vec{n} \dot \vec{m}

  boost_back(1:3) =    &
     mo(1:3) + n_vec*(dot_dot/(cosh_a + 1) - mo(4))
  boost_back(4) = mo(4)*cosh_a - dot_dot

  END FUNCTION BOOST_BACK



  FUNCTION BOOST_RFREC(rec,mo1)   !!boosts mo to rest-frame of rec

  real (kind=prec), intent(in):: rec(4),mo1(4) 
  real (kind=prec)  :: energy, beta(3), gamma, betaSQ, dot_dot, betaABS
  real (kind=prec)  :: boost_rfrec(4)

  energy = rec(4)
  beta = rec(1:3)/energy
  betaSQ = beta(1)**2 + beta(2)**2 + beta(3)**2
  betaABS = sqrt(betaSQ)
  if(betaSQ > zero) then
    beta = beta/betaABS
  else
    beta = (/0._prec,0._prec,0._prec /)
  endif
  gamma = 1./sqrt(1._prec-betaSQ)
  
  dot_dot = sum(beta*mo1(1:3))  ! \vec{beta} \dot \vec{m}

  boost_rfrec(1:3) =    &
     mo1(1:3) + beta*dot_dot*(gamma-1.) - beta*betaABS*gamma*mo1(4)

  boost_rfrec(4) = gamma*(mo1(4)-betaABS*dot_dot)

!!$  boost_rfrec(1:3) =    &
!!$     mo1(1:3) + beta*(dot_dot*(gamma-1.)/betaSQ - gamma*mo1(4))
!!$
!!$  boost_rfrec(4) = gamma*(mo1(4)-dot_dot)

  END FUNCTION BOOST_RFREC



  ! this used to be in phase_space but we need it here as well

  FUNCTION SQ_LAMBDA(ss,m1,m2) !! square root of \lambda(s,m1,m2)

  implicit none
  real (kind=prec) ss,m1,m2,sq_lambda

  sq_lambda = (ss - (m1+m2)**2)*(ss - (m1-m2)**2)
  sq_lambda = sqrt(sq_lambda)

  END FUNCTION SQ_LAMBDA
  FUNCTION CSQ_LAMBDA(ss,m1,m2) !! complex square root of \lambda(s,m1,m2)

  implicit none
  real (kind=prec) ss,m1,m2
  complex (kind=prec) csq_lambda

  csq_lambda = (ss - (m1+m2)**2)*(ss - (m1-m2)**2)
  csq_lambda = sqrt(csq_lambda)

  END FUNCTION CSQ_LAMBDA

!======================    Package X abbr.   ======================!

  FUNCTION DISCBSMN(s,m1,m2)
  ! Package X function DiscB[s,m1,m2]
  real(kind=prec), intent(in) :: s, m1, m2
  complex(kind=prec) :: discbsmn, l
  l = csq_lambda(s,m1,m2)
  discbsmn = l*log( (m1**2+m2**2-s+l)/(2*m1*m2) )/s

  END FUNCTION DISCBSMN

  FUNCTION DISCBSMM(s,m1)
  ! Package X function Re[DiscB[s,m1,m1]]
  real(kind=prec), intent(in) :: s, m1
  real(kind=prec) :: discbsmm, b


  b = sqrt(s*(s-4*m1**2))

  ! Ignore im part
  discbsmm = b*log(abs(1+(b-s)/(2*m1**2)))/s

  END FUNCTION DISCBSMM


  FUNCTION SCALARC0IR6SMM(s,m1)
  ! Package X function Re[ScalarC0IR6[s,m1,m1]]
  real(kind=prec), intent(in) :: s, m1
  real(kind=prec) :: scalarc0ir6smm, b, discS,ImdiscS

  b = sqrt(s*(s-4*m1**2))
  discS = b*log(abs(1+(b-s)/(2*m1**2))) ! Re[DiscB[S,m,m]] S
  imdiscs = b*pi


  scalarc0ir6smm = 0.5*discS*log(s**2*abs(2*m1**4-m1**2*(S+b))/(2*b**4))/b**2 ! Re
  scalarc0ir6smm = scalarc0ir6smm - 0.5*pi*imdiscs/b**2                       ! Im
  scalarc0ir6smm = scalarc0ir6smm - pi**2/(6*b) - 2*Li2((S-b)/(2*m1**2)-1) / b


  END FUNCTION SCALARC0IR6SMM
  FUNCTION SCALARC0IR6SMN(s,m1,m2)
  ! Package X function ScalarC0IR6[s,m1,m2]
  real(kind=prec), intent(in) :: s, m1,m2
  real(kind=prec) :: scalarc0ir6smn
    complex(kind=prec) C0
    external C0
  ! TODO relies on LoopTools

  call setlambda(0.)
  call setmudim(1.)
  scalarc0ir6smn = -0.5*s*DiscB(s,m1,m2) / CSQ_LAMBDA(s,m1,m2)**2
  scalarc0ir6smn = scalarc0ir6smn * ( 2*log(m1/m2) - 4*log(m1) )
  scalarc0ir6smn = scalarc0ir6smn + real(C0(m1**2,S,m2**2,0.,m1**2,m2**2))


  END FUNCTION SCALARC0IR6SMN


  FUNCTION SCALARD0IR16(s2,s3,s12,s23,m1,m2,m3)
  real(kind=prec), intent(in) :: s2,s3,s12,s23,m1,m2,m3
  real(kind=prec) :: scalard0ir16
    integer, parameter :: dd0  = 1
    complex(kind=prec) D0, D0i
    external D0, D0i
  ! TODO relies on LoopTools

  call setlambda(0.)
  call setmudim(1.)
  scalard0ir16 = real(DiscB(s23,m2,m3))*s23 / SQ_LAMBDA(s23,m2,m3)**2 / (m1**2-s12)
  scalard0ir16 = scalard0ir16 * ( log(m2/m3) -2*log(m2) )
  scalard0ir16 = scalard0ir16 + D0i(dd0,m2**2,m3**2,s3,s2,s23,s12,m2**2,0.,m3**2,m1**2)

  END FUNCTION SCALARD0IR16

!=====================  end Package X abbr ====================!
 


!======================    Former user.f95   ======================!

  FUNCTION BOOST_RF(rec,mo)   !!boosts mo to (non-unique) rest frame of rec

  real (kind=prec), intent(in):: rec(4), mo(4)
  real (kind=prec)  :: energy,  mass, dot_dot, boost_rf(4)
  integer :: throw_away

  mass = rec(4)**2 - rec(1)**2 - rec(2)**2 - rec(3)**2

  if(mass < 1.0E-16_prec) then  ! rec must not be massless
    throw_away = 1
    mass = 1.0E-12_prec
  else
    mass = sqrt(mass)
  end if

  energy = rec(4)
  dot_dot = sum(rec(1:3)*mo(1:3))

  boost_rf(4) = (energy*mo(4)-dot_dot)/mass
  boost_rf(1:3) =    &
     mo(1:3) + rec(1:3)*(dot_dot/(energy+mass) - mo(4))/mass

  END FUNCTION BOOST_RF



  FUNCTION ETA(k1)   ! pseudo-rapidity

  real (kind=prec), intent(in) :: k1(4)
  real (kind=prec):: eta
  real (kind=prec) :: kv

  kv = sqrt(k1(1)**2+k1(2)**2+k1(3)**2)

  eta = 0.5*log((kv+k1(3))/(kv-k1(3)))

  END FUNCTION ETA



  FUNCTION RAP(k1)   ! rapidity

  real (kind=prec), intent(in) :: k1(4)
  real (kind=prec):: rap

  rap = 0.5*log((k1(4)+k1(3))/(k1(4)-k1(3)))

  END FUNCTION RAP



  FUNCTION PT(q1)   ! transverse momentum

  real (kind=prec), intent(in) :: q1(4)
  real (kind=prec):: pt

  pt = sqrt(q1(1)**2 + q1(2)**2)

  END FUNCTION PT



  FUNCTION ABSVEC(q1)   ! magnitude of 3-vector

  real (kind=prec), intent(in) :: q1(4)
  real (kind=prec):: absvec

  absvec = sqrt(q1(1)**2 + q1(2)**2 + q1(3)**2)

  END FUNCTION ABSVEC



  FUNCTION RIJ(k1,k2)

  real (kind=prec), intent(in) :: k1(4), k2(4)
  real (kind=prec) ::  rij, yij, phiij, Dres

  Dres = 0.7_prec

  yij = 0.5*log((k1(4)+k1(3))*(k2(4)-k2(3))/(k1(4)-k1(3))/(k2(4)+k2(3)))
  phiij = k1(1)*k2(1) + k1(2)*k2(2)
  phiij = phiij/sqrt(k1(1)**2+k1(2)**2)/sqrt(k2(1)**2+k2(2)**2)
  phiij = acos(phiij)

!  phiij = atan2(k1(2),k1(1)) - atan2(k2(2),k2(1))
!  if(phiij > pi) phiij = phiij - 2*pi
!  if(phiij < -pi) phiij = phiij + 2*pi

  rij = (yij**2 + phiij**2)/Dres**2

  END FUNCTION RIJ



  FUNCTION COS_TH(k1,k2)

  real (kind=prec), intent(in) :: k1(4), k2(4)
  real (kind=prec) ::  cos_th, mag1, mag2

  mag1 = sqrt(k1(1)**2 + k1(2)**2 + k1(3)**2)
  mag2 = sqrt(k2(1)**2 + k2(2)**2 + k2(3)**2)

  cos_th = sum(k1(1:3)*k2(1:3))/mag1/mag2

  END FUNCTION COS_TH


!======================  end user.f95        ======================!


                            !!!!!!!!!!!!!!!!!!!!!!!!!
                              END MODULE FUNCTIONS
                            !!!!!!!!!!!!!!!!!!!!!!!!!



