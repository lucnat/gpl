 MODULE GLOBAL_DEF

 
 implicit none 
  
        !!  ----------
      	!!  parameters
        !!  ----------

 integer, parameter :: prec = selected_real_kind(15,32)
 real (kind=prec), parameter :: cw = 0.876613
 real (kind=prec), parameter :: sw = 0.481196
 real (kind=prec), parameter :: pi = 3.14159265358979323846_prec
 real (kind=prec), parameter :: z3 = 1.20205690315959428540_prec
 real (kind=prec), parameter :: log2 = 0.693147180559945309417_prec
 real (kind=prec), parameter :: conv = 3.893850E+8  ! convert GeV to pb
 real (kind=prec), parameter :: xsc = 0._prec  ! FDH=>1 vs HV=>0
 real (kind=prec), parameter :: Nc = 3._prec
 real (kind=prec), parameter :: Tf = 0.5_prec
 real (kind=prec), parameter :: Cf = (Nc**2-1)/(2*Nc)
 real (kind=prec), parameter :: Nh = 1._prec
 real (kind=prec), parameter :: Nf = 5._prec

 complex (kind=prec), parameter :: imag = (0.0_prec,1.0_prec)
 real (kind=prec), parameter :: zero = 1.0E-50_prec
 
 real (kind=prec), parameter :: alpha_ew = 0.03394_prec   
! real (kind=prec), parameter :: alpha = 1/127.9_prec
! real (kind=prec), parameter :: alpha = 1./137.0359997_prec
! real (kind=prec), parameter :: GF = 1.16637E-11_prec ! MeV^-2
 real (kind=prec), parameter :: GF = 1._prec
 real (kind=prec), parameter :: alpha = 1._prec

 real (kind=prec), parameter :: Mmu = 105.658372_prec   ! MeV
 real (kind=prec), parameter :: Mel = 0.51099893_prec  ! MeV
! real (kind=prec), parameter :: Mel = 10._prec  ! MeV
 real (kind=prec), parameter :: Mtau = 1776.82_prec   ! MeV

 real (kind=prec), parameter :: xi_sep = 1.0E-10_prec
 real (kind=prec), parameter :: del_sep = 1.0E-10_prec
! character (len=3), parameter :: cgamma = "exp"
 character (len=3), parameter :: cgamma = "gam"
 
 integer print_ok, throw_away



        !!  ---------
	!!  variables
        !!  ---------

 integer :: ran_seed = 1

 real (kind=prec) :: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),   &
             p8(4),p9(4), pol1(4)
 real (kind=prec) :: mu, musq, delcut, xinormcut 
 real (kind=prec) :: xinormcut1, xinormcut2
 character (len=8) :: flavour

 real (kind=prec) :: Mm ! MeV
 real (kind=prec) :: Me ! MeV

contains

  SUBROUTINE CRASH(function_name)

  character(len=*) :: function_name

  write(6,*) "Program crashes because of a call to the function ", &
            function_name
  stop

  END SUBROUTINE CRASH


 END MODULE GLOBAL_DEF











