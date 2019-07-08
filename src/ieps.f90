
MODULE ieps
  use globals
  implicit none
  type inum
    complex(kind=prec) :: c
    integer(8) :: i0
  end type inum

  integer(8), parameter :: di0 = +1

  type(inum), parameter :: izero=inum( 0.,di0)
  type(inum), parameter :: imone=inum(-1.,di0)
  type(inum), parameter :: ione=inum(+1.,di0)


  interface operator (*)
    module procedure multinumss, multinumvs
  end interface operator (*)
  interface operator (+)
    module procedure addinumss, addinumvs
  end interface operator (+)
  interface operator (-)
    module procedure subinumss,subinumvs,subinumsv
  end interface operator (-)
  interface operator (**)
    module procedure powinum
  end interface operator (**)
  interface operator (/)
    module procedure divint, divinumss, divinumvs
  end interface operator (/)
  interface abs
    module procedure absinum, absinumv
  end interface abs
  interface log
    module procedure loginum
  end interface log

  interface toinum
    module procedure toinum_cmplx, toinum_real, toinum_int
  end interface toinum
  interface tocmplx
    module procedure tocmplxv, tocmplxs
  end interface tocmplx
CONTAINS


  FUNCTION MULTINUMSS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1, n2
  type(inum) :: multinumss
  multinumss = inum( n1%c*n2%c, int(sign(1._prec,real(n1%c)*n2%i0 + real(n2%c)*n1%i0)) )
  END FUNCTION MULTINUMSS

  FUNCTION MULTINUMVS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1(:), n2
  type(inum) :: multinumvs(size(n1))
  integer i
  do i = 1,size(n1)
    multinumvs(i) = inum( n1(i)%c*n2%c, int(sign(1._prec,real(n1(i)%c)*n2%i0 + real(n2%c)*n1(i)%i0)) )
  enddo
  END FUNCTION MULTINUMVS

  FUNCTION ADDINUMSS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1, n2
  type(inum) :: addinumss
  !TODO: what *is* the sum?
  addinumss = inum(n1%c + n2%c, n1%i0 )
  END FUNCTION ADDINUMSS

  FUNCTION ADDINUMVS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1(:), n2
  type(inum) :: addinumvs(size(n1))
  !TODO: what *is* the sum?
  integer i
  do i = 1,size(n1)
    addinumvs(i) = inum(n1(i)%c + n2%c, n1(i)%i0 )
  enddo
  END FUNCTION ADDINUMVS

  FUNCTION SUBINUMSS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1, n2
  type(inum) :: subinumss
  !TODO: what *is* the sum?
  subinumss = inum(n1%c - n2%c, n1%i0 )
  END FUNCTION SUBINUMSS

  FUNCTION SUBINUMVS(n1, n2)
  implicit none
  type(inum), intent(in) :: n1(:), n2
  type(inum) :: subinumvs(size(n1))
  !TODO: what *is* the sum?
  integer i
  do i = 1,size(n1)
    subinumvs(i) = inum(n1(i)%c - n2%c, n1(i)%i0 )
  enddo
  END FUNCTION SUBINUMvs
  FUNCTION SUBINUMSV(n2, n1)
  implicit none
  type(inum), intent(in) :: n1(:), n2
  type(inum) :: subinumsv(size(n1))
  !TODO: what *is* the sum?
  integer i
  do i = 1,size(n1)
    subinumsv(i) = inum(n2%c - n1(i)%c, n1(i)%i0 )
  enddo
  END FUNCTION SUBINUMSV

  FUNCTION ABSINUM(n1)
  implicit none
  type(inum), intent(in) :: n1
  real(kind=prec) :: absinum
  absinum = sqrt(real(n1%c)**2+aimag(n1%c)**2)
  END FUNCTION ABSINUM

  FUNCTION ABSINUMV(n1)
  implicit none
  type(inum), intent(in) :: n1(:)
  real(kind=prec) :: absinumv(size(n1))
  absinumv = abs(n1%c)
  END FUNCTION ABSINUMV


  FUNCTION POWINUM(n1, m)
  implicit none
  type(inum), intent(in) :: n1
  integer, intent(in) :: m
  type(inum) :: powinum
  if (aimag(n1%c)<zero) then
    powinum = inum( cmplx(real(n1%c)**m,0.), int(sign(1._prec,real(n1%c)**m)) )
  else
    powinum = inum( n1%c**m, n1%i0 )
  endif
  END FUNCTION POWINUM

  FUNCTION DIVINT(n1, m)
  implicit none
  type(inum), intent(in) :: n1
  integer, intent(in) :: m
  type(inum) :: divint
  divint = inum( n1%c/m, n1%i0*sign(1,m))
  END FUNCTION DIVINT

  FUNCTION DIVINUMss(n1, n2)
  implicit none
  type(inum), intent(in) :: n1, n2
  type(inum) :: divinumss
  divinumss = inum( n1%c/n2%c, int(sign(1., real(n2%c)*n1%i0 - real(n1%c)*n2%i0)))
  END FUNCTION DIVINUMss

  FUNCTION DIVINUMvs(n1, n2)
  implicit none
  type(inum), intent(in) :: n1(:), n2
  type(inum) :: divinumvs(size(n1))
  integer i
  do i = 1,size(n1)
    divinumvs(i) = inum( n1(i)%c/n2%c, int(sign(1., real(n2%c)*n1(i)%i0 - real(n1(i)%c)*n2%i0)))
  enddo
  END FUNCTION DIVINUMvs

  FUNCTION LOGINUM(n1)
  implicit none
  type(inum), intent(in) :: n1
  complex(kind=prec) :: loginum
  loginum = log(n1%c)
  END FUNCTION LOGINUM


  FUNCTION TOINUM_cmplx(z, s)
  complex(kind=prec) :: z(:)
  type(inum) :: toinum_cmplx(size(z))
  integer(8),optional :: s
  integer(8) :: ss
  integer i
  if (present(s)) then
    ss = s
  else
    ss = di0
  endif
  do i=1,size(z)
    toinum_cmplx(i) = inum(z(i), ss)
  enddo
  END FUNCTION TOINUM_cmplx

  FUNCTION TOINUM_real(z, s)
  real(kind=prec) :: z(:)
  type(inum) :: toinum_real(size(z))
  integer(8),optional :: s
  integer(8) :: ss
  integer i
  if (present(s)) then
    ss = s
  else
    ss = di0
  endif
  do i=1,size(z)
    toinum_real(i) = inum(z(i), ss)
  enddo
  END FUNCTION TOINUM_real


  FUNCTION TOINUM_int(z, s)
  integer :: z(:)
  type(inum) :: toinum_int(size(z))
  integer(8),optional :: s
  integer(8) :: ss
  integer i
  if (present(s)) then
    ss = s
  else
    ss = di0
  endif
  do i=1,size(z)
    toinum_int(i) = inum(z(i), ss)
  enddo
  END FUNCTION TOINUM_int

  FUNCTION TOCMPLXv(z)
  type(inum) :: z(:)
  complex(kind=prec) tocmplxv(size(z))
  tocmplxv = z%c
  END FUNCTION
  FUNCTION TOCMPLXs(z)
  type(inum) :: z
  complex(kind=prec) tocmplxs
  tocmplxs = z%c
  END FUNCTION


END MODULE IEPS
