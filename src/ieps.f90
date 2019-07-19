
MODULE ieps
  use globals
  implicit none
  type inum
    complex(kind=prec) :: c
    integer(1) :: i0
  end type inum

  integer(1), parameter :: di0 = +1

  type(inum), parameter :: izero=inum( 0.,di0)
  type(inum), parameter :: marker=inum(0.,5)


  interface abs
    module procedure absinum, absinumv
  end interface abs

  interface toinum
    module procedure toinum_cmplx, toinum_real, toinum_reals, toinum_int
  end interface toinum
  interface tocmplx
    module procedure tocmplxv, tocmplxs
  end interface tocmplx
  interface real
    module procedure realis, realiv
  end interface real
  interface aimag
    module procedure imags, imagv
  end interface aimag
CONTAINS
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

  FUNCTION TOINUM_cmplx(z, s)
  complex(kind=prec) :: z(:)
  type(inum) :: toinum_cmplx(size(z))
  integer(1),optional :: s
  integer(1) :: ss
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
  integer(1),optional :: s
  integer(1) :: ss
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
  integer(1),optional :: s
  integer(1) :: ss
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

  FUNCTION TOINUM_reals(z, s)
  real(kind=prec) :: z
  type(inum) :: toinum_reals
  integer(1),optional :: s
  integer(1) :: ss
  if (present(s)) then
    ss = s
  else
    ss = di0
  endif
  toinum_reals = inum(z, ss)
  END FUNCTION TOINUM_reals


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


  FUNCTION REALIV(z)
  type(inum) :: z(:)
  real(kind=prec) realiv(size(z))
  realiv = real(z%c)
  END FUNCTION
  FUNCTION REALIS(z)
  type(inum) :: z
  real(kind=prec) realis
  realis = real(z%c)
  END FUNCTION

  FUNCTION IMAGV(z)
  type(inum) :: z(:)
  real(kind=prec) imagv(size(z))
  imagv = aimag(z%c)
  END FUNCTION
  FUNCTION IMAGS(z)
  type(inum) :: z
  real(kind=prec) imags
  imags = aimag(z%c)
  END FUNCTION

END MODULE IEPS
