
MODULE gpl_module
  use mpl_module

  implicit none

CONTAINS 

  RECURSIVE FUNCTION factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res
    res = merge(1,n*factorial(n-1),n==1)
  END FUNCTION factorial

  FUNCTION GPL_has_convergent_series(m,z,y,k)
    ! tests if GPL has a convergent series representation
    integer :: m(:), k
    complex(kind=prec) :: z(:), y
    logical :: GPL_has_convergent_series

    GPL_has_convergent_series = .false.

    if(all(abs(y) <= abs(z))) then
      if(m(1) == 1) then 
        GPL_has_convergent_series = (y/z(1) /= 1)
      else 
        GPL_has_convergent_series = .true.
      end if
    end if

  END FUNCTION GPL_has_convergent_series

  FUNCTION GPL_zero_zi(l,y)
    ! used to compute the value of GPL when all zi are zero
    integer :: l
    complex(kind=prec) :: y, GPL_zero_zi
    GPL_zero_zi = 1.0d0/factorial(l) * log(y) ** l

  END FUNCTION GPL_zero_zi

  FUNCTION GPL(m,z,y,k)
    ! computes the generalized polylogarithm G_{m1,..mk} (z1,...zk; y)
    ! assumes zero arguments expressed through the m's

    integer :: m(:), k, i
    complex(kind=prec) :: z(:), x(k), y, GPL

    ! are all z_i = 0 ? 
    if(k == 1 .and. z(1) == 0) then
      ! for that we assume that only one argument was passed, the rest through m1
      GPL = GPL_zero_zi(m(1)-1,y)
      return
    end if

    ! do they have convergent series rep?
    if(.not. GPL_has_convergent_series(m,z,y,k)) then
      print*, '  ', 'does not have convergent series representation'
    end if

    do i = 1,k
      x(i) = merge(y/z(1), z(i-1)/z(i),i == 1)
    end do
    GPL = (-1)**k * MPL(m,x)

  END FUNCTION GPL

END MODULE gpl_module

! PROGRAM test
!   ! used to test this module
!   use gpl_module
  
!   integer :: m(2) = (/ 1,1 /)
!   complex(kind=prec) :: z(2) = dcmplx((/ 1.3d0, 1.1d0 /))
!   complex(kind=prec) :: y = 0.4
!   complex(kind=prec) :: res, ref
  
!   res = GPL(m,z,y,2)
!   ref = dcmplx(0.0819393734128676)
!   print*, 'res=',res
!   print*, 'ref=',ref

! END PROGRAM test

