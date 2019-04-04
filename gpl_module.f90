
MODULE gpl_module
  use mpl_module

  implicit none

CONTAINS 

  FUNCTION GPL(m,z,y,k)
    ! computes the generalized polylogarithm G_{m1,..mk} (z1,...zk; y)

    integer :: m(:), k, i
    complex(kind=prec) :: z(:), x(k), y, GPL

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

