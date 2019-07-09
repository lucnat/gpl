 
MODULE mpl_module
  use globals
  use utils
  use maths_functions
  implicit none

CONTAINS 

  FUNCTION polylog_series(m,x,n_passed) result(res)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    integer, optional :: n_passed
    complex(kind=prec) :: x, res
    integer :: i,n
    integer, allocatable :: j(:)
    n = merge(n_passed,GPLInfinity,present(n_passed))  
    j = (/(i, i=1,n,1)/) 
    res = sum(x**j / j**m)
  END FUNCTION polylog_series

  FUNCTION MPL_converges(m,x)
    ! checks if the MPL converges 
    complex(kind=prec) :: x(:)
    integer :: m(:)
    logical :: MPL_converges
    MPL_converges = .false.
    if(abs(product(x)) < 1) then
      if(m(1) /= 1 .or. abs(x(1) - 1) < zero) then
        MPL_converges = .true.
      end if
    end if
  END FUNCTION MPL_converges

  recursive FUNCTION MPL(m, x, n_passed) result(res)
    ! Computes the multiple polylogarithm Li_{m1,...,mk} (x1,...,xk) up to order n
    integer :: m(:)
    complex(kind=prec) :: x(:)
    complex(kind=prec) :: res
    integer, optional :: n_passed
    integer :: i, n

    n = merge(n_passed,GPLInfinity,present(n_passed))  

    if(size(m) /= size(x)) then
      print*, 'Error: m and x must have the same length'
    end if
    
    if(size(m) == 1) then
      ! base case
      res = polylog_series(m(1),x(1), n)
    else 
      ! recursion step
      res = 0
      do i = 1, n    
        res = res + x(1)**i / i**m(1) * MPL(m(2:), x(2:), i - 1)
      end do

      ! a nicer way to do it would be but problem is i
      ! i = (/(j, j=1,n, 1)/)
      ! res = sum( x(1)**i / i**m(1) * MPL(m(2:), x(2:), i(1) - 1) )
      ! we could get around this problem by rewriting MPL to operate on each i and returning an array
      
    end if
  END FUNCTION MPL

END MODULE mpl_module

! PROGRAM test
!   use mpl_module
!   logical :: result
!   result = MPL_converges( dcmplx((/1.0d0,.7d0,.3d0/)), (/ 1,2,1 /) )
!   print*, result
! end PROGRAM test
