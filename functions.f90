MODULE functions
  implicit none

  integer, parameter :: prec = selected_real_kind(15,32)  
  integer :: GPLInfinity = 30   ! the default n if it is not passed

CONTAINS 

  FUNCTION dilog(x,n)
    ! Computes the dilog Li_2(x) using the series representation up to order n
    integer :: n
    complex(kind=prec) :: x, dilog
    integer :: i
    integer :: j(n)

    j = (/(i, i=1,n,1)/) 
    dilog = sum(x**j / j**2)
  END FUNCTION dilog

  FUNCTION polylog(m,x,n)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    complex(kind=prec) :: x, polylog
    integer :: i,n
    integer, allocatable :: j(:)

    j = (/(i, i=1,n,1)/) 
    polylog = sum(x**j / j**m)
    
  END FUNCTION polylog

  recursive FUNCTION multiple_polylog(m, x, n_passed) result(res)
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
      res = polylog(m(1),x(1), n)
    else 
      ! recursion step
      res = 0
      do i = 1, n    
        res = res + x(1)**i / i**m(1) * multiple_polylog(m(2:), x(2:), i - 1)
      end do

      ! a nicer way to do it would be but problem is i
      ! i = (/(j, j=1,n, 1)/)
      ! res = sum( x(1)**i / i**m(1) * multiple_polylog(m(2:), x(2:), i(1) - 1) )
      
    end if
  END FUNCTION multiple_polylog

END MODULE functions

