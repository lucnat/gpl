
MODULE functions
  implicit none

CONTAINS 

  FUNCTION dilog(x,n)
    ! Computes the dilog Li_2(x) using the series representation up to order n
    integer :: n
    complex :: x, dilog
    integer :: i
    integer :: j(n)

    j = (/(i, i=1,n,1)/) 
    dilog = sum(x**j / j**2)
  END FUNCTION dilog

  FUNCTION polylog(m,x,n)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m,n
    complex :: x, polylog
    integer :: i
    integer :: j(n)

    j = (/(i, i=1,n,1)/) 
    polylog = sum(x**j / j**m)
  END FUNCTION polylog

  recursive FUNCTION multiple_polylog(m, x, n) result(res)
    ! Computes the multiple polylogarithm Li_{m1,...,mk} (x1,...,xk) up to order n
    integer :: m(:)
    complex :: x(:)
    complex :: res
    integer :: n, i

    ! print*, 'm = ', m
    ! print*, 'x = ', x
    ! print*, 'n = ', n

    if(size(m) /= size(x)) then
      print*, 'Error: m and x must have the same length'
    end if
    
    if(size(m) == 1) then
      ! base case
      ! print*, 'landed in base case'
      res = polylog(m(1),x(1), n)
    else 
      ! recursion step
      ! print*, 'landed in step'
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


PROGRAM test
  ! Used to test the procedures defined in this module
  use functions
  use utils
  implicit none
  complex :: result

  ! print*, dilog((0.8,0),100)     ! should be 1.07479 + 0i 
  ! print*, dilog((0.2,0.5),100)   ! should be 0.133909 + 0.537628i

  ! print*, polylog(2,(0.2,0.5),100)   ! should be 0.133909 + 0.537628i 
  ! print*, polylog(5,(0.2,0.5),100)   ! should be 0.192872345 + 0.505898833i

  ! result = multiple_polylog((/ 5 /),(/ (0.2,0.5) /), 100)
  ! print*, 'result = ', result
  ! result = multiple_polylog((/ 5, 5 /),(/ (0.8,0),(0.3,0.5) /), 100)
  ! print*, 'result = ', result

END PROGRAM test
