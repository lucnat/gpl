 
MODULE mpl_module
  use globals
  use utils
  use maths_functions
  implicit none

CONTAINS 

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

    res = 0


    if(size(m) == 1) then
      ! base case, we get a polylog
      do i = 1, n
    if(i**m(1) < 0) return ! roll over
    if(abs(x(1)**i) < 1.e-250) return
        res = res + x(1)**i / i**m(1)
      end do
    else 
      ! recursion step
      do i = 1, n    
        res = res + x(1)**i / i**m(1) * MPL(m(2:), x(2:), i - 1)
      end do
      
    end if
  END FUNCTION MPL

END MODULE mpl_module
