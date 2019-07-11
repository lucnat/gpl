 
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

  FUNCTION MPL(m, x) result(res)
    integer :: m(:)
    complex(kind=prec) :: x(:)
    complex(kind=prec) :: res
    complex(kind=prec) :: t(size(x))
    integer :: q, j, k
    j = size(x)


    if(size(m) /= size(x)) then
      print*, 'Error: m and x must have the same length'
    end if
    res=0
    q=0
    t=0
    do while(.true.)
      res = t(1)
      q = q+1
      t(j) = t(j) + x(j)**q/q**m(j)
      do k=1,j-1
        t(j-k) = t(j-k) + t(j-k+1) * x(j-k)**(k+q)/(k+q)**m(j-k)
      enddo

      if (mod(q,2) .eq. 1) then
        if (abs(t(1)-res).lt.MPLdel) exit
      endif
    enddo
    res = t(1)


  END FUNCTION MPL

END MODULE mpl_module
