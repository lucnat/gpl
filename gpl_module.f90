 
MODULE gpl_module
  use mpl_module

  implicit none

CONTAINS 

  FUNCTION GPL()
    ! integer :: m(:)
    ! complex(kind=prec) :: y(:)
    complex(kind=prec) GPL

    print*, 'hello, GPL'
  END FUNCTION GPL

END MODULE gpl_module
