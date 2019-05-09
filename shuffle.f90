
MODULE shuffle
  use globals
  use utils
  implicit none

CONTAINS
  
  FUNCTION append_to_each_row(a, m) result(res)
    ! appends element a to each row of m
    complex(kind=prec) :: a, m(:,:)
    integer :: i
    complex(kind=prec) :: res(size(m,1),size(m,2)+1)
    do i=1,size(m,1)
      res(i,:) = [a,m(i,:)]
    end do
  END FUNCTION append_to_each_row

  FUNCTION stack_matrices_vertically(m1, m2) result(res)
    ! appends to matrix m1 the rows of matrix m2
    complex(kind=prec) :: m1(:,:), m2(:,:)
    complex(kind=prec) :: res(size(m1,1)+size(m2,1), size(m1,2))
    res(1:size(m1,1), :) = m1
    res(size(m1,1)+1:size(res,1),:) = m2 
  END FUNCTION stack_matrices_vertically

  RECURSIVE FUNCTION shuffle_product(v1, v2) result(res)
    complex(kind=prec) :: v1(:), v2(:)
    integer :: i
    complex(kind=prec) :: res(product((/(i,i=1,size(v1)+size(v2))/))/  & 
      (product((/(i,i=1,size(v1))/))*product((/(i,i=1,size(v2))/))), & 
      size(v1) + size(v2))
    complex(kind=prec) :: p1(product((/(i,i=1,size(v1)+size(v2)-1)/))/  & 
      (product((/(i,i=1,size(v1)-1)/))*product((/(i,i=1,size(v2))/))), & 
      size(v1) + size(v2) - 1)
    complex(kind=prec) :: p2(product((/(i,i=1,size(v1)+size(v2)-1)/))/  & 
      (product((/(i,i=1,size(v1))/))*product((/(i,i=1,size(v2)-1)/))), & 
      size(v1) + size(v2) - 1)
    complex(kind=prec) :: alpha, beta, w1(size(v1)-1), w2(size(v2)-1)

    res = 0
    if(size(v1) == 0) then 
      res(1,:) = v2
      return
    else if(size(v2) == 0) then 
      res(1,:) = v1
      return
    end if

    alpha = v1(1)
    beta = v2(1)
    w1 = v1(2:size(v1))
    w2 = v2(2:size(v2))

    res = stack_matrices_vertically( &
      append_to_each_row(alpha, shuffle_product(w1, v2)), & 
      append_to_each_row(beta, shuffle_product(v1, w2)) )

  END FUNCTION shuffle_product

END MODULE shuffle

PROGRAM test
  use utils
  use shuffle
  implicit none

  complex(kind=prec) :: v1(3), v2(2)
  integer :: amount_shuffles

  v1 = cmplx((/1,2,3/))
  v2 = cmplx((/4,5/))

  call print_matrix(shuffle_product(v1,v2))

END PROGRAM test

