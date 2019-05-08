
MODULE shuffle
  use globals
  implicit none

CONTAINS
  
  FUNCTION append_to_each_row(a, m) result(res)
    ! appends element a to each row of m
    integer :: a, m(:,:), i
    integer :: res(size(m,1),size(m,2)+1)
    do i=1,size(m,1)
      res(i,:) = [a,m(i,:)]
    end do
  END FUNCTION append_to_each_row

  FUNCTION stack_matrices_vertically(m1, m2) result(res)
    ! appends to matrix m1 the rows of matrix m2
    integer :: m1(:,:), m2(:,:)
    integer :: res(size(m1,1)+size(m2,1), size(m1,2))
    res(1:size(m1,1), :) = m1
    res(size(m1,1)+1:size(res,1),:) = m2 
  END FUNCTION stack_matrices_vertically
  
  RECURSIVE FUNCTION factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res
    res = merge(1,n*factorial(n-1),n==0)
  END FUNCTION factorial

  RECURSIVE FUNCTION shuffle_product(v1, v2) result(res)
    integer :: v1(:), v2(:)
    integer :: i
    integer :: res(product((/(i,i=1,size(v1)+size(v2))/))/  & 
      (product((/(i,i=1,size(v1))/))*product((/(i,i=1,size(v2))/))), & 
      size(v1) + size(v2))
    integer :: p1(product((/(i,i=1,size(v1)+size(v2)-1)/))/  & 
      (product((/(i,i=1,size(v1)-1)/))*product((/(i,i=1,size(v2))/))), & 
      size(v1) + size(v2) - 1)
    integer :: p2(product((/(i,i=1,size(v1)+size(v2)-1)/))/  & 
      (product((/(i,i=1,size(v1))/))*product((/(i,i=1,size(v2)-1)/))), & 
      size(v1) + size(v2) - 1)
    integer :: alpha, beta, w1(size(v1)-1), w2(size(v2)-1)

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

! PROGRAM test
!   use shuffle_algebra
!   implicit none

!   integer :: v1(3), v2(3)
!   integer :: amount_shuffles
!   integer :: res(3,3)

!   v1 = (/1,2,3/)
!   v2 = (/-1,-2,-3/)

!   call print_as_matrix(shuffle_product(v1, v2))

! END PROGRAM test

