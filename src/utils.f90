

MODULE utils
  use globals
  use ieps
  implicit none

  ! logical :: print_enabled = .true.
  ! logical :: warnings_enabled = .true.

CONTAINS
  
  FUNCTION  get_condensed_m(z) result(m)
    ! returns condensed m where the ones not needed are filled with 0 (returns same size as z)
    type(inum), intent(in) :: z(:)
    integer :: m(size(z)), pos, i 
    m = 1
    pos = 1
    do i = 1,size(z)
      if(abs(z(i)) < zero) then
        if(i == size(z)) then
          pos = pos + 1
        else 
          m(pos) = m(pos) + 1
        end if
      else 
        pos = pos + 1
      end if
    end do
    m(pos:) = 0
  END FUNCTION get_condensed_m

  FUNCTION get_condensed_z(m, z_in) result(z_out)
    ! returns condensed z vector
    integer :: m(:), i, pos
    type(inum) :: z_in(:), z_out(size(m))
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(i) = z_in(pos)
    end do
  END FUNCTION get_condensed_z

  FUNCTION  get_flattened_z(m,z_in) result(z_out)
    ! returns flattened version of z based on m and z
    integer :: m(:), i, pos
    type(inum) :: z_in(:), z_out(sum(m))
    z_out = izero
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(pos) = z_in(i)
    end do
  END FUNCTION get_flattened_z

  FUNCTION find_amount_trailing_zeros(z) result(res)
    type(inum) :: z(:)
    integer :: res, i
    res = 0
    do i = size(z), 1, -1
      if( abs(z(i)) < zero ) then
        res = res + 1
      else
        exit
      end if
    end do
  END FUNCTION find_amount_trailing_zeros

  FUNCTION find_first_zero(v) result(res)
    ! returns index of first zero, or -1 if there is no zero
    integer :: v(:), i, res
    res = -1
    do i = 1,size(v)
      if(v(i) == 0) then
        res = i
        return
      end if
    end do
  END FUNCTION find_first_zero

  FUNCTION find_first_true(v) result(res)
    ! returns index of first element in v that is true
    logical :: v(:)
    integer :: i, res
    do i = 1, size(v)
      if(v(i)) then
        res = i
        return
      end if
    end do
  END FUNCTION find_first_true

  FUNCTION min_index(v)
    ! returns the index of the smallest element in v
    real(kind=prec) :: v(:), minimum
    integer :: min_index, i
    min_index = 1
    minimum = 1e15
    do i = 1,size(v)
      if(v(i) < minimum .and. v(i) > zero) then
        minimum = v(i)
        min_index = i
      end if
    end do
  END FUNCTION min_index

  FUNCTION zeroes(n) result(res)
    integer :: n
    type(inum) :: res(n)
    res = izero
  END FUNCTION zeroes

  RECURSIVE FUNCTION factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res
    res = merge(1,n*factorial(n-1),n==0)
  END FUNCTION factorial

  FUNCTION add_ieps(x) result(res)
    ! adds small imaginary part to x
    complex(kind=prec) :: x, res
    res = x + (0.0,epsilon)
  END FUNCTION add_ieps

  FUNCTION sub_ieps(x) result(res)
    ! subtracts small imaginary part to x
    complex(kind=prec) :: x, res
    res = x - (0.0,epsilon)
  END FUNCTION sub_ieps

  FUNCTION shuffle_with_zero(a) result(res)
    ! rows of result are shuffles of a with 0
    type(inum) :: a(:)
    type(inum) :: res(size(a)+1,size(a)+1)
    integer :: i,j, N
    N = size(a)+1
    do i = 1,N
      ! i is the index of the row
      ! j is the index of the zero
      j  = N+1-i
      res(i,j) = izero
      res(i,1:j-1) = a(1:j-1)
      res(i,j+1:N) = a(j:)
    end do
  END FUNCTION shuffle_with_zero

  SUBROUTINE print_matrix(m) 
    complex(kind=prec) :: m(:,:)
    integer :: s(2), i
    s = shape(m)
    do i = 1,s(1)
      print*, abs(m(i,:))
    end do
  END SUBROUTINE print_matrix

  SUBROUTINE print_logical_matrix(m) 
    logical :: m(:,:)
    integer :: s(2), i
    s = shape(m)
    do i = 1,s(1)
      print*, m(i,:)
    end do
  END SUBROUTINE print_logical_matrix

  ! subroutine print(s1,s2,s3,s4,s5)
  !   character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
  !   if(print_enabled) then
  !     print*, s1, s2, s3, s4, s5
  !   end if
  ! end subroutine print

  ! subroutine warn(s1,s2,s3,s4,s5)
  !   character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
  !   if(warnings_enabled) then
  !     print*, 'Warning: ', s1, s2, s3, s4, s5
  !   end if
  ! end subroutine warn

END MODULE utils


! PROGRAM test
!   use globals
!   use  utils
!   implicit none
  
!   complex(kind=prec) :: z(4)
!   integer :: m_prime(4), condensed_size
!   z = cmplx((/0.0,1.7,0.0,0.0/))

!   ! transform to condensed notation
!   m_prime = get_condensed_m(z)
!   print*, abs(z)
!   m_prime = get_condensed_m(z)
!   print*, abs(z)
!   if(find_first_zero(m_prime) == -1) then
!     condensed_size = size(m_prime)
!   else
!     condensed_size = find_first_zero(m_prime)-1 
!   end if
!   print*, condensed_size


! END  PROGRAM test



! Contains some functions that might be useful later

! Write your own print function with ability to suppress print
! Muss immer alle prints und warnings ausschalten können
! Test Programm schreiben mit exit codes -> gfortran 'test.f90' und dann 'echo $?'
