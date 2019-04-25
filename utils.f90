
! Contains some functions that might be useful later

! Write your own print function with ability to suppress print
! Muss immer alle prints und warnings ausschalten können
! Test Programm schreiben mit exit codes -> gfortran 'test.f90' und dann 'echo $?'
! Define GPL infinity
! Mach n optional
! Kommentar schreiben zu anderer Notation 
! Funktion überprüfen! Tests schreiben!

MODULE utils
  implicit none

  logical :: print_enabled = .true.
  logical :: warnings_enabled = .true.
  integer, parameter :: prec = selected_real_kind(15,32)  

CONTAINS

  FUNCTION  get_condensed_m(z) result(m)
    ! returns condensed m where the ones not needed are filled with 0
    complex(kind=prec) :: z(:), m(size(z))
    integer :: pos = 1, i 
    m = 1
    do i = 1,size(z)
      if(z(i) == 0) then
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

  FUNCTION get_condenced_z(m, z_in) result(z_out)
    ! returns condensed z vector
    integer :: m(:), i, pos
    complex(kind=prec) :: z_in(:), z_out(size(m)) 
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(i) = z_in(pos)
    end do
  END FUNCTION get_condenced_z

  FUNCTION  get_flattened_z(m,z_in) result(z_out)
    ! returns flattened version of z based on m and z
    integer :: m(:), i, pos
    complex(kind=prec) :: z_in(:), z_out(sum(m))
    z_out = 0
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(pos) = z_in(i)
    end do
  END FUNCTION get_flattened_z

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

  SUBROUTINE print_as_matrix(m) 
    ! prints a 2d array as a matrix
    complex :: m(:,:)
    integer :: s(2), i
    s = shape(m)
    do i = 1,s(1)
      print*, m(i,:)
    end do
  END SUBROUTINE print_as_matrix

  FUNCTION shuffle_with_zero(a) result(res)
    ! rows of result are shuffles of a with 0
    complex :: a(:)
    complex :: res(size(a)+1,size(a)+1)
    integer :: i,j, N
    N = size(a)+1
    do i = 1,N
      ! i is the index of the row
      ! j is the index of the zero
      j  = N+1-i
      res(i,j) = 0
      res(i,1:j-1) = a(1:j-1)
      res(i,j+1:N) = a(j:)
    end do
  END FUNCTION shuffle_with_zero

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
!   use  utils
!   implicit none

!   complex(kind=prec), dimension(3) :: a = cmplx((/1,2,3/))
!   complex(kind=prec) :: z_flat(7)
!   complex(kind=prec), allocatable :: z(:)
!   integer :: m_prime(7), condensed_size
!   integer, allocatable :: m(:)
!   complex(kind=prec) :: b(size(a)+1,size(a)+1)

!   ! ! test shuffling
!   ! b = 1
!   ! b = shuffle_with_zero(a)
!   ! call print_as_matrix(b)

!   ! ! test condensing
!   ! z_flat = cmplx((/0,0,0,2,1,0,1/))
!   ! m_prime = get_condensed_m(z_flat)
!   ! condensed_size = find_first_zero(m_prime)-1 
!   ! allocate(m(condensed_size))
!   ! allocate(z(condensed_size))
!   ! m = m_prime(1:condensed_size)
!   ! z = get_condenced_z(m,z_flat)
!   ! z_flat = get_flattened_z(m,z)
!   ! print*, z_flat
!   ! print*, m
!   ! deallocate(m)
!   ! deallocate(z)

! END  PROGRAM test

