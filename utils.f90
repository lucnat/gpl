
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

CONTAINS
  
  subroutine print(s1,s2,s3,s4,s5)
    character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
    if(print_enabled) then
      print*, s1, s2, s3, s4, s5
    end if
  end subroutine print

  subroutine warn(s1,s2,s3,s4,s5)
    character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
    if(warnings_enabled) then
      print*, 'Warning: ', s1, s2, s3, s4, s5
    end if
  end subroutine warn

END MODULE utils

