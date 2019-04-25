
! This is currently a stand alone program which will merely be used as a 
! guide for the implementation of the shuffle algebra for GPL functions

PROGRAM shuffle_algebra
  implicit none

  ! Currently words can be no longer than the following values
  ! Might need to be adjusted

  integer, parameter :: max_word_size = 6
  integer, parameter :: max_word_sum_size = 1000

  type word
    character(len=max_word_size) :: letters
    integer :: length
  endtype word

  type word_sum
    type(word), dimension(max_word_sum_size) :: words
    integer :: length
  end type word_sum

  type(word) :: v1 = word("abc",3)
  type(word) :: v2 = word("123",3)
  type(word) :: w1
  type(word_sum) :: ws, ws1, ws2

  ws = shuffle_product(v1,v2)
  call print_word_sum(ws)
CONTAINS
     
  RECURSIVE FUNCTION shuffle_product(v1, v2) result(res)
    ! takes two words and returns shuffle product as a word sum 
    type(word_sum) :: res, p1, p2, s1, s2
    type(word) :: v1,v2,w1,w2
    type(character) :: alpha, beta

    ! print*, '----------------------'
    ! print*, 'v1 = ', v1
    ! print*, 'v2 = ', v2

    if(v1%length == 0) then
      res = word_sum((/ v2 /),1)
    else if(v2%length == 0) then
      res = word_sum((/ v1 /),1)
    else
      alpha = v1%letters(1:1)
      beta = v2%letters(1:1)
      w1 = word(v1%letters(2:),v1%length-1)
      w2 = word(v2%letters(2:),v2%length-1)
      p1 = shuffle_product(w1,v2)
      p2 = shuffle_product(v1,w2)
      s1 = times(alpha, p1)
      s2 = times(beta, p2)

      res = combined_word_sum(s1,s2)

    end if
  END FUNCTION shuffle_product

  FUNCTION combined_word_sum(s1, s2)
    ! combines word sums s1 and s2 into s
    type(word_sum) :: s1, s2, combined_word_sum
    type(word), dimension(s1%length + s2%length) :: combined_words
    integer :: length
    length = s1%length + s2%length
    combined_words(1:s1%length) = s1%words(1:s1%length)
    combined_words(s1%length+1:length) = s2%words(1:s2%length)
    combined_word_sum = word_sum(combined_words,length)
  END FUNCTION combined_word_sum

  FUNCTION times(l, ws) 
    ! computes word sum from letter times word sum, e.g. a(bc + cd) = abc + acd
    character :: l
    type(word_sum) :: ws
    type(word_sum) :: times
    integer :: i

    times%length = ws%length 

    do i=1,ws%length
      times%words(i)%letters = l // ws%words(i)%letters
    end do
  END FUNCTION times

  SUBROUTINE print_word_sum(ws)
    ! prints a sum of word with plusses for easy readibility
    integer :: i
    type(word_sum) :: ws
    do i = 1,ws%length
      write(*, fmt="(1xai0)", advance="no") ws%words(i)%letters
      if(i /= ws%length) then
        write(*, fmt="(1xai0)", advance="no") " + "
      end if
    end do
  END SUBROUTINE print_word_sum

END PROGRAM shuffle_algebra

