
PROGRAM eval
  use GPL
  implicit none

  complex(kind=prec) :: res
  character(len=20) line
  integer io
  io=0
  do while(io.eq.0)
    line(:)="*"
    read(5,fmt='(A)',iostat=io) line
    if (len(trim(line)).eq.0) exit
    if (io.eq.0) then
      res = parseline(line)
      if (abs(aimag(res)).lt.1e-15) then
        print*,real(res)
      elseif (aimag(res)<0) then
        print*,real(res),aimag(res),"I"
      elseif (aimag(res)>0) then
        print*,real(res),'+',aimag(res),"I"
      endif
    endif
  enddo

CONTAINS

  FUNCTION PARSELINE(line)
  implicit none
  character(len=*) line
  complex(kind=prec) :: parseline
  complex(kind=prec) :: input(20)
  real(kind=prec) :: buf(2)
  character c
  integer i, last, nin


  last=1
  nin = 1
  buf = 0.
  do i=1,len(trim(line))+1
    if (i>len(trim(line))) then
      c=' '
    else
      c = line(i:i)
    endif

    if (last==i) cycle
    select case(c)
      case(' ',',',';')
        read(line(last:i-1), *) buf(1)
        last=i+1
        input(nin) = cmplx(buf(1), buf(2))
        nin = nin + 1
        buf = 0.

      case('+','-')
        read(line(last:i-1), *) buf(1)
        last=i+1
      case('I','i','j')
        read(line(last:i-1), *) buf(2)
        last=i+1
        input(nin) = cmplx(buf(1), buf(2))
        nin = nin + 1
        buf = 0.
    endselect
  enddo

  parseline = G(input(:nin-1))

  END FUNCTION PARSELINE

END PROGRAM eval

