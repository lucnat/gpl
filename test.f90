
PROGRAM TEST
  use functions
  implicit none
  complex(kind=prec) :: result, x(4)

  ! print*, dilog((0.8,0),20)     ! should be 1.07479 + 0i 
  ! print*, dilog((0.2,0.5),20)   ! should be 0.133909 + 0.537628i

  ! print*, polylog(2,(0.2,0.5),20)   ! should be 0.133909 + 0.537628i 
  ! result = polylog(5,(0.2d0,0.5d0),20)   ! should be 0.192872345 + 0.505898833i

  ! result = multiple_polylog((/ 5 /),(/ (0.2,0.5) /),10)
  ! print*, 'result = ', result
  ! result = multiple_polylog((/ 5, 5 /),(/ (0.8,0),(0.3,0.5) /), 20)

  result = multiple_polylog((/ 0, 3, 4, 5 /),cmplx((/ 0.3,0.8,0.3,0.5 /)))   ! 6.929162361968684E-7

  print*, 'result ', result
  print*, 'compare', cmplx(6.9291623623242096179E-7)
  
END PROGRAM TEST
