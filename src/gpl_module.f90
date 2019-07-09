
MODULE gpl_module
  use globals
  use utils
  use maths_functions
  use shuffle
  use mpl_module
  implicit none

  INTERFACE GPL
    module procedure G_flat, G_condensed, G_superflat, G_real, G_int
  END INTERFACE GPL

CONTAINS 

  FUNCTION GPL_has_convergent_series(m,z,y)
    ! tests if GPL has a convergent series representation
    integer :: m(:)
    complex(kind=prec) :: z(:), y
    logical :: GPL_has_convergent_series

    GPL_has_convergent_series = .false.

    if(all(abs(y) <= abs(z))) then
      if(m(1) == 1) then 
        GPL_has_convergent_series = .true. !(abs((y/z(1) - 1)) < zero)
      else 
        GPL_has_convergent_series = .true.
      end if
    end if
  END FUNCTION GPL_has_convergent_series

  FUNCTION GPL_zero_zi(l,y)
    ! used to compute the value of GPL when all zi are zero
    integer :: l
    complex(kind=prec) :: y, GPL_zero_zi
    GPL_zero_zi = 1.0_prec/factorial(l) * log(y) ** l
  END FUNCTION GPL_zero_zi

  FUNCTION is_convergent(z,y)
    ! returns true if G(z,y) convergent, otherwise false
    ! can be used in either notation (flat or condensed)
    complex(kind=prec) :: z(:), y
    logical :: is_convergent
    integer :: i

    is_convergent = .true.
    do i = 1,size(z)
      if(abs(z(i)) < zero) cycle  ! skip zero values
      if(abs(y) > abs(z(i))) is_convergent = .false.
    end do
  END FUNCTION is_convergent

  SUBROUTINE print_G(z_flat, y)
    complex(kind=prec) :: z_flat(:)
    complex(kind=prec), optional :: y
    if(present(y)) print*, 'G(', real(z_flat), real(y), ')'
    if(.not. present(y)) print*, 'G(', abs(z_flat), ')'
  END SUBROUTINE print_G

  RECURSIVE FUNCTION remove_sr_from_last_place_in_PI(a,y2,m,p) result(res)
    ! here what is passed is not the full a vector, only a1, ..., ak without the trailing zeroes
    integer :: m, i, j, n
    complex(kind=prec) :: a(:), y2, s(m), p(:), res
    complex(kind=prec) :: alpha(product((/(i,i=1,size(a)+size(s))/))/  & 
        (product((/(i,i=1,size(a))/))*product((/(i,i=1,size(s))/))), & 
        size(a) + size(s))

    s = [zeroes(m-1),cmplx(42e50)]
    alpha = shuffle_product(a,s)
    if(verb >= 50) then
      print*, 'mapping to '
      call print_G(a,y2)
      print*, 'PI with p=',real(p),'i=',m,'g =',real([zeroes(m-1),y2])
    end if
    res = GPL(a,y2)*pending_integral(p,m,[zeroes(m-1),y2])
    if(verb >= 50) print*, 'also mapping to'
    do j = 2,size(alpha, 1)
      ! find location of s_r
      n = find_first_true(abs(alpha(j,:) - 42e50) < zero)
      if(verb >= 50) print*, 'PI with p=',real(p),'i=',n,'g =',real([alpha(j,1:n-1),alpha(j,n+1:size(alpha,2)),y2])
      res = res - pending_integral(p, n, [alpha(j,1:n-1),alpha(j,n+1:size(alpha,2)),y2])
    end do
  END FUNCTION remove_sr_from_last_place_in_PI

  RECURSIVE FUNCTION pending_integral(p,i,g) result(res)
    ! evaluates a pending integral by reducing it to simpler ones and g functions
    complex(kind=prec) :: p(:), g(:), res
    complex(kind=prec) :: y1, y2, b(size(p)-1), a(size(g)-1)
    integer :: i, m
    res = 0

    if(verb >= 30) print*, 'evaluating PI with p=',real(p),'i=',real(i),'g =',real(g)   

    y1 = p(1)
    b = p(2:size(p))

    ! if integration variable is not in G-function
    if(i == 0 .or. size(g) == 0) then
      if(verb >= 30) print*, 'only integrals in front, we get G-function'
      res = G_flat(b,y1)
      return
    end if

    ! if integration variable at end -> we gat a G function 
    if(i == size(g)+1) then
      if(verb >= 30) print*, 'is just a G-function'
      res = G_flat([p(2:size(p)),g], p(1))
      return
    end if
  

    ! if depth one and m = 1 use my (59)
    if(size(g) == 1) then
      if(verb >= 30) print*, 'case depth one and m = 1'
      res = pending_integral(p,2,[sub_ieps(g(1))]) - pending_integral(p,2,[cmplx(0.0)]) &
        + G_flat(p(2:size(p)), p(1)) * log(-sub_ieps(g(1)))
      return
    end if
  
    a = g(1:size(g)-1)
    y2 = g(size(g)) 
    m = size(g)  ! actually, size(g)-1+1

    ! if depth one and m > 1 use my (60)
    if(all( abs( g(1:size(g)-1) ) < zero)) then       
      if(verb >= 30) print*, 'case depth one and m > 1'
      if(verb >= 50) then 
        print*, 'map to'
        print*, 'zeta(',m,')'
        print*, 'PI with p=',real(p),'i=',0,'g =',zeroes(0)
        print*, 'PI with p=',real([y2,cmplx(0.0)]),'i=',m-1,'g =',[zeroes(m-2),y2]
        print*, 'PI with p=',real(p),'i=',0,'g =',zeroes(0)
        print*, 'PI with p=',[p,cmplx(0.0)],'i=',m-1,'g =',zeroes(0)
      end if
      res = -zeta(m)*pending_integral(p,0,zeroes(0)) &
        + pending_integral([y2,cmplx(0.0)],m-1,[zeroes(m-2),y2])*pending_integral(p,0,zeroes(0)) &
        - pending_integral([p,cmplx(0.0)],m-1,[zeroes(m-2),y2])
      return
    end if
  
    ! case of higher depth, s_r at beginning, use my (68)
    if(i == 1) then
      if(verb >= 30) print*, 'case higher depth, sr at beginning'
      
      if(verb >= 50) then
        print*, 'map to (using 68)'
        print*, 'PI with p=',real(p),'i=',0,'g =',zeroes(0) 
        call print_G([cmplx(0.0),a],y2)
        print*, 'PI with p=',real([p,y2]),'i=',0,'g =',zeroes(0) 
        call print_G(a,y2)
        print*, 'PI with p=',real([p,a(1)]),'i=',1,'g =',real([a(2:size(a)),y2])
        print*, 'PI with p=',real([p,a(1)]),'i=',0,'g =',zeroes(0)
        call print_G(a,y2)
      end if

      res = pending_integral(p,0,zeroes(0)) * G_flat([cmplx(0.0),a],y2) & 
        + pending_integral([p,y2],0,zeroes(0)) * G_flat(a,y2) &
        + pending_integral([p,a(1)],1,[a(2:size(a)),y2]) &
        - pending_integral([p,a(1)],0,zeroes(0)) * G_flat(a,y2)
        return
    end if
    
    ! case higher depth, s_r at the end, use (62)
    if(i == size(g)) then
      if(verb >= 30) print*, 's_r at the end under PI, need to shuffle'
      m = find_amount_trailing_zeros(a) + 1
      res = remove_sr_from_last_place_in_PI(a(1:size(a)-(m-1)), y2, m, p)
      return
    end if
    
    ! case higher depth, s_r in middle, use my (67)
    if(verb >= 30) print*, 's_r in the middle under PI'

    res =  -pending_integral(p,1,zeroes(0)) * G_flat([a(1:i-1),cmplx(0.0),a(i+1:size(a))],y2) &
      + pending_integral([p,a(i-1)],i-1,[a(1:i-2),a(i:size(a)),y2]) &
      + pending_integral([p,a(i-1)],1,zeroes(0)) * G_flat(a,y2) &
      + pending_integral([p,a(i)], i, [a(1:i-1), a(i+1:size(a)),y2]) & 
      - pending_integral([p,a(i)],1,zeroes(0)) * G_flat(a,y2)
  END FUNCTION pending_integral

  RECURSIVE FUNCTION remove_sr_from_last_place_in_G(a,y2,m,sr) result(res)
    complex(kind=prec) :: a(:), sr, res,y2
    integer :: m,i,j
    complex(kind=prec) :: alpha(product((/(i,i=1,size(a)+m)/))/  & 
      (product((/(i,i=1,size(a))/))*product((/(i,i=1,m)/))), & 
      size(a) + m)
    alpha = shuffle_product(a,[zeroes(m-1),sr])
    res = G_flat(a,y2)*G_flat([zeroes(m-1),sr],y2)
    do j = 2,size(alpha,1)
      res = res - G_flat(alpha(j,:),y2)
    end do
  END FUNCTION remove_sr_from_last_place_in_G

  RECURSIVE FUNCTION make_convergent(a,y2) result(res)
    ! goes from G-functions to pending integrals and simpler expressions according to (62),(64),(67) and (68)

    complex(kind=prec) :: a(:), y2, res, sr
    integer :: i, mminus1

    res = 0
    i = min_index(abs(a))
    sr = a(i)

    if(i == size(a)) then
      ! sr at the end, thus shuffle as in (62)
      if(verb >= 30) print*, 'sr at the end'
      mminus1 = find_amount_trailing_zeros(a(1:size(a)-1))
      res = remove_sr_from_last_place_in_G(a(1:size(a)-mminus1-1),y2,mminus1+1,sr)
      return
    end if    

    if(i == 1) then
      !s_r at beginning, thus use (68)

      if(verb >= 30) then 
        print*, '--------------------------------------------------'
        print*, 'sr at beginning, map to: '
        call print_G([cmplx(0), a(i+1:size(a))], y2)
        call print_G([y2], sr)
        call print_G(a(i+1:size(a)), y2)
        print*, 'PI with p=',real([sr, a(i+1)]),'i=',i,'g =', real([a(i+2:size(a)), y2])
        call print_G([a(i+1)], sr)
        call print_G(a(i+1:size(a)), y2)
        print*, '--------------------------------------------------'
      end if

      res = G_flat([cmplx(0), a(i+1:size(a))], y2) &
        + G_flat([y2], sr) * G_flat(a(i+1:size(a)), y2) &
        + pending_integral([sr, a(i+1)], i, [a(i+2:size(a)), y2]) &
        - G_flat([a(i+1)], sr) * G_flat(a(i+1:size(a)), y2)
      return
    end if

    ! so s_r in middle, use (67)
    if(verb >= 30) then 
      print*, '--------------------------------------------------'
      print*, 'sr in the middle, map to: '
      call print_G([a(1:i-1),cmplx(0),a(i+1:size(a))],y2)
      print*, 'PI with p=',real([sr,a(i-1)]),'i=', i-1,'g =', real([a(1:i-2),a(i+1:size(a)),y2])
      call print_G([a(i-1)],sr)
      call print_G([a(1:i-1),a(i+1:size(a))],y2)
      print*, 'and PI with p=',real([sr,a(i+1)]),'i=',i,'g =', real([a(1:i-1),a(i+2:size(a)),y2])
      call print_G([a(i+1)],sr)
      call print_G([a(1:i-1),a(i+1:size(a))],y2)
      print*, '--------------------------------------------------'
    end if

    res = G_flat([a(1:i-1),cmplx(0),a(i+1:size(a))],y2) &
      - pending_integral([sr,a(i-1)], i-1, [a(1:i-2),a(i+1:size(a)),y2])  &
      + G_flat([a(i-1)],sr) * G_flat([a(1:i-1),a(i+1:size(a))],y2)        &
      + pending_integral([sr,a(i+1)], i, [a(1:i-1),a(i+2:size(a)),y2])    &
      - G_flat([a(i+1)],sr) * G_flat([a(1:i-1),a(i+1:size(a))],y2)        
  END FUNCTION make_convergent

  RECURSIVE FUNCTION improve_convergence(z) result(res)
    ! improves the convergence by applying the Hoelder convolution to G(z1,...zk,1)
    complex(kind=prec) :: z(:),oneminusz(size(z)), res
    complex(kind=prec), parameter :: p = 2.0
    integer :: k, j
    if(verb >= 30) print*, 'requires Hoelder convolution'
    oneminusz = 1-z
    k = size(z)

    res = G_flat(z,1.0/p) ! first term of the sum
    res = res + (-1)**k * G_flat(oneminusz(k:1:-1), 1.0-1.0/p)
    do j = 1,k-1
      res = res + (-1)**j * G_flat(oneminusz(j:1:-1),1.0-1.0/p) * G_flat(z(j+1:k),1.0/p)
    end do
  END FUNCTION improve_convergence

  RECURSIVE FUNCTION G_flat(z_flat,y) result(res)
    ! Calls G function with flat arguments, that is, zeroes not passed through the m's. 
    complex(kind=prec) :: z_flat(:), y, res
    complex(kind=prec), allocatable :: z(:), s(:,:)
    integer :: m_prime(size(z_flat)), condensed_size, kminusj, j, k, i, m_1
    integer, allocatable :: m(:)
    logical :: is_depth_one

    if(verb >= 50) call print_G(z_flat,y)

    ! catch G(1,1)
    if(size(z_flat) == 1) then
      if( abs(z_flat(1) - y) <= zero ) then
        print*, 'catch G(1,1)'
        res = 0
        return
      end if
    end if

    ! add small imaginary part if not there
    ! do i = 1,size(z_flat)
    !   if(abs(aimag(z_flat(i))) < 1e-25) z_flat(i) = add_ieps(z_flat(i))
    !   if(abs(aimag(y)) < 1e-25) y = add_ieps(y)
    ! end do

    ! is just a logarithm? 
    if(all(abs(z_flat) < zero)) then
      if(verb >= 70) print*, 'all z are zero'
      res = log(y)**size(z_flat) / factorial(size(z_flat))
      return
    end if
    if(size(z_flat) == 1) then
      if(verb >= 70) print*, 'is just a logarithm'
      if(abs(z_flat(1)) <= zero) then 
        res = log(y)
        return 
      end if
      res = log((z_flat(1) - y)/z_flat(1))
      return
    end if

    ! is it a polylog? 
    m_prime = get_condensed_m(z_flat)
    m_1 = m_prime(1)
    is_depth_one = (count((m_prime>0)) == 1)
    if(is_depth_one) then
      ! case m >= 2, other already handled above
      if(verb >= 70) print*, 'is just a polylog'
      res = -polylog(m_1,y/z_flat(m_1))
      return
    end if

    ! need remove trailing zeroes?
    k = size(z_flat)
    kminusj = find_amount_trailing_zeros(z_flat)
    j = k - kminusj
    if(all(abs(z_flat) < zero)) then 
      res = GPL_zero_zi(k,y)
      return
    else if(kminusj > 0) then
      if(verb >= 30) print*, 'need to remove trailing zeroes'
      allocate(s(j,j))
      s = shuffle_with_zero(z_flat(1:j-1))
      res = log(y)*G_flat(z_flat(1:size(z_flat)-1),y)
      do i = 1,size(s,1)
        res = res - G_flat([s(i,:),z_flat(j),zeroes(kminusj-1)], y)
      end do
      res = res / kminusj
      deallocate(s)
      return
    end if

    ! need make convergent?
    if(.not. is_convergent(z_flat,y)) then
      if(verb >= 10) print*, 'need to make convergent'
      res = make_convergent(z_flat, y)
      return
    end if

    ! requires Hoelder convolution?
    if( any(1.0 <= abs(z_flat/y) .and. abs(z_flat/y) <= 1.1) ) then
      res = improve_convergence(z_flat/y)
      return
    end if



    ! thus it is convergent, and has no trailing zeros
    ! -> evaluate in condensed notation -> which will give series representation
    m_prime = get_condensed_m(z_flat)
    if(find_first_zero(m_prime) == -1) then
      condensed_size = size(m_prime)
    else
      condensed_size = find_first_zero(m_prime)-1 
    end if
    allocate(m(condensed_size))
    allocate(z(condensed_size))
    m = m_prime(1:condensed_size)
    z = get_condensed_z(m,z_flat)
    res = G_condensed(m,z,y,size(m))
    deallocate(m)
    deallocate(z)
  END FUNCTION G_flat

  FUNCTION G_superflat(g) result(res)
    ! simpler notation for flat evaluation
    complex(kind=prec) :: g(:), res
    res = G_flat(g(1:size(g)-1), g(size(g)))
  END FUNCTION G_superflat

  FUNCTION G_real(g) result(res)
    ! simpler notation for flat evaluation
    real(kind=prec) :: g(:)
    complex(kind=prec) :: res
    res = G_flat(cmplx(g(1:size(g)-1)), cmplx(g(size(g))))
  END FUNCTION G_real

  FUNCTION G_int(g) result(res)
    ! simpler notation for flat evaluation
    integer:: g(:)
    complex(kind=prec) :: res
    res = G_flat(cmplx(g(1:size(g)-1)), cmplx(g(size(g))))
  END FUNCTION G_int

  RECURSIVE FUNCTION G_condensed(m,z,y,k) result(res)
    ! computes the generalized polylogarithm G_{m1,..mk} (z1,...zk; y)
    ! assumes zero arguments expressed through the m's
    
    integer :: m(:), k, i
    complex(kind=prec) :: z(:), x(k), y, res, z_flat(sum(m))

    ! print*, 'called G_condensed with args'
    ! print*, 'm = ', m
    ! print*, 'z = ', abs(z)

    ! are all z_i = 0 ? 
    if(k == 1 .and. abs(z(1)) < zero) then
      ! assumes that the zeros at the beginning are passed through m1
      res = GPL_zero_zi(m(1),y)
      return
    end if

    ! has trailing zeroes?
    if(abs(z(k)) < zero ) then
      ! we remove them in flat form
      z_flat = get_flattened_z(m,z)
      res = G_flat(z_flat,y)
      return
    end  if

    ! need make convergent?
    if(.not. GPL_has_convergent_series(m,z,y)) then
      z_flat = get_flattened_z(m,z)
      res = G_flat(z_flat,y)
      return
    end if

    x(1) = y/z(1)
    do i = 2,k
      x(i) = z(i-1)/z(i)
    end do
    ! print*, 'computed using Li with '
    ! print*, 'm = ', m
    ! print*, 'x = ', x
    res = (-1)**k * MPL(m,x)
  END FUNCTION G_condensed



  FUNCTION G_SUPERFLATN(c0,n)
    integer, intent(in) :: n
    complex(kind=prec), intent(in) :: c0(n)
    complex(kind=prec) g_superflatn
    G_superflatn=G_superflat(c0)
  END FUNCTION

END MODULE gpl_module

