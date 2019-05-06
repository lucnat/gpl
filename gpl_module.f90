
MODULE gpl_module
  use mpl_module
  use utils
  implicit none

CONTAINS 

  RECURSIVE FUNCTION factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res
    res = merge(1,n*factorial(n-1),n==1)
  END FUNCTION factorial

  FUNCTION zeta(n) 
    real(kind=prec) :: values(9), zeta
    integer :: n
    values = (/1.6449340668482262, 1.2020569031595942, 1.0823232337111381, &
               1.03692775514337, 1.0173430619844488, 1.008349277381923, & 
               1.0040773561979441, 1.0020083928260821, 1.000994575127818/)
    zeta = values(n-1)
  END FUNCTION zeta

  FUNCTION GPL_has_convergent_series(m,z,y,k)
    ! tests if GPL has a convergent series representation
    integer :: m(:), k
    complex(kind=prec) :: z(:), y
    logical :: GPL_has_convergent_series

    GPL_has_convergent_series = .false.

    if(all(abs(y) <= abs(z))) then
      if(m(1) == 1) then 
        GPL_has_convergent_series = (y/z(1) /= 1)
      else 
        GPL_has_convergent_series = .true.
      end if
    end if
  END FUNCTION GPL_has_convergent_series

  FUNCTION GPL_zero_zi(l,y)
    ! used to compute the value of GPL when all zi are zero
    integer :: l
    complex(kind=prec) :: y, GPL_zero_zi
    print*, 'computed value using zi = 0'
    GPL_zero_zi = 1.0d0/factorial(l) * log(y) ** l
  END FUNCTION GPL_zero_zi

  FUNCTION is_convergent(z,y)
    ! returns true if G(z,y) convergent, otherwise false
    complex(kind=prec) :: z(:), y
    logical :: is_convergent
    integer :: i

    is_convergent = .true.
    do i = 1,size(z)
      if(abs(z(i)) < zero) cycle  ! skip zero values
      if(abs(y) > abs(z(i))) is_convergent = .false.
    end do

  END FUNCTION is_convergent

  RECURSIVE FUNCTION G_flat(z_flat,y) result(res)
    ! Calls G function with flat arguments, that is, zeroes not passed through the m's. 
    complex(kind=prec) :: z_flat(:), y, res
    complex(kind=prec), allocatable :: z(:), s(:,:)
    integer :: m_prime(size(z_flat)), condensed_size, kminusj, j, k, i
    integer, allocatable :: m(:)

    print*, 'G_flat called with args', abs(z_flat)

    ! remove trailing zeroes
    k = size(z_flat)
    kminusj = find_amount_trailing_zeros(z_flat)
    j = k - kminusj
    if(all(abs(z_flat) < zero)) then 
      print*, 'all are zero', abs(z_flat)
      res = GPL_zero_zi(k,y)
      return
    else if(kminusj > 0) then
      print*, 'we have',kminusj,'trailing zeroes'
      allocate(s(j,j))
      s = shuffle_with_zero(z_flat(1:j-1))
      res = log(y)*G_flat(z_flat(1:size(z_flat)-1),y)
      do i = 1,size(s,1)
        res = res - G_flat([s(i,:),z_flat(j),zero_array(kminusj-1)], y)
      end do
      res = res / kminusj
      deallocate(s)
      return
    end if

    ! transform to condensed notation
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
  END  FUNCTION G_flat

  RECURSIVE FUNCTION G_condensed(m,z,y,k) result(res)
    ! computes the generalized polylogarithm G_{m1,..mk} (z1,...zk; y)
    ! assumes zero arguments expressed through the m's
    
    integer :: m(:), k, i
    complex(kind=prec) :: z(:), x(k), y, res, c(sum(m)+1,sum(m)+1), z_flat(sum(m)), a(sum(m)-1)

    print*, 'called G_condensed with args'
    print*, 'm = ', m
    print*, 'z = ', abs(z)

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
    end  if

    ! need make convergent?
    if(.not. GPL_has_convergent_series(m,z,y,k)) then
      print*, 'need to make convergent'
      if(k == 1 .and. m(1) == 1) then
        print*, 'now we use the easy case. ha. ha.'
      else
        print*, '  ', 'does not have convergent series representation'
      end if
      res = 0
      return
    end if

    do i = 1,k
      x(i) = merge(y/z(1), z(i-1)/z(i),i == 1)
    end do
    print*, 'computed using MPL'
    res = (-1)**k * MPL(m,x)
  END FUNCTION G_condensed

END MODULE gpl_module


