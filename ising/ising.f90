program main

!*****************************************************************************80
!
!! MAIN is the main program for ISING.
!
!  Modified:
!
!    30 November 2007
!
  implicit none

  logical, allocatable, dimension ( :, :, : ) :: ising
  integer iterations
  integer n
  real, dimension ( 0:6 ) :: prob

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING'
  write ( *, '(a)' ) '  Monte Carlo simulation of a 3D Ising model.'
!
!  Set data.
!
  n = 25 
  iterations = 1000

! prob = (/ 0.97E+00, 0.95E+00, 0.85E+00, 0.50E+00, 0.15E+00, 0.05E+00, 0.03E+00 /)
  prob = (/ 0.99E+00, 0.99E+00, 0.99E+00, 0.80E+00, 0.30E+00, 0.20E+00, 0.10E+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The linear dimension of the system is N = ', n
  write ( *, '(a,i8)' ) '  The number of iterations taken is ITERATIONS = ', iterations
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The transition probability table, based on the number of'
  write ( *, '(a)' ) '  neighbors with the same spin.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      0         1         2         3         4         5         6'
  write ( *, '(a)' ) ' '
  write ( *, '(7f10.4)' ) prob(0:6)
!
!  Initialize the system.
!
  allocate ( ising(n,n,n) )

  call ising_initialize ( n, ising )
!
!  Do the simulation.
!
  call transition ( n, ising, iterations, prob )

  deallocate ( ising )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ising_initialize ( n, ising )

!*****************************************************************************80
!
!! ISING_INITIALIZE initializes the Ising array.
!
!  Modified:
!
!    01 December 2007
!
!  Parameters:
!
!    Input, integer N, the number of cells in each spatial dimension.
!
!    Output, logical ISING(N,N,N), the initial Ising array.
!
  implicit none

  integer n

  logical ising(n,n,n)
  real r(n,n,n)

  call random_number ( harvest = r )

  ising = ( r <= 0.5E+00 )

  return
end
subroutine transition ( n, ising, iterations, prob )

!*****************************************************************************80
!
!! TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
!
!  Modified:
!
!    30 November 2007
!
!  Reference:
!
!    American National Standard for Programming Language: Fortran - Extended,
!    American National Standards Institute, 1992,
!    pages 296-299.
!
!  Parameters:
!
!    Input, integer N, the number of cells in each spatial dimension.
!
!    Input/output, logical ISING(N,N,N).  On input, the current state of the
!    system.  On output, the state of the system after the iterations.
!
!    Input, integer ITERATIONS, the number of iterations to carry out.
!
!    Input, real PROB(0:6).  PROB(I) represents the probability that the spin
!    of a given cell will be reversed, given that it has I immediate 
!    neighbors with spin the same as its own.
!
  implicit none

  integer n

  integer count(n,n,n)
  integer flip(n,n,n)
  integer flip_count
  integer i
  logical ising(n,n,n)
  integer iterations
  integer ones(n,n,n)
  integer ones_count
  real prob(0:6)
  real r(n,n,n)
  real threshhold(n,n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step     Positives       Flipped'
  write ( *, '(a)' ) ' '

  flip_count = 0

  ones = 0
  where ( ising ) ones = 1

  ones_count = sum ( ones )

  write ( *, '(2x,i4,2x,i12,2x,i12)' ) 0, ones_count

  do i = 1, iterations
!
!  COUNT contains the number of immediate neighbors of (I,J,K) with a value of 1.
!
    count = cshift ( ones, -1,  1 ) &
          + cshift ( ones, +1,  1 ) &
          + cshift ( ones, -1,  2 ) &
          + cshift ( ones, +1,  2 ) &
          + cshift ( ones, -1,  3 ) &
          + cshift ( ones, +1,  3 )
!
!  Now COUNT contains the number of immediate neighbors of (I,J,K) that are the same
!  as the (I,J,K) value.
!
    where ( .not. ising ) count = 6 - count

    where ( count == 0 ) threshhold = prob(0)
    where ( count == 1 ) threshhold = prob(1)
    where ( count == 2 ) threshhold = prob(2)
    where ( count == 3 ) threshhold = prob(3)
    where ( count == 4 ) threshhold = prob(4)
    where ( count == 5 ) threshhold = prob(5)
    where ( count == 6 ) threshhold = prob(6)

    call random_number ( harvest = r )
!
!  "Flip" the value of (I,J,K) with probability PROB(COUNT).
!  
    where ( threshhold < r )
      flip = 1
    elsewhere
      flip = 0
    endwhere

    flip_count = sum ( flip )

    write ( *, '(2x,4x,2x,12x,2x,i12)' ) flip_count

    where ( flip == 1 ) 
      ising = .not. ising
    endwhere

    where ( ising ) 
      ones = 1
    elsewhere
      ones = 0
    endwhere

    ones_count = sum ( ones )

    write ( *, '(2x,i4,2x,i12)' ) i, ones_count

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
