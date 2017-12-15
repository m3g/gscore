module functionpars

  implicit none
  integer :: nbins
  double precision, allocatable :: s(:), sprob(:)

end module functionpars


program callfithist

  use functionpars
  implicit none
  integer :: ibin, ioerr, i
  double precision :: xread, yread
  character(len=200) :: record
  double precision :: fitpars(4), error
  double precision :: sprobfit

  open(10,file="hist.dat")
  nbins = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) xread, yread
    if ( ioerr /= 0 ) cycle
    nbins = nbins + 1
  end do
  rewind(10)

  allocate(s(nbins),sprob(nbins))

  ibin = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) xread, yread
    if ( ioerr /= 0 ) cycle
    ibin = ibin + 1
    s(ibin) = xread
    sprob(ibin) = yread
  end do
  close(10)

  call fithist(fitpars,error)
  write(*,*) fitpars, error

  open(10,file='fit.dat')
  do i = 1, nbins
    write(10,*) s(i), sprob(i), sprobfit(fitpars,s(i))
  end do
  close(10)

end program callfithist

subroutine fithist(x,fbest)

  implicit none
  integer :: i
  integer, parameter :: n = 4 ! Number of variables
  integer :: seed
  integer :: ntrial, itrial, nbest, best_repeat

  double precision :: x(*)

  double precision :: f, fbest
  double precision :: xbest(4)
  double precision :: random

  ! For cgnewton

  integer :: optpars(10), maxcg, maxfunc
  external :: computef, computeg

  ! Parameters for optimization method

  maxfunc = 50
  maxcg = 20
  optpars(1) = maxfunc ! Maximum number of functional evaluations
  optpars(2) = maxcg   ! Maximum number of CG iterations

  ! initialize random number generator

  call seed_from_time(seed)
  !seed = 12345
  call init_random_number(seed)

  ntrial = 1000
  nbest = 5
  best_repeat = 0
  itrial = 0
  fbest = 1.d30
  do while( best_repeat < nbest .and. itrial < ntrial )
    itrial = itrial + 1

    ! Initial link guess

    call random_number(random) ; x(1) = 0.d0 + 1.d0*random
    call random_number(random) ; x(2) = 10.d0 + 90.d0*random
    call random_number(random) ; x(3) = 0.d0 + 1.d0*random
    call random_number(random) ; x(4) = -1.d0 + 2.d0*random

    ! Test analytical gradient (debugging purposes only)
    ! call test_grad(n,x,computef,computeg)
    ! stop

    ! Minimize the energy of the linker

    call callcgnewton(n,x,f,0,computef,computeg,optpars)

    ! Save best fit up to now

    call computef(n,x,f)
    if ( abs((f - fbest)/fbest) < 1.d-3 ) then
      best_repeat = best_repeat + 1
    end if
    if ( f < fbest ) then
      fbest = f 
      do i = 1, n
        xbest(i) = x(i)
      end do
      !write(*,*) ' best f = ', f
    end if

  end do

end subroutine fithist

!
! Function that computes the value of the fit, given the fit parametes and the
! histogram point
!

double precision function sprobfit(x,scurrent)

  use functionpars
  implicit none
  double precision :: x(*)
  double precision :: scurrent, oneovers
 
  oneovers = 1.d0 / scurrent - 1.d0
  sprobfit = x(1)* ( 1.d0/(1.d0+exp(-x(2)*(scurrent-x(3)))) ) * oneovers * exp(-x(4)*oneovers)

end function sprobfit






