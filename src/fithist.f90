module functionpars

  implicit none
  integer :: nbins
  double precision, allocatable :: s(:), probs(:)

end module functionpars

program fithist

  use functionpars
  implicit none
  integer :: n 
  integer :: ioerr, ibin, seed
  integer :: ntrial, itrial, nbest, best_repeat

  double precision :: xread, yread
  double precision, allocatable :: x(:), g(:)
  double precision :: random

  character(len=200) record

  external :: computef, computeg

  ! Number of variables
  n = 4
  allocate(x(n),g(n))

  ! Parameters for optimization method

  !dbond2 = dbond**2
  !maxfunc = 50
  !maxcg = 20
  !optpars(1) = maxfunc ! Maximum number of functional evaluations
  !optpars(2) = maxcg   ! Maximum number of CG iterations
  !seed = 0
  !iguess = 1

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
  ibin = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) xread, yread
    if ( ioerr /= 0 ) cycle
    ibin = ibin + 1
    s(ibin) = xread
    probs(ibin) = yread
  end do
  close(10)
write(*,*) nbins
stop

  ! initialize random number generator

  !call seed_from_time(seed)
  seed = 12345
  call init_random_number(seed)

  ntrial = 1
  nbest = 1
  best_repeat = 0
  itrial = 0
  do while( best_repeat < nbest .and. itrial < ntrial )
    itrial = itrial + 1

    ! Initial link guess

    call random_number(random) ; x(1) = -100.d0 + 200.d0*random
    call random_number(random) ; x(2) = -100.d0 + 200.d0*random
    call random_number(random) ; x(3) = -100.d0 + 200.d0*random
    call random_number(random) ; x(4) = -100.d0 + 200.d0*random

    ! Test analytical gradient (debugging purposes only)
    call test_grad(n,x,g,computef,computeg)
    stop

    ! Minimize the energy of the linker

    !call callcgnewton(n,x,f,0,computef,computeg,optpars)

  end do

end program fithist

