!
! Program CheckConstraints
!
! Reads a protein structure PDB file and the list of constraints
! and evaluates how many constraints are satisfied, and outputs
! the list of constraints
!
! L. Martinez
! Instiute of Chemistry - University of Campinas
! November 22, 2016
! http://leandro.iqm.unicamp.br/gscore
!
program evalconstraints

  use types
  use file_operations
  implicit none
  type constraint_type
    integer :: i, j
    real :: d
  end type constraint_type
  integer :: i, j, iconst, ires, nconstraints, nres
  integer :: nargs, ioerr 
  real :: dij, d
  character(len=200) :: record, constraintsfile
  real, allocatable :: x(:,:)
  logical, allocatable :: cmodel(:)
  type(model_type) :: model
  type(constraint_type), allocatable :: constraint(:)

  ! Print title

  write(*,"(a)") "#" 
  write(*,"(a)") "# G-score calculator"
  call title()
  write(*,"(a)") '# CHECKCONSTRAINTS: Evaluate one model relative to the constraint set. '
  write(*,"(a)") '#'

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 2 ) then
    write(*,*)
    write(*,*) ' Run with: checkconstraints model.pdb constraints.dat'
    write(*,*)
    write(*,*) ' Where: model.pdb is the PDB file of the model.'
    write(*,*) '        constraints.dat is the file containing the set of constraints. '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/gscore '
    write(*,*)
    stop
  end if

  ! Read file names from command line

  call getarg(1,record)
  model%file = record
  call getarg(2,constraintsfile)

  ! Read number of constraints from contraint file

  write(*,"(a)") "# Reading constraint file ... "
  open(10,file=constraintsfile,action='read',status='old')
  nconstraints = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) i, j, d
    if ( ioerr /= 0 ) cycle
    nconstraints = nconstraints + 1 
  end do
  write(*,"(a,i8)") '# Number of constraints: ', nconstraints
  allocate(constraint(nconstraints),cmodel(nconstraints))
 
  ! Read constraints

  rewind(10)
  iconst = 0
  nres = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) i, j, d
    if ( ioerr /= 0 ) cycle
    iconst = iconst + 1
    constraint(iconst)%i = i
    constraint(iconst)%j = j
    constraint(iconst)%d = d
    nres = max(nres,i)
    nres = max(nres,j)
  end do
  close(10)
  write(*,"(a,i8)") '# Greatest residue number associated with constraints: ', nres
  allocate(x(nres,3))

  ! Evaluate the satisfaction of the constraints by the models

  write(*,"(a)") "# Reading model coordinates and computing contacts ... "
  
  open(10,file=model%file,action='read',status='old',iostat=ioerr) 
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open PDB file: ', trim(adjustl(model%file))
    stop
  end if
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( record(1:4) == "ATOM" .and. trim(adjustl(record(13:16))) == "CA" ) then
      read(record(23:26),*,iostat=ioerr) ires
      if ( ires > nres ) exit
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read residue number in file: ',&
                   trim(adjustl(model%file))
        write(*,*) '  Line: ', trim(adjustl(record))
        stop
      end if
      read(record(31:38),*) x(ires,1)
      read(record(39:46),*) x(ires,2)
      read(record(47:54),*) x(ires,3)
    end if
  end do
  close(10)

  ! Compute logical contact vector for this model

  model%ncontacts = 0
  do iconst = 1, nconstraints
    i = constraint(iconst)%i
    j = constraint(iconst)%j
    d = constraint(iconst)%d
    dij = sqrt( ( x(i,1) - x(j,1) )**2 + &
                ( x(i,2) - x(j,2) )**2 + &
                ( x(i,3) - x(j,3) )**2 )
    if ( dij <= d ) then
      cmodel(iconst) = .true.
      model%ncontacts = model%ncontacts + 1
    else
      cmodel(iconst) = .false.
    end if
  end do

  ! Write output

  write(*,"(a,i5)") '# Number of constraints satisfied by this model: ', model%ncontacts

  write(*,"(a)") '#'
  write(*,"(a)") '# END'

end program evalconstraints

