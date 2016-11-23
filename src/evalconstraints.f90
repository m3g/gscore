!
! Program EvalConstraints
!
! Reads a series of structure files, a G-score output file for these
! structures, and a list of contraints, and computes the constraint
! statistics for each structure. Order the structure using the G-score,
! to check the convergence of the number of constraints satisfied as
! a function of the number of constraints considered
! 
! Also outputs a matrix that provides, for each constraint, the correlation
! of the constraints in the model set
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
  end type
  integer :: i, j, iconst, imodel, ires, model_index, nsatisfied, nconstraints, nres
  integer :: nargs, nmodels, ioerr
  real :: dij, gscore, d
  character(len=200) :: record, name, output, gscorefile, pdblist, constraintsfile, &
                        correlationout, format
  logical :: stop = .true.
  integer, allocatable :: satisfied(:)
  real, allocatable :: x(:,:), correlation(:,:)
  type(model_type), allocatable :: model(:)
  type(constraint_type), allocatable :: constraint(:)

  ! Print title

  call title()
  write(*,"(a)") '# EVALCONSTRAINTS: Evaluate the models from the constraint set. '
  write(*,"(a)") '#'

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 4 ) then
    write(*,*)
    write(*,*) ' Run with: evalconstraints pdblist.txt gscore.dat constraints.dat output.dat'
    write(*,*)
    write(*,*) ' Where: pdblist.txt is the file containing the list of PDB models to consider'
    write(*,*) '        gscores.dat it file containing the scores of each structure. '
    write(*,*) '        constraints.dat is the file containing the set of constraints. '
    write(*,*) '        output.dat is the name of the output file to be created. '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/gscore '
    write(*,*)
    stop
  end if

  ! Read file names from command line

  call getarg(1,pdblist)
  call getarg(2,gscorefile)
  call getarg(3,constraintsfile)
  call getarg(4,output)

  ! Read number of models

  open(10,file=pdblist,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open pdb list file: ', trim(adjustl(pdblist))
    stop
  end if
  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    nmodels = nmodels + 1
  end do
  write(*,"(a)") '# Number of models found in pdb list file: ', nmodels
  allocate(model(nmodels))

  ! Read model file names

  rewind(10)
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    imodel = imodel + 1
    model(imodel)%file = record
    model(imodel)%name = basename(record)
  end do
  close(10)

  ! Order models by name

  write(*,"(a)") "# Sorting models by name ... "
  call sort_by_name(model,nmodels)

  ! Read the G-scores of the models

  write(*,"(a)") "# Reading gscore file ... "
  open(10,file=gscorefile,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(gscorefile))
    stop
  end if
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) gscore, name
    imodel = model_index(name,model,nmodels,stop)
    model(imodel)%gscore = gscore
  end do
  close(10)

  ! Read number of constraints from contraint file

  write(*,"(a)") "# Reading constraint file ... "
  open(10,file=constraintsfile,action='read',status='old')
  nconstraints = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) i, j, d
    if ( ioerr /= 0 ) nconstraints = nconstraints + 1 
  end do
  write(*,"(a,i8)") '# Number o constraints: ', nconstraints
  allocate(constraint(nconstraints))
  do imodel = 1, nmodels
    allocate(model(imodel)%constraint(nconstraints))
  end do
 
  ! Read constraints

  rewind(10)
  iconst = 0
  nres = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) i, j, d
    if ( ioerr /= 0 ) iconst = iconst + 1
    constraint(iconst)%i = i
    constraint(iconst)%j = j
    constraint(iconst)%d = d
    d = d**2
    nres = max(nres,i)
    nres = max(nres,j)
  end do
  close(10)
  write(*,"(a,i8)") '# Greatest residue number associated with constraints: ', nres
  allocate(x(nres,3))

  ! Evaluate the satisfaction of the constraints by the models

  write(*,"(a)") "# Reading model coordinates and computing contacts ... "
  do imodel = 1, nmodels
    call progress(i,1,nmodels)
   
    ! Read coordinates of atoms that belong to constraints in this model

    open(10,file=model(imodel)%file,action='read',status='old',iostat=ioerr) 
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open PDB file: ', trim(adjustl(model(imodel)%file))
      stop
    end if
    do
      read(10,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( record(1:4) == "ATOM" .and. record(13:16) == "CA" ) then
        read(record(23:26),*,iostat=ioerr) ires
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Could not read residue number in file: ',&
                     trim(adjustl(model(imodel)%file))
          write(*,*) '  Line: ', trim(adjustl(record))
          stop
        end if
      end if
      read(record(31:38),*) x(ires,1)
      read(record(39:46),*) x(ires,2)
      read(record(47:54),*) x(ires,3)
    end do
    close(10)

    ! Compute logical contact vector for this model

    model(imodel)%ncontacts = 0
    do iconst = 1, nconstraints
      i = constraint(iconst)%i
      j = constraint(iconst)%j
      d = constraint(iconst)%d
      dij = ( x(i,1) - x(j,1) )**2 + &
            ( x(i,2) - x(j,2) )**2 + &
            ( x(i,3) - x(j,3) )**2
      if ( dij <= d ) then
        model(imodel)%constraint(iconst) = .true.
        model(imodel)%ncontacts = model(imodel)%ncontacts + 1
      else
        model(imodel)%constraint(iconst) = .false.
      end if
    end do
  end do

  ! Computing constraint correlation matrix

  allocate(correlation(nconstraints,nconstraints))
  do i = 1, nconstraints
    correlation(i,j) = 1.
  end do
  do i = 1, nconstraints - 1
    do j = i + 1,  nconstraints
      correlation(i,j) = 0.
    end do
  end do
  do imodel = 1, nmodels - 1
    do i = 1, nconstraints - 1
      do j = i + 1, nconstraints
        if ( model(imodel)%constraint(i) .and. model(imodel)%constraint(j) ) then
          correlation(i,j) = correlation(i,j) + 1
        else if( ( model(imodel)%constraint(i) .and. .not. model(imodel)%constraint(j) ) .or. &
                 ( .not. model(imodel)%constraint(i) .and. model(imodel)%constraint(j) ) ) then
          correlation(i,j) = correlation(i,j) - 1
        end if
      end do
    end do
  end do
  do i = 1, nconstraints - 1
    do j = i + 1, nconstraints
      correlation(i,j) = correlation(i,j) / nmodels
      correlation(j,i) = correlation(i,j)
    end do
  end do

  !
  ! Write output files
  ! 

  correlationout = trim(adjustl(remove_extension(output)))//"_correlation."//&
                   &trim(adjustl(file_extension(output)))
  output = trim(adjustl(remove_extension(output)))//"_models."//&
           &trim(adjustl(file_extension(output)))

  ! Write constraint correlation output file
 
  write(*,"(a,a)") "# Writting constraint correlation file: ", trim(adjustl(correlationout))
  open(10,file=correlationout,iostat=ioerr)
  if ( ioerr /= 0 ) then 
    write(*,*) ' ERROR: Could not open output file: ', trim(adjustl(output))
    stop
  end if
  write(format,"(a,i8,a)") "(",nconstraints,"(tr1,f5.2) )"
  do i = 1, nconstraints
    write(10,format) (correlation(i,j),j=1,nconstraints)
  end do
  close(10)

  ! Write file containing list of of models with constraint analysis  

  ! Ordering models by G-score

  call sort_by_gscore(model,nmodels)

  write(*,*) ' Writing constraint analysis file : ', trim(adjustl(output))
  open(10,file=output,iostat=ioerr)
  if ( ioerr /= 0 ) then 
    write(*,*) ' ERROR: Could not open output file: ', trim(adjustl(output))
    stop
  end if
  write(10,"(a)") "# G-score"
  write(10,"(a)") "#"
  write(10,"(a)") "# EvalConstraints output file. "
  write(10,"(a)") "#"
  write(10,"(a,a)") "# PDB list: ", trim(adjustl(pdblist))
  write(10,"(a,a)") "# G-score file: ", trim(adjustl(gscorefile))
  write(10,"(a,a)") "# Constraints file: ", trim(adjustl(constraintsfile))
  write(10,"(a,i8)") "# Number of models ", nmodels
  write(10,"(a,i8)") "# Number of constraints: ", nconstraints
  write(10,"(a)") "#"
  write(10,"(a)") "# Constraint list: "
  do i = 1, nconstraints
    write(10,"(a,i5,2(tr1,i5)tr1,f8.3)") "# ", i, constraint(i)%i, constraint(i)%j, constraint(i)%d
  end do
  write(10,"(a)") "#"
  write(10,"(a)") "# Nmodel: Number of constraints satisfied by this model. "
  write(10,"(a)") "# RelatP: Relative probability of this model (G-score ratio to best model)."
  write(10,"(a)") "# DeltaG: RelatP converted to DeltaG (kcal/mol)."
  write(10,"(a)") "# Ntot: Total number of constraints satisfied by the ensemble up to this model."
  write(10,"(a)") "# Next: constraint indexes according to list above."
  write(10,"(a)") "#"
  write(record,*) "(a,",nconstraints,"(tr1,i3))"
  write(10,record) "#             Model  Nmodel     RelatP       DeltaG  Ntot",(i,i=1,nconstraints)
  write(format,*) "(i8,tr1,a,tr1,i5,2(tr1,f12.5),tr1,i5,",nconstraints,"(tr1,i3))"
  nsatisfied = 0
  allocate(satisfied(nconstraints))
  do i = 1, nconstraints
    satisfied(i) = 0
  end do
  do imodel = 1, nmodels
    do i = 1, nconstraints
      if ( model(imodel)%constraint(i) ) then
        if ( satisfied(i) == 0 ) then
          nsatisfied = nsatisfied + 1
          satisfied(i) = 1
        end if
      end if
    end do
    write(10,format) &
             imodel, &
             trim(adjustl(model(imodel)%name)),&
             model(imodel)%ncontacts,&
             model(imodel)%gscore / model(1)%gscore, &
             -1.987*0.298*dlog(model(imodel)%gscore / model(1)%gscore), &
             nsatisfied, &
             (satisfied(i),i=1,nconstraints)
  end do
  close(10)
  write(*,"(a)") '#'
  write(*,"(a)") '# END'

end program evalconstraints

