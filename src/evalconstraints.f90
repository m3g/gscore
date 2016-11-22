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

  use ioformat 
  use types
  use file_operations
  use progress_bar
  implicit none
  type( constraint_type ) 
    integer :: i, j
    real :: d
  end type
  integer :: i, j, iconst, imodel, ires, model_index
  integer :: nargs, nmodels, ioerr
  real :: dij, gscore
  character(len=200) :: record, name, output, gscorefile, pdblist, constraintsfile 
  logical :: stop = .true.
  real, allocatable :: x(:,:)
  logical, allocatable :: cmodel(:), correlation(:,:)
  type(modeldata), allocatable :: model(:)
  type(constraint_type), allocatable :: constraint(:)

  ! Print title

  call title()
  write(*,"(a)") '# EVALCONSTRAINTS: Evaluate the models from the constraint set. '
  write(*,*)
  write(*,dashes)

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
    write(*,*) ' If the score list is not a LovoAlign log file, -c[int] indicates the column of the '
    write(*,*) ' list containing the score. '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/gscore '
    write(*,*)
    write(*,hashes)
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
  open(10,file=constraintfile,action='read',status='old')
  nconstraints = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) i, j, d
    if ( ioerr /= 0 ) nconstraints = nconstraints + 1 
  end do
  write(*,"(a,i8)") '# Number o constraints: ', nconstraints
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
    if ( ierr /= 0 ) then
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
    do iconst = 1, ncontacts
      i = constraint(iconst)%i
      j = constraint(iconst)%j
      d = constraint(iconst)%d
      dij = ( x(i,1) - x(j,1) )**2 + &
            ( x(i,2) - x(j,2) )**2 + &
            ( x(i,3) - x(j,3) )**2
      if ( dij <= d ) then
        cmodel(iconst) = .true.
        model(imodel)%ncontacts = model(imodel)%ncontacts + 1
      else
        cmodel(iconst) = .false.
      end if
    end do
  end do

  ! Computing constraint correlation matrix

  allocate(correlation(nconstraints,nconstraints))
  do imodel = 1, nmodels - 1
    do j = imodel + 1, nmodels
      correlation(i,j) = 0
      if ( cmodel(i) .and. cmodel(j) ) then
        correlation(i,j) = correlation(i,j) + 1
      else ( ( cmodel(i) .and. .not. cmodel(j) ) .or. &
             ( .not. cmodel(i) .and. cmodel(j) ) ) then
        correlation(i,j) = correlation(i,j) - 1
      end if
    end do
  end do


  !
  ! Write output file
  ! 

  write(*,*) ' Writing output file ... '
  open(10,file=output,iostat=ioerr)
  if ( ioerr /= 0 ) then 
    write(*,*) ' ERROR: Could not open output file: ', trim(adjustl(output))
    stop
  end if
  
  write(10,"(a)") "# TopoLink" 
  write(10,"(a)") "#"
  write(10,"(a)") "# EvalModels output file. " 
  write(10,"(a)") "#"
  write(10,"(a,a)") "# Log file list: ", trim(adjustl(loglist))
  write(10,"(a,a)") "# Score (possibly LovoAlign log) file: ", trim(adjustl(scorelist))
  write(10,"(a,i8)") "# Number of models ", nmodels
  write(10,"(a)") "#"
  write(10,"(a,i5,a)") "# Score: Model quality score, obtained from column ", scorecol,&
                       " of the score file. " 
  write(10,"(a)") "#"
  write(10,"(a)") "# RESULT0: Number of consistent observations. "
  write(10,"(a)") "# RESULT1: Number of topological distances consistent with all observations. "
  write(10,"(a)") "# RESULT2: Number of topological distances NOT consistent with observations. "
  write(10,"(a)") "# RESULT3: Number of missing links in observations. "
  write(10,"(a)") "# RESULT4: Number of distances with min and max bounds that are consistent."
  write(10,"(a)") "# RESULT5: Sum of the scores of observed links in all observations. "
  write(10,"(a)") "# RESULT6: Likelyhood of the structural model, based on observations. "
  write(10,"(a)") "#"
  write(10,"(a)") "# More details at: http://leandro.iqm.unicamp.br/topolink"
  write(10,"(a)") "#"
  write(10,"(a)") "#      Score   RESULT0   RESULT1   RESULT2   RESULT3   RESULT4       RESULT5       RESULT6  MODEL"
  do imodel = 1, nmodels
    call progress(imodel,1,nmodels)
    write(10,"( f12.5,5(tr2,i8),tr2,f12.5,tr2,e12.5,tr2,a )") &
                model(imodel)%score, &
                model(imodel)%nobscons, &
                model(imodel)%ntopcons, &
                model(imodel)%ntopnot, &
                model(imodel)%nmiss, &
                model(imodel)%nminmax, &
                model(imodel)%sumscores, &
                model(imodel)%likely,&
                trim(adjustl(model(imodel)%name))
  end do

  close(10)

  write(*,*) ' Wrote output file: ', trim(adjustl(output))
  write(*,*)
  write(*,hashes)

end program evalconstraints









