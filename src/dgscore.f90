!
! G-score calculation
!
! DG score: Computes the difference between G-scores computed relative
! to different ensembles of modes, to compute the DeltaG score (or bias score)
!
! L. Martinez
! Instiute of Chemistry - University of Campinas
! Sep 19, 2016
! 

program dgscore

  use types
  use file_operations
  integer :: narg, ioerr, nmodels, imodel
  double precision :: gscore
  character(len=200) file1, file2, output, record, name
  type(model_type), allocatable :: model(:)

  write(*,"(a)") "#"
  write(*,"(a)") "# Delta-G-score calculator"
  call title()
  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with dgscore [gscore file 1] [gscore file 2] [output]'
    write(*,*)
    stop
  end if

  call getarg(1,file1)
  call getarg(2,file2)
  call getarg(3,output)

  open(10,file=file1,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open gscore file 1: ', trim(adjustl(file1))
    stop
  end if
  
  open(20,file=file2,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open gscore file 2: ', trim(adjustl(file2))
    stop
  end if

  ! Read number of models from file 1

  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) gscore, name
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read model data from file: ', trim(adjustl(file1))
      write(*,*) '        at line: ', trim(adjustl(record))
      stop
    end if
    nmodels = nmodels + 1
  end do
  write(*,"(a,i8)") '# Number of models in file 1: ', nmodels
  allocate(model(nmodels))

  ! Read data from file 1

  rewind(10)
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) gscore, name
    imodel = imodel + 1
    model(imodel)%name = name
    model(imodel)%gscore = gscore
  end do
  close(10)

  ! Sort models by name

  call sort_by_name(nmodels,model)

  ! Read data from file 2, and compute deltas

  do 
    read(20,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) gscore, name
    if ( ioerr /= 0 ) cycle
    imodel = model_index(name,model,nmodels,error)
    model(imodel)%dgscore = model(imodel)%gscore - gscore
  end do
  close(20)

  ! Sort models by dgscore

  call sort_by_gscore(nmodels,model)

  ! Write output file

  open(10,file=output)
  write(10,"(a)") "# Output of DeltaG-score "
  write(10,"(a,a)") "# First G-score file: ", trim(adjustl(file1))
  write(10,"(a,a)") "# Second G-score file: ", trim(adjustl(file2))
  do imodel = 1, nmodels
    write(10,"(tr2,f12.5,tr2,a)") model(imodel)%dgscore, trim(adjustl(model(imodel)%name))
  end do
  close(10)
  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(output))
  write(*,"(a)") "#"
  write(*,"(a)") "# Finished."
 
end program dgscore







