!
! G-score correlation calculator
!
! Program gcorrelation
!
! This program writes the similarity of the models to a reference,
! as function of their g-scores, to test if the there is a correlation.
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 26, 2016
! http://leandro.iqm.unicamp.br
!
program gcorrelation

  use types
  use file_operations
  implicit none
  integer :: i, j, iref, imodel
  integer :: narg, ioerr, nmodels, model_index
  double precision :: gscore
  double precision, allocatable :: scores(:,:)
  character(len=200) :: compactlog, gscorefile, record, output, reference, name
  logical :: error
  type(model_type), allocatable :: model(:)

  narg = iargc()
  if ( narg /= 4 ) then
    write(*,*) ' ERROR: Run with: ./gcorrelation [reference] [compact align log] [gscore output] [output file]'
    stop
  end if
  call getarg(1,reference)
  call getarg(2,compactlog)
  call getarg(3,gscorefile)
  call getarg(4,output)

  ! Print the input options

  write(*,"(a)") "#" 
  write(*,"(a)") "# G-score correlation calculator " 
  write(*,"(a)") "#" 
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# Reference model: ", trim(adjustl(reference)) 
  write(*,"(a,a)") "# Alignment log file (compact form): ", trim(adjustl(compactlog)) 
  write(*,"(a,a)") "# G-score data file: ", trim(adjustl(gscorefile)) 
  write(*,"(a,a)") "# Output file: ", trim(adjustl(output)) 
  write(*,"(a)") "#" 

  ! Open the compact log file

  open(10,file=compactlog,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open alignment log file: ', trim(adjustl(compactlog))
    stop
  end if

  ! Read model list from log file (first two lines are comments to be ignored)

  read(10,*)
  read(10,*)
  read(10,*) nmodels
  allocate(scores(nmodels,nmodels),model(nmodels))

  ! Read scores and find index of reference model
 
  write(*,"(a)") "# Reading similarity scores from file ... "
  do i = 1, nmodels
    read(10,*,iostat=ioerr) model(i)%index, model(i)%name
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Problem reading model list from align log. '
      stop
    end if
    if ( reference == model(i)%name ) iref = i
  end do
  do i = 1, nmodels-1
    call progress(i,1,nmodels)
    read(10,*,iostat=ioerr) (scores(i,j),j=i+1,nmodels)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Problem reading alignment scores from log. '
      stop
    end if
  end do
  close(10)
  call progress(nmodels,1,nmodels)

  ! Compute G-score for all models

  write(*,"(a)") "# Reading G-scores from file ... "
  open(10,file=gscorefile,status='old',action='read') 
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open G-scores file: ', trim(adjustl(gscorefile))
    stop
  end if
  i = 0
  do 
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    i = i + 1
    call progress(i,1,nmodels)
    read(record,*,iostat=ioerr) gscore, name
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read gscore and model name from line: '
      write(*,*) '       ', trim(adjustl(record))
      stop
    end if
    imodel = model_index(name,model,nmodels,error)
    if ( error ) then
      write(*,*) ' ERROR: Model listed in G-score file is not in align log file. '
      write(*,*) '        Model name: ', trim(adjustl(name))
      stop
    end if
    model(imodel)%gscore = gscore
  end do
  close(10)

  ! Now check the similarity of each model to the reference model
  
  write(*,"(a)") "# Checking the similarity of each model to the reference ... "
  do i = 1, nmodels
    call progress(i,1,nmodels)
    if ( i < iref ) then
      model(i)%similarity = scores(i,iref)
    else if ( i > iref ) then
      model(i)%similarity = scores(iref,i)
    end if
  end do

  !
  ! Sort models from greater to lower similarity
  !

  call sort_by_similarity(nmodels,model)

  !
  ! Write output file 
  !

  write(*,"(a)") "# Writing output file ... "
  open(10,file=output,action='write')
  write(10,"(a)") "# Output of gcorrelation"
  write(10,"(a,a)") "# Input alignment log: ", trim(adjustl(compactlog))
  write(10,"(a,a)") "# Input gscore file: ", trim(adjustl(gscorefile))
  write(10,"(a,a)") "# Reference model: ", trim(adjustl(reference))
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score    Similarity  Model"
  do i = 1, nmodels
    if ( model(i)%index == iref ) cycle
    write(10,"(f12.5,tr2,f12.5,tr2,a)") model(i)%gscore, model(i)%similarity, &
                                        trim(adjustl(model(i)%name))
  end do
  close(10)

  write(*,"(a)") "# Finished. " 

end program gcorrelation

