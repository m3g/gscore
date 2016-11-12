!
! Subroutine that reads the data from the compact alignment log
! file. 
!
! Important: model and scores are allocated here.
!

subroutine read_compactlog(unit)

  use types
  use compactlog_data
  implicit none
  integer, intent(in) :: unit
  integer :: i, j, ioerr
  character(len=200) :: record

  ! Open file

  open(unit,file=compactlog,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open alignment log file: ', trim(adjustl(compactlog))
    stop
  end if   

  ! Read model list from log file 

  write(*,"(a)") "# Reading scores from compactlog file ... "

  read(unit,"(a200)",iostat=ioerr) record
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read data from file: ', trim(adjustl(compactlog))
    stop
  end if
  read(unit,*) ! Ignore line
  read(unit,*) ! Ignore line
  read(unit,"(a200)") record
  if ( index(record,"Score type:") /= 0 ) then
    if ( index(record,"GDT_TS") /= 0 ) then
      write(*,"(a)") '# File contains GDT_TS score information. '
      score_type = 1
    end if
    if ( index(record,"TM-score") /= 0 ) then
      write(*,"(a)") '# File contains TM-score information. '
      score_type = 2
    end if
    if ( index(record,"Contact-correlation") /= 0 ) then
      write(*,"(a)") '# File contains Contact-correlation information. '
      score_type = 3
    end if
  end if
  read(unit,*) nmodels
  allocate(scores(nmodels,nmodels),model(nmodels))

  ! Read scores 
 
  do i = 1, nmodels
    read(unit,*) j, model(i)%name
  end do
  do i = 1, nmodels-1
    call progress(i,1,nmodels)
    read(unit,*) (scores(i,j),j=i+1,nmodels)
  end do
  call progress(nmodels,1,nmodels)
  close(unit)

end subroutine read_compactlog                       

!
! Subroutine that reads the data from the xcompactlog file
!

subroutine read_xcompactlog(unit)

  use types
  use xcompactlog_data
  implicit none
  integer :: i, j, ioerr, unit
  character(len=200) :: record

  open(unit,file=compactlog,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open alignment log file: ', trim(adjustl(compactlog))
    stop
  end if

  ! Read model list from log file

  read(unit,"(a200)",iostat=ioerr) record
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read data from file: ', trim(adjustl(compactlog))
    stop
  end if
  read(unit,*) ! Ignore line
  read(unit,*) ! Ignore line
  read(unit,"(a200)") record
  if ( index(record,"Score type:") /= 0 ) then
    if ( index(record,"GDT_TS") /= 0 ) then
      write(*,"(a)") '# File contains GDT_TS score information. '
      score_type = 1
    end if
    if ( index(record,"TM-score") /= 0 ) then
      write(*,"(a)") '# File contains TM-score information. '
      score_type = 2
    end if
  end if

  ! Number of models of first set
  read(unit,*) nmodels1
  allocate(model1(nmodels1))

  ! Read model1 names
  do i = 1, nmodels1
    read(unit,*) model1(i)%index, model1(i)%name
  end do

  ! Number of models of second set
  read(unit,*) nmodels2
  allocate(model2(nmodels2))

  ! Read model2 names
  do i = 1, nmodels2
    read(unit,*) model2(i)%index, model2(i)%name
  end do

  ! Reading scores
  write(*,"(a)") "# Reading scores from file ... "
  allocate(scores(nmodels1,nmodels2))
  do i = 1, nmodels1
    call progress(i,1,nmodels1)
    read(unit,*) (scores(i,j),j=1,nmodels2)
  end do
  close(unit)                     

end subroutine read_xcompactlog




