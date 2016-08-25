module types

  type model_type
    character(len=200) :: name
    double precision :: pas_gdt
    double precision :: pas_tm
  end type model_type

end module types

program pbias

  use types
  implicit none
  integer :: i1, i2, i, j, model_index, imodel
  integer :: narg, ioerr, nmodels, nlogs, ilog
  double precision :: dummy, tmscore_read, gdt_read, gdt_cut, tm_cut
  double precision, allocatable :: gdt(:,:), tmscore(:,:)
  character(len=200) :: align_list, firstlog, align_log
  character(len=200) :: basename, record, file1, file2
  logical :: comment, repeated
  type(model_type), allocatable :: model(:)

  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with: ./pas [align_list] [GDT cut] [TM-Score cut] '
    stop
  end if

  ! Open first log list

  call getarg(1,align_list)
  open(10,file=align_list,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(align_list))
    stop
  end if

  ! Read gdt and tmscore cutoffs

  call getarg(2,record)
  read(record,*,iostat=ioerr) gdt_cut
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read GDT cutoff from command line. '
    stop
  end if
  call getarg(3,record)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read TM-Score cutoff from command line. '
    stop
  end if
  read(record,*,iostat=ioerr) tm_cut

  ! Print the input options

  write(*,"(a)") "#" 
  write(*,"(a)") "# PAS score calculator " 
  write(*,"(a)") "#" 
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# List of alignment files: ", trim(adjustl(align_list)) 
  write(*,"(a)") "#" 
  write(*,"(a,f12.5)") "# GDT cutoff: ", gdt_cut
  write(*,"(a,f12.5)") "# TM-score cutoff: ", tm_cut
  write(*,"(a)") "#" 
 
  ! Count the number of lovoalign log files and models

  nlogs = 0
  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) align_log
    if ( ioerr /= 0 ) exit
    if ( comment(align_log) ) cycle 
    nlogs = nlogs + 1
    
    ! Read number of alignments in this file (is the number of models
    ! of the target model database)

    if ( nlogs == 1 ) then
      firstlog = align_log
      open(20,file=align_log,action='read',status='old',iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not open file: ', trim(adjustl(align_log))
      end if
      nmodels = 0
      do
        read(20,"(a200)",iostat=ioerr) record
        if ( ioerr /= 0 ) exit
        if ( comment(record) ) cycle
        nmodels = nmodels + 1
      end do
      close(20)
    end if
  end do
  write(*,"(a,i10)") '# Number of lovoalign log files in list: ', nlogs
  write(*,"(a,i10)") '# Number of alignments in each log file: ', nmodels
  write(*,"(a)") "#" 
  if ( nlogs /= nmodels ) then
    write(*,*) ' ERROR: The number of alignment log files is not the same as '
    write(*,*) '        the number of models. This is not what is expected from '
    write(*,*) '        an all-to-all alignment. '
    stop
  end if

  allocate(gdt(nmodels,nmodels),tmscore(nmodels,nmodels),model(nmodels))

  !
  ! Assign an index to each model name
  !

  write(*,"(a)") "# Assigning indexes for each model ... "
  open(20,file=firstlog,action='read',status='old')
  imodel = 0
  do
    read(20,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*) file1, file2
    !
    ! First model will be target alignment model from this file
    !
    if ( imodel == 0 ) then
      imodel = imodel + 1
      model(1)%name = basename(file2)
      call progress(imodel,1,nmodels)
    end if
    !
    ! Index alignment models
    !
    repeated = .false.
    do i = 1, imodel
      if ( basename(file1) == model(i)%name ) then
        repeated = .true.
        exit
      end if
    end do
    if ( .not. repeated ) then
      imodel = imodel + 1
      model(imodel)%name = basename(file1)
      call progress(imodel,1,nmodels)
    end if
  end do
  close(20)

  !
  ! Now, reading all alignment log files and annotating the scores
  ! of the alignment of each pair
  !
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reading all alignment files ... "
  rewind(10)
  ilog = 0
  call progress(ilog,1,nlogs)
  do
    read(10,"(a200)",iostat=ioerr) align_log
    if ( ioerr /= 0 ) exit
    if ( comment(align_log) ) cycle
    open(20,file=align_log,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open alignment log file: ', trim(adjustl(align_log))
      write(*,*) '        Listed in: ', trim(adjustl(align_list))
      stop
    end if
    do
      read(20,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( comment(record) ) cycle
      read(record,*,iostat=ioerr) file1, file2, tmscore_read, (dummy,i=1,4), gdt_read
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read data in alignment log file: ', trim(adjustl(align_log))
        write(*,*) '        Content: ', trim(adjustl(record))
        stop
      end if
      file1 = basename(file1)
      file2 = basename(file2)
      i1 = model_index(file1,model,nmodels)
      i2 = model_index(file2,model,nmodels)
      tmscore(i1,i2) = tmscore_read
      gdt(i1,i2) = gdt_read
      tmscore(i2,i1) = tmscore_read
      gdt(i2,i1) = gdt_read
    end do
    close(20)
    ilog = ilog + 1
    call progress(ilog,1,nlogs)
  end do
  close(10)

  !
  ! End of file reading, now will compute the scores
  !

  ! Compute the PAS score for models of list1

  do i = 1, nmodels
    model(i)%pas_gdt = 0.d0
    model(i)%pas_tm = 0.d0
    do j = 1, nmodels
      if ( i == j ) cycle
      if ( gdt(i,j) > gdt_cut ) then 
        model(i)%pas_gdt = model(i)%pas_gdt + 1.d0
      end if
      if ( tmscore(i,j) > tm_cut ) then 
        model(i)%pas_tm = model(i)%pas_tm + 1.d0
      end if
    end do
    model(i)%pas_gdt = model(i)%pas_gdt / dble(nmodels-1)
    model(i)%pas_tm = model(i)%pas_tm / dble(nmodels-1)
  end do

  !
  ! Order models from best to worst accoding to PAS_GDT
  ! 
!voltar

  ! 
  ! Order models from best to worst according to PAS_TM
  !

  ! 
  ! Write lists of ordered models by each type of PAS score
  !

end program pbias

!
! Function that determines the index of a model from its name
!

function model_index(name,model,n)
 
  use types
  implicit none
  integer :: model_index, n
  character(len=200) :: name
  type(model_type) :: model(n)
      
  model_index = 1
  do while( name /= model(model_index)%name ) 
    model_index = model_index + 1
  end do

end function model_index

!
! Function that determines the basename of a file,
! removing the path and the extension
!

function basename(filename)

  integer :: i
  character(len=200) :: filename, basename

  basename = trim(adjustl(filename))
  i = length(basename)
  idot = i+1
  do while(basename(i:i) /= "/")
    if ( basename(i:i) == "." ) then
      idot = i
    end if
    i = i - 1
    if ( i == 0 ) exit
  end do
  i = i + 1
  basename = basename(i:idot-1)
  do i = idot, 200
    basename(i:i) = achar(32)
  end do

end function basename

!
! Function that determines the length of a string
!

function length(string)

  implicit none
  integer :: length
  character(len=200) :: string
  logical :: empty_char
  length = 200
  do while( empty_char(string(length:length)) ) 
    length = length - 1
  end do

end function length

!
! Function that determines if a character is empty
!

function empty_char(char)

  implicit none
  logical :: empty_char
  character :: char
  empty_char = .false.
  if ( char == achar(9) .or. &
       char == achar(32) .or. &
       char == '' ) then
    empty_char = .true.
  end if 

end function empty_char

!
! Function that checks if a line is a comment line
!

function comment(string)

  implicit none
  integer :: i
  logical :: comment, empty_char
  character(len=200) :: string
  i = 1
  do while( empty_char(string(i:i)) .and. i < 200 ) 
    i = i + 1
  end do
  comment = .false.
  if ( string(i:i) == "#" .or. i == 200 ) comment = .true.

end function comment

!
! Subroutine that writes a beatiful progress bar
!

subroutine progress(current,start,end)

  integer :: current, start, end

  if ( current == start ) then
    write(*,"('# Progress: ',i10,' of ', i10$)") start, end
    return
  end if
  if ( current > start .and. current < end ) then
    write(*,"(24a,$)") (achar(8),i=1,24)
    write(*,"(i10,' of ',i10,$)") current, end
  end if
  if ( current == end ) then
    write(*,"(24a,$)") (achar(8),i=1,24)
    write(*,"(i10,' of ',i10)") current, end
  end if
 
end subroutine progress




