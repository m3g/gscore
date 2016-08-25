!
! PAS score calculator
!
! Program compactlog
!
! This program compacts the output of many multiple lovoalign
! runs into a single file with the data of GDT and TM-scores
! for all pairs, so that the PAS scores can be compued faster
! without having to read all the align logs again.
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 25, 2016
! http://leandro.iqm.unicamp.br
!

module types

  type model_type
    character(len=200) :: name
    double precision :: pas_gdt
    double precision :: pas_tm
  end type model_type

end module types

program compactlog

  use types
  implicit none
  integer :: i1, i2, i, j, model_index, imodel
  integer :: narg, ioerr, nmodels, nlogs, ilog
  double precision :: dummy, tmscore_read, gdt_read
  double precision, allocatable :: gdt(:,:), tmscore(:,:)
  character(len=200) :: align_list, firstlog, align_log, gdt_log, tm_log
  character(len=200) :: basename, record, file1, file2, format
  logical :: comment
  type(model_type), allocatable :: model(:)

  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with: ./compactlog [align_list] [gdt output] [tmscore output]  '
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

  call getarg(2,gdt_log)
  open(30,file=gdt_log,status='new',action='write',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Trying to create file: ', trim(adjustl(gdt_log))
    write(*,*) '        but file already exists, or some other access problem. '
    stop
  end if
  call getarg(3,tm_log)
  open(40,file=tm_log,status='new',action='write',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Trying to create file: ', trim(adjustl(tm_log))
    write(*,*) '        but file already exists, or some other access problem. '
    stop
  end if

  ! Print the input options

  write(*,"(a)") "#" 
  write(*,"(a)") "# PAS score calculator " 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Align log file conversion " 
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
  write(*,"(a,a)") "# Will create compact log for GDTs: ", trim(adjustl(gdt_log))
  write(*,"(a,a)") "# Will create compact log for TM-scores: ", trim(adjustl(tm_log))
  write(*,"(a)") "#" 

  ! Write titles to compact log files

  write(30,*) ' This a compact lovoalign alignment file, with GDT scores '
  write(30,*) ' Alignment files obtained from ', trim(adjustl(align_list))
  write(40,*) ' This a compact lovoalign alignment file, with TM-scores '
  write(40,*) ' Alignment files obtained from ', trim(adjustl(align_list))

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
    read(record,*) file1
    imodel = imodel + 1
    model(imodel)%name = basename(file1)
  end do
  !
  ! Sort model names according to string comparisons
  !
  do i = 1, nmodels-1
    j = i + 1
    do while( j > 1 .and. model(j-1)%name > model(j)%name )
      file1 = model(j-1)%name
      model(j-1)%name = model(j)%name
      model(j)%name = file1
      j = j - 1
    end do
  end do
  write(30,*) nmodels
  write(40,*) nmodels
  do imodel = 1, nmodels
    write(30,*) imodel, trim(adjustl(model(imodel)%name))
  end do

  !
  ! Now, reading all alignment log files and annotating the scores
  ! of the alignment of each pair
  !
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reading all alignment files ... this can take a while. "
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
      i1 = model_index(file1,model,nmodels)
      file2 = basename(file2)
      i2 = model_index(file2,model,nmodels)
      if ( i1 > i2 ) then
        tmscore(i1,i2) = tmscore_read
        gdt(i1,i2) = gdt_read
      else
        tmscore(i2,i1) = tmscore_read
        gdt(i2,i1) = gdt_read
      end if
    end do
    close(20)
    ilog = ilog + 1
    call progress(ilog,1,nlogs)
  end do
  close(10)

  !
  ! Write scores to files
  !

  write(format,*) nmodels
  format = "("//trim(adjustl(format))//"(tr1,f8.3))"
  do imodel = 1, nmodels - 1
    write(30,format) (gdt(imodel,i),i=imodel+1,nmodels)
    write(40,format) (tmscore(imodel,i),i=imodel+1,nmodels)
  end do
  close(30)
  close(40)

  write(*,"(a)") '# Finished. '

end program compactlog

!
! Function that determines the index of a model from its name
!

function model_index(name,model,n)
 
  use types
  implicit none
  integer :: model_index, n, imax, imin, iavg
  character(len=200) :: name
  type(model_type) :: model(n)

  imin = 1
  imax = n
  do
    if ( name == model(imin)%name ) then 
      model_index = imin
      return
    end if
    if ( name == model(imax)%name ) then
      model_index = imax
      return
    end if
    iavg = imin + ( imax - imin ) / 2
    if ( name >= model(iavg)%name ) imin = iavg
    if ( name <= model(iavg)%name ) imax = iavg
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




