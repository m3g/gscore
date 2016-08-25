module types

  type structure 
    character(len=200) :: name
    integer :: index
    double precision :: pas1_gdt, pas2_gdt, pbias_gdt
    double precision :: pas1_tm, pas2_tm, pbias_tm
  end type structure

end module types

program pbias

  use types
  implicit none
  integer :: i1, i2, i, j, structure_index, imodel, itarget
  integer :: narg, ioerr, nlist1, nlist2, ntarget1, ntarget2
  double precision :: dummy, tmscore, gdt
  double precision, allocatable :: gdt1(:,:), gdt2(:,:), tmscore1(:,:), tmscore2(:,:)
  character(len=200) :: align_list1, align_list2, align_log, firstlog1, firstlog2
  character(len=200) :: basename, record, file1, file2
  logical :: comment, repeated
  type(structure), allocatable :: models1(:), models2(:), targets1(:), targets2(:)

  narg = iargc()
  if ( narg /= 4 ) then
    write(*,*) ' ERROR: Run with: ./pas [output_file] [GDT cut] [TM-Score cut] [align_list1] '
    stop
  end if

  ! Open first log list

  call getarg(2,align_list1)
  open(10,file=align_list1,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(align_list1))
    stop
  end if

  ! Count the number of files

  nlist1 = 0
  do
    read(10,"(a200)",iostat=ioerr) align_log
    if ( ioerr /= 0 ) exit
    if ( comment(align_log) ) cycle 
    nlist1 = nlist1 + 1
    
    ! Read number of alignments in this file (is the number of structures
    ! of the target structure database)

    if ( nlist1 == 1 ) then
      firstlog1 = align_log
      open(20,file=align_log,action='read',status='old',iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not open file: ', trim(adjustl(align_log))
      end if
      ntarget1 = 0
      do
        read(20,"(a200)",iostat=ioerr) record
        if ( ioerr /= 0 ) exit
        if ( comment(record) ) cycle
        ntarget1 = ntarget1 + 1
      end do
      close(20)
    end if
  end do
  write(*,"(a,i10)") '# Number of lovoalign log files in list 1: ', nlist1
  write(*,"(a,i10)") '# Number of alignments in files of list 1: ', ntarget1

  allocate(gdt1(nlist1,ntarget1),tmscore1(nlist1,ntarget1),models1(nlist1),targets1(ntarget1))

  !
  ! Assign an index to each structure name
  !

  open(20,file=firstlog1,action='read',status='old')
  imodel = 0
  itarget = 0
  do
    read(20,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*) file1, file2
    !
    ! Index alignment models
    !
    repeated = .false.
    do i = 1, imodel
      if ( basename(file1) == models1(i)%name ) then
        repeated = .true.
        exit
      end if
    end do
    if ( .not. repeated ) then
      imodel = imodel + 1
      models1(imodel)%name = basename(file1)
    end if
    !
    ! Index alignment targets
    !
    repeated = .false.
    do i = 1, itarget
      if ( basename(file2) == targets1(i)%name ) then
        repeated = .true.
        exit
      end if
    end do
    if ( .not. repeated ) then
      itarget = itarget + 1
      targets1(itarget)%name = basename(file2)
    end if
  end do
  close(20)

  !
  ! Now, reading all alignment log files and annotating the scores
  ! of the alignment of each pair
  !
  rewind(10)
  do
    read(10,"(a200)",iostat=ioerr) align_log
    if ( ioerr /= 0 ) exit
    if ( comment(align_log) ) cycle
    open(20,file=align_log,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open alignment log file: ', trim(adjustl(align_log))
      write(*,*) '        Listed in: ', trim(adjustl(align_list1))
      stop
    end if
    do
      read(20,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( comment(record) ) cycle
      read(record,*,iostat=ioerr) file1, file2, tmscore, (dummy,i=1,4), gdt
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read data in alignment log file: ', trim(adjustl(align_log))
        write(*,*) '        Content: ', trim(adjustl(record))
        stop
      end if
      i1 = structure_index(file1,models1,nlist1)
      i2 = structure_index(file2,targets1,ntarget1)
      tmscore1(i1,i2) = tmscore
      gdt1(i1,i2) = gdt
      tmscore1(i2,i1) = tmscore
      gdt1(i2,i1) = gdt
    end do
    close(20)
  end do
  close(10)

  !
  ! End of file reading, now will compute the scores
  !

  ! Compute the PAS score for models of list1

  do i = 1, nlist1
    models1(i)%pas1_gdt = 0.d0
    models1(i)%pas1_tm = 0.d0
    do j = 1, ntargets1
      if ( gdt1(i,j) > gdt_cut ) then 
        models(i)%pas1_gdt = models(i)%pas1_gdt + 1.d0
      end if
      if ( tmscore1(i,j) > tm_cut ) then 
        models(i)%pas1_tm = models(i)%pas1_tm + 1.d0
      end if
    end do
    models1(i)%pas1_gdt = models1(i)%pas1_gdt / dble(nlist)
    models1(i)%pas1_tm = models1(i)%pas1_tm / dble(nlist)
  end do

  !
  ! Order models from best to worst accoding to PAS_GDT
  ! 
voltar

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

function structure_index(name,models,n)
 
  use types
  implicit none
  integer :: structure_index, n
  character(len=200) :: name
  type(structure) :: models(n)
      
  structure_index = 0
  do while( name /= models(structure_index)%name ) 
    structure_index = structure_index + 1
  end do

end function structure_index

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
  do while(basename(i:i) /= "/" .and. i > 0)
    if ( basename(i:i) == "." ) then
      idot = i
    end if
    i = i - 1
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




