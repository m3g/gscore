!
! Program xcompare
!
! Reads two ASCII tables with file model names and some scores,
! and writes a file comparing the scores chosen for each model
!
! L. Martinez
! Institute of Chemistry, University of Campinas
! http://leandro.iqm.unicamp.br
! Nov 16, 2016
!
program xcompare

  use types
  use file_operations
  implicit none
  integer :: i, narg, ioerr, ifile, model_index
  integer :: name(2), col(2), nmodels(2), imodel, jmodel
  real :: value
  character(len=200) :: file(2), record, string, output, modelname
  logical :: stop
  type(model_type), allocatable :: model1(:), model2(:)


  write(*,"(a)") "#" 
  write(*,"(a)") "# Xcompare: Compare scores for lists in different files" 
  call title()
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 

  narg = iargc()
  if ( narg /= 7 ) then
    write(*,*) ' ERROR: Run with: ./xcompare [file1] [name1] [col1] [file2] [name2] [col2] [output] '
    write(*,*) ' Where: [file1] and [file2] are the files containing the data. '
    write(*,*) '        [name1] and [col1] are the columns containing model names '
    write(*,*) '           and scores in file1'
    write(*,*) '        [name2] and [col2] are the columns containing model names '
    write(*,*) '           and scores in file2'
    write(*,*) '        [output] is the name of the output file. '
    stop
  end if
  call getarg(1,file(1))
  call getarg(2,record)
  read(record,*,iostat=ioerr) name(1)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read column containing names from 2nd argument. '
    stop
  end if
  call getarg(3,record)
  read(record,*,iostat=ioerr) col(1)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read column containing scores from 3rd argument. '
    stop
  end if
  call getarg(4,file(2))
  call getarg(5,record)
  read(record,*,iostat=ioerr) name(2)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read column containing names from 5th argument. '
    stop
  end if
  call getarg(6,record)
  read(record,*,iostat=ioerr) col(2)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read column containing scores from 6th argument. '
    stop
  end if
  call getarg(7,output)

  ! Read number of models listed in the file1

  do ifile = 1, 2
    open(10,file=file(ifile),status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file(ifile)))
    end if
    nmodels(ifile) = 0
    do
      read(10,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( comment(record) ) cycle
      read(record,*,iostat=ioerr) (modelname, i = 1, name(ifile))
      if ( ioerr /= 0 ) cycle
      read(record,*,iostat=ioerr) (string, i = 1, col(ifile))
      if ( ioerr /= 0 ) cycle
      read(string,*,iostat=ioerr) value 
      if ( ioerr /= 0 ) cycle
      nmodels(ifile) = nmodels(ifile) + 1
    end do
    close(10)
    write(*,"(a,i2,a,i10)") '# Number of models read in the file ', ifile,':', nmodels(ifile)
  end do

  ! Read file names and scores

  allocate(model1(nmodels(1)),model2(nmodels(2)))

  open(10,file=file(1),status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file(1)))
  end if
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (modelname, i = 1, name(1))
    if ( ioerr /= 0 ) cycle
    read(record,*,iostat=ioerr) (string, i = 1, col(1))
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) value 
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model1(imodel)%name = basename(modelname)
    model1(imodel)%gscore = value
  end do
  close(10)

  open(10,file=file(2),status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file(2)))
  end if
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (modelname, i = 1, name(2))
    if ( ioerr /= 0 ) cycle
    read(record,*,iostat=ioerr) (string, i = 1, col(2))
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) value 
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model2(imodel)%name = basename(modelname)
    model2(imodel)%gscore = value
  end do
  close(10)

  ! Sort models by name in both sets
  write(*,"(a)") "# Sorting models of first list by name ... "
  call sort_by_name(nmodels(1),model1)
  write(*,"(a)") "# Sorting models of second list by name ... "
  call sort_by_name(nmodels(2),model2)

  ! Write output

  write(*,"(a,a)") "# Writting output file: ", trim(adjustl(output))
  open(10,file=output)
  do imodel = 1, nmodels(1)
    stop = .false.
    jmodel = model_index(model1(imodel)%name,model2,nmodels(2),stop,.true.)
    if ( stop ) cycle
    write(10,*) model1(imodel)%gscore, model2(jmodel)%gscore, trim(adjustl(model1(imodel)%name)) 
  end do
  close(10)
  write(*,"(a)") "#"
  write(*,"(a)") "# END"

end program xcompare

