!
! G-score calculator
!
! Program xcompactlog
!
! This program compacts the output of many multiple lovoalign
! runs performed for sets of disjoint sets of models 
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 25, 2016
! http://leandro.iqm.unicamp.br
!

program xcompactlog

  use types
  use file_operations
  implicit none
  integer :: i1, i2, i, model_index
  integer :: imodel1, imodel2, nmodels1, nmodels2
  integer :: narg, ioerr, nlogs, ilog
  double precision :: dummy, tmscore_read, gdt_read
  double precision, allocatable :: gdt(:,:), tmscore(:,:)
  character(len=200) :: pdb_list1, pdb_list2
  character(len=200) :: align_list, align_log, gdt_log, tm_log, output
  character(len=200) :: record, file1, file2, format
  logical :: error
  type(model_type), allocatable :: model1(:), model2(:)

  narg = iargc()
  if ( narg /= 4 ) then
    write(*,*) ' ERROR: Run with: ./xcompactlog [pdb list 1] [pdb list 2] [align list] [output]  '
    stop
  end if

  ! Read PDB file lists from comand line

  call getarg(1,pdb_list1)
  call getarg(2,pdb_list2)

  ! Read align log file list from command line

  call getarg(3,align_list)

  ! Names of output files for GDT_TS and TM-score files.

  call getarg(4,output)
  gdt_log = trim(adjustl(basename(output)))//"-GDT_TS"&
            //trim(adjustl(output(length(basename(output))+1:length(output))))
  tm_log = trim(adjustl(basename(output)))//"-TMscore"&
           //trim(adjustl(output(length(basename(output))+1:length(output))))

  call checkfile(gdt_log)
  call checkfile(tm_log)

  ! Print the input options

  write(*,"(a)") "#" 
  write(*,"(a)") "# G-score calculator " 
  write(*,"(a)") "#" 
  write(*,"(a)") "# XCOMPACTLOG: Align log file conversion to compact form, from " 
  write(*,"(a)") "#              alignments of disjoint sets of structures. "
  write(*,"(a)") "#" 
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# First list of model PDB files: ", trim(adjustl(pdb_list1)) 
  write(*,"(a,a)") "# Second list of model PDB files: ", trim(adjustl(pdb_list2)) 
  write(*,"(a,a)") "# List of alignment files: ", trim(adjustl(align_list)) 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# Will create compact log for GDT_TS: ", trim(adjustl(gdt_log))
  write(*,"(a,a)") "# Will create compact log for TM-score: ", trim(adjustl(tm_log))
  write(*,"(a)") "#" 

  !
  ! Count the number of models in the first PDB list
  !

  open(10,file=pdb_list1,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(pdb_list1))
    stop
  end if
  nmodels1 = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle 
    nmodels1 = nmodels1 + 1
  end do

  write(*,"(a,i10)") '# Number of PDB files in first list: ', nmodels1
  allocate(model1(nmodels1))
  rewind(10)

  !
  ! Assign an index to each model name
  !

  write(*,"(a)") "# Reading list of model files ... "
  imodel1 = 0
  do
    read(10,"(a200)",iostat=ioerr) file1
    if ( ioerr /= 0 ) exit
    if ( comment(file1) ) cycle
    imodel1 = imodel1 + 1
    model1(imodel1)%name = basename(file1)
  end do
  close(10)

  !
  ! Count the number of models in the second PDB list
  !

  open(10,file=pdb_list2,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(pdb_list2))
    stop
  end if
  nmodels2 = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle 
    nmodels2 = nmodels2 + 1
  end do

  write(*,"(a,i10)") '# Number of PDB files in second list: ', nmodels2
  allocate(model2(nmodels2))
  rewind(10)

  !
  ! Assign an index to each model name
  !

  write(*,"(a)") "# Reading second list of PDB files ... "
  imodel2 = 0
  do
    read(10,"(a200)",iostat=ioerr) file1
    if ( ioerr /= 0 ) exit
    if ( comment(file1) ) cycle
    imodel2 = imodel2 + 1
    model2(imodel2)%name = basename(file1)
  end do
  close(10)

  !
  ! Sort model names according to string comparisons
  !

  call sort_by_name(nmodels1,model1)
  call sort_by_name(nmodels2,model2)

  ! Allocate arrays that will contain the scores read

  allocate(gdt(nmodels1,nmodels2),tmscore(nmodels1,nmodels2))
  
  ! Count the number of lovoalign log files

  open(10,file=align_list,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(align_list))
    stop
  end if
  nlogs = 0
  do
    read(10,"(a200)",iostat=ioerr) align_log
    if ( ioerr /= 0 ) exit
    if ( comment(align_log) ) cycle 
    nlogs = nlogs + 1
  end do
  rewind(10)
  write(*,"(a,i10)") '# Number of lovoalign log files in list: ', nlogs
  write(*,"(a)") "#" 

  !
  ! Now, reading all alignment log files and annotating the scores
  ! of the alignment of each pair
  !
  write(*,"(a)") "# Reading all alignment files ... this can take a while. "
  ilog = 0
  do
    call progress(ilog,1,nlogs)
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
      i1 = model_index(file1,model1,nmodels1,error)
      if ( error ) cycle
      file2 = basename(file2)
      i2 = model_index(file2,model2,nmodels2,error)
      if ( error ) cycle
      tmscore(i1,i2) = tmscore_read
      gdt(i1,i2) = gdt_read
    end do
    close(20)
    ilog = ilog + 1
  end do
  close(10)

  ! Create output files

  write(*,"(a)") "#"
  write(*,"(a)") "# Writing output files ... "

  ! Write titles to compact log files

  open(10,file=gdt_log,action='write',iostat=ioerr)
  open(20,file=tm_log,action='write',iostat=ioerr)

  write(10,"(a)") '# This a compact lovoalign alignment file, with GDT_TS scores '
  write(10,"(a,a)") '# Alignment files obtained from ', trim(adjustl(align_list))
  write(10,"(a,2(tr1,a))") '# With PDB lists: ', trim(adjustl(pdb_list1)), trim(adjustl(pdb_list2))
  write(20,"(a)") ' This a compact lovoalign alignment file, with TM-scores '
  write(20,"(a,a)") ' Alignment files obtained from ', trim(adjustl(align_list))
  write(20,"(a,2(tr1,a))") ' With PDB lists: ', trim(adjustl(pdb_list1)), trim(adjustl(pdb_list2))

  ! Write list of models to output files

  write(10,*) nmodels1
  write(20,*) nmodels1
  do imodel1 = 1, nmodels1
    write(10,*) imodel1, trim(adjustl(model1(imodel1)%name))
    write(20,*) imodel1, trim(adjustl(model1(imodel1)%name))
  end do
  write(10,*) nmodels2
  write(20,*) nmodels2
  do imodel2  = 1, nmodels2
    write(10,*) imodel2, trim(adjustl(model2(imodel2)%name))
    write(20,*) imodel2, trim(adjustl(model2(imodel2)%name))
  end do

  ! Write scores to files

  write(format,*) nmodels2
  format = "("//trim(adjustl(format))//"(tr1,f8.3))"
  do imodel1 = 1, nmodels1
    write(10,format) (gdt(imodel1,imodel2),imodel2=1,nmodels2)
    write(20,format) (tmscore(imodel1,imodel2),imodel2=1,nmodels2)
  end do
  close(10)
  close(20)

  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(gdt_log))
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(tm_log))
  write(*,"(a)") "#"
  write(*,"(a)") '# Finished. '

end program xcompactlog
