!
! G-score correlation calculator
!
! Program xgcorrelation
!
! This program writes the similarity of the models to a reference,
! as function of their g-scores, the reference being a model that
! does not belong to the original modeling set (for instance, the
! crystallographic model). Therefore, the alignment of the models
! to the reference is provided as one of the input files.
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 26, 2016
! http://leandro.iqm.unicamp.br
!
program xgcorrelation

  use types
  use file_operations
  implicit none
  integer :: i, i1, imodel
  integer :: narg, ioerr, nmodels, model_index
  double precision :: gscore, dummy, gdt_read, tmscore_read
  double precision, allocatable :: gdt(:), tmscore(:)
  character(len=200) :: alignlog, gscorefile, record, name
  character(len=200) :: gdt_output, tm_output, file1, file2, pdblist
  logical :: error
  type(model_type), allocatable :: model(:)

  narg = iargc()
  if ( narg /= 4 ) then
    write(*,*) ' ERROR: Run with: ./gcorrelation [reference align log] [pdb list] [gscore output] [output]'
    stop
  end if
!voltar:
!
! Needs that gscore and xgscore inform on header which kind of similarity was used (GDT or TM-score)
! Then, this program will read that from that file and use gdt or tmscore accordingly.
! In this case, there is no need for gdt(:) and tmscore(:), just read the model%similarity
! on the correct column
!
  call getarg(1,alignlog)
  call getarg(2,gscorefile)
  call getarg(3,pdblist)
  call getarg(4,gdt_output)
  call getarg(5,tm_output)

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
  write(*,"(a,a)") "# Log of alignment to reference: ", trim(adjustl(alignlog)) 
  write(*,"(a,a)") "# G-score data file: ", trim(adjustl(gscorefile)) 
  write(*,"(a,a)") "# Output file for GDT scores: ", trim(adjustl(gdt_output)) 
  write(*,"(a,a)") "# Output file for TM-scores: ", trim(adjustl(tm_output)) 
  write(*,"(a)") "#" 

  !
  ! Count the number of models in the PDB list
  !

  open(10,file=pdblist,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(pdblist))
    stop
  end if
  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    nmodels = nmodels + 1
  end do

  write(*,"(a,i10)") '# Number of PDB files in list: ', nmodels
  allocate(gdt(nmodels),tmscore(nmodels),model(nmodels))
  rewind(10)

  !
  ! Assign an index to each model name
  !

  write(*,"(a)") "# Reading list of model files ... "
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) file1
    if ( ioerr /= 0 ) exit
    if ( comment(file1) ) cycle
    imodel = imodel + 1
    model(imodel)%name = basename(file1)
  end do
  close(10)

  !
  ! Sort model names according to string comparisons
  !

  call sort_by_name(nmodels,model)

  ! Open the align log file

  open(10,file=alignlog,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open alignment log file: ', trim(adjustl(alignlog))
    stop
  end if
  i = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) file1, file2, tmscore_read, (dummy,i=1,4), gdt_read
    i = i + 1
    call progress(i,1,nmodels)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read data in alignment log file: ', trim(adjustl(alignlog))
      write(*,*) '        Content: ', trim(adjustl(record))
      stop
    end if
    file1 = basename(file1)
    i1 = model_index(file1,model,nmodels,error)
    if ( error ) cycle
    tmscore(i1) = tmscore_read
    gdt(i1) = gdt_read
  end do
  close(10)

  ! Open the gscore output and read the scores

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

  !
  ! Sort models from greater to lower GDT
  !

  do i = 1, nmodels
    model(i)%similarity = gdt(i)
  end do
  call sort_by_similarity(nmodels,model)
  
  ! Write GDT output file

  write(*,"(a)") "# Writing GDT output file ... "
  open(10,file=gdt_output,action='write')
  write(10,"(a)") "# Output of xgcorrelation (GDT scores)"
  write(10,"(a,a)") "# Input alignment log: ", trim(adjustl(alignlog))
  write(10,"(a,a)") "# Input gscore file: ", trim(adjustl(gscorefile))
  write(10,"(a,a)") "# List of PDB models: ", trim(adjustl(pdblist))
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score    Similarity  Model"
  do i = 1, nmodels
    write(10,"(f12.5,tr2,f12.5,tr2,a)") model(i)%gscore, model(i)%similarity, &
                                        trim(adjustl(model(i)%name))
  end do
  close(10)
  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(gdt_output))
  write(*,"(a)") "#"

  !
  ! Sort models from greater to lower TM-score
  !

  do i = 1, nmodels
    model(i)%similarity = tmscore(i)
  end do
  call sort_by_similarity(nmodels,model)
  
  ! Write TM-score output file

  write(*,"(a)") "# Writing TM-score output file ... "
  open(10,file=tm_output,action='write')
  write(10,"(a)") "# Output of xgcorrelation (TM-scores)"
  write(10,"(a,a)") "# Input alignment log: ", trim(adjustl(alignlog))
  write(10,"(a,a)") "# Input gscore file: ", trim(adjustl(gscorefile))
  write(10,"(a,a)") "# List of PDB models: ", trim(adjustl(pdblist))
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score    Similarity  Model"
  do i = 1, nmodels
    write(10,"(f12.5,tr2,f12.5,tr2,a)") model(i)%gscore, model(i)%similarity, &
                                        trim(adjustl(model(i)%name))
  end do
  close(10)

  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(tm_output))
  write(*,"(a)") "#"
  write(*,"(a)") "# Finished. " 

end program xgcorrelation

