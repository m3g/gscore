!
! G-score calculator
!
! Program xgscore
!
! This program computes the G score given a compact alignment log
! generated with xcompactlog, for the gscore computed for a set
! of models in relation to a set of models generated by another 
! alignment
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 26, 2016
! http://leandro.iqm.unicamp.br
!
program xgscore

  use types
  use file_operations
  implicit none
  integer :: i, j 
  integer :: narg, ioerr, nmodels1, nmodels2
  double precision :: scorecut
  double precision, allocatable :: scores(:,:)
  character(len=200) :: compactlog, record, output
  type(model_type), allocatable :: model1(:), model2(:)

  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with: ./xgscore [compact align log] [score cut] [output file]'
    stop
  end if
  call getarg(1,compactlog)
  call getarg(2,record)
  read(record,*,iostat=ioerr) scorecut
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read score cut from second command line argument. '
    stop
  end if
  call getarg(3,output)

  ! Print the input options

  write(*,"(a)") "#" 
  write(*,"(a)") "# G-score calculator " 
  write(*,"(a)") "# XGscore: Compute gscore of a set of models relative to another set." 
  write(*,"(a)") "#" 
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 
  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# Alignment log file (compact form): ", trim(adjustl(compactlog)) 
  write(*,"(a)") "#" 
  write(*,"(a,f12.5)") "# Score similarity cutoff: ", scorecut
  write(*,"(a)") "#" 

  ! Open the compact log file

  open(10,file=compactlog,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not find or open alignment log file: ', trim(adjustl(compactlog))
    stop
  end if

  ! Read model list from log file (first three lines are comments to be ignored)

  read(10,*)
  read(10,*)
  read(10,*)

  ! Number of models of first set
  read(10,*) nmodels1
  allocate(model1(nmodels1))

  ! Read model1 names
  do i = 1, nmodels1
    read(10,*) model1(i)%index, model1(i)%name
  end do

  ! Number of models of second set
  read(10,*) nmodels2
  allocate(model2(nmodels2))

  ! Read model2 names
  do i = 1, nmodels2
    read(10,*) model2(i)%index, model2(i)%name
  end do

  ! Reading scores
  write(*,"(a)") "# Reading scores from file ... "
  allocate(scores(nmodels1,nmodels2))
  do i = 1, nmodels1
    call progress(i,1,nmodels1)
    read(10,*) (scores(i,j),j=1,nmodels2)
  end do
  close(10)

  ! Compute G-score for all models of second list

  write(*,"(a)") "# Computing the G-scores ... "
  do i = 1, nmodels2
    model2(i)%gscore = 0.d0
  end do
  do i = 1, nmodels1
    call progress(i,1,nmodels1)
    do j = 1, nmodels2
      if ( scores(i,j) > scorecut ) then
        model2(j)%gscore = model2(j)%gscore + 1.d0
      end if
    end do
  end do
  call progress(nmodels2,1,nmodels2)
  do i = 1, nmodels2
    model2(i)%gscore = model2(i)%gscore / dble(nmodels2)
  end do

  ! Order models from greater to lower G-scores

  call sort_by_gscore(nmodels2,model2)

  !
  ! Write output file 
  !

  write(*,"(a)") "# Writing output file ... "
  open(10,file=output,action='write')
  write(10,"(a)") "# Output of G-score"
  write(10,"(a,a)") "# Input alignment log: ", trim(adjustl(compactlog))
  write(10,"(a,f12.5)") "# Score cutoff: ", scorecut
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score  Model"
  do i = 1, nmodels2
    write(10,"(f12.5,tr2,a)") model2(i)%gscore, trim(adjustl(model2(i)%name))
  end do
  close(10)

  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(output))
  write(*,"(a)") "#"
  write(*,"(a)") "# Finished. " 

end program xgscore



