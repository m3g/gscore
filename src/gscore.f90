!
! G-score calculator
!
! Program gscore
!
! This program computes the G score given a compact alignment log
! generated with compactlog.
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 26, 2016
! http://leandro.iqm.unicamp.br
!
program gscore

  use types
  use file_operations
  implicit none
  integer :: i, j 
  integer :: narg, ioerr, nmodels
  double precision :: scorecut
  double precision, allocatable :: scores(:,:)
  character(len=200) :: compactlog, record, output
  type(model_type), allocatable :: model(:)

  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with: ./gscore [compact align log] [score cut] [output file]'
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

  ! Read model list from log file (first two lines are comments to be ignored)

  read(10,*)
  read(10,*)
  read(10,*) nmodels
  allocate(scores(nmodels,nmodels),model(nmodels))

  ! Read scores 
 
  write(*,"(a)") "# Reading scores from file ... "
  do i = 1, nmodels
    read(10,*) j, model(i)%name
  end do
  do i = 1, nmodels-1
    call progress(i,1,nmodels)
    read(10,*) (scores(i,j),j=i+1,nmodels)
  end do
  close(10)
  call progress(nmodels,1,nmodels)

  ! Compute G-score for all models

  write(*,"(a)") "# Computing the G-scores ... "
  do i = 1, nmodels
    model(i)%gscore = 0.d0
  end do
  do i = 1, nmodels-1
    call progress(i,1,nmodels)
    do j = i+1, nmodels
      if ( scores(i,j) > scorecut ) then
        model(i)%gscore = model(i)%gscore + 1.d0
        model(j)%gscore = model(j)%gscore + 1.d0
      end if
    end do
  end do
  call progress(nmodels,1,nmodels)
  do i = 1, nmodels
    model(i)%gscore = model(i)%gscore / dble(nmodels)
  end do

  ! Order models from greater to lower G-scores

  call sort_by_gscore(nmodels,model)

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
  do i = 1, nmodels
    write(10,"(f12.5,tr2,a)") model(i)%gscore, trim(adjustl(model(i)%name))
  end do
  close(10)

  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(output))
  write(*,"(a)") "#"
  write(*,"(a)") "# Finished. " 

end program gscore




