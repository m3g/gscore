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

module wparameters

  double precision :: wpars(4)

end module wparameters

program gscore

  use types
  use file_operations
  use compactlog_data
  use wparameters
  implicit none
  integer :: i, j 
  integer :: narg, ioerr
  double precision :: cutoff, pcontact
  character(len=200) :: record, output, normtype
  double precision :: sumwdegree

  write(*,"(a)") "#" 
  write(*,"(a)") "# G-score calculator " 
  call title()
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 

  narg = iargc()
  if ( narg < 3 ) then
    write(*,*) ' Run with: ./gscore [compact log] [score cut] [output file] [norm type]'
    write(*,*) '    Where: [compact log] is the output of compactlog or qmatrix '
    write(*,*) '           [score cut] is the similarity cutoff for the score '
    write(*,*) '           [output file] is the name of the output file. '
    write(*,*) '           [norm type] is the normalization to be used for the contact score.'
    write(*,*) '                       options: none, maxcontacts, ncontacts, ijmax (default: maxcontacts)'
    stop
  end if
  call getarg(1,compactlog)
  call getarg(2,record)
  read(record,*,iostat=ioerr) cutoff
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read score cut from second command line argument. '
    stop
  end if
  call getarg(3,output)
  normtype = "maxcontacts"
  if ( narg == 4 ) then
    call getarg(4,normtype)
  end if
  if ( narg > 4 ) then
    call getarg(5,record)
    read(record,*) wpars(1)
    call getarg(6,record)
    read(record,*) wpars(2)
    call getarg(7,record)
    read(record,*) wpars(3)
    call getarg(8,record)
    read(record,*) wpars(4)
  end if

  ! Print the input options

  write(*,"(a)") "# Reference:" 
  write(*,"(a)") "# L. Martinez, A. Ferrari, F. C. Gozzo," 
  write(*,"(a)") "# A model evaluation score for ... 2016" 
  write(*,"(a)") "#" 
  write(*,"(a,a)") "# Alignment log file (compact form): ", trim(adjustl(compactlog)) 
  write(*,"(a)") "#" 
  write(*,"(a,f12.5)") "# Score similarity cutoff: ", cutoff
  write(*,"(a)") "#" 

  ! Read compact log file

  call read_compactlog(10)

  ! Norm type must only be used for contact score

  if ( score_type == 1 ) then
    write(*,"(a)") "# Score type: GDT_TS"  
  end if
  if ( score_type == 2 ) then
    write(*,"(a)") "# Score type: TM-score"  
  end if
  if ( score_type == 3 ) then
    write(*,"(a)") "# Score type: Contact-correlation"  
    write(*,"(a,a)") "# Normalization used: ", trim(adjustl(normtype))
  end if
  if ( score_type /= 3 .and. narg == 4 ) then
    write(*,*) ' ERROR: Normalization options (4th argument) is for contact scores only. '
    stop
  end if

  ! Compute G-score for all models

  write(*,"(a)") "# Computing the G-scores ... "
  do i = 1, nmodels
    model(i)%gscore = 0.d0
    model(i)%degree = 0.d0
    model(i)%wdegree = 0.d0
  end do

  ! If the score is a contact score and normalization is maxcontact or ijmax, normalize

  if ( score_type == 3 .and. ( normtype == "maxcontacts" .or. normtype == "ijmax" ) ) then
    if ( normtype == "maxcontacts" ) then
      write(*,"(a,i10)") "# Normalizing contact score by maxcontacts: ", maxcontacts
    else
      write(*,"(a)") "# Normalizing contact score by ijmax."
    end if
    do i = 1, nmodels-1
      call progress(i,1,nmodels)
      do j = i+1, nmodels
        if ( normtype == "maxcontacts" ) then 
          scores(i,j) = scores(i,j) / maxcontacts
        else if ( normtype == "ijmax" ) then 
          if ( model(i)%ncontacts >= model(j)%ncontacts ) then
            scores(i,j) = scores(i,j) / model(i)%ncontacts
          else
            scores(i,j) = scores(i,j) / model(j)%ncontacts
          end if
        end if
      end do
    end do
    call progress(nmodels,1,nmodels)
  end if

  ! For alignment scores or contact score without model-dependent norm

  if ( score_type == 1 .or. &
       score_type == 2 .or. &
       ( score_type == 3 .and. normtype /= "ncontacts" ) ) then
    do i = 1, nmodels-1
      call progress(i,1,nmodels)
      do j = i+1, nmodels
        model(i)%wdegree = model(i)%wdegree + sumwdegree(scores(i,j))
        model(j)%wdegree = model(j)%wdegree + sumwdegree(scores(i,j)) 
        if ( scores(i,j) >= cutoff ) then
          model(i)%degree = model(i)%degree + 1.d0
          model(j)%degree = model(i)%degree + 1.d0
        end if
      end do
    end do
    call progress(nmodels,1,nmodels)
  end if

  ! If the normalization of contact scores is the number of
  ! contacts of each model, it is not symmetric 

  if ( score_type == 3 .and. normtype == "ncontacts" ) then 
    do i = 1, nmodels-1
      call progress(i,1,nmodels)
      do j = i + 1, nmodels
        model(i)%wdegree = model(i)%wdegree + sumwdegree(scores(i,j))
        model(j)%wdegree = model(j)%wdegree + sumwdegree(scores(i,j)) 
        pcontact = scores(i,j) / model(i)%ncontacts
        if ( pcontact >= cutoff ) then
          model(i)%degree = model(i)%degree + 1.d0
        end if
        pcontact = scores(i,j) / model(j)%ncontacts
        if ( pcontact >= cutoff ) then
          model(j)%degree = model(j)%degree + 1.d0
        end if
      end do
    end do
    call progress(nmodels,1,nmodels)
  end if

  do i = 1, nmodels
    model(i)%degree = model(i)%degree / dble(nmodels-1)
    model(i)%gscore = -0.593d0*dlog(model(i)%degree+1.d-30)
    model(i)%wdegree = model(i)%wdegree/dble(nmodels-1)
  end do

  ! Order models from greater to lower G-scores

  !call sort_by_gscore(nmodels,model)
  call sort_by_wdegree(nmodels,model)

  !
  ! Write output file 
  !

  write(*,"(a)") "# Writing output file ... "
  open(10,file=output,action='write')
  write(10,"(a)") "# Output of G-score"
  write(10,"(a,a)") "# Input alignment log: ", trim(adjustl(compactlog))
  if ( score_type == 1 ) then
    write(10,"(a)") "# Score type: GDT_TS"
  end if
  if ( score_type == 2 ) then
    write(10,"(a)") "# Score type: TM-score"
  end if
  if ( score_type == 3 ) then
    write(10,"(a)") "# Score type: Contact-correlation"
  end if
  write(10,"(a,f12.5)") "# Score cutoff: ", cutoff
  if ( score_type == 3 ) then
    write(10,"(a,a)") "# Contact score normalization: ", trim(adjustl(normtype))
  end if
  write(10,"(a)") "#"
  write(10,"(a)") "# G-score is -RTln(P) for RT=0.593 kcal/mol (T=298.15)"
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score     Degree(P)       WDegree  Model"
  do i = 1, nmodels
    write(10,"(3(f12.5,tr2),a)") model(i)%gscore, model(i)%degree, model(i)%wdegree, trim(adjustl(model(i)%name))
  end do
  close(10)

  write(*,"(a)") "#"
  write(*,"(a,a)") "# Wrote file: ", trim(adjustl(output))
  write(*,"(a)") "#"
  write(*,"(a)") "# Finished. " 

end program gscore

function sumwdegree(x)

  use wparameters
  implicit none
  double precision :: sumwdegree, x
  double precision :: a0, a1, a2

  !sumwdegree = wpars(1)*(1/(x**wpars(2))-1.d0) + wpars(3)*(1.d0-x**(wpars(4)))

  !sumwdegree = 1.d0 - 1.d0 / (1.d0+dexp(-1.d0*wpars(1)*(x-wpars(2))))
  !sumwdegree = 1.d0/x -1.d0
  !sumwdegree = x**4
  !sumwdegree = (1/x - 1) + 8*(1-x**12)
  !sumwdegree = (1/x - 1)
  !sumwdegree = 1-x
  !sumwdegree = 1.d0 / ( 1.d0/min(x,0.99d0) - 1.d0 )
  sumwdegree = x

end function sumwdegree



