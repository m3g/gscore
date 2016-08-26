!
! Subroutines to sort models
!

!
! Sort by names
!

subroutine sort_by_name(n,model)

  use types
  implicit none
  integer :: i, j, n
  type(model_type) :: model(n), modeltemp

  write(*,"(a)") "# Sorting models by file name ... "
  do i = 1, n-1
    call progress(i,1,n)
    j = i + 1
    do while( model(j-1)%name > model(j)%name )
      modeltemp = model(j-1)
      model(j-1) = model(j) 
      model(j) = modeltemp
      j = j - 1
      if ( j == 1 ) exit
    end do
  end do
  call progress(n,1,n)

end subroutine sort_by_name

!
! Order models from greater to lower G-scores
!

subroutine sort_by_gscore(n,model)

  use types
  implicit none
  integer :: i, j, n
  type(model_type) :: model(n), modeltemp

  write(*,"(a)") "# Sorting models by G-score ... "
  do i = 1, n-1
    call progress(i,1,n)
    j = i + 1
    do while( model(j-1)%gscore < model(j)%gscore )
      modeltemp = model(j-1)
      model(j-1) = model(j)
      model(j) = modeltemp
      j = j - 1
      if ( j == 1 ) exit
    end do
  end do
  call progress(n,1,n)

end subroutine sort_by_gscore

!
! Sort models from greater to lower similarity
!

subroutine sort_by_similarity(n,model)

  use types
  implicit none
  integer :: i, j, n
  type(model_type) :: model(n), modeltemp

  write(*,"(a)") "# Sorting models by similarity to reference ... "
  do i = 1, n-1
    call progress(i,1,n)
    j = i + 1
    do while( model(j-1)%similarity < model(j)%similarity )
      modeltemp = model(j-1)
      model(j-1) = model(j)
      model(j) = modeltemp
      j = j - 1
      if ( j == 1 ) exit
    end do
  end do
  call progress(n,1,n)

end subroutine sort_by_similarity




