!
! Subroutine that writes a beatiful progress indicator
!

subroutine progress(current,start,end)

  implicit none
  integer :: i, factor
  integer :: current, start, end
  real :: p

  ! If first step, print title

  if ( current == start ) then
    write(*,"('# Progress:   0.0%'$)")
  end if

  ! If last step, print and return

  if ( current == end ) then
    write(*,"(6a,'100.0%')") (achar(8),i=1,6)
    return
  end if

  ! Only display steps of 0.1%

  factor = max((end - start + 1)/1000,1)
  if ( ( mod(current-start,factor) /= 0 ) ) return

  p = 100.*current/(end-start+1)
  write(*,"(6a,f5.1,'%'$)") (achar(8),i=1,6), p
 
end subroutine progress

!
! Progress bar using number instead of percentages
!

subroutine progress2(current,start,end)

  implicit none
  integer :: i, factor
  integer :: current, start, end

  ! If first step, print title

  if ( current == start ) then
    write(*,"('# Progress: ',i10,' of ', i10$)") start, end
  end if

  ! If last step, print and return

  if ( current == end ) then
    write(*,"(24a,$)") (achar(8),i=1,24)
    write(*,"(i10,' of ',i10)") current, end
    return
  end if

  ! Only display steps of 0.01%

  factor = max((end - start + 1)/1000,1)
  if ( ( mod(current-start,factor) /= 0 ) ) return

  write(*,"(24a,$)") (achar(8),i=1,24)
  write(*,"(i10,' of ',i10,$)") current, end
 
end subroutine progress2


