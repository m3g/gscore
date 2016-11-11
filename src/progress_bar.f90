!
! Subroutine that writes a beatiful progress bar
!

subroutine progress(current,start,end)

  implicit none
  integer :: i, factor
  integer :: current, start, end

  ! Only display steps of 0.01%

  factor = (end - start)/10000
  if ( ( mod(current-start,factor) /= 0 ) ) return

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




