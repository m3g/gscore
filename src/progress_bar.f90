!
! Subroutine that writes a beatiful progress bar
!

subroutine progress(current,start,end)

  implicit none
  integer :: i, factor
  integer :: current, start, end

  ! If last step, print and return

  if ( current == end ) then
    write(*,"(24a,$)") (achar(8),i=1,24)
    write(*,"(i10,' of ',i10)") current, end
    return
  end if

  ! If first step, print

  if ( current == start ) then
    write(*,"('# Progress: ',i10,' of ', i10$)") start, end
    return
  end if

  ! Only display steps of 0.01%

  factor = max((end - start)/10000,1)
  if ( ( mod(current-start,factor) /= 0 ) ) return

  write(*,"(24a,$)") (achar(8),i=1,24)
  write(*,"(i10,' of ',i10,$)") current, end
 
end subroutine progress




