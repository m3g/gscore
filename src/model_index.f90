!
! Function that determines the index of a model from its name
!

function model_index(name,model,n,stop)
 
  use types
  implicit none
  integer :: model_index, n, imax, imin, iavg
  character(len=200) :: name
  logical :: stop, error
  type(model_type) :: model(n)

  imin = 1
  imax = n
  error = .false.
  if ( name < model(1)%name ) error = .true.
  if ( name > model(n)%name ) error = .true.
  if ( .not. error ) then
    do
      if ( name == model(imin)%name ) then 
        model_index = imin
        return
      end if
      if ( name == model(imax)%name ) then
        model_index = imax
        return
      end if
      if ( (imax - imin) > 1 ) then
        iavg = imin + ( imax - imin ) / 2
        if ( name == model(iavg)%name ) then
          model_index = iavg
          return
        end if
      else
        error = .true.
        exit
      end if
      if ( name > model(iavg)%name ) imin = iavg + 1
      if ( name < model(iavg)%name ) imax = iavg - 1
    end do
  end if
  if ( error ) then
    if ( stop ) then
      write(*,*)
      write(*,*) ' ERROR: A model is listed in a log file but was not found in list: '
      write(*,*) '        Model: ', trim(adjustl(name))
      stop
    else
      ! Return stop true to outer program, to report error
      stop = .true.
    end if
  end if

end function model_index
