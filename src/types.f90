module types

  type model_type
    character(len=200) :: name
    character(len=200) :: file
    integer :: index
    integer :: ncontacts
    double precision :: gscore
    double precision :: similarity
    double precision :: dgscore
    double precision :: degree
    double precision :: wdegree
  end type model_type

end module types 
