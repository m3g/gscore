module types

  type model_type
    character(len=200) :: name
    integer :: index
    double precision :: gscore
    double precision :: similarity
    double precision :: dgscore
  end type model_type

end module types 
