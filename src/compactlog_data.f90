module compactlog_data

  use types
  character(len=200) :: compactlog
  integer :: nmodels, score_type
  double precision, allocatable :: scores(:,:)
  type(model_type), allocatable :: model(:)

end module compactlog_data

module xcompactlog_data

  use types
  character(len=200) :: compactlog

  integer :: score_type
  double precision, allocatable :: scores(:,:)

  integer :: nmodels1
  type(model_type), allocatable :: model1(:)

  integer :: nmodels2
  type(model_type), allocatable :: model2(:)

end module xcompactlog_data

