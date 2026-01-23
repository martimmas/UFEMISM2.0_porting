module SMB_model_state

  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_SMB_model_state

  type :: type_SMB_model_state

    real(dp), dimension(:), contiguous, pointer :: SMB => null()
    type(MPI_WIN) :: wSMB

  end type type_SMB_model_state

end module SMB_model_state