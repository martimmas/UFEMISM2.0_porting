module basic_program_info

  private

  public :: program_name, routine_path, n_MPI_windows_used

  character(len=1024) :: program_name = 'GIVE ME A NAME'
  character(len=1024) :: routine_path
  integer             :: n_MPI_windows_used

end module basic_program_info