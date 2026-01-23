module SMB_model

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use SMB_model_basic, only: atype_SMB_model

  implicit none

  private

  public :: atype_SMB_model, create_SMB_model

contains

  subroutine create_SMB_model( SMB, choice_SMB_model)
    !< Allocate a concrete implementation of a SMB_model

    ! In/output variables:
    class(atype_SMB_model), allocatable, intent(inout) :: SMB
    character(len=*),                    intent(in   ) :: choice_SMB_model

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (choice_SMB_model)
    case default
      call crash('invalid choice_SMB_model "' // trim( choice_SMB_model) // '"')
    ! case ('idealised')
    !   allocate( type_SMB_model_idealised :: SMB)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_SMB_model

end module SMB_model