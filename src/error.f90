MODULE error
    !
    ! ERROR
    !
    ! Contains subroutines for dealing with program errors
    !
    IMPLICIT NONE

    PUBLIC

    CONTAINS

    SUBROUTINE print_error(routine, desc)
        !
        ! This subroutine prints an error message
        !
        ! INPUT :
        !
        ! routine [char]
        !   Name of the calling routine
        ! desc [char]
        !   Description of error
        !
        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: routine
        CHARACTER(LEN=*), INTENT(IN) :: desc

        WRITE (*,"(60('*'))")
        WRITE (*,*) "    ERROR in ", trim(routine)
        WRITE (*,*) "> ", trim(desc)
        WRITE (*,"(60('*'))")

    END SUBROUTINE print_error

    SUBROUTINE print_warning(routine, desc)
        !
        ! This subroutine prints a warning message
        !
        ! INPUT :
        !
        ! routine [char]
        !   Name of the calling routine
        ! desc [char]
        !   Description of warning
        !
        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: routine
        CHARACTER(LEN=*), INTENT(IN) :: desc

        WRITE (*,"(60('*'))")
        WRITE (*,*) "    WARNING from ", trim(routine)
        WRITE (*,*) "> ", trim(desc)
        WRITE (*,"(60('*'))")

    END SUBROUTINE print_warning

END MODULE error
