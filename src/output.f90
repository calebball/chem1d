MODULE output
    !
    ! OUTPUT
    !
    ! Contains subroutines for printing output
    !
    USE input,          ONLY : nuclear_position
    USE hartree_fock,   ONLY : hf_energy, &
                             & curr_hf_cycle
    USE moller_plesset, ONLY : mp2_correction, &
                             & mp3_correction
    IMPLICIT NONE

    PUBLIC

    CONTAINS

    SUBROUTINE print_header
        !
        ! This subroutine prints the program header for the
        ! beginning of the file
        !
        IMPLICIT NONE

        WRITE (*,"(60('='))")
        WRITE (*,*)
        WRITE (*,*) "  ####                      #  ######  "
        WRITE (*,*) " ##   # #   # ##### #   #  ##  ##   ## "
        WRITE (*,*) " ##     #   # #     ## ##   #  ##   ## "
        WRITE (*,*) " ##     ##### ###   # # #   #  ##   ## "
        WRITE (*,*) " ##   # #   # #     #   #   #  ##   ## "
        WRITE (*,*) "  ####  #   # ##### #   #  ### ######  "

!       WRITE (*,"(60('-'))")
        WRITE (*,*)

    END SUBROUTINE print_header

    SUBROUTINE print_footer
        !
        ! This subroutine prints the program header for the
        ! beginning of the file
        !
        IMPLICIT NONE

        WRITE (*,"(60('='))")

    END SUBROUTINE print_footer

    SUBROUTINE print_scan_start
        IMPLICIT NONE

        WRITE (*,"(60('='))")
        WRITE (*,*)
    END SUBROUTINE print_scan_start

    SUBROUTINE print_scan_step(step)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: step

!        WRITE (*,"(2(' '), F6.3, 4(' '), F22.15, 4(' '), F22.15, 4(' '), F22.15, 4(' '), i3)") &
!            &nuclear_position(2), hf_energy(0), mp2_correction, mp3_correction, curr_hf_cycle

        WRITE (*,"(2(' '), F6.3, 4(' '), F6.3, 4(' '), F6.3, 4(' '), F22.15, 4(' '), i3)") &
            &nuclear_position(2), nuclear_position(3), &
            &nuclear_position(3) - nuclear_position(2), hf_energy(0), curr_hf_cycle

    END SUBROUTINE print_scan_step

    SUBROUTINE print_scan_end
        IMPLICIT NONE

        WRITE (*,*)
        WRITE (*,"(60('='))")
    END SUBROUTINE print_scan_end

END MODULE output
