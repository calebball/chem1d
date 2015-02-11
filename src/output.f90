MODULE output
    !
    ! OUTPUT
    !
    ! Contains subroutines for printing output
    !
    USE input,          ONLY : nuclear_position, &
                             & scan_step
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

        WRITE (*,*)

    END SUBROUTINE print_header

    SUBROUTINE print_footer
        !
        ! This subroutine prints the end of the output
        !
        IMPLICIT NONE

        WRITE (*,"(60('='))")

    END SUBROUTINE print_footer

    SUBROUTINE print_scan_start
        !
        ! This subroutine prints the start of the scan
        ! section in the output
        !
        IMPLICIT NONE

        WRITE (*,"(60('='))")
        WRITE (*,*)
        WRITE (*,"(6(' '), '--- ', a, ' ---')") "START OF SCAN"
        WRITE (*,"(2(' '), a, 2(' '), a, 10(' '), a, 7(' '), a)") &
            & "Step", "Displacement", "HF Energy", "SCF Cycles"

    END SUBROUTINE print_scan_start

    SUBROUTINE print_scan_step(step)
        !
        ! This subroutine prints the output for an 
        ! individual step during the scan process
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: step

        WRITE (*,"(2(' '), i4, 5(' '), f6.3, 5(' '), f22.15, 5(' '), i3)") &
            & step, step * scan_step, hf_energy(0), curr_hf_cycle

    END SUBROUTINE print_scan_step

    SUBROUTINE print_scan_end
        !
        ! This subroutine prints the end of the scan
        ! section in the program output
        !
        IMPLICIT NONE

        WRITE (*,*)
        WRITE (*,"(60('='))")
    END SUBROUTINE print_scan_end

END MODULE output
