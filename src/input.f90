MODULE input
    !
    ! INPUT
    ! The INPUT module contains subroutines that collect data
    ! from supplied input files and variables that store 
    ! information about the system
    !
    USE constants
    USE error

    IMPLICIT NONE
    INTEGER :: n_nuclei, n_domains, &
             & max_scf_cycles, &
             & conv_check_type, &
             & diis_length, &
             & scan_length, scan_type, &
             & matrix_vec_length
    REAL(dp) :: alpha, scf_threshold, basis_threshold, scan_step
    LOGICAL  :: input_set, default_input, print_orbitals, scan_job
    LOGICAL  :: debug = .FALSE.

    INTEGER, TARGET, ALLOCATABLE :: nuclear_charge(:), &
                                  & electrons_in_domain(:), &
                                  & functions_in_domain(:), &
                                  & orbitals_in_domain(:), &
                                  & evens_in_domain(:), odds_in_domain(:), &
                                  & evens_kept(:), odds_kept(:), &
                                  & start_of_domain(:), end_of_domain(:)
    REAL(dp), ALLOCATABLE :: nuclear_position(:)
    CHARACTER(LEN=2), ALLOCATABLE :: name_of_nucleus(:)
    LOGICAL, ALLOCATABLE :: scanning_nucleus(:)

    PUBLIC
    PRIVATE :: element_to_charge

    CONTAINS

    SUBROUTINE read_input(input_channel, exit_state)
        !
        ! This subroutine reads an input file that has
        ! already been attached to a unit
        !
        ! INPUT :
        !   [int ] record
        ! The unit to which the input file has been
        ! attached
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Fatal input file read error
        !   3 = Fatal memory allocation error
        !       (in linked lists)
        !   4 = Fatal memory allocation error
        !       (in static storage)
        !   5 = Insufficient temporary input
        !       to fill static storage
        !   6 = Fatal option set
        !
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)  :: input_channel
        INTEGER, INTENT(OUT) :: exit_state

        INTEGER  :: stat, line, i, d
        INTEGER  :: basis, exps, evens, odds
        REAL(dp) :: distance, value
        CHARACTER(LEN=2)  :: nucleus_type
        CHARACTER(LEN=10) :: option
        CHARACTER(LEN=58) :: error_message
        CHARACTER(LEN=80) :: buffer
        LOGICAL :: next_line_is_domain = .TRUE., molecule_read = .FALSE.

        TYPE domain_node
            TYPE(domain_node), POINTER :: next_node
            INTEGER :: electrons = 0
        END TYPE domain_node

        TYPE nucleus_node
            TYPE(nucleus_node), POINTER :: next_node
            INTEGER :: charge = 0
            REAL(dp) :: pos
            CHARACTER(LEN=2) :: symbol
        END TYPE nucleus_node

        TYPE(domain_node), TARGET   :: domain_list_head
        TYPE(domain_node), POINTER  :: dll_last, dll_curr
        TYPE(nucleus_node), TARGET  :: nucleus_list_head
        TYPE(nucleus_node), POINTER :: nll_last, nll_curr

        IF (input_set) THEN
            WRITE (error_message,'(a)') &
                & "Reading input a second time"
            CALL print_warning("read_input", error_message)
        ENDIF

        ! Default option values
        scf_threshold = 1.0e-6
        basis_threshold = 1.0e-7
        max_scf_cycles = 20
        print_orbitals = .FALSE.
        basis = 6
        exps = 0
        evens = 0
        odds = 0
        alpha = 0.d0
        conv_check_type = 1
        diis_length = 4
        scan_length = 0
        scan_step = 0.d0
        scan_type = 0
        scan_job = .FALSE.

        ! Setup linked lists for temporarily storing
        ! molecule specification
        NULLIFY(domain_list_head%next_node)
        NULLIFY(nucleus_list_head%next_node)

        nucleus_list_head%pos = 0.d0

        dll_curr => domain_list_head
        nll_curr => nucleus_list_head

        n_nuclei = 0

        line = 0
        ! Read the molecule specification
        DO
            READ (input_channel,'(a80)',IOSTAT=stat) buffer
            line = line + 1

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a, i4, a)') &
                    & "Could not read line ", line, " to buffer"
                CALL print_error("read_input", error_message)
                exit_state = 2
                RETURN
            ENDIF

            ! Check to see if we've reached the end of the molecule
            IF (buffer.EQ."") THEN
                molecule_read = .TRUE.
                EXIT
            ENDIF

            IF (next_line_is_domain) THEN
                READ (buffer,*,IOSTAT=stat) dll_curr%electrons

                IF (stat.NE.0) THEN
                    WRITE (error_message,'(a, i4)') &
                        & "Could not unpack the buffer on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 2
                    RETURN
                ENDIF

                ! Extend the temporary storage
                ALLOCATE(dll_curr%next_node, STAT=stat)
                dll_curr => dll_curr%next_node
                NULLIFY(dll_curr%next_node)

                IF (stat.NE.0) THEN
                    WRITE (error_message,'(a)') &
                        & "Could not allocate temporary memory"
                    CALL print_error("read_input", error_message)
                    exit_state = 3
                    RETURN
                ENDIF

                next_line_is_domain = .FALSE.
            ELSE ! the next line is a nucleus
                IF (line.EQ.2) THEN
                    READ (buffer,*,IOSTAT=stat) nucleus_type
                ELSE
                    READ (buffer,*,IOSTAT=stat) nucleus_type, distance
                ENDIF

                IF (stat.NE.0) THEN
                    WRITE (error_message,'(a, i4)') &
                        & "Could not unpack the buffer on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 2
                    RETURN
                ENDIF

                ! Pack the input data
                IF (line.GT.2) nll_curr%pos = nll_last%pos + distance
                nll_curr%charge = element_to_charge(nucleus_type, stat)
                nll_curr%symbol = nucleus_type
                n_nuclei = n_nuclei + 1

                IF (stat.NE.0) THEN
                    WRITE (error_message,'(a, i4)') &
                        & "Unrecognised element name on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 2
                    RETURN
                ENDIF

                ! Extend the temporary storage
                ALLOCATE(nll_curr%next_node, STAT=stat)
                nll_last => nll_curr
                nll_curr => nll_curr%next_node
                NULLIFY(nll_curr%next_node)

                IF (stat.NE.0) THEN
                    WRITE (error_message,'(a)') &
                        & "Could not allocate temporary memory"
                    CALL print_error("read_input", error_message)
                    exit_state = 3
                    RETURN
                ENDIF

                next_line_is_domain = .TRUE.
            ENDIF

        ENDDO

        ! Read any options that are left
        ALLOCATE(scanning_nucleus(n_nuclei))
        FORALL (i = 1:n_nuclei) scanning_nucleus(i) = .FALSE.
        DO
            READ (input_channel,'(a80)',IOSTAT=stat) buffer
            line = line + 1

            ! Check for end of file
            IF (stat.NE.0) THEN
                EXIT
            ENDIF

            READ (buffer,*,IOSTAT=stat) option, value

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a, i4, a)') &
                    & "Could not read line ", line, " to the buffer"
                CALL print_error("read_input", error_message)
                exit_state = 2
                RETURN
            ENDIF

            IF (option.EQ."debug") THEN
                IF (value.NE.0.d0) debug = .TRUE.

            ELSEIF (option.EQ."basis") THEN
                basis = int(value)
                IF (basis.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set basis to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."exps") THEN
                exps = int(value)
                IF (exps.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set exps to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."evens") THEN
                evens = int(value)
                IF (evens.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set evens to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."odds") THEN
                odds = int(value)
                IF (odds.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set odds to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."alpha") THEN
                alpha = value

            ELSEIF (option.EQ."maxscf") THEN
                max_scf_cycles = int(value)
                IF (max_scf_cycles.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set maxscf to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."p_orb") THEN
                IF (int(value).NE.0) print_orbitals = .TRUE.

            ELSEIF (option.EQ."scf_th") THEN
                scf_threshold = 10.d0**(-value)

            ELSEIF (option.EQ."bas_th") THEN
                basis_threshold = 10.d0**(-value)

            ELSEIF (option.EQ."check") THEN
                conv_check_type = int(value)
                IF ((conv_check_type.NE.0).AND.(conv_check_type.NE.1)) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "check must be set to 0 or 1", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."diis") THEN
                diis_length = int(value)
                IF (diis_length.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set diis to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."scanlength") THEN
                scan_length = int(value)
                IF (scan_length.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Attempting to set scanlength to a non-positive number", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."scanstep") THEN
                scan_step = value

            ELSEIF (option.EQ."scantype") THEN
                scan_type = int(value)
                IF (scan_type.LE.0) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "scantype must be set to 0 or 1", &
                        & " on line ", line
                    CALL print_error("read_input", error_message)
                    exit_state = 6
                    RETURN
                ENDIF

            ELSEIF (option.EQ."scannuc") THEN
                IF (int(value).GT.0.AND.int(value).LE.n_nuclei) THEN
                    scanning_nucleus(int(value)) = .TRUE.
                ELSE
                    WRITE (error_message,'(a, i4)') &
                        & "Attempting to scan non existent nucleus on line ", &
                        & line
                    CALL print_warning("read_input", error_message)
                ENDIF
                IF (exit_state.EQ.0) exit_state = 1

            ELSE
                WRITE (error_message,'(a, i4)') &
                    & "Unrecognised option on line ", line
                CALL print_warning("read_input", error_message)
                IF (exit_state.EQ.0) exit_state = 1
            ENDIF

        ENDDO

        ! WARNING : There are a number of error catches that
        !           should be added here
        scan_job = .FALSE.
        DO i = 1, n_nuclei
            scan_job = scan_job.OR.scanning_nucleus(i)
        ENDDO
        IF (scan_job) THEN
            IF (scan_length.LE.0) THEN
                WRITE (error_message,'(a)') &
                    & "scannuc has been set but scanlength is 0"
                CALL print_error("read_input", error_message)
                exit_state = 6
                RETURN
            ENDIF

            IF (scan_step.EQ.0.d0) THEN
                WRITE (error_message,'(a)') &
                    & "scannuc has been set but scanstep is 0"
                CALL print_error("read_input", error_message)
                exit_state = 6
                RETURN
            ENDIF
        ELSE
            IF (scan_length.GT.0.OR.scan_step.NE.0.d0) THEN
                WRITE (error_message,'(a)') &
                    & "scannuc must be set to perform a scan"
                CALL print_error("read_input", error_message)
                exit_state = 6
                RETURN
            ENDIF
        ENDIF

        ! Allocate static memory
        n_domains = n_nuclei + 1
        ALLOCATE( nuclear_position(n_nuclei),     &
                & nuclear_charge(n_nuclei),       &
                & name_of_nucleus(n_nuclei),      &
                & electrons_in_domain(n_domains), &
                & functions_in_domain(n_domains), &
                & evens_in_domain(n_domains),     &
                & odds_in_domain(n_domains),      &
                & orbitals_in_domain(n_domains),      &
                & evens_kept(n_domains),          &
                & odds_kept(n_domains),           &
                & start_of_domain(n_domains),     &
                & end_of_domain(n_domains),       &
                & STAT=stat                       )

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Could not allocate memory"
            CALL print_error("read_input", error_message)
            exit_state = 4
            RETURN
        ENDIF

        ! Unpack the domain information from the linked list
        dll_curr => domain_list_head
        i = 1

        DO WHILE (associated(dll_curr%next_node))

            ! Do we have too much temporary data?
            IF (i.GT.n_domains) THEN
                WRITE (error_message,'(a)') &
                    & "Temporary storage contains too much domain data"
                CALL print_warning("read_input", error_message)
                ! This probably isn't fatal, since all the static
                ! data must be full
                IF (exit_state.EQ.0) exit_state = 1
            ENDIF

            electrons_in_domain(i) = dll_curr%electrons
            IF (electrons_in_domain(i).GT.0) THEN

                IF (i.GT.1.AND.i.LT.n_domains) THEN
                    IF (evens.GT.0) THEN
                        evens_in_domain(i) = evens
                    ELSE
                        evens_in_domain(i) = basis
                    ENDIF

                    IF (odds.GT.0) THEN
                        odds_in_domain(i) = odds
                    ELSE
                        odds_in_domain(i) = basis
                    ENDIF

                    functions_in_domain(i) = evens_in_domain(i) + &
                                           & odds_in_domain(i)
                    matrix_vec_length = matrix_vec_length + &
                                      & functions_in_domain(i)**2
                ELSE
                    IF (exps.GT.0) THEN
                        functions_in_domain(i) = exps
                    ELSE
                        functions_in_domain(i) = basis
                    ENDIF

                    evens_in_domain(i) = 0
                    odds_in_domain(i) = 0
                    matrix_vec_length = matrix_vec_length + &
                                      & functions_in_domain(i)**2
                ENDIF

                ! Check that there are sufficient functions to house
                ! the electrons
                IF (functions_in_domain(i).LT.electrons_in_domain(i)) THEN
                    WRITE (error_message,'(a, a, i2)') &
                        & "Not enough basis functions to house electrons ", &
                        & "in domain ", i
                    CALL print_error("read_input", error_message)
                    exit_state = 5
                    RETURN
                ENDIF

            ELSE
                functions_in_domain(i) = 0
                evens_in_domain(i) = 0
                odds_in_domain(i) = 0
            ENDIF

            i = i + 1
            dll_curr => dll_curr%next_node
        ENDDO

        ! Do we have not enough temporary data?
        IF (i.LT.n_domains) THEN
            WRITE (error_message,'(a)') &
                & "Temporary storage contains too little domain data"
            CALL print_error("read_input", error_message)
            ! This is definitely fatal though
            exit_state = 5
            RETURN
        ENDIF

        ! Now unpack the nucleus information
        nll_curr => nucleus_list_head
        i = 1

        DO WHILE (associated(nll_curr%next_node))

            ! Do we have too much temporary data?
            IF (i.GT.n_nuclei) THEN
                WRITE (error_message,'(a)') &
                    & "Temporary storage contains too much nucleus data"
                CALL print_warning("read_input", error_message)
                ! This probably isn't fatal, since all the static
                ! data must be full
                IF (exit_state.EQ.0) exit_state = 1
            ENDIF

            nuclear_charge(i) = nll_curr%charge
            nuclear_position(i) = nll_curr%pos
            name_of_nucleus(i) = nll_curr%symbol

            i = i + 1
            nll_curr => nll_curr%next_node
        ENDDO

        ! Do we have not enough temporary data?
        IF (i.LT.n_nuclei) THEN
            WRITE (error_message,'(a)') &
                & "Temporary storage contains too little nucleus data"
            CALL print_error("read_input", error_message)
            ! This is definitely fatal though
            exit_state = 5
            RETURN
        ENDIF

        ! Compute the remaining input
        i = 0
        DO d = 1, n_domains
            orbitals_in_domain(d) = functions_in_domain(d)
            start_of_domain(d) = i + 1
            end_of_domain(d) = i + functions_in_domain(d)
            i = i + functions_in_domain(d)
        ENDDO

        IF (alpha.EQ.0.d0) THEN
            alpha = dble(sum(nuclear_charge)) / &
                & dble(max(functions_in_domain(1), &
                & functions_in_domain(n_domains))**2)
        ENDIF

        ! Deallocate the linked lists
        nll_curr => nucleus_list_head
        DO WHILE (associated(nll_curr%next_node))
            nll_last => nll_curr
            nll_curr => nll_curr%next_node
            DEALLOCATE(nll_last%next_node)
        ENDDO
        NULLIFY(nll_curr, nll_last)

        dll_curr => domain_list_head
        DO WHILE (associated(dll_curr%next_node))
            dll_last => dll_curr
            dll_curr => dll_curr%next_node
            DEALLOCATE(dll_last%next_node)
        ENDDO
        NULLIFY(dll_curr, dll_last)

        input_set = .TRUE.
        default_input = .FALSE.

    END SUBROUTINE read_input

    INTEGER FUNCTION element_to_charge(element, exit_state)
        !
        ! This function returns a nuclear charge given
        ! the two letter symbol for an element
        !
        ! INPUT :
        !
        !   element [char]
        !   Two letter symbol for an element
        !
        ! OUTPUT :
        !
        !   exit_state [int ]
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Input string unrecognised
        !
        IMPLICIT NONE

        CHARACTER(LEN=2), INTENT(IN) :: element
        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: i, Z = 0
        CHARACTER(LEN=2) :: periodic_table(18)
        DATA periodic_table / "H", "He", &
                            & "Li", "Be",  "B",  "C",  "N",  "O",  "F", "Ne", &
                            & "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar" /

        DO i = 1, 18
            IF (element.EQ.periodic_table(i)) THEN
                Z = i
                EXIT
            ENDIF
        ENDDO

        IF (Z.EQ.0) READ (element,*,IOSTAT=exit_state) Z

        element_to_charge = Z

    END FUNCTION element_to_charge

    SUBROUTINE clean_up_input(exit_state)
        !
        ! This subroutine deallocates the input storage
        ! declared by read_input/set_default_input
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Memory deallocation error
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: stat
        CHARACTER(LEN=58) :: error_message

        IF (.NOT.input_set) THEN
            WRITE (error_message,'(a)') &
                & "Input is not currently set"
            CALL print_warning("clean_up_input", error_message)
            IF (exit_state.EQ.0) exit_state = 1
            RETURN
        ENDIF

        DEALLOCATE( nuclear_position,    &
                  & nuclear_charge,      &
                  & name_of_nucleus,     &
                  & electrons_in_domain, &
                  & functions_in_domain, &
                  & evens_in_domain,     &
                  & odds_in_domain,      &
                  & start_of_domain,     &
                  & end_of_domain,       &
                  & STAT=stat            )

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Memory could not be deallocated"
            CALL print_error("clean_up_input", error_message)
            exit_state = 2
            RETURN
        ENDIF

        input_set = .FALSE.
        default_input = .FALSE.

    END SUBROUTINE clean_up_input

    SUBROUTINE print_input(file_name)
        !
        ! This subroutine writes the current input data
        ! to standard out
        !
        ! INPUT :
        !   [char] file_name
        !       Name of the input file
        !
        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: file_name

        INTEGER :: i

        WRITE (*,"(6('='), (' '), a, (' '), 47('='))") "INPUT"
        WRITE (*,*)
        WRITE (*,"(2(' '), ('Reading input from : '), a)") file_name

        WRITE (*,*)

        ! Molecule specification
        WRITE (*,"(5(' '), ('--- '), a, (' ---'))") "MOLECULE"
        WRITE (*,"(2(' '), a)") "Electrons"
        WRITE (*,"(3(' '), a, 6(' '), a)") "Nucleus", "Position"

        DO i = 1, n_nuclei
            WRITE (*,"(6(' '), i2)") electrons_in_domain(i)
            WRITE (*,"(6(' '), a, 8(' '), F8.3)") name_of_nucleus(i), &
                                                & nuclear_position(i)
        ENDDO
        WRITE (*,"(6(' '), i2)") electrons_in_domain(n_domains)

        WRITE (*,*)
        WRITE (*,"(5(' '), ('--- '), a, (' ---'))") "OPTIONS"

        ! Basis set options
        WRITE (*,"(2(' '), a, (' : '), i4)") &
            & "Exponential functions", functions_in_domain(1)
        IF (n_domains.GT.2) THEN
            WRITE (*,"(2(' '), a, ('       : '), i4)") &
                & "Evens functions", evens_in_domain(2)
            WRITE (*,"(2(' '), a, ('        : '), i4)") &
                & "Odds functions", odds_in_domain(2)
        ENDIF
        WRITE (*,"(2(' '), a, ('                 : '), F10.6)") &
            & "Alpha", alpha

        WRITE (*,*)

        ! Geometry scan options
        IF (scan_job) THEN
            WRITE (*,"(2(' '), a)") &
                & "Running a geometry scan"
            WRITE (*,"(2(' '), a)") &
                & "Moving the following nuclei"
            DO i = 1, n_nuclei
                IF (scanning_nucleus(i)) WRITE (*,"(6(' '), i2)") &
                    & i
            ENDDO
            WRITE (*,"(2(' '), a, ' : ', i6)") &
                & "Number of steps", scan_length
            WRITE (*,"(2(' '), a, ' : ', F7.4)") &
                & "Step length", scan_step
            WRITE (*,*)
            IF (scan_type.EQ.1) THEN
                WRITE (*,"(2(' '), a)") &
                    & "Resetting orbitals for each step"
            ENDIF
        ENDIF

        ! Other variables
        WRITE (*,"(2(' '), a, ('          : '), i4)") &
            & "Maximum SCF cycles", max_scf_cycles
        WRITE (*,"(2(' '), a, ('          : '), i4)") &
            & "DIIS memory length", diis_length
        WRITE (*,"(2(' '), a, ('   : '), es7.1)") &
            & "SCF convergence threshold", scf_threshold
        WRITE (*,"(2(' '), a, (' : '), es7.1)") &
            & "Linear dependence threshold", basis_threshold

        WRITE (*,*)

        ! Discrete options
        IF (conv_check_type.EQ.0) THEN
            WRITE (*,"(2(' '), a, a, a)") &
                & "Checking ", &
                & "max element of error matrix", &
                & " for SCF convergence"
        ELSEIF (conv_check_type.EQ.1) THEN
            WRITE (*,"(2(' '), a, a, a)") &
                & "Checking ", &
                & "RMS of error matrix", &
                & " for SCF convergence"
        ENDIF

        IF (print_orbitals) THEN
            WRITE (*,"(2(' '), a)") &
                & "Printing occupied HF orbitals"
        ENDIF

        WRITE (*,*)
    END SUBROUTINE print_input

END MODULE input
