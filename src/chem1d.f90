PROGRAM chem1d
    USE constants
    USE input
    USE output
    USE storage,        ONLY : init_storage, &
                             & populate_one_e_matrices, &
                             & clean_up_storage, &
                             & hf_in_domain
    USE eri,            ONLY : build_eri,            &
                             & mo_eri_transform
    USE hartree_fock,   ONLY : init_hartree_fock,    &
                             & clean_up_hf,          &
                             & hf_one_electron,      &
                             & hf_scf_cycle,         &
                             & print_hf_start,       &
                             & print_hf_end,        &
                             & curr_hf_cycle
    USE moller_plesset, ONLY : compute_mp2, &
                             & compute_mp3, &
                             & print_mp
    IMPLICIT NONE

    CHARACTER(LEN=256) :: in_file, out_file
    INTEGER :: in_channel, out_channel, stat
    REAL(sp) :: cpu_time_start, cpu_time_finish
    INTEGER :: sys_clock_start, sys_clock_finish, sys_clock_rate
    INTEGER :: step, nucleus
    INTEGER :: d, i, j

    CALL cpu_time(cpu_time_start)
    CALL system_clock(sys_clock_start, sys_clock_rate)

    ! Print a nice little header
    CALL print_header

    ! Start by opening and reading the input and output files
    in_channel  = 10
    IF (command_argument_count().GT.0) THEN
        CALL get_command_argument(1, in_file)
        OPEN (in_channel, FILE=trim(in_file), ACTION="read", STATUS="old", IOSTAT=stat)

        IF (stat.NE.0) THEN
            CALL print_error("chem1d", "Could not open supplied input file")
            STOP
        ENDIF

        CALL read_input(in_channel, stat)
        IF (stat.GT.1) STOP

        CLOSE (in_channel, IOSTAT=stat)
        IF (stat.NE.0) THEN
            CALL print_warning("chem1d", "Could not close input channel")
        ENDIF

        CALL print_input(trim(in_file))

    ELSE

        CALL print_error("chem1d", "No input file supplied")
        STOP

    ENDIF

    ! Initialise storage for globally held data
    CALL init_storage(stat)
    IF (stat.GT.1) STOP

    ! Are we scanning some geometry parameter?
    IF (scan_job) THEN
        CALL print_scan_start

        CALL init_hartree_fock(stat)
        IF (stat.GT.1) STOP

        DO step = 0, scan_length

            ! Move geometry
            IF (step.NE.0) THEN
                DO nucleus = 2, n_nuclei
                    IF (scanning_nucleus(nucleus)) THEN
                        DO i = nucleus, n_nuclei
                            nuclear_position(i) = nuclear_position(i) + &
                                & scan_step
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
            curr_hf_cycle = 0

            IF (scan_type.eq.1) THEN ! Start from initial orbital guess
                DO d = 1, n_domains
                    DO i = 1, functions_in_domain(d)
                    DO j = 1, functions_in_domain(d)
                        hf_in_domain(d)%orbitals(i,j) = 0.d0
                        hf_in_domain(d)%density(i,j) = 0.d0
                    ENDDO
                    ENDDO
                ENDDO
            ENDIF

            CALL populate_one_e_matrices

            ! Build ERI arrays
            CALL build_eri(stat)
            IF (stat.GT.1) STOP

            ! Setup and run a HF calculation
            CALL hf_one_electron(stat)
            IF (stat.GT.1) STOP
            CALL hf_scf_cycle(stat)
            IF (stat.GT.1) STOP

            CALL print_scan_step(step)

        ENDDO

    ELSE ! We're just doing a normal calculation

        CALL populate_one_e_matrices

        ! Build ERI arrays
        CALL build_eri(stat)
        IF (stat.GT.1) STOP

        ! Setup and run a HF calculation
        CALL print_hf_start
        CALL init_hartree_fock(stat)
        IF (stat.GT.1) STOP
        CALL hf_one_electron(stat)
        IF (stat.GT.1) STOP
        CALL hf_scf_cycle(stat)
        IF (stat.GT.1) STOP

        CALL print_hf_end

!DO d = 1, n_domains
!DO i = 1, functions_in_domain(d)
!DO j = functions_kept(d) + 1, functions_in_domain(d)
!    hf_in_domain(d)%orbitals(i,j) = 0.d0
!ENDDO
!ENDDO
!ENDDO

        ! Transform ERIs to MO basis
        CALL mo_eri_transform(stat)
        IF (stat.GT.1) STOP

        ! Run MP2 and MP3 correlation calculations
        CALL compute_mp2
        CALL compute_mp3
        CALL print_mp

    ENDIF

    ! Clean up memory
    CALL clean_up_hf(stat)
    CALL clean_up_storage(stat)
    CALL clean_up_input(stat)

    CALL cpu_time(cpu_time_finish)
    CALL system_clock(sys_clock_finish)
    WRITE (*,"(6('='), (' '), a, (' '), 46('='))") "TIMING"
    WRITE (*,*)
    WRITE (*,"(2(' '), a, F10.3)") "cpu_time     : ", &
        cpu_time_finish - cpu_time_start
    WRITE (*,"(2(' '), a, F10.3)") "system_clock : ", &
        dble(sys_clock_finish - sys_clock_start) / dble(sys_clock_rate)
    write (*,*)

    CALL print_footer

END PROGRAM chem1d
