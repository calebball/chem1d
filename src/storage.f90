MODULE storage
    !
    ! STORAGE
    !
    ! This module contains type definitions, declarations and
    ! initialisation routines of various data structures used
    ! throughout the program
    !
    USE constants
    USE error
    USE input
    USE one_e_integrals, ONLY : overlap_int_in_domain, &
                              & kinetic_int_in_domain, &
                              & potential_int_in_domain

    IMPLICIT NONE
    PUBLIC

    TYPE one_e_matrices
        REAL(dp), ALLOCATABLE :: overlap(:,:)
        REAL(dp), ALLOCATABLE :: kinetic(:,:)
        REAL(dp), ALLOCATABLE :: potential(:,:)
    END TYPE one_e_matrices

    TYPE hf_properties
        REAL(dp), ALLOCATABLE :: orbitals(:,:)
        REAL(dp), ALLOCATABLE :: density(:,:)
        REAL(dp), ALLOCATABLE :: orb_energies(:)
    END TYPE hf_properties

    TYPE dense_eri_storage
        REAL(dp), ALLOCATABLE :: integral(:,:,:,:)
    END TYPE dense_eri_storage
    TYPE eri_storage
        TYPE(dense_eri_storage), ALLOCATABLE :: with(:)
    END TYPE eri_storage


    TYPE(one_e_matrices), TARGET, ALLOCATABLE :: one_e_in_domain(:)
    TYPE(hf_properties),  TARGET, ALLOCATABLE :: hf_in_domain(:)
    TYPE(eri_storage), ALLOCATABLE :: eri_of_domain(:)
    TYPE(eri_storage), ALLOCATABLE :: mo_eri_domain(:)
    REAL(dp), TARGET, ALLOCATABLE :: double_bar(:,:,:,:)
    REAL(dp), ALLOCATABLE :: mo_double_bar(:,:,:,:)


    CONTAINS

    SUBROUTINE init_storage(exit_state)
        !
        ! This subroutine initialises storage that is to be
        ! used throughout the runtime of the program
        !
        ! OUTPUT :
        !
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Memory not allocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, d_two, functions, stat
        CHARACTER(LEN=58) :: error_message
        INTEGER, POINTER :: m, m_two

        ALLOCATE(one_e_in_domain(n_domains), &
                &hf_in_domain(n_domains),    &
                &eri_of_domain(n_domains),   &
                &mo_eri_domain(n_domains),   &
                &STAT=stat                   )

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Derived type variables not allocated"
            CALL print_error("init_storage", error_message)
            exit_state = 2
            RETURN
        ENDIF

        functions = 0
        DO d = 1, n_domains
            m => functions_in_domain(d)
            functions = functions + m

            ! Allocate memory for one electron matrices
            ALLOCATE(one_e_in_domain(d)%overlap(m,m),   &
                    &one_e_in_domain(d)%kinetic(m,m),   &
                    &one_e_in_domain(d)%potential(m,m), &
                    &STAT=stat)

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a)') &
                    & "One electrons matrices not allocated"
                CALL print_error("init_storage", error_message)
                exit_state = 2
                RETURN
            ENDIF

            ! Allocate memory for storing the results of
            ! Hartree-Fock calculations
            ALLOCATE(hf_in_domain(d)%orbitals(m,m),   &
                    &hf_in_domain(d)%density(m,m),    &
                    &hf_in_domain(d)%orb_energies(m), &
                    &STAT=stat)

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a)') &
                    & "Hartree-Fock storage not allocated"
                CALL print_error("init_storage", error_message)
                exit_state = 2
                RETURN
            ENDIF

            ! Allocate memory for storing electron
            ! repulsion integrals
            ALLOCATE(eri_of_domain(d)%with(n_domains), &
                    &mo_eri_domain(d)%with(n_domains), &
                    &STAT=stat)
            DO d_two = 1, n_domains
                m_two => functions_in_domain(d_two)
                ALLOCATE(eri_of_domain(d)%with(d_two)%integral(m,m,m_two,m_two), &
                        &mo_eri_domain(d)%with(d_two)%integral(m,m,m_two,m_two), &
                        &STAT=stat)
            ENDDO

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a)') &
                    & "Electron repulsion integral storage not allocated"
                CALL print_error("init_storage", error_message)
                exit_state = 2
                RETURN
            ENDIF

        ENDDO

        ALLOCATE(double_bar(functions,functions,functions,functions),    &
                &mo_double_bar(functions,functions,functions,functions), &
                &STAT=stat)

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Double bar integral storage not allocated"
            CALL print_error("init_storage", error_message)
            exit_state = 2
            RETURN
        ENDIF

    END SUBROUTINE

    SUBROUTINE populate_one_e_matrices
        !
        ! This subroutine fills the overlap, kinetic energy
        ! and potential energy matrices
        !
        IMPLICIT NONE

        INTEGER :: d

        ! WARNING : Does this procedure need error checks?
        DO d = 1, n_domains
            CALL overlap_int_in_domain(d, one_e_in_domain(d)%overlap)
            CALL kinetic_int_in_domain(d, one_e_in_domain(d)%kinetic)
            CALL potential_int_in_domain(d, one_e_in_domain(d)%overlap, &
                                        &one_e_in_domain(d)%potential)
        ENDDO
    END SUBROUTINE populate_one_e_matrices

    SUBROUTINE clean_up_storage(exit_state)
        !
        ! This subroutine deallocates all the storage 
        ! that was set by init_storage
        !
        ! OUTPUT :
        !
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Memory not deallocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, d_two, stat
        CHARACTER(LEN=58) :: error_message

        DO d = 1, n_domains

            DO d_two = 1, n_domains
                DEALLOCATE(eri_of_domain(d)%with(d_two)%integral, &
                          &mo_eri_domain(d)%with(d_two)%integral, &
                          &STAT=stat)
            ENDDO

            DEALLOCATE(one_e_in_domain(d)%overlap,   &
                      &one_e_in_domain(d)%kinetic,   &
                      &one_e_in_domain(d)%potential, &
                      &hf_in_domain(d)%orbitals,     &
                      &hf_in_domain(d)%density,      &
                      &hf_in_domain(d)%orb_energies, &
                      &eri_of_domain(d)%with,        &
                      &mo_eri_domain(d)%with,        &
                      &STAT=stat)

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a, i4, a)') &
                    & "Storage for domain ", d, " not deallocated"
                CALL print_error("clean_up_storage", error_message)
                exit_state = 2
            ENDIF

        ENDDO

        DEALLOCATE(one_e_in_domain, &
                  &hf_in_domain,    &
                  &eri_of_domain,   &
                  &mo_eri_domain,   &
                  &double_bar,      &
                  &mo_double_bar,   &
                  &STAT=stat)

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a, i4, a)') &
                & "Top level storage not deallocated"
            CALL print_error("clean_up_storage", error_message)
            exit_state = 2
        ENDIF

    END SUBROUTINE clean_up_storage

END MODULE storage
