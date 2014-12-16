MODULE hartree_fock
    !
    ! HARTREE_FOCK
    !
    ! Contains the subroutines and data structures
    ! for running a Hartree-Fock self-consistent
    ! field calculation
    !
    USE constants
    USE input
    USE error
    USE storage,         ONLY : one_e_in_domain,      &
                              & hf_in_domain,         &
                              & double_bar
    USE one_e_integrals, ONLY : overlap_int_in_domain

    IMPLICIT NONE

    TYPE hf_one_electron_block
        REAL(dp), ALLOCATABLE :: s_evals(:)
        REAL(dp), ALLOCATABLE :: u_matrix(:,:)
        REAL(dp), ALLOCATABLE :: x_matrix(:,:)
    END TYPE hf_one_electron_block

    TYPE hf_two_electron_block
        REAL(dp), ALLOCATABLE :: g(:,:)
        REAL(dp), ALLOCATABLE :: fock(:,:)
        REAL(dp), ALLOCATABLE :: fock_p(:,:)
    END TYPE hf_two_electron_block

    TYPE hf_diis_block
        REAL(dp), ALLOCATABLE :: stored_fock(:,:,:)
        REAL(dp), ALLOCATABLE :: stored_error(:,:,:)
        REAL(dp), ALLOCATABLE :: diis_matrix(:,:)
        REAL(dp), ALLOCATABLE :: diis_vector(:)
        REAL(dp), ALLOCATABLE :: spf(:,:)
        REAL(dp), ALLOCATABLE :: fps(:,:)
    END TYPE hf_diis_block

    TYPE hf_work_block
        REAL(dp), ALLOCATABLE :: matrix(:,:,:)
        INTEGER, ALLOCATABLE  :: ISUPPZ(:)
        REAL(dp), ALLOCATABLE :: WORK(:)
        INTEGER, ALLOCATABLE  :: IWORK(:)
        INTEGER, ALLOCATABLE  :: IPIV(:)
    END TYPE hf_work_block

    TYPE(hf_one_electron_block), ALLOCATABLE, TARGET :: hf_one_e(:)
    TYPE(hf_two_electron_block), ALLOCATABLE :: hf_two_e(:)
    TYPE(hf_diis_block), ALLOCATABLE :: hf_diis(:)
    TYPE(hf_work_block), ALLOCATABLE :: hf_work(:)

    PRIVATE

    PUBLIC :: init_hartree_fock
    PUBLIC :: clean_up_hf
    PUBLIC :: hf_one_electron
    PUBLIC :: hf_scf_cycle
    PUBLIC :: hf_energy

    PUBLIC :: print_hf_start
    PUBLIC :: print_hf_end

    INTEGER, PUBLIC :: curr_hf_cycle
    INTEGER :: curr_diis, curr_basis_th
    LOGICAL, PUBLIC :: hf_converged

    CONTAINS

!====== COMPUTATION HANDLERS =================================================!

    SUBROUTINE hf_one_electron(exit_state)
        !
        ! This subroutine computes the one electron matrices
        ! specific to a Hartree-Fock calculation
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = DSYEVR failed
        !   3 = DSYEVR called incorrectly
        !   4 = Negative eigenvalues found
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, i, j, M, INFO, throw_count
        REAL(dp) :: condition, temp_thresh
        REAL(dp), POINTER :: overlap_e(:,:), u_matrix_e(:,:), evals_e(:)
        REAL(dp), POINTER :: overlap_o(:,:), u_matrix_o(:,:), evals_o(:)
        INTEGER, POINTER :: evens, odds
        CHARACTER(LEN=58) :: error_message

        DO d = 1, n_domains
            IF (functions_in_domain(d).EQ.0) CYCLE

            ! Start by finding the eigenvalues
            CALL DSYEV('N',                   & ! JOBZ
                &'U',                         & ! UPLO
                &functions_in_domain(d),      & ! N
                &one_e_in_domain(d)%overlap,  & ! A
                &functions_in_domain(d),      & ! LDA
                &hf_one_e(d)%s_evals,         & ! W
                &hf_work(d)%WORK,             & ! WORK
                &26 * functions_in_domain(d), & ! LWORK
                &INFO                         ) ! INFO

!           CALL DSYEVR('N',                  & ! JOBZ
!               &'A',                         & ! RANGE
!               &'U',                         & ! UPLO
!               &functions_in_domain(d),      & ! N
!               &one_e_in_domain(d)%overlap,  & ! A
!               &functions_in_domain(d),      & ! LDA
!               &0.d0, 0.d0,                  & ! VL, VU
!               &0, 0,                        & ! IL, IU
!               &0,                           & ! ABSTOL
!               &M,                           & ! M
!               &hf_one_e(d)%s_evals,         & ! W
!               &hf_one_e(d)%u_matrix,        & ! Z
!               &functions_in_domain(d),      & ! LDZ
!               &hf_work(d)%ISUPPZ,           & ! ISUPPZ
!               &hf_work(d)%WORK,             & ! WORK
!               &26 * functions_in_domain(d), & ! LWORK
!               &hf_work(d)%IWORK,            & ! IWORK
!               &10 * functions_in_domain(d), & ! LIWORK
!               &INFO                         ) ! INFO

            ! Choose how many functions to keep
            IF (primary_basis_th.GT.0.AND.&
               &secondary_basis_th.GT.primary_basis_th) THEN

                functions_kept(d) = functions_in_domain(d)
                throw_count = 1
                temp_thresh = 10.d0**(-primary_basis_th)
                DO WHILE (hf_one_e(d)%s_evals(throw_count).LT.basis_threshold)
                    functions_kept(d) = functions_kept(d) - 1
                    throw_count = throw_count + 1
                ENDDO
                curr_basis_th = primary_basis_th

            ELSE

                functions_kept(d) = functions_in_domain(d)
                throw_count = 1
                DO WHILE (hf_one_e(d)%s_evals(throw_count).LT.basis_threshold)
                    functions_kept(d) = functions_kept(d) - 1
                    throw_count = throw_count + 1
                ENDDO

            ENDIF

!do i = 1, functions_in_domain(d)
!if (i.EQ.throw_count) write (*,"(20('-'), a, 20('-'))") " keep the following "
!write (*,"(F30.20, E30.20)") hf_one_e(d)%s_evals(i), hf_one_e(d)%s_evals(i)
!enddo
!write (*,*)

            ! DSYEVR destroys the overlap matrix, so we have to repair it
            FORALL (i=1:functions_in_domain(d), j=1:functions_in_domain(d))
                one_e_in_domain(d)%overlap(i,j) = 0.d0
            ENDFORALL
            CALL overlap_int_in_domain(d, one_e_in_domain(d)%overlap)

            ! Now find the SVD
            CALL DGESVD('S',                    & ! JOBU
                & 'O',                          & ! JOBVT
                & functions_in_domain(d),       & ! M
                & functions_in_domain(d),       & ! N
                & one_e_in_domain(d)%overlap,   & ! A
                & functions_in_domain(d),       & ! LDA
                & hf_one_e(d)%s_evals,          & ! S
                & hf_one_e(d)%u_matrix,         & ! U
                & functions_in_domain(d),       & ! LDU
                & hf_work(d)%matrix(:,:,1),     & ! VT
                & functions_in_domain(d),       & ! LDVT
                & hf_work(d)%WORK,              & ! WORK
                & 26 * functions_in_domain(d),  & ! LWORK
                & INFO                          ) ! INFO

            IF (INFO.GT.0) THEN
                WRITE (error_message,'(a, i4)') &
                    & "DSYEVR failed in domain ", 1
                CALL print_error("hf_one_electron", error_message)
                exit_state = 2
                RETURN
            ELSEIF (INFO.LT.0) THEN
                WRITE (error_message,'(a, i4, a, i2, a)') &
                    & "DSYEVR in domain ", 1, "received incorrect ", &
                    & abs(INFO), "th argument"
                CALL print_error("hf_one_electron", error_message)
                exit_state = 3
                RETURN
            ENDIF

            ! Compute the transformation matrix
            FORALL (i = 1:functions_in_domain(d), &
                   &j = 1:functions_in_domain(d)  )
                hf_work(d)%matrix(i,j,1) = 0.d0
            END FORALL

            FORALL (i = 1:functions_in_domain(d))
                hf_work(d)%matrix(i,i,1) = &
                    1.d0 / sqrt(hf_one_e(d)%s_evals(i))
            END FORALL

            CALL DGEMM('N',                   & ! TRANSA
                &'N',                         & ! TRANSB
                &functions_in_domain(d),      & ! M
                &functions_in_domain(d),      & ! N
                &functions_in_domain(d),      & ! K
                &1.d0,                        & ! ALPHA
                &hf_one_e(d)%u_matrix,        & ! A
                &functions_in_domain(d),      & ! LDA
                &hf_work(d)%matrix(:,:,1),    & ! B
                &functions_in_domain(d),      & ! LDB
                &0.d0,                        & ! BETA
                &hf_one_e(d)%x_matrix,        & ! C
                &functions_in_domain(d)       ) ! LDC

            ! DGESVD destroys the overlap matrix, 
            ! so we have to repair it again
            FORALL (i=1:functions_in_domain(d), j=1:functions_in_domain(d))
                one_e_in_domain(d)%overlap(i,j) = 0.d0
            ENDFORALL
            CALL overlap_int_in_domain(d, one_e_in_domain(d)%overlap)

            ! X * Xt orbital guess
!           hf_work(d)%matrix(:,:,1) = hf_one_e(d)%x_matrix
!           CALL DGEMM('N', 'T',&
!               &functions_in_domain(d), &
!               &functions_in_domain(d), &
!               &functions_in_domain(d), &
!               &1.d0, hf_one_e(d)%x_matrix, functions_in_domain(d), &
!               &hf_work(d)%matrix(:,:,1), functions_in_domain(d), 0.d0, &
!               &hf_in_domain(d)%orbitals, functions_in_domain(d))

!write (*,*) "KINETIC", d
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') one_e_in_domain(d)%kinetic(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)

!write (*,*) "POTENTIAL", d
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') one_e_in_domain(d)%potential(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)

        ENDDO

if (.NOT.scan_job) then
write (*,*) "functions kept"
do d = 1, n_domains
write (*,*) functions_kept(d), "of", functions_in_domain(d)
enddo
write (*,*)
endif

    END SUBROUTINE hf_one_electron

    RECURSIVE SUBROUTINE hf_scf_cycle(exit_state)
        !
        ! This subroutine performs one Hartree_Fock 
        ! self-consistent field cycle and checks for
        ! convergence
        ! Calls itself recursively until either
        ! convergence is achieved or the number of
        ! cycles exceed max_scf_cycles
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Fock transform failed
        !   3 = Fock transform called incorrectly
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, i, j, M, INFO, throw_count
        CHARACTER(LEN=80) :: error_message
        REAL(dp) :: error, last_energy, temp_thresh

        ! Assume we converge, check later if a domain hasn't
        hf_converged = .TRUE.

        ! Update the density matrix for all the domains
        DO d = 1, n_domains
            if (curr_hf_cycle.GT.0) CALL build_density_matrix(d)

if (debug) then
    ploop: do i = 1, functions_in_domain(d)
    do j = i, functions_in_domain(d)
        if (hf_in_domain(d)%density(i,j).ne.hf_in_domain(d)%density(i,j)) then
            write (*,*) "density in domain", d, " not symmetric"
            exit ploop
        endif
    enddo
    enddo ploop
endif

        ENDDO

        ! Build two-electron matrix
        CALL build_g_matrix

!d = 3
!write (*,*)
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') hf_in_domain(d)%density(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)

!do d = 1, n_domains
!write (*,*) "G", d
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') hf_two_e(d)%g(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)
!enddo

        ! Solve each domain individually
        DO d = 1, n_domains
            IF (functions_in_domain(d).EQ.0) CYCLE

            ! Build fock matrix
            hf_two_e(d)%fock = one_e_in_domain(d)%kinetic - &
                             & one_e_in_domain(d)%potential + &
                             & hf_two_e(d)%g

if (debug) then
    floop: do i = 1, functions_in_domain(d)
    do j = i, functions_in_domain(d)
        if (hf_two_e(d)%fock(i,j).ne.hf_two_e(d)%fock(i,j)) then
            write (*,*) "fock in domain", d, " not symmetric"
            exit floop
        endif
    enddo
    enddo floop
endif

!write (*,*) "FOCK", d
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') hf_two_e(d)%fock(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)

            ! Shuffle the stored DIIS matrices
            curr_diis = min(curr_hf_cycle + 1, diis_length)
            DO i = curr_diis - 1, 1, -1
                hf_diis(d)%stored_fock(:,:,i+1) = &
                    & hf_diis(d)%stored_fock(:,:,i)
                hf_diis(d)%stored_error(:,:,i+1) = &
                    & hf_diis(d)%stored_error(:,:,i)
            ENDDO
            hf_diis(d)%stored_fock(:,:,1) = hf_two_e(d)%fock

            ! Build DIIS error matrix
            CALL build_diis_error(d, hf_diis(d)%stored_error(:,:,1))
            CALL orthogonalize(d, hf_diis(d)%stored_error(:,:,1))

            ! Build the new Fock prime matrix
            IF (curr_hf_cycle.LT.diis_length.OR.diis_length.EQ.0&
                &.OR.check_convergence(d, .TRUE.)) THEN
!           IF (curr_hf_cycle.LT.diis_length) THEN
                hf_two_e(d)%fock_p = hf_two_e(d)%fock
            ELSE

!write (*,'(i1)',advance='no') d
                FORALL (i=1:functions_in_domain(d), j=1:functions_in_domain(d))
                    hf_two_e(d)%fock_p(i,j) = 0.d0
                END FORALL
                ! Construct and solve the DIIS matrix equation
                CALL build_diis_matrix(d, hf_diis(d)%diis_matrix)

                DO i = 1, curr_diis
                    hf_diis(d)%diis_vector(i) = 0.d0
                ENDDO
                hf_diis(d)%diis_vector(diis_length + 1) = -1.d0

                CALL DGESV(curr_diis + 1,           & ! N
                          &1,                       & ! NRHS
                          &hf_diis(d)%diis_matrix,  & ! A
                          &diis_length + 1,         & ! LDA
                          &hf_work(d)%IPIV,         & ! IPIV
                          &hf_diis(d)%diis_vector,  & ! B
                          &diis_length + 1,         & ! LDB
                          &INFO                     ) ! INFO

                DO i = 1, curr_diis
                    hf_two_e(d)%fock_p = hf_two_e(d)%fock_p + &
                        & hf_diis(d)%diis_vector(i) * hf_diis(d)%stored_fock(:,:,i)
                ENDDO

            ENDIF
            CALL orthogonalize(d, hf_two_e(d)%fock_p)

if (debug) then
    fploop: do i = 1, functions_in_domain(d)
    do j = i, functions_in_domain(d)
        if (hf_two_e(d)%fock_p(i,j).ne.hf_two_e(d)%fock_p(i,j)) then
            write (*,*) "fock in domain", d, " not symmetric"
            exit fploop
        endif
    enddo
    enddo fploop
endif

!write (*,*) "FOCK_P", d
!do i = 1, functions_in_domain(d)
!do j = 1, functions_in_domain(d)
!    write (*,"(F10.5)",advance='no') hf_two_e(d)%fock_p(i,j)
!enddo
!write (*,*)
!enddo
!write (*,*)

            ! Diagonalise transformed fock matrix
            hf_work(d)%matrix(:,:,1) = hf_two_e(d)%fock_p
            CALL DSYEV('V', 'U',         & ! JOBZ, UPLO
                &functions_kept(d),            & ! N
                &hf_work(d)%matrix(:,:,1),     & ! A
                &functions_in_domain(d),       & ! LDA
                &hf_in_domain(d)%orb_energies, & ! W
                &hf_work(d)%WORK,              & ! WORK
                &26 * functions_in_domain(d),  & ! LWORK
                &INFO                          ) ! INFO

!           CALL DSYEVR('V', 'A', 'U',         & ! JOBZ, RANGE, UPLO
!               &functions_kept(d),            & ! N
!               &hf_two_e(d)%fock_p,           & ! A
!               &functions_in_domain(d),       & ! LDA
!               &0.d0, 0.d0, 0, 0,             & ! VL, VU, IL, IU
!               &0,                            & ! ABSTOL
!               &M,                            & ! M
!               &hf_in_domain(d)%orb_energies, & ! W
!               &hf_work(d)%matrix(:,:,1),     & ! Z
!               &functions_in_domain(d),       & ! LDZ
!               &hf_work(d)%ISUPPZ,            & ! ISUPPZ
!               &hf_work(d)%WORK,              & ! WORK
!               &26 * functions_in_domain(d),  & ! LWORK
!               &hf_work(d)%IWORK,             & ! IWORK
!               &10 * functions_in_domain(d),  & ! LIWORK
!               &INFO                          ) ! INFO

            IF (INFO.GT.0) THEN
                WRITE (error_message,'(a, i4, a)') &
                    & "Fock matrix transform in domain ", d, " failed"
                CALL print_error("hf_scf_cycle", error_message)
                exit_state = 2
                RETURN
            ELSEIF (INFO.LT.0) THEN
                WRITE (error_message,'(a, i4, a, i2, a)') &
                    & "Fock matrix transform in domain ", d, &
                    & "received incorrect ", abs(INFO), "th argument"
                CALL print_error("hf_scf_cycle", error_message)
                exit_state = 3
                RETURN
            ENDIF

            ! Back transform the eigenvectors
            CALL DGEMM('N', 'N',           & ! TRANSA, TRANSB
                &functions_in_domain(d),   & ! M
                &functions_kept(d),        & ! N
                &functions_kept(d),        & ! K
                &1.d0,                     & ! ALPHA
                &hf_one_e(d)%x_matrix(:,1:functions_kept(d)),   & ! A
                &functions_in_domain(d),   & ! LDA
                &hf_work(d)%matrix(:,:,1), & ! B
                &functions_in_domain(d),   & ! LDB
                &0.d0,                     & ! BETA
                &hf_in_domain(d)%orbitals, & ! C
                &functions_in_domain(d)    ) ! LDC

            ! Check convergence
            IF (curr_hf_cycle.GT.0) THEN
                hf_converged = hf_converged.AND.&
!                   &maxval(abs(hf_diis(d)%stored_error(:,:,1))).LT.scf_threshold
                    &check_convergence(d, .FALSE.)
            ELSE
                hf_converged = .FALSE.
            ENDIF

        ENDDO

        curr_hf_cycle = curr_hf_cycle + 1
        IF (.NOT.scan_job) CALL print_hf_scf_cycle

        ! If not converged, recur
        IF (.NOT.hf_converged.AND.(curr_hf_cycle.LT.max_scf_cycles)) &
            & CALL hf_scf_cycle(exit_state)

        ! Check if we're slowly increasing our basis and
        ! restart the SCF if we are
        IF (primary_basis_th.GT.0.AND.&
           &secondary_basis_th.GT.curr_basis_th) THEN

!write (*,*) "KICK!"
            curr_basis_th = curr_basis_th + 1
            temp_thresh = 10.d0**(-curr_basis_th)

            DO d = 1, n_domains
                functions_kept(d) = functions_in_domain(d)
                throw_count = 1
                DO WHILE (hf_one_e(d)%s_evals(throw_count).LT.basis_threshold)
                    functions_kept(d) = functions_kept(d) - 1
                    throw_count = throw_count + 1
                ENDDO
            ENDDO

            curr_hf_cycle = 0
            IF (curr_hf_cycle.LT.max_scf_cycles) CALL hf_scf_cycle(exit_state)

        ENDIF

        ! We've stopped iterating, so we build a new density
        ! matrix from the final orbitals
!       DO d = 1, n_domains
!           CALL build_density_matrix(d)
!       ENDDO

        IF (.NOT.hf_converged) exit_state = 1
    END SUBROUTINE hf_scf_cycle

    REAL(dp) FUNCTION hf_energy(id, exit_state)
        !
        ! This function calculates the current Hartree-Fock
        ! expectation value for the energy
        !
        ! INPUT:
        !   [int ] id
        ! Requests a specific energy component
        !   0 = Total energy
        !   1 = Kinetic energy
        !   2 = Nuclear attraction energy
        !   3 = Nuclear repulsion energy
        !   4 = Electron repulsion energy
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: id
        INTEGER, INTENT(OUT), OPTIONAL :: exit_state

        INTEGER :: d, i, n1, n2
        REAL(dp) :: energy
        CHARACTER(LEN=58) :: error_message
        LOGICAL :: electronic, halve_electronic, nuclear_repulsion

        electronic = .FALSE.
        halve_electronic = .FALSE.
        nuclear_repulsion = .FALSE.
        energy = 0.d0

        SELECT CASE (id)
            CASE (0) ! Total energy
                electronic = .TRUE.
                halve_electronic = .TRUE.
                nuclear_repulsion = .TRUE.

                DO d = 1, n_domains
                    hf_work(d)%matrix(:,:,1) = &
                        & one_e_in_domain(d)%kinetic - &
                        & one_e_in_domain(d)%potential + &
                        & hf_two_e(d)%fock
                ENDDO

            CASE (1) ! Kinetic energy
                electronic = .TRUE.

                DO d = 1, n_domains
                    hf_work(d)%matrix(:,:,1) = &
                        & one_e_in_domain(d)%kinetic
                ENDDO

            CASE (2) ! Nuclear attraction energy
                electronic = .TRUE.

                DO d = 1, n_domains
                    hf_work(d)%matrix(:,:,1) = &
                        & - one_e_in_domain(d)%potential
                ENDDO

            CASE (3) ! Nuclear repulsion energy
                nuclear_repulsion = .TRUE.

            CASE (4) ! Electron repulsion energy
                electronic = .TRUE.
                halve_electronic = .TRUE.

                DO d = 1, n_domains
                    hf_work(d)%matrix(:,:,1) = &
                        & hf_two_e(d)%g
                ENDDO

            CASE DEFAULT
                WRITE (error_message,'(a)') &
                    & "Requested energy type not recognised"
                CALL print_error("hf_energy", error_message)
                IF (present(exit_state)) exit_state = 1

        END SELECT

        IF (electronic) THEN
            DO d = 1, n_domains
                IF (electrons_in_domain(d).EQ.0) CYCLE

                CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
                    &functions_in_domain(d),    & ! M
                    &functions_in_domain(d),    & ! N
                    &functions_in_domain(d),    & ! K
                    &1.d0,                      & ! ALPHA
                    &hf_in_domain(d)%density,   & ! A
                    &functions_in_domain(d),    & ! LDA
                    &hf_work(d)%matrix(:,:,1),  & ! B
                    &functions_in_domain(d),    & ! LDB
                    &0.d0,                      & ! BETA
                    &hf_work(d)%matrix(:,:,2),  & ! C
                    &functions_in_domain(d)     ) ! LDC

                DO i = 1, functions_in_domain(d)
                    energy = energy + hf_work(d)%matrix(i,i,2)
                ENDDO
            ENDDO
        ENDIF

        IF (halve_electronic) energy = energy / 2.d0

        IF (nuclear_repulsion) THEN
            DO n1 = 1, n_nuclei
                DO n2 = n1 + 1, n_nuclei
                    energy = energy + (nuclear_charge(n1) * nuclear_charge(n2)) / &
                        & (nuclear_position(n2) - nuclear_position(n1))
                ENDDO
            ENDDO
        ENDIF

        hf_energy = energy
    END FUNCTION hf_energy

!====== COMPUTATION FUNCTIONS ================================================!

    SUBROUTINE build_density_matrix(d)
        !
        ! This routine constructs a density matrix within a
        ! domain from the current Hartree-Fock orbital 
        ! coefficients
        !
        ! INPUT :
        !   [int ]  d
        ! Index of the domain
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: d
        REAL(dp), POINTER :: occupied_orbitals(:,:)

        IF (functions_in_domain(d).NE.0) THEN
            occupied_orbitals => &
                &hf_in_domain(d)%orbitals(:,1:electrons_in_domain(d))

            CALL DGEMM('N',               & ! TRANSA
                &'T',                     & ! TRANSB
                &functions_in_domain(d),  & ! M
                &functions_in_domain(d),  & ! N
                &electrons_in_domain(d),  & ! K
                &1.d0,                    & ! ALPHA
                &occupied_orbitals,       & ! A
                &functions_in_domain(d),  & ! LDA
                &occupied_orbitals,       & ! B
                &functions_in_domain(d),  & ! LDB
                &0.d0,                    & ! BETA
                &hf_in_domain(d)%density, & ! C
                &functions_in_domain(d)   ) ! LDC
        ENDIF

    END SUBROUTINE build_density_matrix

    SUBROUTINE build_g_matrix
        !
        ! Constructs the matrix representing the mean-field
        ! two electron repulsion energy within the molecule
        ! Referred to as the G matrix in Szabo & Ostlund
        !

        INTEGER :: d, int_d, mu, nu, i, j
        INTEGER, POINTER :: length

        DO d = 1, n_domains
            FORALL (mu = 1:functions_in_domain(d), &
                   &nu = 1:functions_in_domain(d))
                hf_two_e(d)%g(mu,nu) = 0.d0
            END FORALL
            IF (electrons_in_domain(d).LE.0) CYCLE

            ! int_d is short for interacting_domain
            DO int_d = 1, n_domains
                IF (electrons_in_domain(int_d).LE.0) CYCLE
                length => functions_in_domain(int_d)

                DO mu = 1, functions_in_domain(d)
                    DO nu = 1, functions_in_domain(d)

                        CALL DSYMM('L', 'U', length, length, 1.d0, &
                            & hf_in_domain(int_d)%density, length, &
                            & double_bar(start_of_domain(d) + mu - 1, &
                                        &start_of_domain(int_d):end_of_domain(int_d), &
                                        &start_of_domain(d) + nu - 1, &
                                        &start_of_domain(int_d):end_of_domain(int_d)), &
                            & length, 0.d0, hf_work(int_d)%matrix(:,:,1), length)

                        DO i = 1, length
                            hf_two_e(d)%g(mu,nu) = hf_two_e(d)%g(mu,nu) + &
                                &hf_work(int_d)%matrix(i,i,1)
                        ENDDO

                    ENDDO
                ENDDO

!mu = 7
!nu = 1
!write (*,*) "DOMAIN ", d, " with ", int_d
!do i = 1, functions_in_domain(int_d)
!do j = 1, functions_in_domain(int_d)
!write (*,"(F10.6)",advance='no') double_bar(start_of_domain(d) + mu - 1, &
!    &start_of_domain(int_d) + i - 1, &
!    &start_of_domain(d) + nu - 1, &
!    &start_of_domain(int_d) + j - 1)
!enddo
!write (*,*)
!enddo
!write (*,*)

            ENDDO

        ENDDO

    END SUBROUTINE build_g_matrix

    SUBROUTINE build_diis_error(d, error)
        !
        ! Computes the error matrix used in the
        ! DIIS method for domain d
        !
        ! INPUT :
        !   [int ] d
        ! Working domain
        !
        ! OUTPUT :
        !   [real] error(:,:)
        ! DIIS error matrix
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN)   :: d
        REAL(dp), INTENT(OUT) :: error(:,:)

        ! Compute FPS
        CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
            &functions_in_domain(d),    & ! M
            &functions_in_domain(d),    & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &hf_two_e(d)%fock,          & ! A
            &functions_in_domain(d),    & ! LDA
            &hf_in_domain(d)%density,   & ! B
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &hf_work(d)%matrix(:,:,1),  & ! C
            &functions_in_domain(d)     ) ! LDC

        CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
            &functions_in_domain(d),    & ! M
            &functions_in_domain(d),    & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &hf_work(d)%matrix(:,:,1),  & ! A
            &functions_in_domain(d),    & ! LDA
            &one_e_in_domain(d)%overlap,& ! B
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &hf_diis(d)%fps,            & ! C
            &functions_in_domain(d)     ) ! LDC

        ! Compute SPF
        CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
            &functions_in_domain(d),    & ! M
            &functions_in_domain(d),    & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &one_e_in_domain(d)%overlap,& ! B
            &functions_in_domain(d),    & ! LDA
            &hf_in_domain(d)%density,   & ! A
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &hf_work(d)%matrix(:,:,1),  & ! C
            &functions_in_domain(d)     ) ! LDC

        CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
            &functions_in_domain(d),    & ! M
            &functions_in_domain(d),    & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &hf_work(d)%matrix(:,:,1),  & ! A
            &functions_in_domain(d),    & ! LDA
            &hf_two_e(d)%fock,          & ! B
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &hf_diis(d)%spf,            & ! C
            &functions_in_domain(d)     ) ! LDC

        ! Pack into the output
        error = hf_diis(d)%fps - hf_diis(d)%spf

    END SUBROUTINE build_diis_error

    SUBROUTINE build_diis_matrix(d, diis_matrix)
        !
        ! Computes the DIIS matrix for the given domain
        !
        ! INPUT :
        !   [int ] d
        ! Working domain
        !
        ! OUTPUT :
        !   [real] diis_matrix(:,:)
        ! Contains the DIIS matrix for the given domain
        ! on exit
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: d
        REAL(dp), INTENT(OUT) :: diis_matrix(:,:)

        INTEGER :: i, j, k

        DO i = 1, curr_diis
            DO j = 1, curr_diis

                diis_matrix(i,j) = 0.d0

                CALL DGEMM('N', 'T',                & ! TRANSA, TRANSB
                    &functions_kept(d),             & ! M
                    &functions_kept(d),             & ! N
                    &functions_kept(d),             & ! K
                    &1.d0,                          & ! ALPHA
                    &hf_diis(d)%stored_error(:,:,i),& ! A
                    &functions_in_domain(d),        & ! LDA
                    &hf_diis(d)%stored_error(:,:,j),& ! B
                    &functions_in_domain(d),        & ! LDB
                    &0.d0,                          & ! BETA
                    &hf_work(d)%matrix(:,:,1),      & ! C
                    &functions_in_domain(d)         ) ! LDC

                DO k = 1, functions_in_domain(d)
                    diis_matrix(i,j) = diis_matrix(i,j) + &
                                     & hf_work(d)%matrix(k,k,1)
                ENDDO

            ENDDO

            j = curr_diis + 1
            diis_matrix(i,j) = -1.d0
            diis_matrix(j,i) = -1.d0
        ENDDO

        diis_matrix(curr_diis + 1, curr_diis + 1) = 0.d0

    END SUBROUTINE build_diis_matrix

    SUBROUTINE orthogonalize(d, matrix)
        !
        ! Transforms the supplied matrix into an
        ! orthogonalized basis using the X matrix
        ! computed in hf_one_electron
        !
        ! INPUT :
        !   [int ] d
        ! Working domain
        !
        ! IN/OUT :
        !   [real] matrix
        ! The matrix to be orthogonalized
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: d
        REAL(dp), INTENT(INOUT) :: matrix(:,:)

        REAL(dp), POINTER :: x_truncated(:,:)

        x_truncated => hf_one_e(d)%x_matrix(:,1:functions_kept(d))

        CALL DGEMM('T', 'N',            & ! TRANSA, TRANSB
            &functions_kept(d),         & ! M
            &functions_in_domain(d),    & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &x_truncated,               & ! A
            &functions_in_domain(d),    & ! LDA
            &matrix,                    & ! B
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &hf_work(d)%matrix(:,:,1),  & ! C
            &functions_in_domain(d)     ) ! LDC

        CALL DGEMM('N', 'N',            & ! TRANSA, TRANSB
            &functions_in_domain(d),    & ! M
            &functions_kept(d),         & ! N
            &functions_in_domain(d),    & ! K
            &1.d0,                      & ! ALPHA
            &hf_work(d)%matrix(:,:,1),  & ! A
            &functions_in_domain(d),    & ! LDA
            &x_truncated,               & ! B
            &functions_in_domain(d),    & ! LDB
            &0.d0,                      & ! BETA
            &matrix,                    & ! C
            &functions_in_domain(d)     ) ! LDC

    END SUBROUTINE orthogonalize

    LOGICAL FUNCTION check_convergence(d, near)
        !
        ! Determines whether an SCF cycle has converged 
        ! within a domain
        !
        ! INPUT :
        !   [int ] d
        ! Working domain
        !   [bool] near
        ! If .TRUE. then check if we're close to convergence
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: d
        LOGICAL, INTENT(IN) :: near

        INTEGER :: i, j
        REAL(dp) :: max_error, rms

        SELECT CASE (conv_check_type)
            CASE (0) ! Max value in the current error matrix
                max_error = maxval(abs(hf_diis(d)%stored_error(&
                    &:functions_kept(d),:functions_kept(d),1)))
                IF (near) THEN
                    check_convergence = max_error.LT.(10*scf_threshold)
                ELSE
                    check_convergence = max_error.LT.scf_threshold
                ENDIF

            CASE (1) ! RMS of the current error matrix
                rms = 0.d0
                DO i = 1, functions_kept(d)
                    DO j = 1, functions_kept(d)
                        rms = rms + hf_diis(d)%stored_error(i,j,1)**2
                    ENDDO
                ENDDO

                rms = sqrt(rms / functions_kept(d)**2)
                IF (near) THEN
                    check_convergence = rms.LT.(10*scf_threshold)
                ELSE
                    check_convergence = rms.LT.scf_threshold
                ENDIF

        END SELECT

    END FUNCTION check_convergence

!====== UTILITIES ============================================================!

    SUBROUTINE init_hartree_fock(exit_state)
        !
        ! This subroutine initialises storage space for
        ! quantities required throughout the SCF procedure
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Top level storage could not be allocated
        !   3 = Domain level storage could not be allocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, i, j, stat
        INTEGER, POINTER :: m
        CHARACTER(LEN=58) :: error_message

        ALLOCATE ( hf_one_e(n_domains), &
                 & hf_two_e(n_domains), &
                 & hf_diis(n_domains),  &
                 & hf_work(n_domains),  &
                 & STAT=stat)

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Could not allocate derived types"
            CALL print_error("init_hartree_fock", error_message)
            exit_state = 2
            RETURN
        ENDIF

        DO d = 1, n_domains
            m => functions_in_domain(d)

            ! Allocate memory for one electron objects
            ALLOCATE ( hf_one_e(d)%s_evals(m),      &
                     & hf_one_e(d)%u_matrix(m,m),   &
                     & hf_one_e(d)%x_matrix(m,m),   &
                     & STAT=stat)

            ! Allocate memory for two electron objects
            ALLOCATE ( hf_two_e(d)%g(m,m),          &
                     & hf_two_e(d)%fock(m,m),       &
                     & hf_two_e(d)%fock_p(m,m),     &
                     & STAT=stat)

            ! Allocate memory for DIIS procedure
            ALLOCATE ( hf_diis(d)%stored_fock(m,m,diis_length), &
                     & hf_diis(d)%stored_error(m,m,diis_length), &
                     & hf_diis(d)%diis_matrix(diis_length+1,diis_length+1), &
                     & hf_diis(d)%diis_vector(diis_length+1), &
                     & hf_diis(d)%spf(m,m),         &
                     & hf_diis(d)%fps(m,m),         &
                     & STAT=stat)

            ! Allocate working space
            ALLOCATE ( hf_work(d)%matrix(m,m,2),        &
                     & hf_work(d)%ISUPPZ(2 * m),        &
                     & hf_work(d)%WORK(26 * m),         &
                     & hf_work(d)%IWORK(10 * m),        &
                     & hf_work(d)%IPIV(diis_length + 1),&
                     & STAT=stat)

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a, i4)') &
                    & "Could not allocate storage for domain ", d
                CALL print_error("init_hartree_fock", error_message)
                exit_state = 3
                RETURN
            ENDIF

            ! Generate orbital guess

            ! Null guess
            FORALL (i=1:m, j=1:m)
                hf_in_domain(d)%orbitals(i,j) = 0.d0
            END FORALL

            ! Construct initial density
            CALL build_density_matrix(d)
        ENDDO

        curr_hf_cycle = 0
    END SUBROUTINE init_hartree_fock

    SUBROUTINE clean_up_hf(exit_state)
        !
        ! This subroutine cleans up the scratch space
        ! allocated by init_hartree_fock
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Not all memory deallocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, stat
        CHARACTER(LEN=80) :: error_message

        DO d = 1, n_domains

            ! Deallocate memory for one electron objects
            DEALLOCATE ( hf_one_e(d)%s_evals,   &
                       & hf_one_e(d)%u_matrix,  &
                       & hf_one_e(d)%x_matrix,  &
                       & STAT=stat)

            ! Deallocate memory for two electron objects
            DEALLOCATE ( hf_two_e(d)%g,         &
                       & hf_two_e(d)%fock,      &
                       & hf_two_e(d)%fock_p,    &
                       & STAT=stat)

            ! Deallocate memory for DIIS procedure
            DEALLOCATE ( hf_diis(d)%stored_fock,    &
                       & hf_diis(d)%stored_error,   &
                       & hf_diis(d)%diis_matrix,    &
                       & hf_diis(d)%diis_vector,    &
                       & hf_diis(d)%spf,            &
                       & hf_diis(d)%fps,            &
                       & STAT=stat)

            ! Deallocate working space
            DEALLOCATE ( hf_work(d)%matrix, &
                       & hf_work(d)%ISUPPZ, &
                       & hf_work(d)%WORK,   &
                       & hf_work(d)%IWORK,  &
                       & hf_work(d)%IPIV,   &
                       & STAT=stat)


            IF (stat.NE.0) THEN
                WRITE (error_message,'(a, i4)') &
                    & "Could not deallocate storage for domain ", d
                CALL print_error("clean_up_hf", error_message)
                exit_state = 2
                RETURN
            ENDIF

        ENDDO

        DEALLOCATE ( hf_one_e, hf_two_e, hf_diis, hf_work, STAT=stat)

        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Could not deallocate derived types"
            CALL print_error("clean_up_hf", error_message)
            exit_state = 2
            RETURN
        ENDIF

    END SUBROUTINE clean_up_hf

    SUBROUTINE print_hf_start
        IMPLICIT NONE

        WRITE (*,"(6('='), (' '), a, (' '), 40('='))") "HARTREE-FOCK"
        WRITE (*,*)
        WRITE (*,"(5(' '), ('--- '), a, (' ---'))") "START OF SCF"
        WRITE (*,"(2(' '), a, 5(' '), a)") "Cycle", "Energy"

    END SUBROUTINE

    SUBROUTINE print_hf_scf_cycle
        IMPLICIT NONE

        INTEGER :: d, i, j
        REAL :: error, rms

        SELECT CASE (conv_check_type)
            CASE (0) ! Max value in the current error matrix
                error = 0.d0
                DO d = 1, n_domains
                    IF (electrons_in_domain(d).EQ.0) CYCLE
                    error = max(error, maxval(abs(hf_diis(d)%stored_error(&
                        &:functions_kept(d),:functions_kept(d),1))))
                ENDDO

            CASE (1) ! RMS of the current error matrix
                error = 0.d0
                DO d = 1, n_domains
                    IF (electrons_in_domain(d).EQ.0) CYCLE
                    rms = 0.d0
                    DO i = 1, functions_kept(d)
                        DO j = 1, functions_kept(d)
                            rms = rms + hf_diis(d)%stored_error(i,j,1)**2
                        ENDDO
                    ENDDO

                    rms = sqrt(rms / functions_kept(d)**2)
                    error = max(error, rms)
                ENDDO

        END SELECT

        WRITE (*,"(3(' '), i3, 6(' '), F22.15, ('  '), F22.15)") &
            & curr_hf_cycle, hf_energy(0), error

    END SUBROUTINE

    SUBROUTINE print_hf_end
        IMPLICIT NONE

        INTEGER :: d, e, f

        WRITE (*,"(5(' '), ('---  '), a, ('  ---'))") "END OF SCF"

        IF (hf_converged) THEN
            WRITE (*,"(2(' '), a, i3, a)") "Calculation has CONVERGED after ", &
                      & curr_hf_cycle, " iterations"
        ELSE
            WRITE (*,"(2(' '), a)") "Calculation has NOT CONVERGED"
        ENDIF

        WRITE (*,*)
        WRITE (*,"(5(' '), ('--- '), a, (' ---'))") "RESULTS"
        WRITE (*,"(2(' '), a, F22.15)") "Kinetic            = ", hf_energy(1)
        WRITE (*,"(2(' '), a, F22.15)") "Nuclear attraction = ", hf_energy(2)
        WRITE (*,"(2(' '), a, F22.15)") "Nuclear repulsion  = ", hf_energy(3)
        WRITE (*,"(2(' '), a, F22.15)") "Electron repulsion = ", hf_energy(4)
        WRITE (*,*)
        WRITE (*,"(2(' '), a, F22.15)") "Total energy       = ", hf_energy(0)
        WRITE (*,*)

        IF (print_orbitals) THEN
            WRITE (*,"(5(' '), ('--- '), a, (' ---'))") "ORBITALS"
            DO d = 1, n_domains
                IF (electrons_in_domain(d).EQ.0) CYCLE
                WRITE (*,"(2(' '), a, i3)") "Domain ", d
                DO f = 1, functions_in_domain(d)
                    DO e = 1, electrons_in_domain(d)

                        WRITE (*,"(2(' '), F22.15)",ADVANCE="NO") &
                            & hf_in_domain(d)%orbitals(f, e)

                    ENDDO
                    WRITE (*,*)
                ENDDO
                WRITE (*,*)
            ENDDO
        ENDIF

    END SUBROUTINE

END MODULE hartree_fock
