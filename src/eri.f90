MODULE eri
    !
    ! ERI
    !
    ! This module manages the construction of all
    ! electron repulstion integrals
    !
    USE constants
    USE input
    USE error
    USE storage,       ONLY : one_e_in_domain, &
                            & eri_of_domain,   &
                            & double_bar,      &
                            & hf_in_domain,    &
                            & mo_eri_domain,   &
                            & mo_double_bar
    USE true_eri_exp,  ONLY : true_llrr,       &
                            & true_llee
    USE true_eri_poly, ONLY : true_eeee
    USE quasi_eri,     ONLY : quasi_llll,      &
                            & quasi_poly

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE :: scratch_matrix(:,:)
    REAL(dp), ALLOCATABLE :: scratch_vector(:)

    PRIVATE
    PUBLIC :: build_eri
    PUBLIC :: build_double_bar_eri
    PUBLIC :: mo_eri_transform
    PUBLIC :: build_mo_double_bar

    CONTAINS

    SUBROUTINE build_eri(exit_state)
        !
        ! This subroutine constructs the full set of electron
        ! repulsion integrals for all domains
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Scratch memory space could not be allocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: d, d1, d2, max_m, stat
        INTEGER :: m, n, l, s
        REAL(dp) :: A, B, R, max_integral
        CHARACTER(LEN=58) :: error_message

        IF (integral_check) THEN
            WRITE (*,*) "Checking integrals"

            DO d1 = 1, n_domains
            DO d2 = 1, n_domains
                FORALL (m=1:functions_in_domain(d1),n=1:functions_in_domain(d1),&
                       &l=1:functions_in_domain(d2),s=1:functions_in_domain(d2))
                    eri_of_domain(d1)%with(d2)%integral(m,n,l,s) = 0.d0
                END FORALL
            ENDDO
            ENDDO
        ENDIF

        max_m = 0
        DO d = 1, n_domains
            max_m = max(max_m, 2 * evens_in_domain(d), &
                    &functions_in_domain(d), 2 * odds_in_domain(d))
        ENDDO
        max_m = max_m + 1

        IF (.NOT.allocated(scratch_matrix)) THEN
            ALLOCATE(scratch_matrix(-1:max_m,-1:max_m),STAT=stat)
            FORALL (m = -1:max_m, n = -1:max_m)
                scratch_matrix(m,n) = 0.d0
            ENDFORALL

            IF (stat.NE.0) THEN
                WRITE (error_message,'(a)') &
                    & "Could not allocate scratch space"
                CALL print_error("build_eri", error_message)
                exit_state = 2
                RETURN
            ENDIF
        ELSE

            FORALL (m = -1:max_m, n = -1:max_m)
                scratch_matrix(m,n) = 0.d0
            ENDFORALL

        ENDIF

        ! Handle the exterior domains seperately to avoid
        ! IF blocks within loops
        d1 = 1
        ! Quasi-integrals (LLLL)
        IF (integral_check) THEN
            DO m = 1, functions_in_domain(d1)
            DO n = 1, functions_in_domain(d1)
            DO l = 1, functions_in_domain(d1)
            DO s = 1, functions_in_domain(d1)
                IF (eri_of_domain(d1)%with(d1)%integral(m,n,l,s).NE.0.d0) THEN
                    CALL print_error("eri", "overlapping integral in LLLL")
                ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
        ENDIF

        CALL quasi_llll(functions_in_domain(1), alpha, &
                       &one_e_in_domain(1)%overlap, &
                       &eri_of_domain(1)%with(1)%integral)

        ! True integrals (LLEE / LLEO / LLOO)
        DO d2 = 2, n_domains - 1
            A = (nuclear_position(d2) - nuclear_position(d2 - 1)) / 2
            R = nuclear_position(d2 - 1) + A - nuclear_position(d1)

            ! Call the routine which will compute and pack up
            ! all the integrals
            IF (integral_check) THEN
                DO m = 1, functions_in_domain(d1)
                DO n = 1, functions_in_domain(d1)
                DO l = 1, functions_in_domain(d2)
                DO s = 1, functions_in_domain(d2)
                    IF (eri_of_domain(d1)%with(d2)%integral(m,n,l,s).NE.0.d0) THEN
                        CALL print_error("eri", "overlapping integral in LLPP")
                    ENDIF
                    IF (eri_of_domain(d2)%with(d1)%integral(l,s,m,n).NE.0.d0) THEN
                        CALL print_error("eri", "overlapping integral in PPLL")
                    ENDIF
                ENDDO
                ENDDO
                ENDDO
                ENDDO
            ENDIF

            CALL true_llee(.TRUE., A, R, evens_in_domain(d2), odds_in_domain(d2), &
                          &functions_in_domain(d1), &
                          &one_e_in_domain(d1)%overlap, &
                          &one_e_in_domain(d2)%overlap, &
                          &eri_of_domain(d1)%with(d2)%integral, &
                          &eri_of_domain(d2)%with(d1)%integral)
        ENDDO

        ! LLRR
        d2 = n_domains
        A = nuclear_position(d1)
        B = nuclear_position(d2 - 1)
        R = B - A

        ! Call the routine which will compute and pack up
        ! all the integrals
        IF (integral_check) THEN
            DO m = 1, functions_in_domain(d1)
            DO n = 1, functions_in_domain(d1)
            DO l = 1, functions_in_domain(d2)
            DO s = 1, functions_in_domain(d2)
                IF (eri_of_domain(d1)%with(d2)%integral(m,n,l,s).NE.0.d0) THEN
                    CALL print_error("eri", "overlapping integral in LLRR")
                ENDIF
                IF (eri_of_domain(d2)%with(d1)%integral(l,s,m,n).NE.0.d0) THEN
                    CALL print_error("eri", "overlapping integral in RRLL")
                ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
        ENDIF

        CALL true_llrr(R, functions_in_domain(d1), functions_in_domain(d2), &
                      &one_e_in_domain(d1)%overlap, &
                      &one_e_in_domain(d2)%overlap, &
                      &scratch_matrix, &
                      &eri_of_domain(d1)%with(d2)%integral, &
                      &eri_of_domain(d2)%with(d1)%integral)

        ! Quasi-integrals (RRRR)
        d1 = n_domains
        IF (integral_check) THEN
            DO m = 1, functions_in_domain(d1)
            DO n = 1, functions_in_domain(d1)
            DO l = 1, functions_in_domain(d1)
            DO s = 1, functions_in_domain(d1)
                IF (eri_of_domain(d1)%with(d1)%integral(m,n,l,s).NE.0.d0) THEN
                    CALL print_error("eri", "overlapping integral in RRRR")
                ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
        ENDIF

        CALL quasi_llll(functions_in_domain(d1), alpha, &
                       &one_e_in_domain(d1)%overlap, &
                       &eri_of_domain(d1)%with(d1)%integral)

        ! True integrals (RREE / RREO / RROO)
        DO d2 = 2, n_domains - 1
            A = (nuclear_position(d2) - nuclear_position(d2 - 1)) / 2
            R = nuclear_position(n_nuclei) - nuclear_position(d2 - 1) - A

            ! Call the routine which will compute and pack up
            ! all the integrals
            IF (integral_check) THEN
                DO m = 1, functions_in_domain(d1)
                DO n = 1, functions_in_domain(d1)
                DO l = 1, functions_in_domain(d2)
                DO s = 1, functions_in_domain(d2)
                    IF (eri_of_domain(d1)%with(d2)%integral(m,n,l,s).NE.0.d0) THEN
                        CALL print_error("eri", "overlapping integral in RRPP")
                    ENDIF
                    IF (eri_of_domain(d2)%with(d1)%integral(l,s,m,n).NE.0.d0) THEN
                        CALL print_error("eri", "overlapping integral in PPRR")
                    ENDIF
                ENDDO
                ENDDO
                ENDDO
                ENDDO
            ENDIF

            CALL true_llee(.FALSE., A, R, evens_in_domain(d2), odds_in_domain(d2), &
                          &functions_in_domain(d1), &
                          &one_e_in_domain(d1)%overlap, &
                          &one_e_in_domain(d2)%overlap, &
                          &eri_of_domain(d1)%with(d2)%integral, &
                          &eri_of_domain(d2)%with(d1)%integral)
        ENDDO



        ! Now compute the interactions between finite
        ! domains
        DO d1 = 2, n_domains - 1
            ! Quasi-integrals
            A = (nuclear_position(d1) - nuclear_position(d1 - 1)) / 2
            IF (integral_check) THEN
                DO m = 1, functions_in_domain(d1)
                DO n = 1, functions_in_domain(d1)
                DO l = 1, functions_in_domain(d1)
                DO s = 1, functions_in_domain(d1)
                    IF (eri_of_domain(d1)%with(d1)%integral(m,n,l,s).NE.0.d0) THEN
                        CALL print_error("eri", "overlapping integral in quasi PPPP")
                    ENDIF
                ENDDO
                ENDDO
                ENDDO
                ENDDO
            ENDIF


            CALL quasi_poly(A, evens_in_domain(d1), odds_in_domain(d1), &
                           &one_e_in_domain(d1)%overlap, &
                           &eri_of_domain(d1)%with(d1)%integral)

            ! True integrals
            DO d2 = d1 + 1, n_domains - 1
                B = (nuclear_position(d2) - nuclear_position(d2 - 1)) / 2
                R = (nuclear_position(d2 - 1) + B) - &
                  & (nuclear_position(d1 - 1) + A)

                max_m = max(evens_in_domain(d1) + evens_in_domain(d1), &
                           &evens_in_domain(d1) + odds_in_domain(d1),  &
                           &odds_in_domain(d1) + odds_in_domain(d1),   &
                           &evens_in_domain(d2) + evens_in_domain(d2), &
                           &evens_in_domain(d2) + odds_in_domain(d2),  &
                           &odds_in_domain(d2) + odds_in_domain(d2)    )

                ! Call the routine which will compute and pack up
                ! all the integrals
                IF (integral_check) THEN
                    DO m = 1, functions_in_domain(d1)
                    DO n = 1, functions_in_domain(d1)
                    DO l = 1, functions_in_domain(d2)
                    DO s = 1, functions_in_domain(d2)
                        IF (eri_of_domain(d1)%with(d2)%integral(m,n,l,s).NE.0.d0) THEN
                            CALL print_error("eri", "overlapping integral in true PPPP")
                        ENDIF
                        IF (eri_of_domain(d2)%with(d1)%integral(l,s,m,n).NE.0.d0) THEN
                            CALL print_error("eri", "overlapping integral in true PPPP")
                        ENDIF
                    ENDDO
                    ENDDO
                    ENDDO
                    ENDDO
                ENDIF

                CALL true_eeee(A, B, R, &
                              &evens_in_domain(d1), odds_in_domain(d1),     &
                              &evens_in_domain(d2), odds_in_domain(d2),     &
                              &max_m, one_e_in_domain(d1)%overlap,          &
                              &one_e_in_domain(d2)%overlap, scratch_matrix, &
                              &eri_of_domain(d1)%with(d2)%integral,         &
                              &eri_of_domain(d2)%with(d1)%integral          )
            ENDDO
        ENDDO

        IF (allocated(scratch_matrix)) THEN
            DEALLOCATE(scratch_matrix)
            IF (stat.NE.0) THEN
                WRITE (error_message,'(a)') &
                    & "Could not deallocate scratch space"
                CALL print_warning("build_eri", error_message)
                exit_state = 1
            ENDIF
        ENDIF

        IF (debug) THEN
            DO d1 = 1, n_domains
                DO d2 = d1, n_domains
                    max_integral = 0.d0

                    DO m = 1, functions_in_domain(d1)
                    DO n = 1, functions_in_domain(d1)
                    DO l = 1, functions_in_domain(d2)
                    DO s = 1, functions_in_domain(d2)
                        max_integral = max(max_integral, &
                            &abs(eri_of_domain(d1)%with(d2)%integral(m,n,l,s)))
                    ENDDO
                    ENDDO
                    ENDDO
                    ENDDO

                    WRITE (*,*) "DOMAINS", d1, " and ", d2, " : ", max_integral
                ENDDO
            ENDDO
        ENDIF

    END SUBROUTINE build_eri

    SUBROUTINE build_double_bar_eri
        !
        ! This subroutine builds the full set of double
        ! bar integrals for all domains
        !
        ! NOTE: Now in physicists notation!
        !
        IMPLICIT NONE

        INTEGER :: di, dj, dk, dl, i, j, k, l, shifti, shiftj, shiftk, shiftl

        FORALL (i=1:end_of_domain(n_domains), &
               &j=1:end_of_domain(n_domains), &
               &k=1:end_of_domain(n_domains), &
               &l=1:end_of_domain(n_domains))
            double_bar(i,j,k,l) = 0.d0
        END FORALL

        DO di = 1, n_domains
        DO dj = 1, n_domains
        DO dk = 1, n_domains
        DO dl = 1, n_domains
            IF (di.EQ.dk.AND.dj.EQ.dl.AND.di.EQ.dj) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    IF (integral_check) THEN
                        IF (double_bar(i,j,k,l).NE.0.d0) THEN
                            CALL print_error("build_double_bar_eri", &
                                &"overlapping integrals encountered (section 1)")
write (*,*) di, "<", i, j, "|", k, l, ">"
                        ENDIF
                    ENDIF

                    double_bar(i, j, k, l) = &
                        &eri_of_domain(di)%with(di)%integral(i-shifti,k-shiftk,j-shiftj,l-shiftl) - &
                        &eri_of_domain(di)%with(di)%integral(i-shifti,l-shiftl,j-shiftj,k-shiftk)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ELSEIF (di.EQ.dk.AND.dj.EQ.dl.AND.di.NE.dj) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    IF (integral_check) THEN
                        IF (double_bar(i,j,k,l).NE.0.d0) THEN
                            CALL print_error("build_double_bar_eri", &
                                &"overlapping integrals encountered (section 1)")
write (*,*) di, "<", i, j, "|", k, l, ">"
                        ENDIF
                    ENDIF

                    double_bar(i, j, k, l) = &
                        &eri_of_domain(di)%with(dj)%integral(i-shifti,k-shiftk,j-shiftj,l-shiftl)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ELSEIF (di.EQ.dl.AND.dj.EQ.dk) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    IF (integral_check) THEN
                        IF (double_bar(i,j,k,l).NE.0.d0) THEN
                            CALL print_error("build_double_bar_eri", &
                                &"overlapping integrals encountered (section 1)")
write (*,*) di, "<", i, j, "|", k, l, ">"
                        ENDIF
                    ENDIF

                    double_bar(i, j, k, l) = &
                        & - eri_of_domain(di)%with(dj)%integral(i-shifti,l-shiftl,j-shiftj,k-shiftk)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ENDIF
        ENDDO
        ENDDO
        ENDDO
        ENDDO

    END SUBROUTINE build_double_bar_eri

    SUBROUTINE mo_eri_transform(exit_state)
        !
        ! This subroutine transforms previously calculated
        ! electron repulsion integrals into the molecular
        ! orbital basis.
        !
        ! Currently only transforms to the Hartree-Fock
        ! molecular orbitals
        !
        ! OUTPUT :
        !   [int ] exit_state
        ! Describes the exit condition of the routine
        !   0 = Procedure terminated successfully
        !   1 = Non-fatal error encountered
        !   2 = Scratch memory space could not be allocated
        !
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: exit_state

        INTEGER :: max_f, d1, d2, m, n, l, s, stat
        REAL(dp), ALLOCATABLE :: temp_matrix(:,:)
        CHARACTER(LEN=58) :: error_message

        max_f = maxval(functions_in_domain)

        ALLOCATE(temp_matrix(max_f, max_f), STAT=stat)
        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Could not allocate scratch space"
            CALL print_error("mo_eri_transform", error_message)
            exit_state = 2
            RETURN
        ENDIF

        DO d1 = 1, n_domains
            IF (functions_in_domain(d1).EQ.0) CYCLE
            DO d2 = 1, n_domains
                IF (functions_in_domain(d2).EQ.0) CYCLE

                DO l = 1, functions_in_domain(d2)
                DO s = 1, functions_in_domain(d2)

                    CALL DGEMM('T',                                 & ! TRANSA
                        &'N',                                       & ! TRANSB
                        &functions_in_domain(d1),                   & ! M
                        &functions_in_domain(d1),                   & ! N
                        &functions_in_domain(d1),                   & ! K
                        &1.d0,                                      & ! ALPHA
                        &hf_in_domain(d1)%orbitals,                 & ! A
                        &functions_in_domain(d1),                   & ! LDA
                        &eri_of_domain(d1)%with(d2)%integral(:,:,l,s), & ! B
                        &functions_in_domain(d1),                   & ! LDB
                        &0.d0,                                      & ! BETA
                        &temp_matrix,                               & ! C
                        &max_f                                      ) ! LDC

                    CALL DGEMM('N',                                 & ! TRANSA
                        &'N',                                       & ! TRANSB
                        &functions_in_domain(d1),                   & ! M
                        &functions_in_domain(d1),                   & ! N
                        &functions_in_domain(d1),                   & ! K
                        &1.d0,                                      & ! ALPHA
                        &temp_matrix,                               & ! A
                        &max_f,                                     & ! LDA
                        &hf_in_domain(d1)%orbitals,                 & ! B
                        &functions_in_domain(d1),                   & ! LDB
                        &0.d0,                                      & ! BETA
                        &mo_eri_domain(d1)%with(d2)%integral(:,:,l,s), & ! C
                        &functions_in_domain(d1)                    ) ! LDC

                ENDDO
                ENDDO

                DO m = 1, functions_in_domain(d1)
                DO n = 1, functions_in_domain(d1)

                    CALL DGEMM('T',                                 & ! TRANSA
                        &'N',                                       & ! TRANSB
                        &functions_in_domain(d2),                   & ! M
                        &functions_in_domain(d2),                   & ! N
                        &functions_in_domain(d2),                   & ! K
                        &1.d0,                                      & ! ALPHA
                        &hf_in_domain(d2)%orbitals,                 & ! A
                        &functions_in_domain(d2),                   & ! LDA
                        &mo_eri_domain(d1)%with(d2)%integral(m,n,:,:), & ! B
                        &functions_in_domain(d2),                   & ! LDB
                        &0.d0,                                      & ! BETA
                        &temp_matrix,                               & ! C
                        &maxval(functions_in_domain)                ) ! LDC

                    CALL DGEMM('N',                                 & ! TRANSA
                        &'N',                                       & ! TRANSB
                        &functions_in_domain(d2),                   & ! M
                        &functions_in_domain(d2),                   & ! N
                        &functions_in_domain(d2),                   & ! K
                        &1.d0,                                      & ! ALPHA
                        &temp_matrix,                               & ! A
                        &maxval(functions_in_domain),               & ! LDA
                        &hf_in_domain(d2)%orbitals,                 & ! B
                        &functions_in_domain(d2),                   & ! LDB
                        &0.d0,                                      & ! BETA
                        &mo_eri_domain(d1)%with(d2)%integral(m,n,:,:), & ! C
                        &functions_in_domain(d2)                    ) ! LDC

                ENDDO
                ENDDO

            ENDDO
        ENDDO

        DEALLOCATE(temp_matrix, STAT=stat)
        IF (stat.NE.0) THEN
            WRITE (error_message,'(a)') &
                & "Could not deallocate scratch space"
            CALL print_warning("mo_eri_transform", error_message)
            exit_state = 1
        ENDIF

    END SUBROUTINE mo_eri_transform

    SUBROUTINE build_mo_double_bar
        !
        ! This subroutine builds the full set of double
        ! bar integrals in the MO basis for all domains
        !
        ! NOTE: Unlike the single bar integrals, these
        !       integrals use physicist's notation
        !
        IMPLICIT NONE

        INTEGER :: di, dj, dk, dl, i, j, k, l, shifti, shiftj, shiftk, shiftl

        FORALL (i=1:end_of_domain(n_domains), &
               &j=1:end_of_domain(n_domains), &
               &k=1:end_of_domain(n_domains), &
               &l=1:end_of_domain(n_domains))
            mo_double_bar(i,j,k,l) = 0.d0
        END FORALL

        DO di = 1, n_domains
        DO dj = 1, n_domains
        DO dk = 1, n_domains
        DO dl = 1, n_domains
            IF (di.EQ.dk.AND.dj.EQ.dl.AND.di.EQ.dj) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    mo_double_bar(i, j, k, l) = &
                        &mo_eri_domain(di)%with(di)%integral(i-shifti,k-shiftk,j-shiftj,l-shiftl) - &
                        &mo_eri_domain(di)%with(di)%integral(i-shifti,l-shiftl,j-shiftj,k-shiftk)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ELSEIF (di.EQ.dk.AND.dj.EQ.dl.AND.di.NE.dj) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    mo_double_bar(i, j, k, l) = &
                        &mo_eri_domain(di)%with(dj)%integral(i-shifti,k-shiftk,j-shiftj,l-shiftl)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ELSEIF (di.EQ.dl.AND.dj.EQ.dk) THEN
                shifti = start_of_domain(di) - 1
                shiftj = start_of_domain(dj) - 1
                shiftk = start_of_domain(dk) - 1
                shiftl = start_of_domain(dl) - 1

                DO i = start_of_domain(di), end_of_domain(di)
                DO j = start_of_domain(dj), end_of_domain(dj)
                DO k = start_of_domain(dk), end_of_domain(dk)
                DO l = start_of_domain(dl), end_of_domain(dl)

                    mo_double_bar(i, j, k, l) = &
                        & - mo_eri_domain(di)%with(dj)%integral(i-shifti,l-shiftl,j-shiftj,k-shiftk)

                ENDDO
                ENDDO
                ENDDO
                ENDDO

            ENDIF
        ENDDO
        ENDDO
        ENDDO
        ENDDO

    END SUBROUTINE build_mo_double_bar

END MODULE eri
