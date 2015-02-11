MODULE quasi_eri
    !
    ! QUASI_ERI
    !
    ! This module contains subroutines for computing electron
    ! repulsion integrals when the two electrons inhabit 
    ! the same domain
    !
    USE constants
    USE error
    USE one_e_integrals, ONLY : overlap_integral
    USE quasi_data
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: quasi_llll
    PUBLIC :: quasi_poly

    CONTAINS

    SUBROUTINE quasi_llll(functions, alpha, overlap, integrals)
        !
        ! This subroutine computes fake electron repulsion 
        ! integrals for electrons occupying the same infinite 
        ! domain. These values produce the correct double bar
        ! integral when appropriately combined
        !
        ! INPUT :
        !   [int ] functions
        ! Maximum index of the basis functions in this domain
        !   [real] alpha
        ! Base exponent for the exponential basis functions
        ! in this domain
        !   [real] overlap(:,:)
        ! Overlap matrix for this domain
        !
        ! OUTPUT :
        !   [real] integrals(:,:,:,:)
        ! Memory space for the resulting electron repulsion
        ! integrals
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN)   :: functions
        REAL(dp), INTENT(IN)  :: alpha, overlap(functions, functions)
        REAL(dp), INTENT(OUT) :: integrals(functions, functions, &
                                         &functions, functions)
        INTEGER :: m, n, l, s
        INTEGER :: zeta, eta

        DO m = 1, functions
        DO n = 1, functions
            zeta = m**2 + n**2

            DO l = 1, functions
            DO s = 1, functions
                eta = l**2 + s**2

                integrals(m,n,l,s) = overlap(m,n) * overlap(l,s) * &
                    & alpha * zeta * eta * &
                    & ((0.5d0 * dble(zeta**2 + 6 * zeta * eta + eta**2) / &
                    & dble(zeta + eta)**3) - (6.d0 * dble(zeta * eta)**2 * &
                    & dlog(dble(zeta * eta))) / dble(zeta + eta)**5)

            ENDDO
            ENDDO
        ENDDO
        ENDDO
    END SUBROUTINE quasi_llll

    SUBROUTINE quasi_poly(A, evens, odds, overlap, integrals)
        !
        ! This subroutine computes fake electron repulsion 
        ! integrals for electrons occupying the same finite 
        ! domain. These values produce the correct double bar
        ! integral when appropriately combined
        !
        ! WARNING : Precomputed integrals are scaled for the
        !           approriate domain length
        !
        ! INPUT :
        !   [real] A
        ! Half-width of the domain
        !   [int ] evens, odds
        ! Number of even and odd polynomial basis functions
        ! used in the domain
        !   [real] overlap(:,:)
        ! Overlap for this domain
        !
        ! OUTPUT :
        !   [real] integrals(:,:,:,:)
        ! Memory space for the resulting integrals
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN)   :: evens, odds
        REAL(dp), TARGET, INTENT(IN)  :: A, overlap(evens + odds, evens + odds)
        REAL(dp), TARGET, INTENT(OUT) :: integrals(evens + odds, evens + odds, &
                                          &evens + odds, evens + odds)

        INTEGER :: m, n, l, s
        REAL(dp), POINTER :: overlap_ee(:,:), overlap_oo(:,:), &
                           & eri_EEEE(:,:,:,:), eri_EEOO(:,:,:,:), &
                           & eri_OOEE(:,:,:,:), eri_OOOO(:,:,:,:), &
                           & eri_EOEO(:,:,:,:), eri_EOOE(:,:,:,:), &
                           & eri_OEEO(:,:,:,:), eri_OEOE(:,:,:,:)

        ! Associate pointers
        overlap_ee => overlap(1:evens, 1:evens)
        overlap_oo => overlap(evens+1:evens+odds, evens+1:evens+odds)

        ! WARNING : Hacky bug fix
        FORALL (m = 1:evens+odds, n = 1:evens+odds, &
               &l = 1:evens+odds, s = 1:evens+odds)
            integrals(m,n,l,s) = 0.d0
        END FORALL

        eri_EEEE => integrals(1:evens, 1:evens, &
                             &1:evens, 1:evens)
        eri_EEOO => integrals(1:evens, 1:evens, &
                             &evens+1:evens+odds, evens+1:evens+odds)
        eri_OOEE => integrals(evens+1:evens+odds, evens+1:evens+odds, &
                             &1:evens, 1:evens)
        eri_OOOO => integrals(evens+1:evens+odds, evens+1:evens+odds, &
                             &evens+1:evens+odds, evens+1:evens+odds)

        eri_EOEO => integrals(1:evens, evens+1:evens+odds, &
                             &1:evens, evens+1:evens+odds)
        eri_EOOE => integrals(1:evens, evens+1:evens+odds, &
                             &evens+1:evens+odds, 1:evens)
        eri_OEEO => integrals(evens+1:evens+odds, 1:evens, &
                             &1:evens, evens+1:evens+odds)
        eri_OEOE => integrals(evens+1:evens+odds, 1:evens, &
                             &evens+1:evens+odds, 1:evens)

        ! (EE|EE)
        DO m = 1, evens
        DO n = 1, evens
            DO l = 1, evens
            DO s = 1, evens

                eri_EEEE(m,n,l,s) = overlap_ee(m,n) * overlap_ee(l,s) * &
                                  & quasi_eeee_data(m + n, l + s) / A

            ENDDO
            ENDDO

            ! (EE|OO)
            DO l = 1, odds
            DO s = 1, odds

                eri_EEOO(m,n,l,s) = overlap_ee(m,n) * overlap_oo(l,s) * &
                    & ((3 + 2 * (l + s)) * quasi_eeee_data(m + n, l + s) - &
                    & 2 * (l + s + 1) * quasi_eeee_data(m + n, l + s + 1)) / A
                eri_OOEE(l,s,m,n) = eri_EEOO(m,n,l,s)

            ENDDO
            ENDDO
        ENDDO

        ! (EO|EO)
        DO n = 1, odds
            DO l = 1, evens
            DO s = 1, odds

                eri_EOEO(m,n,l,s) = overlap_integral(m,n,2) * &
                                  & overlap_integral(l,s,2) * &
                                  & quasi_eoeo_data(m + n, l + s) / &
                                  & (A * sqrt(dble((3 + 4 * m) * (3 + 4 * l))))
                eri_EOOE(m,n,s,l) = eri_EOEO(m,n,l,s)
                eri_OEEO(n,m,l,s) = eri_EOEO(m,n,l,s)
                eri_OEOE(n,m,s,l) = eri_EOEO(m,n,l,s)

            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! (OO|OO)
        DO m = 1, odds
        DO n = 1, odds
            DO l = 1, odds
            DO s = 1, odds

                eri_OOOO(m,n,l,s) = overlap_oo(m,n) * overlap_oo(l,s) * &
                    & ((3 + 2 * (m + n)) * (3 + 2 * (l + s)) * quasi_eeee_data(m + n, l + s) - &
                    & (3 + 2 * (m + n)) * (2 + 2 * (l + s)) * quasi_eeee_data(m + n, l + s + 1) - &
                    & (2 + 2 * (m + n)) * (3 + 2 * (l + s)) * quasi_eeee_data(m + n + 1, l + s) + &
                    & (2 + 2 * (m + n)) * (2 + 2 * (l + s)) * quasi_eeee_data(m + n + 1, l + s + 1)) / A

            ENDDO
            ENDDO
        ENDDO
        ENDDO

    END SUBROUTINE quasi_poly

END MODULE quasi_eri
