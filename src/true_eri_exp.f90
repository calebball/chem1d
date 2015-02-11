MODULE true_eri_exp
    !
    ! TRUE_ERI_EXP
    !
    ! This module contains subroutines for computing electron
    ! repulsion integrals when the two electrons inhabit 
    ! seperate domains and at least one of them occupies an
    ! exterior domain.
    !
    USE constants
    USE input, ONLY : alpha
    USE error
    USE special_functions, ONLY : exp_incomplete_gamma
    USE one_e_integrals, ONLY : overlap_integral

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: true_llrr
    PUBLIC :: true_llee

    CONTAINS

    SUBROUTINE true_llrr(R, func_l, func_r, overlap_ll, overlap_rr, &
                        &scratch, eri_llrr, eri_rrll)
        !
        ! This subroutine correctly packs integrals of 
        ! exponential basis functions into a four index 
        ! ERI table. It makes use of other subroutines
        ! to compute the necessary integrals
        ! WARNING : This routine lacks the ability to use
        !           different alpha exponents for the two
        !           domains
        !
        ! INPUT :
        !   [real] R
        ! Distance between the two domains
        !   [int ] func_l, func_r
        ! Number of basis functions in domain one and two
        !   [real] overlap_ll(:,:), overlap_rr(:,:)
        ! Overlap matrices for domain one and two
        !
        ! OUTPUT :
        !   [real] scratch(:,:)
        ! Scratch space for storing intermediate values.
        !   [real] eri_llrr(:,:,:,:), eri_rrll(:,:,:,:)
        ! Four index array for storing the resulting integrals
        !
        IMPLICIT NONE

        REAL(dp), INTENT(IN) :: R, overlap_ll(:,:), overlap_rr(:,:)
        INTEGER, INTENT(IN)  :: func_l, func_r

        REAL(dp), INTENT(OUT) :: scratch(:,:), &
                               & eri_llrr(:,:,:,:), &
                               & eri_rrll(:,:,:,:)

        REAL(dp) :: integral
        INTEGER  :: m, n, l, s, zeta, eta, max_m

        integral = 0.d0

        ! Check whether the input is an atom
        IF (R.EQ.0d0) THEN

            ! Iterate over function pairs in domain one
            DO m = 1, func_l
            DO n = 1, func_l
                zeta = m**2 + n**2

                ! Iterate over function pairs in domain two
                DO l = 1, func_r
                DO s = 1, func_r
                    eta = l**2 + s**2

                    ! Compute integral value
                    ! NOTE: This IF block unfortunately remains
                    ! The eta == zeta condition will trigger rarely
                    ! when (m,n) != (l,s)
                    IF (eta.NE.zeta) THEN
                        integral = overlap_ll(m,n) * overlap_rr(l,s) * &
                            & ((alpha * zeta * eta) / &
                            &  (2.d0 * dble(zeta - eta)**5)) * &
                            & ((zeta**2 - eta**2) * &
                            &  (zeta**2 - 8.d0 * zeta * eta + eta**2) + &
                            &  12.d0 * zeta**2 * eta**2 * log(dble(zeta) / dble(eta)))
                    ELSE
                        integral = overlap_ll(m,n) * overlap_rr(l,s) * &
                            & alpha * zeta / 5.d0
                    ENDIF

                    eri_llrr(m,n,l,s) = integral
                    eri_rrll(l,s,m,n) = integral

                    integral = 0.d0

                ENDDO
                ENDDO
            ENDDO
            ENDDO

        ELSE

            max_m = max(func_l, func_r)

            ! Start by computing incomplete gamma functions
            DO m = 1, max_m
            DO n = 1, max_m
                scratch(m,n) = dble(m**2 + n**2)**2 * &
                    & exp_incomplete_gamma(0, R * alpha * dble(m**2 + n**2))
            ENDDO
            ENDDO

            ! Iterate over function pairs in the left domain
            DO m = 1, func_l
            DO n = 1, func_l
                zeta = m**2 + n**2

                ! Iterate over function pairs in the right domain
                DO l = 1, func_r
                DO s = 1, func_r
                    eta = l**2 + s**2

                    ! Compute integral value
                    ! Again, the IF block is here since the eta == zeta 
                    ! condition can trigger rarely when (m,n) != (l,s)
                    IF (eta.NE.zeta) THEN
                        integral = overlap_ll(m,n) * overlap_rr(l,s) * &
                            & ((alpha * zeta * eta) / (2.d0 * dble(zeta - eta)**5)) * ( &
                            & (zeta - eta) * (eta**3 - &
                            & zeta**2 * eta * (7.d0 - 2.d0 * R * alpha * eta) - &
                            & zeta**3 * (R * alpha * eta - 1.d0) - &
                            & zeta * eta**2 * (7.d0 + R * alpha * eta) ) - &
                            & eta**2 * (12.d0 - R * alpha * (6.d0 - R * alpha * &
                            & (zeta - eta)) * (zeta - eta)) * scratch(m,n) + &
                            & zeta**2 * (12.d0 - R * alpha * (- 6.d0 - R * alpha * &
                            & (zeta - eta)) * (zeta - eta)) * scratch(l,s) )

                    ELSE
                        integral = overlap_ll(m,n) * overlap_rr(l,s) * &
                            & (alpha * zeta / 120.d0) * &
                            & (24.d0 - R * alpha * zeta * &
                            & (6.d0 - R * alpha * zeta * &
                            & (2.d0 - R * alpha * zeta * &
                            & (1.d0 - R * alpha * zeta ))) - &
                            & (R * alpha)**5 * zeta**3 * scratch(m,n))

                    ENDIF

                    eri_llrr(m,n,l,s) = integral
                    eri_rrll(l,s,m,n) = integral

                    integral = 0.d0

                ENDDO
                ENDDO
            ENDDO
            ENDDO

        ENDIF

    END SUBROUTINE true_llrr

    SUBROUTINE true_llee(left, A, R, e, o, max_m, &
                        &overlap_ll, overlap_poly, &
                        &eri_exp_poly, eri_poly_exp)
        !
        ! Computes the integrals between exponential type pairs
        ! and polynomial type pairs using numerical quadrature
        ! to integrate the potential of the LL pair against the
        ! EE pair (also EO and OO)
        !
        ! WARNING : This routine currently uses the trapezoid
        !           rule with a fixed number of points
        ! WARNING : This routine lacks the ability to use
        !           different alpha exponents for the two
        !           external domains
        !
        ! INPUT :
        !   [bool] left
        ! Logical value stating whether the exponential domain
        ! is to the left of the molecule
        !   [real] A
        ! Half-width of the polynomial domain
        !   [real] R
        ! Distance of the exponential nucleus from the centroid
        ! of the polynomial domain
        !   [int ] e, o
        ! Maximum order of even and odd polynomials respectively
        !   [int ] max_m
        ! Maximum order of exponential basis functions
        !   [real] overlap_ll, overlap_poly
        ! Overlap matrices for the exponential and polynomial
        ! functions respectively
        !
        ! OUTPUT :
        !   [real] eri_exp_poly, eri_poly_exp
        ! Four index arrays for storing the resulting integrals
        ! where the first two indexes refer to the infinite 
        ! domain and finite domain respectively
        ! 
        IMPLICIT NONE

        LOGICAL,  INTENT(IN)    :: left
        REAL(dp), INTENT(IN)    :: A, R, overlap_ll(:,:)
        INTEGER,  INTENT(IN)    :: e, o, max_m
        REAL(dp), TARGET, INTENT(IN)  :: overlap_poly(:,:)
        REAL(dp), TARGET, INTENT(OUT) :: eri_exp_poly(:,:,:,:), &
                                       & eri_poly_exp(:,:,:,:)

        INTEGER  :: m, n, l, s, i, eo_sign
        INTEGER, PARAMETER :: quad_points = 500
        REAL(dp) :: rma, rpa, zeta, x, h, u(quad_points)

        REAL(dp), POINTER :: overlap_ee(:,:), overlap_oo(:,:),     &
                           & eri_LLEE(:,:,:,:), eri_EELL(:,:,:,:), &
                           & eri_LLEO(:,:,:,:), eri_EOLL(:,:,:,:), &
                           & eri_LLOE(:,:,:,:), eri_OELL(:,:,:,:), &
                           & eri_LLOO(:,:,:,:), eri_OOLL(:,:,:,:)

        ! Associate pointers
        overlap_ee => overlap_poly(1:e, 1:e)
        overlap_oo => overlap_poly(e+1:e+o, e+1:e+o)

        eri_LLEE => eri_exp_poly(1:max_m, 1:max_m, 1:e, 1:e)
        eri_LLEO => eri_exp_poly(1:max_m, 1:max_m, 1:e, e+1:e+o)
        eri_LLOE => eri_exp_poly(1:max_m, 1:max_m, e+1:e+o, 1:e)
        eri_LLOO => eri_exp_poly(1:max_m, 1:max_m, e+1:e+o, e+1:e+o)

        eri_EELL => eri_poly_exp(1:e, 1:e, 1:max_m, 1:max_m)
        eri_EOLL => eri_poly_exp(1:e, e+1:e+o, 1:max_m, 1:max_m)
        eri_OELL => eri_poly_exp(e+1:e+o, 1:e, 1:max_m, 1:max_m)
        eri_OOLL => eri_poly_exp(e+1:e+o, e+1:e+o, 1:max_m, 1:max_m)

        ! Useful values
        rma = R - A
        rpa = R + A
        h = 2.d0 * A / dble(quad_points + 1)
        IF (left) THEN
            eo_sign = 1
        ELSE 
            eo_sign = -1
        ENDIF

        ! Loop over LL shell pairs
        DO m = 1, max_m
        DO n = 1, max_m

            zeta = alpha * (m**2 + n**2)

            ! Compute hypergeometric functions at quadrature points
            DO i = 1, quad_points
                x = (rma + i * h)
                u(i) = x**2 * &
                     & exp_incomplete_gamma(-2, zeta * x)
            ENDDO

            ! Loop over EE shell pairs
            DO l = 1, e
            DO s = 1, e

                ! Sum internal integrand points
                eri_LLEE(m,n,l,s) = 0.d0
                DO i = 1, quad_points
                    x = (rma + i * h)
                    eri_LLEE(m,n,l,s) = eri_LLEE(m,n,l,s) + &
                        & ((x - rma) * (rpa - x))**(l + s) * u(i)
                ENDDO
                eri_LLEE(m,n,l,s) = eri_LLEE(m,n,l,s) * 2.d0 * A / dble(quad_points + 1)

                ! Including the prefactor on the sum
                eri_LLEE(m,n,l,s) = eri_LLEE(m,n,l,s) * zeta**3 / &
                    & (A**(2 * (l + s) + 1) * ROOT_PI)
                eri_LLEE(m,n,l,s) = eri_LLEE(m,n,l,s) * ROOT_PI / 2.d0
                DO i = 1, l + s
                    eri_LLEE(m,n,l,s) = eri_LLEE(m,n,l,s) * dble(1 + 2 * i) / dble(2 * i)
                ENDDO
                eri_LLEE(m,n,l,s) = overlap_ll(m,n) * overlap_ee(l,s) * eri_LLEE(m,n,l,s)

                ! Symmetry
                eri_EELL(l,s,m,n) = eri_LLEE(m,n,l,s)

            ENDDO
            ENDDO

            ! Loop over EO shell pairs
            DO l = 1, e
            DO s = 1, o

                ! Sum internal integrand points
                eri_LLEO(m,n,l,s) = 0.d0
                DO i = 1, quad_points
                    x = (rma + i * h)
                    eri_LLEO(m,n,l,s) = eri_LLEO(m,n,l,s) + &
                        & ((x - rma) * (rpa - x))**(l + s) * (R - x) * u(i)
                ENDDO
                eri_LLEO(m,n,l,s) = eri_LLEO(m,n,l,s) * 2.d0 * A / dble(quad_points + 1)

                ! Including the prefactor on the sum
                eri_LLEO(m,n,l,s) = - eri_LLEO(m,n,l,s) * 2.d0 * zeta**3 / &
                    & (A**(2 * (l + s) + 2) * ROOT_PI * sqrt(dble(3 + 4 * l)))
                eri_LLEO(m,n,l,s) = eri_LLEO(m,n,l,s) * 3.d0 * ROOT_PI / 4.d0
                DO i = 1, l + s
                    eri_LLEO(m,n,l,s) = eri_LLEO(m,n,l,s) * dble(3 + 2 * i) / dble(2 * i)
                ENDDO
                eri_LLEO(m,n,l,s) = eo_sign * overlap_ll(m,n) * &
                                  & overlap_integral(l,s,2) * eri_LLEO(m,n,l,s)

                ! Symmetry
                eri_LLOE(m,n,s,l) = eri_LLEO(m,n,l,s)
                eri_EOLL(l,s,m,n) = eri_LLEO(m,n,l,s)
                eri_OELL(s,l,m,n) = eri_LLEO(m,n,l,s)

            ENDDO
            ENDDO

            ! Loop over OO shell pairs
            DO l = 1, o
            DO s = 1, o

                ! Sum internal integrand points
                eri_LLOO(m,n,l,s) = 0.d0
                DO i = 1, quad_points
                    x = (rma + i * h)
                    eri_LLOO(m,n,l,s) = eri_LLOO(m,n,l,s) + &
                        & ((x - rma) * (rpa - x))**(l + s) * (R - x)**2 * u(i)
                ENDDO
                eri_LLOO(m,n,l,s) = eri_LLOO(m,n,l,s) * 2.d0 * A / dble(quad_points + 1)

                ! Including the prefactor on the sum
                eri_LLOO(m,n,l,s) = eri_LLOO(m,n,l,s) * 2.d0 * zeta**3 / &
                    & (A**(2 * (l + s) + 3) * ROOT_PI)
                eri_LLOO(m,n,l,s) = eri_LLOO(m,n,l,s) * 3.d0 * ROOT_PI / 4.d0
                DO i = 1, l + s
                    eri_LLOO(m,n,l,s) = eri_LLOO(m,n,l,s) * dble(3 + 2 * i) / dble(2 * i)
                ENDDO
                eri_LLOO(m,n,l,s) = overlap_ll(m,n) * overlap_oo(l,s) * eri_LLOO(m,n,l,s)

                ! Symmetry
                eri_OOLL(l,s,m,n) = eri_LLOO(m,n,l,s)

            ENDDO
            ENDDO

        ENDDO
        ENDDO

    END SUBROUTINE true_llee

END MODULE true_eri_exp
