MODULE one_e_integrals
    !
    ! ONE_E_INTEGRALS
    !
    ! The one_e_inetgrals module contains subroutines for generating
    ! integrals involving only one electron, i.e. overlap, kinetic 
    ! energy and potential energy integrals
    !
    USE constants
    USE input
    USE error
    USE special_functions, ONLY : incomplete_gamma, exp_incomplete_gamma, two_f_one

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: overlap_int_in_domain
    PUBLIC :: overlap_integral
    PUBLIC :: kinetic_int_in_domain
    PUBLIC :: potential_int_in_domain

    CONTAINS

    SUBROUTINE overlap_int_in_domain(domain, integrals)
        !
        ! Computes the overlap integrals for all pairs
        ! of basis functions in a given domain
        !
        ! INPUT :
        !   [int ] domain
        ! Index of the domain
        !
        ! OUTPUT :
        !   [real] integrals(:,:)
        ! Matrix for storing the resulting integrals
        !
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: domain
        REAL(dp), INTENT(OUT) :: integrals(:,:)

        INTEGER :: i, j
        INTEGER :: n_evens, n_odds
        REAL(dp) :: coeff
        CHARACTER(LEN=58) :: error_message

        IF (domain.EQ.1.OR.domain.EQ.n_domains) THEN

            ! We compute these integrals from an easy to evaluate
            ! analytic expression
            FORALL(i = 1:functions_in_domain(domain), &
                  &j = 1:functions_in_domain(domain))
                integrals(i,j) = (dble(2 * i * j) / dble(i**2 + j**2))**3
            END FORALL

        ELSE IF (domain.GT.1.AND.domain.LT.n_domains) THEN
            ! NOTE : A better way to do this would be to build an
            !        identity matrix and then use recur from the
            !        diagonal to build the rest

            n_evens = evens_in_domain(domain)
            n_odds  = odds_in_domain(domain)

            ! The analytic expressions for these integrals are in
            ! terms of ratios of gamma functions.
            ! We're going to compute them recursively.
            IF (functions_in_domain(domain).GT.0) integrals(1,1) = 1.d0
            IF (n_odds.GT.0) integrals(n_evens+1,n_evens+1) = 1.d0

            DO i = 2, n_evens
                coeff = dble(i + 1) / dble(6 + 4*i)
                coeff = coeff * dsqrt( dble(2 - 32 * i**2) / &
                      &                dble( i - 2 * i**2) )
                integrals(i,1) = coeff * integrals(i-1,1)
            ENDDO

            DO i = 2, n_odds
                coeff = dble(2 * i + 2) / dble(5 + 2 * i)
                coeff = coeff * dsqrt( dble(3 + 16 * (i + i**2)) / &
                      &                dble(16 * i**2 - 8 * i  ) )
                integrals(i+n_evens,n_evens+1) = &
                                & coeff * integrals(i-1+n_evens,n_evens+1)
            ENDDO

            DO i = 1, n_evens
                DO j = 2, n_evens
                    coeff = dble(j + i) / dble(2 + 4 * (i + j))
                    coeff = coeff * dsqrt( dble(2 - 32 * j**2) / &
                          &                dble( j - 2 * j**2) )
                    integrals(i,j) = coeff * integrals(i,j-1)
                ENDDO
            ENDDO

            DO i = 1, n_odds
                DO j = 2, n_odds
                    coeff = dble(2 * (i + j)) / dble(3 + 2 * (i + j))
                    coeff = coeff * dsqrt( dble(3 + 16 * (j + j**2)) / &
                          &                dble(16 * j**2 - 8 * j  ) )
                    integrals(i+n_evens,j+n_evens) = &
                                    & coeff * integrals(i+n_evens,j-1+n_evens)
                ENDDO
            ENDDO

        ELSE
            WRITE (error_message,'(a)') &
                & "Supplied domain index does not exist"
            CALL print_warning("overlap_int_in_domain", error_message)
        ENDIF
    END SUBROUTINE overlap_int_in_domain

    REAL(dp) PURE FUNCTION overlap_integral(mu, nu, id)
        !
        ! Computes the specific overlap integral between
        ! the basis functions mu and nu
        !
        ! INPUT :
        !   [int ] mu, nu
        ! Indices of the two basis functions
        !   [int ] id
        ! Type of basis functions
        !       0 -> Exponentials
        !       1 -> Even polynomials
        !       2 -> Odd polynomials
        !
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: mu, nu, id
        INTEGER  :: m
        REAL(dp) :: mu_p_nu, only_mu, only_nu

        SELECT CASE (id)
        CASE (0)
            overlap_integral = (2.d0 * mu * nu / (mu**2 + nu**2))**3
        CASE (1)
            mu_p_nu = 16.d0 / 15.d0
            DO m = 3, mu + nu
                mu_p_nu = mu_p_nu * dble(2 * m) / dble(1 + 2 * m)
            ENDDO

            only_mu = 15.d0 / 16.d0
            DO m = 2, mu
                only_mu = only_mu * dble(16 * m**2 - 1) / &
                    & dble(8 * m * (2 * m - 1))
            ENDDO

            only_nu = 15.d0 / 16.d0
            DO m = 2, nu
                only_nu = only_nu * dble(16 * m**2 - 1) / &
                    & dble(8 * m * (2 * m - 1))
            ENDDO

            overlap_integral = mu_p_nu * sqrt(only_mu * only_nu)
        CASE (2)
            mu_p_nu = 1.d0
            DO m = 3, mu + nu
                mu_p_nu = mu_p_nu * dble(2 * m) / dble(3 + 2 * m)
            ENDDO

            only_mu = 1.d0
            DO m = 2, mu
                only_mu = only_mu * dble(3 + 16 * m * (m + 1)) / &
                    & dble(8 * m * (2 * m - 1))
            ENDDO

            only_nu = 1.d0
            DO m = 2, nu
                only_nu = only_nu * dble(3 + 16 * m * (m + 1)) / &
                    & dble(8 * m * (2 * m - 1))
            ENDDO

            overlap_integral = mu_p_nu * sqrt(only_mu * only_nu)
        CASE DEFAULT

            ! NOTE : It would be nice to catch an error at this point.
            !        Using an error routine is not pure, so we need
            !        something better.
            overlap_integral = 0.d0

        END SELECT

    END FUNCTION overlap_integral

    SUBROUTINE kinetic_int_in_domain(domain, integrals)
        !
        ! Computes the overlap integrals for all pairs
        ! of basis functions in a given domain
        !
        ! INPUT :
        !   [int ] domain
        ! Index of the domain
        !
        ! OUTPUT :
        !   [real] integrals(:,:)
        ! Matrix for storing the resulting integrals
        !
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: domain
        REAL(dp), INTENT(OUT) :: integrals(:,:)

        INTEGER :: i, j
        INTEGER :: n_exps, n_evens, n_odds
        REAL(dp) :: A, B
        CHARACTER(LEN=58) :: error_message

        CALL overlap_int_in_domain(domain, integrals)

        IF (domain.EQ.1.OR.domain.EQ.n_domains) THEN

            n_exps = functions_in_domain(domain)

            FORALL(i = 1:n_exps, j = 1:n_exps)
                integrals(i,j) = integrals(i,j) * &
                               & dble(alpha * i * j)**2 / 2.d0
            END FORALL

        ELSE IF (domain.GT.1.AND.domain.LT.n_domains) THEN

            n_evens = evens_in_domain(domain)
            n_odds  = odds_in_domain(domain)

            A = nuclear_position(domain - 1)
            B = nuclear_position(domain)

            FORALL (i = 1:n_evens, j = 1:n_evens)
                integrals(i,j) = integrals(i,j) * &
                               & dble(4.d0 * i * j * (i + j + 0.5d0)) / &
                               & dble((i + j)*(i + j - 1)*(A - B)**2)
            END FORALL

            FORALL (i = 1:n_odds, j = 1:n_odds)
                integrals(i+n_evens,j+n_evens) = &
                               & integrals(i+n_evens,j+n_evens) * &
                               & dble(12.d0 * i * j * (i + j + 1.5d0)) / &
                               & dble((i + j)*(i + j - 1)*(A - B)**2)
            END FORALL

        ELSE

            WRITE (error_message,'(a)') &
                & "Supplied domain index does not exist"
            CALL print_warning("kinetic_int_in_domain", error_message)

        ENDIF

    END SUBROUTINE kinetic_int_in_domain

    SUBROUTINE potential_int_in_domain(domain, overlap, integrals)
        ! WARNING! Requires a minimum of 2 bfs in each finite domain
        ! WARNING! Needs error handling
        !
        ! Computes the overlap integrals for all pairs
        ! of basis functions in a given domain
        !
        ! INPUT :
        !   [int ] domain
        ! Index of the domain
        !   [real] overlap
        ! Matrix containing the overlap integrals in the
        ! given domain
        !
        ! OUTPUT :
        !   [real] integrals(:,:)
        ! Matrix for storing the resulting integrals
        !
        IMPLICIT NONE
        INTEGER,  INTENT(IN)  :: domain
        REAL(dp), INTENT(IN)  :: overlap(:,:)
        REAL(dp), INTENT(OUT) :: integrals(:,:)

        INTEGER  :: i, j, k, X, Z
        INTEGER  :: n_exps, n_evens, n_odds
        INTEGER  :: dom_nuc, other_nuc, s
        REAL(dp) :: A, R, arg, F0, F1, F2, c1, c2
        CHARACTER(LEN=58) :: error_message

        IF (domain.EQ.1.OR.domain.EQ.n_domains) THEN

            n_exps = functions_in_domain(domain)

            ! Check if the system is an atom
            IF (n_nuclei.EQ.1) THEN
                Z = nuclear_charge(1)
                FORALL (i=1:n_exps, j=1:n_exps)
                    integrals(i,j) = overlap(i,j) * Z * &
                        & dble(i**2 + j**2) * alpha / 2.d0
                END FORALL

            ELSE

                IF (domain.EQ.1) THEN
                    dom_nuc = 1
                    other_nuc = n_nuclei
                ELSE
                    dom_nuc = n_nuclei
                    other_nuc = 1
                ENDIF

                ! Find the interaction with the outermost nuclei
                Z = nuclear_charge(dom_nuc)
                FORALL (i=1:n_exps, j=1:n_exps)
                    integrals(i,j) = Z * alpha * dble(i**2 + j**2) / 2.d0
                END FORALL

                Z = nuclear_charge(other_nuc)
                R = nuclear_position(n_nuclei) - nuclear_position(1) 
                DO i = 1, n_exps
                DO j = 1, n_exps
                    integrals(i,j) = integrals(i,j) + &
                        & Z * R**2 * (alpha * dble(i**2 + j**2))**3 * &
                        & exp_incomplete_gamma(-2, R * alpha * dble(i**2 + j**2))
                ENDDO
                ENDDO

                ! Find the interaction with the interior nuclei
                DO X = 2, n_nuclei - 1
                    Z = nuclear_charge(X)
                    R = abs(nuclear_position(X) - nuclear_position(dom_nuc))

                    DO i = 1, n_exps
                    DO j = 1, n_exps
                        integrals(i,j) = integrals(i,j) + &
                            & Z * R**2 * (alpha * dble(i**2 + j**2))**3 * &
                            & exp_incomplete_gamma(-2, R * alpha * dble(i**2 + j**2))
                    ENDDO
                    ENDDO
                END DO

                ! Multiply by the overlap
                FORALL (i=1:n_exps, j=1:n_exps)
                    integrals(i,j) = overlap(i,j) * integrals(i,j)
                END FORALL

            ENDIF

        ELSE IF (domain.GT.1.AND.domain.LT.n_domains) THEN

            n_evens = evens_in_domain(domain)
            n_odds  = odds_in_domain(domain)

            A = (nuclear_position(domain) - &
              &  nuclear_position(domain - 1)) / 2.d0

            ! Nuclei adjacent to the domain
            Z = nuclear_charge(domain - 1) + nuclear_charge(domain)
            ! Evens
            FORALL (i=1:n_evens, j=1:n_evens)
                integrals(i,j) = Z * (1.d0 + 1.d0 / dble(2 * (i + j))) / A
            END FORALL

            ! Odds
            FORALL (i=1:n_odds, j=1:n_odds)
                integrals(i + n_evens,j + n_evens) = &
                                & Z * (1.d0 + 3.d0 / dble(2 * (i + j))) / A
            END FORALL

            ! Even/Odds
            Z = - nuclear_charge(domain - 1) + nuclear_charge(domain)
            FORALL (i=1:n_evens, j=1:n_odds)
                integrals(i, j + n_evens) = & 
                                & Z * (1.d0 + 3.d0 / dble(2 * (i + j))) / &
                                & dble(A * sqrt(4.d0 * i + 3.d0))
                integrals(j + n_evens, i) = & 
                                & Z * (1.d0 + 3.d0 / dble(2 * (i + j))) / &
                                & dble(A * sqrt(4.d0 * i + 3.d0))
            END FORALL

            ! Other nuclei
            DO X = 1, n_nuclei
                IF (X.EQ.(domain - 1).OR.X.EQ.domain) CYCLE
                IF (X.GT.domain) THEN
                    R = nuclear_position(X) - nuclear_position(domain) + A
                    s = 1
                ELSE
                    R = nuclear_position(domain-1) - nuclear_position(X) + A
                    s = -1
                ENDIF
                Z = nuclear_charge(X)
                arg = dble(A)**2 / dble(R)**2

                ! EVENS
                ! Bottom right corner
                F2 = two_f_one(1.d0, 0.5d0, 2 * n_evens + 1.5d0, &
                   & arg) / R
                integrals(n_evens, n_evens) = integrals(n_evens, n_evens) + Z * F2

                ! Next to the bottom right corner
                F1 = two_f_one(1.d0, 0.5d0, 2 * n_evens + 0.5d0, &
                   & arg) / R
                integrals(n_evens - 1, n_evens) = &
                                    & integrals(n_evens - 1, n_evens) + Z * F1
                integrals(n_evens, n_evens - 1) = &
                                    & integrals(n_evens, n_evens - 1) + Z * F1

                ! Recur backwards to find the rest
                DO k = 2 * n_evens - 2, 2, -1
                    c1 = dble(3 + 2 * k * (1.d0 - 2 * arg) - 5 * arg) / &
                       & dble((2 * k + 3) * (1.d0 - arg))
                    c2 = dble(2.d0 * (k + 2.d0) * arg) / &
                       & dble((2 * k + 5) * (1.d0 - arg))
                    F0 = c1 * F1 + c2 * F2

                    FORALL (i=1:n_evens, j=1:n_evens, (i + j).EQ.k)
                        integrals(i,j) = integrals(i,j) + Z * F0
                    END FORALL

                    F2 = F1
                    F1 = F0
                ENDDO

                ! ODDS
                ! Bottom right corner
                F2 = two_f_one(1.d0, 1.5d0, 2 * n_odds + 2.5d0, &
                   & arg) / R
                integrals(n_evens+n_odds,n_evens+n_odds) = &
                            & integrals(n_evens+n_odds,n_evens+n_odds) + Z * F2

                ! Next to the bottom right corner
                F1 = two_f_one(1.d0, 1.5d0, 2 * n_odds + 1.5d0, &
                   & arg) / R
                integrals(n_evens+n_odds - 1,n_evens+n_odds) = &
                        & integrals(n_evens+n_odds - 1,n_evens+n_odds) + Z * F1
                integrals(n_evens+n_odds,n_evens+n_odds - 1) = &
                        & integrals(n_evens+n_odds,n_evens+n_odds - 1) + Z * F1

                ! Recur backwards to find the rest
                DO k = 2 * n_odds - 2, 2, -1
                    c1 = dble(5 + 2 * k * (1.d0 - 2 * arg) - 7 * arg) / &
                       & dble((2 * k + 5) * (1.d0 - arg))
                    c2 = dble(2.d0 * (k + 2.d0) * arg) / &
                       & dble((2 * k + 7) * (1.d0 - arg))
                    F0 = c1 * F1 + c2 * F2

                    FORALL (i=1:n_odds, j=1:n_odds, (i + j).EQ.k)
                        integrals(n_evens+i,n_evens+j) = &
                                    & integrals(n_evens+i,n_evens+j) + Z * F0
                    END FORALL

                    F2 = F1
                    F1 = F0
                ENDDO

                ! EVEN/ODDS
                ! Bottom right corner
                F2 = two_f_one(1.d0, 1.5d0, n_evens + n_odds + 2.5d0, &
                   & arg) / R
                integrals(n_evens,n_evens+n_odds) = &
                        & integrals(n_evens,n_evens+n_odds) + &
                        & s * Z * F2 * A / dble(R * sqrt(4.d0 * n_evens + 3.d0))
                integrals(n_evens+n_odds,n_evens) = &
                        & integrals(n_evens+n_odds,n_evens) + &
                        & s * Z * F2 * A / dble(R * sqrt(4.d0 * n_evens + 3.d0))

                ! Next to the bottom right corner
                F1 = two_f_one(1.d0, 1.5d0, n_evens + n_odds + 1.5d0, &
                   & arg) / R
                integrals(n_evens - 1,n_evens+n_odds) = &
                        & integrals(n_evens - 1,n_evens+n_odds) + &
                        & s * Z * F1 * A / dble(R * sqrt(4.d0 * n_evens - 1.d0))
                integrals(n_evens,n_evens+n_odds - 1) = &
                        & integrals(n_evens,n_evens+n_odds - 1) + &
                        & s * Z * F1 * A / dble(R * sqrt(4.d0 * n_evens + 3.d0))
                integrals(n_evens+n_odds - 1,n_evens) = &
                        & integrals(n_evens+n_odds - 1,n_evens) + &
                        & s * Z * F1 * A / dble(R * sqrt(4.d0 * n_evens + 3.d0))
                integrals(n_evens+n_odds,n_evens - 1) = &
                        & integrals(n_evens+n_odds,n_evens - 1) + &
                        & s * Z * F1 * A / dble(R * sqrt(4.d0 * n_evens - 1.d0))

                ! Recur backwards to find the rest
                DO k = n_evens + n_odds - 2, 2, -1
                    c1 = dble(5 + 2 * k * (1.d0 - 2 * arg) - 7 * arg) / &
                       & dble((2 * k + 5) * (1.d0 - arg))
                    c2 = dble(2.d0 * (k + 2.d0) * arg) / &
                       & dble((2 * k + 7) * (1.d0 - arg))
                    F0 = c1 * F1 + c2 * F2

                    FORALL (i=1:n_evens, j=1:n_odds, (i + j).EQ.k)
                        integrals(i,n_evens+j) = integrals(i,n_evens+j) + &
                            & s * Z * F0 * A / dble(R * sqrt(4.d0 * i + 3.d0))
                        integrals(n_evens+j,i) = integrals(n_evens+j,i) + &
                            & s * Z * F0 * A / dble(R * sqrt(4.d0 * i + 3.d0))
                    END FORALL

                    F2 = F1
                    F1 = F0
                ENDDO
            ENDDO

            ! Multiply by the overlap
            FORALL (i=1:n_evens, j=1:n_evens)
                integrals(i,j) = integrals(i,j) * overlap(i,j)
            END FORALL

            FORALL (i=1:n_odds, j=1:n_odds)
                integrals(i + n_evens, j + n_evens) = &
                                    & integrals(i + n_evens, j + n_evens) * &
                                    &  overlap(i + n_evens, j + n_evens)
            END FORALL

            ! NOTE : There must be a way to use the more
            ! intensive overlap calculation less...

            FORALL (i=1:n_evens, j=1:n_odds)
                integrals(i, j + n_evens) = integrals(i, j + n_evens) * &
                                        & overlap_integral(i, j, 2)
                integrals(j + n_evens, i) = integrals(j + n_evens, i) * &
                                        & overlap_integral(i, j, 2)
            END FORALL

        ELSE
            WRITE (error_message,'(a)') &
                & "Supplied domain index does not exist"
            CALL print_warning("potential_int_in_domain", error_message)
        END IF
    END SUBROUTINE potential_int_in_domain

END MODULE one_e_integrals
