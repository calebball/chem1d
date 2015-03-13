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
    USE special_functions, ONLY : incomplete_gamma, &
    & exp_incomplete_gamma, two_f_one, cin, si, log_im
    USE fgsl

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

        !Test Variables
        REAL(dp) :: input
        REAL(dp) :: result


        !Initializing the matrix to zero
        FORALL (i = 1 : functions_in_domain(domain), &
                & j = 1 : functions_in_domain(domain))

            integrals(i,j) = 0

        END FORALL

        !Filling the identity matrix
        FORALL (i = 1 : functions_in_domain(domain))

            integrals(i,i) = 1

        END FORALL

        !TESTING SECTION
        input = 0.d0
        

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

    REAL(dp) PURE FUNCTION outer_kinetic(i, j)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: i , j

        IF (i == j) THEN
            outer_kinetic = alpha ** 2 * ((dble(i * (i + 1) * &
                            & (2 * i + 1)) / &
                            & (dble(3 * sqrt(dble(i * (i + 1) * &
                            & i * (i + 1)))))) - 0.5d0)
        ELSE IF (i > j) THEN
            outer_kinetic = alpha ** 2 * ((dble(j * (j + 1) * &
                            & (2 * j + 1)) / &
                            & (dble(3 * sqrt(dble(i * (i + 1) * &
                            & j * (j + 1)))))))
        ELSE
            outer_kinetic = alpha ** 2 * ((dble(i * (i + 1) * &
                            & (2 * i + 1)) / &
                            & (dble(3 * sqrt(dble(i * (i + 1) * &
                            & j * (j + 1)))))))
        END IF

    END FUNCTION outer_kinetic

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

        !Outside left or right domains
        IF (domain.EQ.1.OR.domain.EQ.n_domains) THEN

            DO i = 1, functions_in_domain(domain)
                DO j = 1, functions_in_domain(domain)
                    integrals(i,j) = outer_kinetic(i,j)
                END DO
            END DO

            !WRITE (*,*) "Alpha = " , alpha
            !DO i = 1, functions_in_domain(domain)
            !DO j = 1, functions_in_domain(domain)
            !    WRITE (*,"(F10.6)",advance='no') integrals(i,j)
            !ENDDO
            !WRITE (*,*)
            !ENDDO

        !Inner domains
        ELSE IF (domain.GT.1.AND.domain.LT.n_domains) THEN

            n_evens = evens_in_domain(domain)
            n_odds  = odds_in_domain(domain)

            A = nuclear_position(domain - 1)
            B = nuclear_position(domain)

            FORALL (i = 1 : functions_in_domain(domain))

                integrals(i,i) = dble(0.5d0 * &
                                 & dble(i ** 2 * PI) / &
                                 & dble( (A-B) **2 ) )

            END FORALL

        ELSE

            WRITE (error_message,'(a)') &
                & "Supplied domain index does not exist"
            CALL print_warning("kinetic_int_in_domain", error_message)

        ENDIF

    END SUBROUTINE kinetic_int_in_domain

    REAL(dp) FUNCTION middle_potential(i, j, A, B, C)
    !Computes the one electron potential in middle domains
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: i , j
        REAL(dp), INTENT(IN) :: A , B , C
        REAL(dp) :: S1, s2, T1, t2

        S1 = ((A - C) / (B - A)) * (i + j) * PI
        s2 = ((A - C) / (B - A)) * (i - j) * PI
        T1 = ((B - C) / (B - A)) * (i + j) * PI
        t2 = ((B - C) / (B - A)) * (i - j) * PI

        !Debug, "Potential (1,2,2,3,-1) = "
        !Working
        middle_potential = ((COS(S1) * (cin(T1) - cin(S1)) - SIN(S1) &
        & * (si(T1) - si(S1))+(1d0 - COS(S1)) * (log_im(T1) - log_im(S1))) / &
        & (B - A)) - (((COS(s2) * (cin(t2) - cin(s2)) - SIN(s2) &
        & * (si(t2) - si(s2))+(1d0 - COS(s2)) * (log_im(t2) - log_im(s2))) / &
        & (B - A)))


    END FUNCTION middle_potential

    REAL(dp) FUNCTION exterior_potential_other(i, j, R)
    !Computes the VLL expression
    !NEEDS TO BE UPDATED
    !If clauses and evaluation of VLL expression in both cases is very inefficient

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: i , j
        REAL(dp), INTENT(IN) :: R

        INTEGER :: k
        REAL(dp) :: factor

        factor = ((2 * alpha * fgsl_sf_fact((i + 1)) * &
                  & fgsl_sf_fact((j + 1))) / sqrt( dble(i*(i + 1d0)*j*(j + 1d0))))

        !Top of Matrix
        IF (i >= j) THEN

            !Sum
            DO k = 1, j

                exterior_potential_other = exterior_potential_other & 
                                           & + (factor * & 
                                           & ((fgsl_sf_fact(i - k + j - k) * R ** (2 * k) * &
                                           & fgsl_sf_hyperg_U_int(i + j + 1, 2 * k + 1, R)) / &
                                           & (fgsl_sf_fact(i - k) * fgsl_sf_fact(j - k) * &
                                           & fgsl_sf_fact(k - 1) * fgsl_sf_fact(k + 1))))

            END DO
                                
        !Bottom Matrix
        ELSE

            !Sum
            DO k = 1, i
                    
                exterior_potential_other = exterior_potential_other & 
                                           & + (factor * & 
                                           & ((fgsl_sf_fact(i - k + j - k) * R ** (2 * k) * &
                                           & fgsl_sf_hyperg_U_int(i + j + 1, 2 * k + 1, R)) / &
                                           & (fgsl_sf_fact(i - k) * fgsl_sf_fact(j - k) * &
                                           & fgsl_sf_fact(k - 1) * fgsl_sf_fact(k + 1))))

            END DO

        END IF

    END FUNCTION exterior_potential_other

    SUBROUTINE potential_int_in_domain(domain, integrals)
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
        REAL(dp), INTENT(OUT) :: integrals(:,:)

        INTEGER  :: i, j, k, X, Z
        INTEGER  :: dom_nuc, other_nuc, left_nuc, right_nuc
        REAL(dp) :: A, B, R
        CHARACTER(LEN=58) :: error_message

        !WRITE (*,*) "Middle Potential (1,2,2,3,-1) = " , middle_potential(1,2,2d0,3d0,-1d0)
        !WRITE (*,*) "Exterior Potential Other (4,6,-2,alpha,1) = " , &
        !            & exterior_potential_other(4,6,2d0 * alpha * (1d0 + 2d0))

        !Do I have to initializa the integral matrix?

        !Exterior Domains
        IF (domain.EQ.1.OR.domain.EQ.n_domains) THEN

            !Adjacent Nuclei \[Zeta] = 0 Case

                !Distinguis between Left or Right Exeterior Case
                
                !Left Case
                IF (domain.EQ.1) THEN
                    dom_nuc = 1
                    other_nuc = n_nuclei
                !Right Case
                ELSE
                    dom_nuc = n_nuclei
                    other_nuc = 1
                ENDIF

                    WRITE(*,*) "Alpha = " , alpha
                    WRITE(*,*) "Functions = " , functions_in_domain
                    WRITE(*,*) "Other Nuc = " , other_nuc
                    WRITE(*,*) "Dom Nuc = " , dom_nuc

                !Calculate potential of adjacent nucleus
                Z = nuclear_charge(dom_nuc)
               
                !Top of Matrix
                FORALL (i=1:functions_in_domain(domain), &
                    & j=1:functions_in_domain(domain), i >= j)
                    
                    integrals(i,j) = integrals(i,j) + (alpha * j * (j + 1d0)) / &
                    & sqrt(i * (i + 1d0) * j * (j + 1d0))

                END FORALL

                !Bottom of Matrix
                FORALL (i=1:functions_in_domain(domain), &
                    & j=1:functions_in_domain(domain), j > i)
                   
                    integrals(i,j) = integrals(i,j) + (alpha * i * (i + 1d0)) / &
                    & sqrt(i * (i + 1d0) * j * (j + 1d0))

                END FORALL

                WRITE (*,*) "PASSED ADJACENT"

                !Calculate potential of all other nuclei
                Z = nuclear_charge(other_nuc)
                !Differentiate different sign for left and right case
                !Left Case
                IF (domain.EQ.1) THEN

                    WRITE(*,*) "Left Exterior"
                    WRITE(*,*) "other_nuc = " , other_nuc
                    WRITE(*,*) "dom_nuc = " , dom_nuc
                    
                    !Iterate over nuclei from the right
                    !WARNING NOT CERRTAIN ABOUT DO FALL THROUGH
                    DO k = other_nuc, dom_nuc + 1, -1 

                        WRITE(*,*) "Left Exterior Iterator"

                        !Set new nuclei position argument
                        R = 2d0 * alpha * (nuclear_position(other_nuc) - nuclear_position(dom_nuc))
                        
                        DO i = 1 , functions_in_domain(domain)
                            DO j = 1 , functions_in_domain(domain)
                                integrals(i,j) = integrals(i,j) + exterior_potential_other(i,j,R)
                            END DO
                        END DO

                    END DO

                !Right Case
                ELSE
                    
                    WRITE(*,*) "Right Exterior"

                    !Iterate over nuclei from the left
                    !WARNING NOT CERRTAIN ABOUT DO FALL THROUGH
                    DO k = other_nuc, dom_nuc - 1, 1 

                        WRITE(*,*) "Right Exterior Iterator"

                        !Set new nuclei position argument
                        R = 2d0 * alpha * (nuclear_position(dom_nuc) - nuclear_position(other_nuc))

                        DO i = 1 , functions_in_domain(domain)
                            DO j = 1 , functions_in_domain(domain)
                                integrals(i,j) = integrals(i,j) + exterior_potential_other(i,j,R)
                            END DO
                        END DO

                    END DO

                END IF



        !Interior Domains
        ELSE IF (domain.GT.1.AND.domain.LT.n_domains) THEN
            
            !Number of Nuclei to Left/Right of specified domain NOT including adjacent
            left_nuc = domain - 2
            right_nuc = n_domains - domain - 1

            !Adjacent Case
            !Set Adjacent Nuclei Positions
            A = nuclear_position(domain - 1)
            B = nuclear_position(domain)

            DO i = 1, functions_in_domain(domain)
                DO j = 1, functions_in_domain(domain)

                
                    integrals(i,j) = (cin(PI * (i + j)) - cin(PI * (i - j))) &
                                   & / (B - A)
                
                END DO
            END DO

            !Leftward Case
            DO k = left_nuc, 1, -1
                WRITE(*,*) "Left Interior"
                integrals(i,j) = integrals(i,j) + &
                               & middle_potential(i,j,A,B,nuclear_position(k))
            END DO

            !Rightward Case
            DO k = right_nuc, n_nuclei, 1
                WRITE(*,*) "Right Interior"
                integrals(i,j) = integrals(i,j) + &
                               & middle_potential(i,j,A,B,nuclear_position(k))
            END DO

        ELSE
            WRITE (error_message,'(a)') &
                & "Supplied domain index does not exist"
            CALL print_warning("potential_int_in_domain", error_message)
        END IF
    END SUBROUTINE potential_int_in_domain

END MODULE one_e_integrals
