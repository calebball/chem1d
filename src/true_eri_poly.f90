MODULE true_eri_poly
    !
    ! TRUE_ERI_POLY
    !
    ! This module contains subroutines for computing electron
    ! repulsion integrals when the two electrons inhabit 
    ! seperate finite domains
    !
    USE constants
    USE error
    USE one_e_integrals, ONLY : overlap_integral
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: true_eeee

    CONTAINS

    SUBROUTINE true_eeee(A, B, R, e1, o1, e2, o2, max_m, &
                        &overlap_one, overlap_two, scratch, &
                        &eri_one_two, eri_two_one)
        !
        ! This subroutine computes a four index array
        ! containing electron repulsion integrals between 
        ! two electrons in disjoint finite domains.
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of each domain
        !   [real] R
        ! Distance between the centroids of each domain
        !   [int ] e1, o1
        ! Number of even and odd polynomials in the first 
        ! domain
        !   [int ] e2, o2
        ! Number of even and odd polynomials in the second 
        ! domain
        !   [int ] max_m
        ! Maximum value of intermediates (mu + nu) required
        !   [real] overlap_one, overlap_two
        ! Overlap matrices for the two domains
        ! 
        ! OUTPUT :
        !   [real] scratch(:,:)
        ! Scratch for storing intermediates
        ! Each dimension must be at least max_m + 2
        !   [real] eri_one_two(:,:,:,:), eri_two_one(:,:,:,:)
        ! Four index arrays for storing integrals
        !
        IMPLICIT NONE

        REAL(dp), INTENT(IN) :: A, B, R
        INTEGER,  INTENT(IN) :: e1, o1, e2, o2, max_m
        REAL(dp), INTENT(IN), TARGET :: overlap_one(:,:), overlap_two(:,:)

        REAL(dp), INTENT(OUT) :: scratch(-1:,-1:)
        REAL(dp), TARGET, INTENT(OUT) :: eri_one_two(:,:,:,:), &
                                       & eri_two_one(:,:,:,:)

        REAL(dp) :: frac
        INTEGER  :: m, n, l, s
        REAL(dp), POINTER :: ov_EE_one(:,:), ov_EE_two(:,:), &
                           & ov_OO_one(:,:), ov_OO_two(:,:)
        REAL(dp), POINTER :: eri_EEEE_one(:,:,:,:), eri_EEEE_two(:,:,:,:), &
                           & eri_EEOO_one(:,:,:,:), eri_EEOO_two(:,:,:,:), &
                           & eri_OOEE_one(:,:,:,:), eri_OOEE_two(:,:,:,:), &
                           & eri_OOOO_one(:,:,:,:), eri_OOOO_two(:,:,:,:)
        REAL(dp), POINTER :: eri_EOEO_one(:,:,:,:), eri_EOEO_two(:,:,:,:), &
                           & eri_EOOE_one(:,:,:,:), eri_EOOE_two(:,:,:,:), &
                           & eri_OEEO_one(:,:,:,:), eri_OEEO_two(:,:,:,:), &
                           & eri_OEOE_one(:,:,:,:), eri_OEOE_two(:,:,:,:)
        REAL(dp), POINTER :: eri_EEEO_one(:,:,:,:), eri_EEEO_two(:,:,:,:), &
                           & eri_EEOE_one(:,:,:,:), eri_EEOE_two(:,:,:,:), &
                           & eri_EOEE_one(:,:,:,:), eri_EOEE_two(:,:,:,:), &
                           & eri_OEEE_one(:,:,:,:), eri_OEEE_two(:,:,:,:)
        REAL(dp), POINTER :: eri_EOOO_one(:,:,:,:), eri_EOOO_two(:,:,:,:), &
                           & eri_OEOO_one(:,:,:,:), eri_OEOO_two(:,:,:,:), &
                           & eri_OOEO_one(:,:,:,:), eri_OOEO_two(:,:,:,:), &
                           & eri_OOOE_one(:,:,:,:), eri_OOOE_two(:,:,:,:)

        ! Associate pointers to chunks of the overlap matrices
        ov_EE_one => overlap_one(1:e1, 1:e1)
        ov_EE_two => overlap_two(1:e2, 1:e2)
        ov_OO_one => overlap_one(e1+1:e1+o1, e1+1:e1+o1)
        ov_OO_two => overlap_two(e2+1:e2+o2, e2+1:e2+o2)

        ! Associate pointers to chunks of the eri arrays
        eri_EEEE_one => eri_one_two(1:e1, 1:e1, 1:e2, 1:e2)
        eri_EEEE_two => eri_two_one(1:e2, 1:e2, 1:e1, 1:e1)

        eri_EEOO_one => eri_one_two(1:e1, 1:e1, e2+1:e2+o2, e2+1:e2+o2)
        eri_EEOO_two => eri_two_one(1:e2, 1:e2, e1+1:e1+o1, e1+1:e1+o1)
        eri_OOEE_one => eri_one_two(e1+1:e1+o1, e1+1:e1+o1, 1:e2, 1:e2)
        eri_OOEE_two => eri_two_one(e2+1:e2+o2, e2+1:e2+o2, 1:e1, 1:e1)

        eri_OOOO_one => eri_one_two(e1+1:e1+o1, e1+1:e1+o1, e2+1:e2+o2, &
                                     &e2+1:e2+o2)
        eri_OOOO_two => eri_two_one(e2+1:e2+o2, e2+1:e2+o2, e1+1:e1+o1, &
                                     &e1+1:e1+o1)

        eri_EOEO_one => eri_one_two(1:e1, e1+1:e1+o1, 1:e2, e2+1:e2+o2)
        eri_EOEO_two => eri_two_one(1:e2, e2+1:e2+o2, 1:e1, e1+1:e1+o1)

        eri_EOOE_one => eri_one_two(1:e1, e1+1:e1+o1, e2+1:e2+o2, 1:e2)
        eri_EOOE_two => eri_two_one(1:e2, e2+1:e2+o2, e1+1:e1+o1, 1:e1)

        eri_OEEO_one => eri_one_two(e1+1:e1+o1, 1:e1, 1:e2, e2+1:e2+o2)
        eri_OEEO_two => eri_two_one(e2+1:e2+o2, 1:e2, 1:e1, e1+1:e1+o1)

        eri_OEOE_one => eri_one_two(e1+1:e1+o1, 1:e1, e2+1:e2+o2, 1:e2)
        eri_OEOE_two => eri_two_one(e2+1:e2+o2, 1:e2, e1+1:e1+o1, 1:e1)

        eri_EEEO_one => eri_one_two(1:e1, 1:e1, 1:e2, e2+1:e2+o2)
        eri_EEEO_two => eri_two_one(1:e2, 1:e2, 1:e1, e1+1:e1+o1)
        eri_EEOE_one => eri_one_two(1:e1, 1:e1, e2+1:e2+o2, 1:e2)
        eri_EEOE_two => eri_two_one(1:e2, 1:e2, e1+1:e1+o1, 1:e1)
        eri_EOEE_one => eri_one_two(1:e1, e1+1:e1+o1, 1:e2, 1:e2)
        eri_EOEE_two => eri_two_one(1:e2, e2+1:e2+o2, 1:e1, 1:e1)
        eri_OEEE_one => eri_one_two(e1+1:e1+o1, 1:e1, 1:e2, 1:e2)
        eri_OEEE_two => eri_two_one(e2+1:e2+o2, 1:e2, 1:e1, 1:e1)

        eri_EOOO_one => eri_one_two(1:e1, e1+1:e1+o1, e2+1:e2+o2, e2+1:e2+o2)
        eri_EOOO_two => eri_two_one(1:e2, e2+1:e2+o2, e1+1:e1+o1, e1+1:e1+o1)
        eri_OEOO_one => eri_one_two(e1+1:e1+o1, 1:e1, e2+1:e2+o2, e2+1:e2+o2)
        eri_OEOO_two => eri_two_one(e2+1:e2+o2, 1:e2, e1+1:e1+o1, e1+1:e1+o1)
        eri_OOEO_one => eri_one_two(e1+1:e1+o1, e1+1:e1+o1, 1:e2, e2+1:e2+o2)
        eri_OOEO_two => eri_two_one(e2+1:e2+o2, e2+1:e2+o2, 1:e1, e1+1:e1+o1)
        eri_OOOE_one => eri_one_two(e1+1:e1+o1, e1+1:e1+o1, e2+1:e2+o2, 1:e2)
        eri_OOOE_two => eri_two_one(e2+1:e2+o2, e2+1:e2+o2, e1+1:e1+o1, 1:e1)

        ! Construct EEEE prototype integrals
        CALL compute_eeee_mn_integrals(A, B, R, max_m + 1, scratch)

        ! Pack the EEEE type integrals
        FORALL (m = 1:e1, n = 1:e1, l = 1:e2, s = 1:e2)
            eri_EEEE_one(m, n, l, s) = ov_EE_one(m,n) * ov_EE_two(l,s) * &
                & scratch(m + n, l + s)
            eri_EEEE_two(l, s, m, n) = ov_EE_one(m,n) * ov_EE_two(l,s) * &
                & scratch(m + n, l + s)
        END FORALL

        ! Pack the EEOO type integrals

        FORALL (m = 1:e1, n = 1:e1, l = 1:o2, s = 1:o2)
            eri_EEOO_one(m, n, l, s) = ov_EE_one(m,n) * ov_OO_two(l,s) * &
                & (dble(3 + 2 * (l + s)) * scratch(m + n, l + s) - &
                &  dble(2 + 2 * (l + s)) * scratch(m + n, l + s + 1))
            eri_OOEE_two(l, s, m, n) = eri_EEOO_one(m, n, l, s)
        END FORALL

        FORALL (m = 1:o1, n = 1:o1, l = 1:e2, s = 1:e2)
            eri_OOEE_one(m, n, l, s) = ov_OO_one(m,n) * ov_EE_two(l,s) * &
                & (dble(3 + 2 * (m + n)) * scratch(m + n, l + s) - &
                &  dble(2 + 2 * (m + n)) * scratch(m + n + 1, l + s))
            eri_EEOO_two(l, s, m, n) = eri_OOEE_one(m, n, l, s)
        END FORALL

        ! Pack the OOOO type integrals
        FORALL (m = 1:o1, n = 1:o1, l = 1:o2, s = 1:o2)
            eri_OOOO_one(m, n, l, s) = ov_OO_one(m,n) * ov_OO_two(l,s) * &
                & (dble((3 + 2 * (m + n)) * (3 + 2 * (l + s))) * &
                &  scratch(m + n, l + s) - &
                &  dble((3 + 2 * (m + n)) * (2 + 2 * (l + s))) * &
                &  scratch(m + n, l + s + 1) - &
                &  dble((2 + 2 * (m + n)) * (3 + 2 * (l + s))) * &
                &  scratch(m + n + 1, l + s) + &
                &  dble((2 + 2 * (m + n)) * (2 + 2 * (l + s))) * &
                &  scratch(m + n + 1, l + s + 1) )
            eri_OOOO_two(l, s, m, n) = eri_OOOO_one(m, n, l, s)
        END FORALL

        ! Pack the EOEO type integrals
        frac = - 4.d0 * B / A
        FORALL (m = 1:e1, n = 1:o1, l = 1:e2, s = 1:o2)
            eri_EOEO_one(m, n, l, s) = overlap_integral(m, n, 2) * &
                & overlap_integral(l, s, 2) * &
                & (frac * (0.5d0 + m + n) * (1.5d0 + m + n) / &
                & sqrt(dble((3 + 4 * m) * (3 + 4 * l))) * &
                & (scratch(l + s + 1, m + n - 1) - scratch(l + s + 1, m + n)))

            eri_EOOE_one(m, n, s, l) = eri_EOEO_one(m, n, l, s)
            eri_OEOE_one(n, m, s, l) = eri_EOEO_one(m, n, l, s)
            eri_OEEO_one(n, m, l, s) = eri_EOEO_one(m, n, l, s)
        END FORALL

        frac = - 4.d0 * B / A
        FORALL (m = 1:e1, n = 1:o1, l = 1:e2, s = 1:o2)
            eri_EOEO_two(l, s, m, n) = overlap_integral(m, n, 2) * &
                & overlap_integral(l, s, 2) * &
                & (frac * (0.5d0 + m + n) * (1.5d0 + m + n) / &
                & sqrt(dble((3 + 4 * m) * (3 + 4 * l))) * &
                & (scratch(l + s + 1, m + n - 1) - scratch(l + s + 1, m + n)))

            eri_EOOE_two(l, s, n, m) = eri_EOEO_two(l, s, m, n)
            eri_OEEO_two(s, l, m, n) = eri_EOEO_two(l, s, m, n)
            eri_OEOE_two(s, l, n, m) = eri_EOEO_two(l, s, m, n)
        END FORALL

        ! Reset scratch space
        FORALL (m=-1:max_m+1, n=-1:max_m+1)
            scratch(m,n) = 0.d0
        END FORALL

        ! Construct EEEO prototype integrals
        CALL compute_eeeo_mn_integrals(A, B, R, max_m + 1, scratch)

        ! Pack the EEEO type integrals

        FORALL (m = 1:e1, n = 1:e1, l = 1:e2, s = 1:o2)
            eri_EEEO_one(m, n, l, s) = - ov_EE_one(m,n) * overlap_integral(l, s, 2) * &
                & B / sqrt(dble(3 + 4 * l)) * scratch(m + n, l + s)

            eri_EEOE_one(m, n, s, l) = eri_EEEO_one(m, n, l, s)
            eri_EOEE_two(l, s, m, n) = eri_EEEO_one(m, n, l, s)
            eri_OEEE_two(s, l, m, n) = eri_EOEE_two(l, s, m, n)
        END FORALL

        ! Pack the EOOO type integrals

        FORALL (m = 1:e1, n = 1:o1, l = 1:o2, s = 1:o2)
            eri_EOOO_two(l, s, m, n) = - overlap_integral(l, s, 2) * ov_OO_one(m,n) * &
                & (B / sqrt(dble(3 + 4 * l))) * &
                & (dble(3 + 2 * (m + n)) * scratch(m + n, l + s) - &
                &  dble(2 + 2 * (m + n)) * scratch(m + n + 1, l + s))

            eri_OEOO_two(s, l, m, n) = eri_EOOO_two(l, s, m, n)
            eri_OOEO_one(m, n, l, s) = eri_EOOO_two(l, s, m, n)
            eri_OOOE_one(m, n, s, l) = eri_EOOO_two(l, s, m, n)
        END FORALL

        ! Construct EEEO prototype integrals
        CALL compute_eeeo_mn_integrals(B, A, R, max_m + 1, scratch)

        ! Pack the EEEO type integrals

        FORALL (m = 1:e1, n = 1:o1, l = 1:e2, s = 1:e2)
            eri_EEEO_two(l, s, m, n) = overlap_integral(m, n, 2) * ov_EE_two(l,s) * &
                & A / sqrt(dble(3 + 4 * m)) * scratch(l + s, m + n)

            eri_EEOE_two(l, s, n, m) = eri_EEEO_two(l, s, m, n)
            eri_EOEE_one(m, n, l, s) = eri_EEEO_two(l, s, m, n)
            eri_OEEE_one(n, m, l, s) = eri_EOEE_one(m, n, l, s)
        END FORALL

        ! Pack the EOOO type integrals

        FORALL (m = 1:o1, n = 1:o1, l = 1:e2, s = 1:o2)
            eri_EOOO_one(m, n, l, s) = ov_OO_two(l,s) * overlap_integral(m, n, 2) * & 
                & (A / sqrt(dble(3 + 4 * m))) * &
                & (dble(3 + 2 * (l + s)) * scratch(l + s, m + n) - &
                &  dble(2 + 2 * (l + s)) * scratch(l + s + 1, m + n))

            eri_OEOO_one(n, m, l, s) = eri_EOOO_one(m, n, l, s)
            eri_OOEO_two(l, s, m, n) = eri_EOOO_one(m, n, l, s)
            eri_OOOE_two(l, s, n, m) = eri_EOOO_one(m, n, l, s)
        END FORALL


    END SUBROUTINE true_eeee


    SUBROUTINE compute_eeee_mn_integrals(A, B, R, max_m, output)
        !
        ! This subroutine computes the integrals between an
        ! electron in an EE shell-pair and one in another
        ! EE shell-pair. This is done by recursion in m, n
        ! space, where m and n are compound descriptions of
        ! the orbital indices of the two shell-pairs.
        !
        ! These integrals may also be used to obtain the 
        ! EEOO, OOOO and EOEO classes of integrals
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of each domain
        !   [real] R
        ! Distance between the centroids of the domains
        !   [int ] max_m
        ! Maximum order of orbitals
        !
        ! OUTPUT :
        !   [real] output(:,:)
        ! Array containing the final integrals
        ! Each dimensions must be at least max_m + 2
        !
        IMPLICIT NONE

        REAL(dp), INTENT(IN) :: A, B, R
        INTEGER, INTENT(IN)  :: max_m
        REAL(dp), INTENT(INOUT) :: output(-1:,-1:)

        INTEGER  :: m, n
        REAL(dp) :: frac

        ! Start by generating the beginning values for the recurrence
        CALL eeee_m_minus_one(A, B, R, 1, max_m, output(:,1), &
                             &output(:,2), output(1:,-1))

        CALL eeee_m_zero(A, B, R, max_m, output(0:,1), &
                        &output(0:,2), output(0:,0))

        FORALL (m = max_m-1:max_m, n = 1:max_m)
            output(m, n) = p_series_eeee_int(A, B, R, n, m)
        END FORALL

        ! Perform the recurrence
        DO n = 1, max_m
            DO m = max_m - 2, 1, -1
                frac = (B**2 * (n - 0.5d0) * (n + 0.5d0)) / &
                     & (A**2 * (m + 1.5d0) * (m + 2.5d0))
                output(m, n) = output(m + 1, n) + frac * &
                             & (output(m + 2, n - 2) - &
                             &  output(m + 2, n - 1))
            ENDDO
        ENDDO

    END SUBROUTINE compute_eeee_mn_integrals

    SUBROUTINE compute_eeeo_mn_integrals(A, B, R, max_m, output)
        !
        ! This subroutine computes the integrals between an
        ! electron in an EE shell-pair and one in another
        ! EO shell-pair. This is done by recursion in m, n
        ! space, where m and n are compound descriptions of
        ! the orbital indices of the two shell-pairs.
        !
        ! These integrals may also be used to obtain the 
        ! EOOO class of integrals
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of each domain
        !   [real] R
        ! Distance between the centroids of the domains
        !   [int ] max_m
        ! Maximum order of orbitals
        !
        ! OUTPUT :
        !   [real] output(:,:)
        ! Array containing the final integrals
        ! Each dimensions must be at least max_m + 2
        !
        IMPLICIT NONE

        REAL(dp), INTENT(IN) :: A, B, R
        INTEGER, INTENT(IN)  :: max_m
        REAL(dp), INTENT(INOUT) :: output(-1:,-1:)

        INTEGER  :: m, n
        REAL(dp) :: frac

        ! Start by generating the beginning values for the recurrence
!           CALL eeeo_m_minus_one(A, B, R, 1, max_m, output(:,1), &
!                                &output(:,2), output(1:,-1))
!           CALL eeeo_m_zero(A, B, R, max_m, output(:,1), &
!                           &output(:,2), output(0:,0))
!           FORALL (m = max_m-1:max_m, n = 1:max_m)
!               output(m, n) = p_series_eeeo_int(A, B, R, m, n)
!           END FORALL

        ! NOTE : This is old code that treats the arrays as row major
        !        rather than column major, resulting in array temporaries
        !        being created. This code is correct, the new code may not be

        ! Start by generating the beginning values for the recurrence
        CALL eeeo_m_minus_one(A, B, R, 1, max_m, output(1,:), &
                             &output(2,:), output(-1,1:))
        CALL eeeo_m_zero(A, B, R, max_m, output(1,:), &
                        &output(2,:), output(0,0:))
        FORALL (m = 1:max_m, n = max_m-1:max_m)
            output(m, n) = p_series_eeeo_int(A, B, R, m, n)
        END FORALL

        ! Perform the recurrence
        DO n = max_m - 2, 1, -1
            DO m = 1, max_m
                frac = (B**2 * (m - 0.5d0) * (m + 0.5d0)) / &
                     & (A**2 * (n + 2.5d0) * (n + 3.5d0))
                output(m, n) = output(m, n + 1) + frac * &
                             & (output(m - 2, n + 2) - &
                             &  output(m - 1, n + 2))
            ENDDO
        ENDDO

    END SUBROUTINE compute_eeeo_mn_integrals

    SUBROUTINE eeee_m_minus_one(A, B, R, min_n, max_n, rma_array, &
                               &rpa_array, output)
        !
        ! This subroutine computes the integrals between 
        ! two EE basis functions pairs in disjoint regions
        ! of space when for one pair m = mu + nu = -1
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of each domain
        !   [real] R
        ! Distance between the centroids of the domains
        !   [int ] min_n, max_n
        ! Minimum and maximum order of orbitals
        !
        ! OUTPUT :
        !   [real] rma_array(:), rpa_array(:)
        ! Scratch arrays for holding intermediate hypergeometrics
        !   [real] output(:)
        ! Array containing the final integrals
        !

        INTEGER, INTENT(IN) :: min_n, max_n
        REAL(dp), INTENT(IN) :: A, B, R
        REAL(dp), INTENT(OUT), DIMENSION(min_n:max_n) :: output, &
                                                       & rma_array, &
                                                       & rpa_array

        INTEGER :: term, n
        REAL(dp) :: rma, rpa, summand_big, summand_less_big, p1_coeff, p2_coeff

        ! Compute the two arguments
        rma = (B / (R - A))**2
        rpa = (B / (R + A))**2

        ! We will compute the necessary hypergeometrics recursively
        ! Compute the starting values from power series
        summand_big = 1d0
        summand_less_big = 1d0
        term = 0

        ! Check for adjacent domains, which requires a simpler calculation
        ! We assume no miniscule domains exist; these will break things
        ! We also assume no ridulously large domains; which also break things
        IF (abs(1-rma).gt.1e-14) THEN
            rma_array(max_n) = 1d0
            rpa_array(max_n) = 1d0
            rma_array(max_n - 1) = 1d0
            rpa_array(max_n - 1) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_less_big.gt.1e-14)
                summand_big = summand_big * (0.5d0 + term) * (1d0 + term) &
                            & / ((1.5d0 + max_n + term) * (term + 1d0))
                summand_less_big = &
                    & summand_less_big * (0.5d0 + term) * (1d0 + term) &
                    & / ((1.5d0 + max_n - 1d0 + term) * (term + 1d0))
                term = term + 1

                rma_array(max_n) = rma_array(max_n) + &
                                 & summand_big * rma**term
                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rma_array(max_n - 1) = rma_array(max_n - 1) + &
                                     & summand_less_big * rma**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_less_big * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, min_n, -1
                p1_coeff = 2d0 / (3d0 + 2d0 * n)
                p2_coeff = (4d0 + 2d0 * n) / (5d0 + 2d0 * n)

                rma_array(n) = p1_coeff * rma_array(n + 1) * &
                    & (((2.5d0 + 2d0 * n) * rma - 1.5d0 - n) / (rma - 1d0)) &
                    & - p2_coeff * (rma / (rma - 1d0)) * rma_array(n + 2)
                rpa_array(n) = p1_coeff * rpa_array(n + 1) * &
                    & (((2.5d0 + 2d0 * n) * rpa - 1.5d0 - n) / (rpa - 1d0)) &
                    & - p2_coeff * (rpa / (rpa - 1d0)) * rpa_array(n + 2)
            ENDDO
        ELSE
            rpa_array(max_n) = 1d0
            rpa_array(max_n - 1) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_less_big.gt.1e-14)
                summand_big = summand_big * (0.5d0 + term) * (1d0 + term) &
                            & / ((1.5d0 + max_n + term) * (term + 1d0))
                summand_less_big = summand_less_big * &
                    & (0.5d0 + term) * (1d0 + term) / &
                    & ((1.5d0 + max_n - 1d0 + term) * (term + 1d0))
                term = term + 1

                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_less_big * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, min_n, -1
                p1_coeff = 2d0 / (3d0 + 2d0 * n)
                p2_coeff = (4d0 + 2d0 * n) / (5d0 + 2d0 * n)

                rpa_array(n) = p1_coeff * rpa_array(n + 1) * &
                    & (((2.5d0 + 2d0 * n) * rpa - 1.5d0 - n) / (rpa - 1d0)) &
                    & - p2_coeff * (rpa / (rpa - 1d0)) * rpa_array(n + 2)
            ENDDO

            ! Use simple expression to compute the rma array
            DO n = min_n, max_n
                rma_array(n) = 1d0 + 1d0 / (2d0 * n)
            ENDDO
        ENDIF

        ! Package the generated hypergeometrics into the output
        output = 0.5d0 * (rma_array / (R - A) + rpa_array / (R + A))
    END SUBROUTINE eeee_m_minus_one

    SUBROUTINE eeee_m_zero(A, B, R, max_n, rma_array, rpa_array, &
                          &output)
        !
        ! This subroutine computes the integrals between 
        ! two EE basis functions pairs in disjoint regions 
        ! of space when for one pair m = mu + nu = 0
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of the two domains
        !   [real] R
        ! Distance between the midpoint of the two domains
        !   [int ] max_n
        ! Maximum order n of the non-zero order EE pair
        !
        ! OUTPUT :
        !   [real] rma_array(:), rpa_array(:)
        ! Scratch arrays for holding intermediate
        ! hypergeometrics 
        ! Vectors of minimum length max_n + 1
        !   [real] output
        ! Array of length max_n, contains final integrals
        ! of order (0, n)
        !

        INTEGER, INTENT(IN) :: max_n
        REAL(dp), INTENT(IN) :: A, B, R
        REAL(dp), INTENT(OUT), DIMENSION(0:max_n) :: output, &
                                                   & rma_array, &
                                                   & rpa_array

        INTEGER :: term, n
        REAL(dp) :: rma, rpa, summand_big, summand_m1, summand_m2, &
                  & p1_coeff, p2_coeff, p3_coeff

        ! Compute the two arguments
        rma = (B / (R - A))**2
        rpa = (B / (R + A))**2

        ! We will compute the necessary hypergeometrics recursively
        ! Compute the starting values from power series
        summand_big = 1d0
        summand_m1  = 1d0
        summand_m2  = 1d0
        term = 0

        ! Check for adjacent domains, which requires a simpler calculation
        ! We assume no miniscule domains exist; these will break things
        ! We also assume no ridulously large domains; which also break things
        IF (abs(1-rma).gt.1e-14) THEN ! This is a nasty hack, it would be much
                                      ! better to pass the domain indices in
            rma_array(max_n) = 1d0
            rpa_array(max_n) = 1d0
            rma_array(max_n - 1) = 1d0
            rpa_array(max_n - 1) = 1d0
            rma_array(max_n - 2) = 1d0
            rpa_array(max_n - 2) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_m2.gt.1e-14) 
                summand_big = summand_big * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n + term) * (term + 1d0))
                summand_m1  = summand_m1 * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n - 1d0 + term) * (term + 1d0))
                summand_m2  = summand_m2 * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n - 2d0 + term) * (term + 1d0))
                term = term + 1

                rma_array(max_n) = rma_array(max_n) + &
                                 & summand_big * rma**term
                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rma_array(max_n - 1) = rma_array(max_n - 1) + &
                                     & summand_m1 * rma**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_m1 * rpa**term
                rma_array(max_n - 2) = rma_array(max_n - 2) + &
                                     & summand_m2 * rma**term
                rpa_array(max_n - 2) = rpa_array(max_n - 2) + &
                                     & summand_m2 * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 3, 0, -1
         
                p1_coeff = ((7d0 + 3d0 * n) * rma / (2.5d0 + n) - 2d0) / &
                         & (rma - 1d0)
                p2_coeff = (1d0 - (24.25d0 + 17d0 * n + 3d0 * n**2) / &
                         & ((2.5d0 + n) * (3.5d0 + n)) * rma) / &
                         & (rma - 1d0)
                p3_coeff = rma * (3d0 + n) * (3.5d0 + n) / &
                         & ((2.5d0 + n) * (4.5d0 + n) * (rma - 1.d0))
                rma_array(n) = p1_coeff * rma_array(n + 1) + &
                             & p2_coeff * rma_array(n + 2) + &
                             & p3_coeff * rma_array(n + 3)

                p1_coeff = ((7d0 + 3d0 * n) * rpa / (2.5d0 + n) - 2d0) / &
                         & (rpa - 1d0)
                p2_coeff = (1d0 - (24.25d0 + 17d0 * n + 3d0 * n**2) / &
                         & ((2.5d0 + n) * (3.5d0 + n)) * rpa) / &
                         & (rpa - 1d0)
                p3_coeff = rpa * (3d0 + n) * (3.5d0 + n) / &
                         & ((2.5d0 + n) * (4.5d0 + n) * (rpa - 1.d0))
                rpa_array(n) = p1_coeff * rpa_array(n + 1) + &
                             & p2_coeff * rpa_array(n + 2) + &
                             & p3_coeff * rpa_array(n + 3)
            ENDDO
        ELSE
            rpa_array(max_n) = 1d0
            rpa_array(max_n - 1) = 1d0
            rpa_array(max_n - 2) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_m2.gt.1e-14)
                summand_big = summand_big * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n + term) * (term + 1d0))
                summand_m1  = summand_m1 * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n - 1d0 + term) * (term + 1d0))
                summand_m2  = summand_m2 * (1d0 + term)**2 * &
                            & (1.5d0 + term) / ((2d0 + term) * &
                            & (2.5d0 + max_n - 2d0 + term) * (term + 1d0))
                term = term + 1

                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_m1 * rpa**term
                rpa_array(max_n - 2) = rpa_array(max_n - 2) + &
                                     & summand_m2 * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 3, 0, -1
                p1_coeff = ((7d0 + 3d0 * n) * rpa / (2.5d0 + n) - 2d0) / &
                         & (rpa - 1d0)
                p2_coeff = (1d0 - (24.25d0 + 17d0 * n + 3d0 * n**2) / &
                         & ((2.5d0 + n) * (3.5d0 + n)) * rpa) / &
                         & (rpa - 1d0)
                p3_coeff = rpa * ((3d0 + n) * (3.5d0 + n)) / &
                         & ((2.5d0 + n) * (4.5d0 + n) * (rpa - 1.d0))

                rpa_array(n) = p1_coeff * rpa_array(n + 1) + &
                             & p2_coeff * rpa_array(n + 2) + &
                             & p3_coeff * rpa_array(n + 3)
            ENDDO

            ! Compute rma array from simplified expression
            rma_array(0) = log(4d0) - 2d0
            DO n = 1, max_n
                rma_array(n) = rma_array(n-1) + 1d0 / (n + 2d0 * n**2)
                rma_array(n-1) = - (1d0 + 2d0 * n) * rma_array(n-1)
            ENDDO
            rma_array(max_n) = - (3d0 + 2d0 * max_n) * rma_array(max_n)
        ENDIF

        DO n = 0, max_n
            ! Package the generated hypergeometrics into the output
            output(n) = (rma * rma_array(n) / (3d0 + 2d0 * n) - &
                      &  rpa * rpa_array(n) / (3d0 + 2d0 * n) - &
                      &  2d0 * (log(R - A) - log(R + A))) / &
                      & (4d0 * A)
        ENDDO
    END SUBROUTINE eeee_m_zero

    REAL(dp) PURE FUNCTION p_series_eeee_int(A, B, R, m, n)
        ! WARNING: Requires variable length working arrays
        !
        ! Generates the integral between two ee pairs from
        ! a power series definition. This converges quickly
        ! when either m or n is large
        !
        ! INPUT :
        !   [real] A, B
        ! Half-length of each domain
        !   [real] R
        ! Distance between the centroids of each domain
        !   [int ] m, n
        ! Compund order (mu + nu) of each function pair
        !
        INTEGER, INTENT(IN) :: m, n
        REAL(dp), INTENT(IN) :: A, B, R

        INTEGER :: k, tk, max_k, l, tl, max_l
        REAL(dp) :: Asq_Rsq, Bsq_Rsq, term, num, den, output
        REAL(dp), DIMENSION(0:1000) :: k_zero_array, l_zero_array

        Asq_Rsq = A**2 / R**2
        Bsq_Rsq = B**2 / R**2

        term = 1d0 / R
        k_zero_array(0) = term
        l = 1
        DO WHILE (term.gt.1e-15)
            term = term * Bsq_Rsq * (l - 0.5d0) / (0.5d0 + l + n)
            k_zero_array(l) = term
            l = l + 1
        ENDDO
        max_l = l - 1

        term = 1d0 / R
        l_zero_array(0) = term
        k = 1
        DO WHILE (term.gt.1e-15)
            term = term * Asq_Rsq * (k - 0.5d0) / (0.5d0 + k + m)
            l_zero_array(k) = term
            k = k + 1
        ENDDO
        max_k = k - 1

        output = SUM(l_zero_array(:max_k)) + SUM(k_zero_array(1:max_l))

        DO k = 0, max_k
            term = l_zero_array(k)
            tk = k
            tl = 0
            DO WHILE (term.gt.1e-15)
                num = (0.5d0 + tk + tl) * (1.d0 + tk + tl) * &
                    & (1.5d0 + tk + tl) * (2.d0 + tk + tl)
                den = (1.d0 + tk) * (1.d0 + tl) * &
                    & (1.5d0 + tk + m) * (1.5d0 + tl + n)
                term = term * Asq_Rsq * Bsq_Rsq * num / den

                output = output + term
                tk = tk + 1
                tl = tl + 1
            ENDDO
        ENDDO

        DO l = 1, max_l
            term = k_zero_array(l)
            tk = 0
            tl = l
            DO WHILE (term.gt.1e-15)
                num = (0.5d0 + tk + tl) * (1.d0 + tk + tl) * &
                    & (1.5d0 + tk + tl) * (2.d0 + tk + tl)
                den = (1.d0 + tk) * (1.d0 + tl) * &
                    & (1.5d0 + tk + m) * (1.5d0 + tl + n)
                term = term * Asq_Rsq * Bsq_Rsq * num / den

                output = output + term
                tk = tk + 1
                tl = tl + 1
            ENDDO
        ENDDO

        p_series_eeee_int = output
    END FUNCTION p_series_eeee_int



    SUBROUTINE eeeo_m_minus_one (A, B, R, min_n, max_n, rma_array, &
                                &rpa_array, output)
        !
        ! This subroutine computes the integrals between 
        ! one EE basis function pair and one EO basis function
        ! pair in disjoint regions of space when for one pair
        ! m = mu + nu = -1
        !
        ! INPUT :
        !   [real] A, B  = Half the length of each domain
        !   [real] R     = Distance between the centroids of 
        !                  the domains
        !   [int ] min_n = Minumum order of orbitals
        !                  Must be at least 1
        !   [int ] max_n = Maximum order of orbitals
        !
        ! OUTPUT :
        !   [real] rma_array(:)
        !   [real] rpa_array(:) = Scratch arrays for holding
        !                         intermediate hypergeometrics
        !   [real] output(:)    = Array containing the final 
        !                         integrals
        !

        INTEGER, INTENT(IN) :: min_n, max_n
        REAL(dp), INTENT(IN) :: A, B, R
        REAL(dp), INTENT(OUT), DIMENSION(min_n:max_n) :: output, &
                                                       & rma_array, &
                                                       & rpa_array

        INTEGER :: term, n
        REAL(dp) :: rma, rpa, summand_big, summand_less_big, p2_coeff

        ! Compute the two arguments
        rma = (B / (R - A))**2
        rpa = (B / (R + A))**2

        ! We will compute the necessary hypergeometrics recursively
        ! Compute the starting values from power series
        summand_big = 1d0
        summand_less_big = 1d0
        term = 0

        ! Check for adjacent domains, which requires a simpler calculation
        ! We assume no miniscule domains exist; these will break things
        ! We also assume no ridulously large domains; which also break things
        IF (abs(1-rma).gt.1e-14) THEN
            rma_array(max_n) = 1.d0
            rpa_array(max_n) = 1.d0
            rma_array(max_n - 1) = 1.d0
            rpa_array(max_n - 1) = 1.d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_less_big.gt.1e-14)
                summand_big = summand_big * (term - 0.5d0) * (1d0 + term) &
                            & / ((2.5d0 + max_n + term) * (term + 1d0))
                summand_less_big = &
                    & summand_less_big * (term - 0.5d0) * (1d0 + term) &
                    & / ((2.5d0 + max_n - 1d0 + term) * (term + 1d0))
                term = term + 1

                rma_array(max_n) = rma_array(max_n) + &
                                 & summand_big * rma**term
                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rma_array(max_n - 1) = rma_array(max_n - 1) + &
                                     & summand_less_big * rma**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_less_big * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, min_n, -1
                !p1_coeff = 2d0 / (3d0 + 2d0 * n)
                p2_coeff = dble(2 * (4 + n)) / dble(7 + 2 * n)

                rma_array(n) = rma_array(n + 1) * &
                    & ((dble(11 + 4 * n) * rma - dble(5 + 2 * n)) / &
                    & (dble(5 + 2 * n) * (rma - 1.d0))) - &
                    & p2_coeff * (rma / (rma - 1d0)) * rma_array(n + 2)
                rpa_array(n) = rpa_array(n + 1) * &
                    & ((dble(11 + 4 * n) * rpa - dble(5 + 2 * n)) / &
                    & (dble(5 + 2 * n) * (rpa - 1.d0))) - &
                    & p2_coeff * (rpa / (rpa - 1d0)) * rpa_array(n + 2)
            ENDDO

            ! Apply hypergeometric coefficient and add the second easy term
            FORALL (n = min_n:max_n)
                rma_array(n) = rma_array(n) * &
                    & dble(4 * (1 + n) * (2 + n)) * (R - A)**2
                rpa_array(n) = rpa_array(n) * &
                    & dble(4 * (1 + n) * (2 + n)) * (R + A)**2

                rma_array(n) = - rma_array(n) + &
                    & dble(3 + 2 * n) * (dble(3 + 2 * n) * (R - A)**2 - B**2)
                rpa_array(n) = - rpa_array(n) + &
                    & dble(3 + 2 * n) * (dble(3 + 2 * n) * (R + A)**2 - B**2)

                rma_array(n) = rma_array(n) / ((R - A)**2 - B**2)**2
                rpa_array(n) = rpa_array(n) / ((R + A)**2 - B**2)**2
            END FORALL

        ELSE
            rpa_array(max_n) = 1d0
            rpa_array(max_n - 1) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_less_big.gt.1e-14)
                summand_big = summand_big * (term + 1.5d0) * (1d0 + term) &
                            & / ((2.5d0 + max_n + term) * (term + 1d0))
                summand_less_big = &
                    & summand_less_big * (term + 1.5d0) * (1d0 + term) &
                    & / ((2.5d0 + max_n - 1d0 + term) * (term + 1d0))
                term = term + 1

                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_less_big * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, min_n, -1
                p2_coeff = dble(2 * (2 + n)) / dble(7 + 2 * n)

                rpa_array(n) = rpa_array(n + 1) * &
                    & ((dble(7 + 4 * n) * rpa - dble(5 + 2 * n)) / &
                    & (dble(5 + 2 * n) * (rpa - 1.d0))) - &
                    & p2_coeff * (rpa / (rpa - 1d0)) * rpa_array(n + 2)
            ENDDO

            ! Apply hypergeometric coefficient
            FORALL (n = min_n:max_n)
                rpa_array(n) = rpa_array(n) / (R + A)**2
            END FORALL

            ! Use simple expression to compute the rma array
            FORALL (n = min_n:max_n)
                rma_array(n) = dble(2 * n + 3) / (dble(2 * n) * B**2)
            END FORALL

        ENDIF

        ! Package the generated hypergeometrics into the output
        output = (rma_array + rpa_array) / 2.d0
    END SUBROUTINE eeeo_m_minus_one

    SUBROUTINE eeeo_m_zero(A, B, R, max_n, rma_array, rpa_array, &
                          &output)
        !
        ! This subroutine computes the integrals between 
        ! one EE basis function pair and one EO basis function 
        ! pair in disjoint regions of space when for one pair
        ! m = mu + nu = 0
        !
        ! INPUT :
        !   [real] A, B   = Half the length of the two domains
        !   [real] R      = Distance between the midpoint of 
        !                   the two domains
        !   [int ] max_n  = Maximum order n for the other EE pair
        !
        ! OUTPUT :
        !   [real] rma_array(:)
        !   [real] rpa_array(:) = Scratch arrays for holding
        !                         intermediate hypergeometrics.
        !                         Must be indexed from 0 to max_n
        !   [real] output = Array of length max_n, contains final 
        !                   integrals of order (0, n)
        !

        INTEGER, INTENT(IN) :: max_n
        REAL(dp), INTENT(IN) :: A, B, R
        REAL(dp), INTENT(OUT), DIMENSION(0:max_n) :: output, &
                                                   & rma_array, &
                                                   & rpa_array

        INTEGER :: term, n
        REAL(dp) :: rma, rpa, summand_big, summand_m1, p1_coeff, p2_coeff

        ! Compute the two arguments
        rma = B**2 / (R - A)**2
        rpa = B**2 / (R + A)**2

        ! We will compute the necessary hypergeometrics recursively
        ! Compute the starting values from power series
        summand_big = 1d0
        summand_m1  = 1d0
        term = 0

        ! Check for adjacent domains, which requires a simpler calculation
        ! We assume no miniscule domains exist; these will break things
        ! We also assume no ridulously large domains; which also break things
        IF (abs(1-rma).gt.1e-14) THEN
            rma_array(max_n) = 1d0
            rpa_array(max_n) = 1d0
            rma_array(max_n - 1) = 1d0
            rpa_array(max_n - 1) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_m1.gt.1e-14) 
                summand_big = summand_big * (0.5d0 + term) * &
                            & dble(1 + term) / ((2.5d0 + max_n + term) * &
                            & (term + 1d0))
                summand_m1  = summand_m1 * (0.5d0 + term) * &
                            & dble(1 + term) / ((term + 1d0) * &
                            & (2.5d0 + max_n - 1d0 + term)) 
                term = term + 1

                rma_array(max_n) = rma_array(max_n) + &
                                 & summand_big * rma**term
                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rma_array(max_n - 1) = rma_array(max_n - 1) + &
                                     & summand_m1 * rma**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_m1 * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, 0, -1
                p1_coeff = (dble(9 + 4 * n) * rma - dble(5 + 2 * n)) / &
                         & (dble(5 + 2 * n) * (rma - 1d0))
                p2_coeff = dble(2 * (3 + n)) / &
                         & dble(7 + 2 * n)
                rma_array(n) = p1_coeff * rma_array(n + 1) - &
                             & p2_coeff * (rma / (rma - 1.d0)) * &
                             & rma_array(n + 2)

                p1_coeff = (dble(9 + 4 * n) * rpa - dble(5 + 2 * n)) / &
                         & (dble(5 + 2 * n) * (rpa - 1d0))
                rpa_array(n) = p1_coeff * rpa_array(n + 1) - &
                             & p2_coeff * (rpa / (rpa - 1.d0)) * &
                             & rpa_array(n + 2)
            ENDDO
        ELSE
            rpa_array(max_n) = 1d0
            rpa_array(max_n - 1) = 1d0

            ! WARNING: This condition may not be satisfactory
            DO WHILE (summand_m1.gt.1e-14)
                summand_big = summand_big * (0.5d0 + term) * &
                            & dble(1 + term) / ((2.5d0 + max_n + term) * &
                            & (term + 1d0))
                summand_m1  = summand_m1 * (0.5d0 + term) * &
                            & dble(1 + term) / ((term + 1d0) * &
                            & (2.5d0 + max_n - 1d0 + term)) 
                term = term + 1

                rpa_array(max_n) = rpa_array(max_n) + &
                                 & summand_big * rpa**term
                rpa_array(max_n - 1) = rpa_array(max_n - 1) + &
                                     & summand_m1 * rpa**term
            ENDDO

            ! Recur backwards to obtain the rest
            DO n = max_n - 2, 0, -1
                p1_coeff = (dble(9 + 4 * n) * rpa - dble(5 + 2 * n)) / &
                         & (dble(5 + 2 * n) * (rpa - 1d0))
                p2_coeff = rpa * dble(2 * (3 + n)) / &
                         & (dble(7 + 2 * n) * (rpa - 1.d0))

                rpa_array(n) = p1_coeff * rpa_array(n + 1) - &
                             & p2_coeff * rpa_array(n + 2)
            ENDDO

            ! Compute rma array from simplified expression
            FORALL (n = 0:max_n)
                rma_array(n) = dble(3 + 2 * n) / dble(2 + 2 * n)
            END FORALL
        ENDIF

        DO n = 0, max_n
            ! Package the generated hypergeometrics into the output
            output(n) = (rma_array(n) / (R - A) - &
                      &  rpa_array(n) / (R + A)) / (2.d0 * A)
        ENDDO
    END SUBROUTINE eeeo_m_zero

    REAL(dp) PURE FUNCTION p_series_eeeo_int(A, B, R, m, n)
        ! WARNING: Requires variable length working arrays
        !
        ! generates an integral from a power series definition?
        !
        INTEGER, INTENT(IN) :: m, n
        REAL(dp), INTENT(IN) :: A, B, R

        INTEGER :: k, tk, max_k, l, tl, max_l
        REAL(dp) :: Asq_Rsq, Bsq_Rsq, term, num, den, output
        REAL(dp), DIMENSION(0:1000) :: k_zero_array, l_zero_array

        Asq_Rsq = A**2 / R**2
        Bsq_Rsq = B**2 / R**2

        term = 1d0 / R**2
        k_zero_array(0) = term
        l = 1
        DO WHILE (term.gt.1e-15)
            term = term * Bsq_Rsq * (l + 0.5d0) / (1.5d0 + l + n)
            k_zero_array(l) = term
            l = l + 1
        ENDDO
        max_l = l - 1

        term = 1d0 / R**2
        l_zero_array(0) = term
        k = 1
        DO WHILE (term.gt.1e-15)
            term = term * Asq_Rsq * (k + 0.5d0) / (0.5d0 + k + m)
            l_zero_array(k) = term
            k = k + 1
        ENDDO
        max_k = k - 1

        output = SUM(l_zero_array(:max_k)) + SUM(k_zero_array(1:max_l))

        DO k = 0, max_k
            term = l_zero_array(k)
            tk = k
            tl = 0
            DO WHILE (term.gt.1e-15)
                num = (1.d0 + tk + tl) * (1.5d0 + tk + tl) * &
                    & (2.d0 + tk + tl) * (2.5d0 + tk + tl)
                den = (1.d0 + tk) * (1.d0 + tl) * &
                    & (1.5d0 + tk + m) * (2.5d0 + tl + n)
                term = term * Asq_Rsq * Bsq_Rsq * num / den

                output = output + term
                tk = tk + 1
                tl = tl + 1
            ENDDO
        ENDDO

        DO l = 1, max_l
            term = k_zero_array(l)
            tk = 0
            tl = l
            DO WHILE (term.gt.1e-15)
                num = (1.d0 + tk + tl) * (1.5d0 + tk + tl) * &
                    & (2.d0 + tk + tl) * (2.5d0 + tk + tl)
                den = (1.d0 + tk) * (1.d0 + tl) * &
                    & (1.5d0 + tk + m) * (2.5d0 + tl + n)
                term = term * Asq_Rsq * Bsq_Rsq * num / den

                output = output + term
                tk = tk + 1
                tl = tl + 1
            ENDDO
        ENDDO

        p_series_eeeo_int = output
    END FUNCTION p_series_eeeo_int

END MODULE true_eri_poly
