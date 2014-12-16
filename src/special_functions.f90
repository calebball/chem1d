MODULE special_functions
    !
    ! SPECIAL_FUNCTIONS
    ! This module contains function definitions for various special
    ! functions required throughout the program. These are general
    ! purpose routines
    !
    USE constants

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: incomplete_gamma
    PUBLIC :: exp_incomplete_gamma
    PUBLIC :: two_f_one

    CONTAINS

        REAL(dp) PURE FUNCTION incomplete_gamma(a, x)
            ! WARNING! Needs error handling
            !          Needs better convergence tracking
            !          Needs better algorithm for small x
            !
            ! Computes the upper incomplete gamma function
            ! for negative integer orders using the continued
            ! fraction method from Numerical Recipes
            !
            ! INPUT :
            ! a ( int) = order
            ! x (real) = argument
            !
            IMPLICIT NONE
            INTEGER, INTENT(IN)  :: a
            REAL(dp), INTENT(IN) :: x

            REAL(dp), PARAMETER :: accuracy = 1.e-15, fpmin = epsilon(1.d0)
            INTEGER, PARAMETER  :: max_iterations = 1000
            INTEGER  :: i
            REAL(dp) :: an, b, c, d, delta, h

            b = x + 1.d0 - dble(a)
            c = 1.d0 / fpmin
            d = 1.d0 / b
            h = d

            DO i = 1, max_iterations
                an = dble(-i * (i - a))
                b = b + 2.d0

                d = an * d + b
                IF (abs(d).LT.fpmin) d = fpmin
                c = b + an / c
                IF (abs(c).LT.fpmin) c = fpmin
                d = 1.d0 / d

                delta = d * c
                h = h * delta
                IF (abs(delta - 1.d0).LT.accuracy) EXIT
            END DO
            incomplete_gamma = exp(-x + a * log(x)) * h
        END FUNCTION incomplete_gamma

        REAL(dp) PURE FUNCTION exp_incomplete_gamma(a, x)
            ! WARNING! Needs error handling
            !          Needs better convergence tracking
            !          Needs better algorithm for small x
            !
            ! Computes the upper incomplete gamma function
            ! for negative integer orders using the continued
            ! fraction method from Numerical Recipes
            !
            ! INPUT :
            ! a ( int) = order
            ! x (real) = argument
            !
            IMPLICIT NONE
            INTEGER, INTENT(IN)  :: a
            REAL(dp), INTENT(IN) :: x

            REAL(dp), PARAMETER :: accuracy = 1.e-15, fpmin = epsilon(1.d0)
            INTEGER, PARAMETER  :: max_iterations = 1000
            INTEGER  :: i
            REAL(dp) :: an, b, c, d, delta, h

            b = x + 1.d0 - dble(a)
            c = 1.d0 / fpmin
            d = 1.d0 / b
            h = d

            DO i = 1, max_iterations
                an = dble(-i * (i - a))
                b = b + 2.d0

                d = an * d + b
                IF (abs(d).LT.fpmin) d = fpmin
                c = b + an / c
                IF (abs(c).LT.fpmin) c = fpmin
                d = 1.d0 / d

                delta = d * c
                h = h * delta
                IF (abs(delta - 1.d0).LT.accuracy) EXIT
            END DO
            exp_incomplete_gamma = x**a * h
        END FUNCTION exp_incomplete_gamma

        REAL(dp) PURE FUNCTION two_f_one(a, b, c, z)
            ! WARNING! Needs error handling
            !
            ! Computes the Gauss hypergeometric function, 2F1,
            ! from its Taylor series expansion. Converges well
            ! when c is large
            !
            ! INPUT :
            ! a, b (real) = numerator parameters
            ! c    (real) = denominator parameter
            ! z    (real) = argument
            !
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: a, b, c, z

            INTEGER, PARAMETER :: max_iterations = 100
            INTEGER :: k
            REAL(dp) :: summand, total, convergence

            total = 1.d0
            summand = 1.d0
            convergence = 1.d0
            k = 1

            DO WHILE (convergence.GT.1.e-15.AND.k.LT.max_iterations)
                summand = summand * (a + k - 1.d0) * (b + k - 1.d0) * z / &
                        & dble(k * (c + k - 1.d0))
                convergence = total
                total = total + summand
                convergence = abs(convergence - total)
                k = k + 1
            END DO

            two_f_one = total
        END FUNCTION two_f_one

        REAL(dp) PURE FUNCTION two_f_two(a1, a2, b1, b2, z)
            !
            ! This function computes the 2F2 Gauss
            ! hypergeometric function
            !
            ! INPUT :
            !   [real] a1, a2
            ! Parameters appearing in the numerator
            !   [real] b1, b2
            ! Parameters appearing in the denominator
            !   [real] z
            ! Argument of the function
            !
            IMPLICIT NONE

            REAL(dp), INTENT(IN) :: a1, a2, b1, b2, z

            INTEGER  :: k
            REAL(dp) :: term, last_term, total

            ! Starting with the first k = 0 term
            total     = 1.d0
            k         = 1
            last_term = 1.d0

            ! WARNING : Hard coded threshold!
            DO WHILE (abs(last_term/total).GT.1.e-12)

                term = last_term
                term = term * ((a1 + k - 1) * (a2 + k - 1) * z) / &
                     & dble((b1 + k - 1) * (b2 + k - 1) * k)

                total     = total + term
                k         = k + 1
                last_term = term

            ENDDO

            two_f_two = total
        END FUNCTION two_f_two

END MODULE special_functions
