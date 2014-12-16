MODULE moller_plesset
    !
    ! MOLLER_PLESSET
    !
    ! Contains subroutines and data structures
    ! for running a Moller-Plesset calculation
    ! using the results from a previously 
    ! completed Hartree-Fock calculation
    !
    USE constants
    USE input
    USE storage, ONLY : hf_in_domain,  &
                      & mo_eri_domain, &
                      & mo_double_bar
    USE hartree_fock
    IMPLICIT NONE

    REAL(dp) :: mp2_correction, mp3_correction

    PRIVATE

!   PUBLIC :: mp2_correction
!   PUBLIC :: mp3_correction
    PUBLIC :: compute_mp2
    PUBLIC :: compute_mp3

    PUBLIC :: mp2_correction
    PUBLIC :: mp3_correction

    PUBLIC :: print_mp

    CONTAINS

!   SUBROUTINE mp2_correction(correction)
!       !
!       ! This subroutine generates the second-order
!       ! correction term in the perturbation series
!       !
!       ! The algorithm used here exploits the domain
!       ! structure
!       ! WARNING : This doesn't work when more than
!       !           one electron is in a domain
!       ! WARNING : No error handling
!       !
!       IMPLICIT NONE

!       REAL(dp), INTENT(INOUT) :: correction

!       INTEGER  :: d1, d2, a, b, r, s
!       REAL(dp) :: denominator

!       correction = 0.d0

!       DO d1 = 1, n_domains
!           IF (electrons_in_domain(d1).EQ.0) CYCLE

!           ! All double excitations within the domain
!           DO a = 1, electrons_in_domain(d1)
!           DO b = a + 1, electrons_in_domain(d1)
!               DO r = electrons_in_domain(d1) + 1, functions_in_domain(d1)
!               DO s = r + 1, functions_in_domain(d1)
!                   denominator = hf_in_domain(d1)%orb_energies(a) + &
!                               & hf_in_domain(d1)%orb_energies(b) - &
!                               & hf_in_domain(d1)%orb_energies(r) - &
!                               & hf_in_domain(d1)%orb_energies(s)

!                   correction = correction + &
!                              & mo_eri_domain(d1)%with(d1)%integral(a, r, b, s)**2 / denominator
!               ENDDO
!               ENDDO
!           ENDDO
!           ENDDO

!           ! All other double excitations
!           DO d2 = d1 + 1, n_domains
!               IF (electrons_in_domain(d2).EQ.0) CYCLE

!               ! Both excitations occuring within the domains
!               ! yields non-zero coulomb integrals
!               DO a = 1, electrons_in_domain(d1)
!               DO b = 1, electrons_in_domain(d2)
!                   DO r = electrons_in_domain(d1) + 1, functions_in_domain(d1)
!                   DO s = electrons_in_domain(d2) + 1, functions_in_domain(d2)
!                       denominator = hf_in_domain(d1)%orb_energies(a) + &
!                                   & hf_in_domain(d2)%orb_energies(b) - &
!                                   & hf_in_domain(d1)%orb_energies(r) - &
!                                   & hf_in_domain(d2)%orb_energies(s)

!                       correction = correction + &
!                                  & mo_eri_domain(d1)%with(d2)%integral(a, r, b, s)**2 / denominator
!                   ENDDO
!                   ENDDO
!               ENDDO
!               ENDDO

!               ! Excitations occuring from one domain to another
!               ! yields non-zero exchange integrals
!               DO a = 1, electrons_in_domain(d1)
!               DO b = a + 1, electrons_in_domain(d1)
!                   DO r = electrons_in_domain(d2) + 1, functions_in_domain(d2)
!                   DO s = r + 1, functions_in_domain(d2)
!                       denominator = hf_in_domain(d1)%orb_energies(a) + &
!                                   & hf_in_domain(d1)%orb_energies(b) - &
!                                   & hf_in_domain(d2)%orb_energies(r) - &
!                                   & hf_in_domain(d2)%orb_energies(s)

!                       ! No sign flip of the integral due to squaring
!                       correction = correction + &
!                                  & mo_eri_domain(d1)%with(d2)%integral(a, b, r, s)**2 / denominator
!                   ENDDO
!                   ENDDO
!               ENDDO
!               ENDDO

!               DO a = 1, electrons_in_domain(d2)
!               DO b = a + 1, electrons_in_domain(d2)
!                   DO r = electrons_in_domain(d1) + 1, functions_in_domain(d1)
!                   DO s = r + 1, functions_in_domain(d1)
!                       denominator = hf_in_domain(d2)%orb_energies(a) + &
!                                   & hf_in_domain(d2)%orb_energies(b) - &
!                                   & hf_in_domain(d1)%orb_energies(r) - &
!                                   & hf_in_domain(d1)%orb_energies(s)

!                       ! No sign flip of the integral due to squaring
!                       correction = correction + &
!                                  & mo_eri_domain(d1)%with(d2)%integral(r, s, a, b)**2 / denominator
!                   ENDDO
!                   ENDDO
!               ENDDO
!               ENDDO
!           ENDDO
!       ENDDO

!   END SUBROUTINE mp2_correction

!   SUBROUTINE mp3_correction(correction)
!       !
!       ! This subroutine generates the third-order
!       ! correction term in the perturbation series
!       !
!       ! The algorithm used here exploits the domain
!       ! structure
!       !
!       ! WARNING : This code was intended to exploit
!       !           the domain structure. Unfortunately
!       !           the initial attempt missed out some
!       !           integrals. The code has been 
!       !           removed pending a fresh attempt
!       !
!       IMPLICIT NONE

!       REAL(dp), INTENT(INOUT) :: correction

!       INTEGER  :: d1, d2, d3, a, b, c, d, r, s, t, u
!       REAL(dp) :: suma, sumb, sumc, numerator, denominator

!       correction = 0.d0

!       suma = 0.d0
!       sumb = 0.d0
!       sumc = 0.d0
!       correction = sumc + (suma + sumb) / 8.d0

!   END SUBROUTINE mp3_correction

    SUBROUTINE compute_mp2
        !
        ! This subroutine generates the second-order
        ! correction term in the perturbation series
        !
        ! A 'brute-force' algorithm is used here that
        ! makes no attempt to exploit the domain
        ! separation of the system
        !
        IMPLICIT NONE

        INTEGER  :: a, b, r, s
        INTEGER  :: da, db, dr, ds
        INTEGER  :: shifta, shiftb, shiftr, shifts
        REAL(dp) :: numerator, denominator

        mp2_correction = 0.d0

        ! choose 2 occupied orbitals
        DO da = 1, n_domains
        DO a = start_of_domain(da), start_of_domain(da) + electrons_in_domain(da) - 1
        DO db = 1, n_domains
        DO b = start_of_domain(db), start_of_domain(db) + electrons_in_domain(db) - 1

!           ! choose 2 virtual orbitals
!           DO dr = 1, n_domains
!           DO r = start_of_domain(dr) + electrons_in_domain(dr), end_of_domain(dr)
!           DO ds = 1, n_domains
!           DO s = start_of_domain(ds) + electrons_in_domain(ds), end_of_domain(ds)

            ! choose 2 virtual orbitals
            DO dr = 1, n_domains
            DO r = start_of_domain(dr) + electrons_in_domain(dr), &
                 & start_of_domain(dr) + functions_kept(dr) - 1
            DO ds = 1, n_domains
            DO s = start_of_domain(ds) + electrons_in_domain(ds), &
                 & start_of_domain(ds) + functions_kept(ds) - 1

                ! The mo_double_bar array indices refer to the
                ! basis functions in all the domains
                ! It is necessary to shift the indices when
                ! when using arrays that index individual
                ! domains
                shifta = start_of_domain(da) - 1
                shiftb = start_of_domain(db) - 1
                shiftr = start_of_domain(dr) - 1
                shifts = start_of_domain(ds) - 1

                numerator = mo_double_bar(a,b,r,s)**2

                denominator = hf_in_domain(da)%orb_energies(a-shifta) + &
                            & hf_in_domain(db)%orb_energies(b-shiftb) - &
                            & hf_in_domain(dr)%orb_energies(r-shiftr) - &
                            & hf_in_domain(ds)%orb_energies(s-shifts)

                mp2_correction = mp2_correction + numerator / denominator

            ENDDO
            ENDDO
            ENDDO
            ENDDO

        ENDDO
        ENDDO
        ENDDO
        ENDDO

        mp2_correction = mp2_correction / 4.d0

    END SUBROUTINE compute_mp2

    SUBROUTINE compute_mp3
        !
        ! This subroutine generates the third-order
        ! correction term in the perturbation series
        !
        ! A 'brute-force' algorithm is used here that
        ! makes no attempt to exploit the domain
        ! separation of the system
        !
        IMPLICIT NONE

        INTEGER  :: a, b, c, d, r, s, t, u
        INTEGER  :: da, db, dc, dd, dr, ds, dt, du
        INTEGER  :: shifta, shiftb, shiftc, shiftd, shiftr, shifts, shiftt, shiftu
        REAL(dp) :: suma, sumb, sumc, numerator, denominator

        mp3_correction = 0.d0
        suma = 0.d0
        sumb = 0.d0
        sumc = 0.d0

        ! suma
        ! choose 4 occupied orbitals
        DO da = 1, n_domains
        DO a = start_of_domain(da), start_of_domain(da) + electrons_in_domain(da) - 1
        DO db = 1, n_domains
        DO b = start_of_domain(db), start_of_domain(db) + electrons_in_domain(db) - 1
        DO dc = 1, n_domains
        DO c = start_of_domain(dc), start_of_domain(dc) + electrons_in_domain(dc) - 1
        DO dd = 1, n_domains
        DO d = start_of_domain(dd), start_of_domain(dd) + electrons_in_domain(dd) - 1

!           ! choose 2 virtual orbitals
!           DO dr = 1, n_domains
!           DO r = start_of_domain(dr) + electrons_in_domain(dr), end_of_domain(dr)
!           DO ds = 1, n_domains
!           DO s = start_of_domain(ds) + electrons_in_domain(ds), end_of_domain(ds)

            ! choose 2 virtual orbitals
            DO dr = 1, n_domains
            DO r = start_of_domain(dr) + electrons_in_domain(dr), &
                 & start_of_domain(dr) + functions_kept(dr) - 1
            DO ds = 1, n_domains
            DO s = start_of_domain(ds) + electrons_in_domain(ds), &
                 & start_of_domain(ds) + functions_kept(ds) - 1

                ! The mo_double_bar array indices refer to the
                ! basis functions in all the domains
                ! It is necessary to shift the indices when
                ! when using arrays that index individual
                ! domains
                shifta = start_of_domain(da) - 1
                shiftb = start_of_domain(db) - 1
                shiftc = start_of_domain(dc) - 1
                shiftd = start_of_domain(dd) - 1
                shiftr = start_of_domain(dr) - 1
                shifts = start_of_domain(ds) - 1

                numerator = mo_double_bar(a,b,r,s) * &
                          & mo_double_bar(c,d,a,b) * &
                          & mo_double_bar(r,s,c,d)

                denominator = (hf_in_domain(da)%orb_energies(a-shifta) +  &
                            &  hf_in_domain(db)%orb_energies(b-shiftb) -  &
                            &  hf_in_domain(dr)%orb_energies(r-shiftr) -  &
                            &  hf_in_domain(ds)%orb_energies(s-shifts)) * &
                            & (hf_in_domain(dc)%orb_energies(c-shiftc) +  &
                            &  hf_in_domain(dd)%orb_energies(d-shiftd) -  &
                            &  hf_in_domain(dr)%orb_energies(r-shiftr) -  &
                            &  hf_in_domain(ds)%orb_energies(s-shifts))

                suma = suma + numerator / denominator

            ENDDO
            ENDDO
            ENDDO
            ENDDO

        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        ! sumb
        ! choose 2 occupied orbitals
        DO da = 1, n_domains
        DO a = start_of_domain(da), start_of_domain(da) + electrons_in_domain(da) - 1
        DO db = 1, n_domains
        DO b = start_of_domain(db), start_of_domain(db) + electrons_in_domain(db) - 1

            ! choose 4 virtual orbitals
!           DO dr = 1, n_domains
!           DO r = start_of_domain(dr) + electrons_in_domain(dr), end_of_domain(dr)
!           DO ds = 1, n_domains
!           DO s = start_of_domain(ds) + electrons_in_domain(ds), end_of_domain(ds)
!           DO dt = 1, n_domains
!           DO t = start_of_domain(dt) + electrons_in_domain(dt), end_of_domain(dt)
!           DO du = 1, n_domains
!           DO u = start_of_domain(du) + electrons_in_domain(du), end_of_domain(du)

            ! choose 4 virtual orbitals
            DO dr = 1, n_domains
            DO r = start_of_domain(dr) + electrons_in_domain(dr), &
                 & start_of_domain(dr) + functions_kept(dr) - 1
            DO ds = 1, n_domains
            DO s = start_of_domain(ds) + electrons_in_domain(ds), &
                 & start_of_domain(ds) + functions_kept(ds) - 1
            DO dt = 1, n_domains
            DO t = start_of_domain(dt) + electrons_in_domain(dt), &
                 & start_of_domain(dt) + functions_kept(dt) - 1
            DO du = 1, n_domains
            DO u = start_of_domain(du) + electrons_in_domain(du), &
                 & start_of_domain(du) + functions_kept(du) - 1

                shifta = start_of_domain(da) - 1
                shiftb = start_of_domain(db) - 1
                shiftr = start_of_domain(dr) - 1
                shifts = start_of_domain(ds) - 1
                shiftt = start_of_domain(dt) - 1
                shiftu = start_of_domain(du) - 1

                numerator = mo_double_bar(a,b,r,s) * &
                          & mo_double_bar(r,s,t,u) * &
                          & mo_double_bar(t,u,a,b)

                denominator = (hf_in_domain(da)%orb_energies(a-shifta) +  &
                            &  hf_in_domain(db)%orb_energies(b-shiftb) -  &
                            &  hf_in_domain(dr)%orb_energies(r-shiftr) -  &
                            &  hf_in_domain(ds)%orb_energies(s-shifts)) * &
                            & (hf_in_domain(da)%orb_energies(a-shifta) +  &
                            &  hf_in_domain(db)%orb_energies(b-shiftb) -  &
                            &  hf_in_domain(dt)%orb_energies(t-shiftt) -  &
                            &  hf_in_domain(du)%orb_energies(u-shiftu))

                sumb = sumb + numerator / denominator

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

        ENDDO
        ENDDO
        ENDDO
        ENDDO

        ! sumc
        ! choose 3 occupied orbitals
        DO da = 1, n_domains
        DO a = start_of_domain(da), start_of_domain(da) + electrons_in_domain(da) - 1
        DO db = 1, n_domains
        DO b = start_of_domain(db), start_of_domain(db) + electrons_in_domain(db) - 1
        DO dc = 1, n_domains
        DO c = start_of_domain(dc), start_of_domain(dc) + electrons_in_domain(dc) - 1

!           ! choose 3 virtual orbitals
!           DO dr = 1, n_domains
!           DO r = start_of_domain(dr) + electrons_in_domain(dr), end_of_domain(dr)
!           DO ds = 1, n_domains
!           DO s = start_of_domain(ds) + electrons_in_domain(ds), end_of_domain(ds)
!           DO dt = 1, n_domains
!           DO t = start_of_domain(dt) + electrons_in_domain(dt), end_of_domain(dt)

            ! choose 3 virtual orbitals
            DO dr = 1, n_domains
            DO r = start_of_domain(dr) + electrons_in_domain(dr), &
                 & start_of_domain(dr) + functions_kept(dr) - 1
            DO ds = 1, n_domains
            DO s = start_of_domain(ds) + electrons_in_domain(ds), &
                 & start_of_domain(ds) + functions_kept(ds) - 1
            DO dt = 1, n_domains
            DO t = start_of_domain(dt) + electrons_in_domain(dt), &
                 & start_of_domain(dt) + functions_kept(dt) - 1

                shifta = start_of_domain(da) - 1
                shiftb = start_of_domain(db) - 1
                shiftc = start_of_domain(dc) - 1
                shiftr = start_of_domain(dr) - 1
                shifts = start_of_domain(ds) - 1
                shiftt = start_of_domain(dt) - 1

                numerator = mo_double_bar(a,b,r,s) * &
                          & mo_double_bar(c,s,t,b) * &
                          & mo_double_bar(r,t,a,c)

                denominator = (hf_in_domain(da)%orb_energies(a-shifta) +  &
                            &  hf_in_domain(db)%orb_energies(b-shiftb) -  &
                            &  hf_in_domain(dr)%orb_energies(r-shiftr) -  &
                            &  hf_in_domain(ds)%orb_energies(s-shifts)) * &
                            & (hf_in_domain(da)%orb_energies(a-shifta) +  &
                            &  hf_in_domain(dc)%orb_energies(c-shiftc) -  &
                            &  hf_in_domain(dr)%orb_energies(r-shiftr) -  &
                            &  hf_in_domain(dt)%orb_energies(t-shiftt))

                sumc = sumc + numerator / denominator

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        mp3_correction = sumc + (suma + sumb) / 8.d0

    END SUBROUTINE compute_mp3

    SUBROUTINE print_mp
        IMPLICIT NONE

        WRITE (*,"(6('='), (' '), a, (' '), 38('='))") "MOLLER-PLESSET"
        WRITE (*,*)
        WRITE (*,"(2(' '), a, F22.15)") &
            & "Second order correction : ", mp2_correction
        WRITE (*,"(2(' '), a, F22.15)") &
            & "Third  order correction : ", mp3_correction
        WRITE (*,*)
        WRITE (*,"(2(' '), a, F22.15)") "Total MP2 energy : ", &
            & hf_energy(0) + mp2_correction
        WRITE (*,"(2(' '), a, F22.15)") "Total MP3 energy : ", &
            & hf_energy(0) + mp2_correction + mp3_correction
        WRITE (*,*)

    END SUBROUTINE print_mp

END MODULE
