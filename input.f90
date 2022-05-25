!==============================================================================
        subroutine input
!==============================================================================
        use stuff
        use material
        implicit none

        integer :: itmp
!==============================================================================

#define VERBOSE

!.... Constant mu or Sutherland's law

 1      write (*,"(/,'(0) for constant Mu, (1) for Sutherland ==> ',$)")
        read (*,*) mattyp
        if (mattyp .ne. 0 .and. mattyp .ne. 1) goto 1

        if (mattyp .eq. 1) then         ! for Sutherland's law (AIR)

!.... get the freestream stagnation temperature

          write (*,"('Enter the freestream T0 (K) ==> ',$)")
          read (*,*) T0
          te    = t0 / ( one + pt5 * gamma1 * Ma**2 )

!.... set fluid properties

          datmat(1) = 1.715336725523065e-05     ! 1.716e-5
          datmat(2) = 273.0
          datmat(3) = 110.4                     ! 111.0
        else
          datmat(1) = one
          datmat(2) = zero
          datmat(3) = zero
        end if

        write (*,"('Enter Ma, Re, Pr ==> ',$)")
        read (*,*) Ma, Re, Pr

!.... get the fluid properties at the reference state

        call getmat(te,   rmue, rlme, cone, dmue, d2mue, &
                    dlme, d2lme, dcone, d2cone)

        write (*,"('Enter ny, Yi, Ymax ==> ',$)")
        read (*,*) ny, yi, ymax
        write (*,"('Compute eigenfuntions (0/1) ==> ',$)")
        read (*,*) ievec
        write (*,"('Use internal derivatives (0/1) ==> ',$)")
        read (*,*) itmp
        if (itmp.eq.0) then
          ider = .false.
        else
          ider = .true.
        end if
        write(*,"('top, wall, wallt, curve ==> ',$)")
        read (*,*) top, wall, wallt, curve

        ymin = zero

 10     continue
#ifdef VERBOSE
        write (*,"('(1) Temporal,(2) Spatial,(3) FD_temporal,',/, &
      &  '(4) FD_spatial (5) Stokes,(6) Bump,',/, &
      &  '(7) Multi Temporal (8) Multi Spatial (9) Multi Bump ==> ')")
#else
        write (*,"('Enter itype ==> ',$)")
#endif
        read (*,*) itype

        if (itype .eq. 1 .or. itype .eq. 3) then
          write (*,"('Enter alphar, alphai ==> ',$)")
          read (*,*) alphar, alphai
          write (*,"('Enter betar, betai ==> ',$)")
          read (*,*) betar, betai
        elseif (itype .eq. 2 .or. itype .eq. 4) then
          write (*,"('Enter omegar, omegai ==> ',$)")
          read (*,*) omegar, omegai
          write (*,"('Enter betar, betai ==> ',$)")
          read (*,*) betar, betai
        elseif (itype .eq. 5) then
          write (*,"('Enter omegar ==> ',$)")
          read (*,*) omegar
          omegai = zero
          write (*,"('Enter alphar, alphai ==> ',$)")
          read (*,*) alphar, alphai
!         alphar = Ma * omegar / ( one + Ma )
!         if (Re .ne. zero) then
!           alphai = Ma * alphar**2 * pt5 / (Re * (one + Ma)) * &
!                    ( 4.0/3.0 + gamma1 / Pr )
!         else
!           alphai = zero
!         end if
          write (*,"('Enter betar, betai ==> ',$)")
          read (*,*) betar, betai
        elseif (itype .eq. 6) then
          omegar = zero
          omegai = zero
          write (*,"('Enter alphar, alphai ==> ',$)")
          read (*,*) alphar, alphai
          write (*,"('Enter betar, betai ==> ',$)")
          read (*,*) betar, betai
        elseif (itype .eq. 7) then
          alphar = zero
          alphai = zero
          betar  = zero
          betai  = zero
        elseif (itype .eq. 8) then
          omegar = zero
          omegai = zero
          alphar = zero
          alphai = zero
          betar  = zero
          betai  = zero
        elseif (itype .eq. 9) then
          omegar = zero
          omegai = zero
          alphar = zero         ! read from a file
          alphai = zero
          write (*,"('Enter betar, betai ==> ',$)")
          read (*,*) betar, betai
        else
          goto 10
        end if

!.... set alpha and beta and omega

        alpha = alphar + im * alphai
        beta  = betar  + im * betai
        omega = omegar + im * omegai

!.... Echo the run parameters (it would be best to write to a file)

!       write (*,*)
!       write (*,*) 'R U N   P A R A M E T E R S'
!       write (*,*)
!       write (*,*) 'Mach   = ',Ma
!       write (*,*) 'Re     = ',Re
!       write (*,*) 'Pr     = ',Pr
!       write (*,*) 'T0     = ',T0
!       write (*,*) 'Te     = ',Te
!       write (*,*) 'mu     = ',rmue
!       write (*,*) 'lambda = ',rlme
!       write (*,*) 'con    = ',cone
!       write (*,*) 'alpha  = ',alpha
!       write (*,*) 'beta   = ',beta
!       write (*,*) 'omega  = ',omega
!       write (*,*)

        return
        end
