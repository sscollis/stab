!=============================================================================!
        module stuff
!
!  Variables for STAB
!
!  Revised:  9-2-96
!=============================================================================!

!.... Run parameters

          real    :: Ma, Re, Pr, T0
          complex :: alpha, beta, omega
          real    :: alphar, alphai, betar, betai, omegar, omegai
          integer :: itype = 1, ievec = 1
          logical :: ider = .true., Navier = .true.
          integer :: top = 0, wall = 0, wallt = 0
          integer :: curve = 0

!.... parameters for receptivity calculations

          real    :: dalpha_r, dalpha_i, lambda = 0.0
          complex :: dalpha

!.... Edge conditions

          real :: Te, rmue, dmue, d2mue, rlme, dlme, d2lme, cone
          real :: dcone, d2cone

!.... Useful constants

          real, parameter :: zero    = 0.0000000000000000000d+0
          real, parameter :: pt25    = 2.5000000000000000000d-1
          real, parameter :: pt33    = 3.3333333333333333333d-1
          real, parameter :: pt5     = 5.0000000000000000000d-1
          real, parameter :: pt66    = 6.6666666666666666666d-1
          real, parameter :: one     = 1.0000000000000000000d+0
          real, parameter :: onept25 = 1.2500000000000000000d+0
          real, parameter :: onept5  = 1.5000000000000000000d+0
          real, parameter :: two     = 2.0000000000000000000d+0
          real, parameter :: four    = 4.0000000000000000000d+0
          real, parameter :: pi      = 3.1415926535897932385d+0
          complex, parameter :: im   = (0.0,1.0)

!.... fluid properties

          real :: gamma = 1.4, gamma1 = 0.4, cv = 716.5, cp = 1003.1
          real :: Rgas  = 286.6

          integer :: mattyp
          real    :: datmat(3)
        
!.... problem dimensions
          
          integer :: ny = 1, nx = 1, nz = 1, nsd = 3, ndof = 5
          real    :: x, ymin, ymax, dy, aa, bb, yi
                  
        end module stuff

!=============================================================================!
        module material
!       
!  Interfaces for material routines
!
!  Revised:  4-16-96
!=============================================================================!
          interface getmat
            subroutine getmat(t, rmu, rlm, con, &
                              dmu, d2mu, dlm, d2lm, dcon, d2con)
              real :: t(:), rmu(:), rlm(:), con(:), dmu(:), d2mu(:)
              real :: dlm(:), d2lm(:), dcon(:), d2con(:)
            end
            subroutine sgetmat(t, rmu, rlm, con, &
                                dmu, d2mu, dlm, d2lm, dcon, d2con)
              real :: t, rmu, rlm, con, dmu, d2mu, dlm, d2lm
              real :: dcon, d2con
            end
          end interface
        
        end module material
