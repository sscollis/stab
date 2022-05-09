!==============================================================================
        program stab 
!  
!  Compressible linear stability solver with surface curvature
!  which uses either Chebyshev collocation or fourth-order finite 
!  difference to solve the stability equations. Both the spatial and 
!  temporal problems are supported as well as parallel flow Finite-Reynolds 
!  number receptivity theory for surface roughness (bump.f90)
!
!  For the Chebyschev collocation version:
!
!  A mapping parameter of Yi=Lmap=0.05 seems to work well for the R=1000, 
!  BL with delta*=0.1 which is the sample profile in the TStest directory.
!  
!  Both IMSL and LAPACK routines can be used although LAPACK is currently
!  turned on.
!
!  Use the -DCRAY flag to tell the compiler that the code is on a CRAY
!
!  To include nonparallel effects, see the code "shoot"
!
!  Author:     S. Scott Collis
!
!  Revised:    6-13-1997
!  Revised:    1-17-2020
!
!  Copyright:  S. Scott Collis
!              Department of Mechanical Engineering and Materials Science
!              Rice University, MS 321
!              Houston, TX 77005-1892
!              (713) 527-8101 x3617
!==============================================================================
        use stuff
        implicit none

        character*80 :: name = 'evec.dat'
        integer :: ind = 0
!==============================================================================
        
!.... input parameters

        call input

!.... select which version to use

        if (itype .eq. 1) then
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          call temporal(name, ind)      ! Chebyschev parallel eigensolver
        elseif (itype .eq. 2) then
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          write(*,"('Enter x ==> ',$)")
          read(*,*) x
          call spatial(name, ind)       ! Chebyschev parallel eigensolver
        elseif (itype .eq. 3) then
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          itype = 1
          call fd_temporal(name, ind)   ! FD parallel eigensolver
        elseif (itype .eq. 4) then
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          write(*,"('Enter x ==> ',$)")
          read(*,*) x
          itype = 2
          call fd_spatial(name, ind)    ! FD parallel eigensolver
        elseif (itype .eq. 5) then
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          call stokes(ind)              ! Stokes wave solver
        elseif (itype .eq. 6) then
          name = 'bump.dat'
          write(*,"('Enter index ==> ',$)")
          read(*,*) ind
          write(*,"('Enter x ==> ',$)")
          read(*,*) x
          write(*,"('Enter dalpha_r, dalpha_i ==> ',$)") 
          read(*,*) dalpha_r, dalpha_i
          dalpha = cmplx(dalpha_r,dalpha_i)
          lambda = zero                 ! Could set to the local sweep angle
          lambda = lambda * pi / 180.0
          call bump(name, ind, 0)       ! Roughness receptivity solver
        elseif (itype .eq. 7) then
          itype = 1
          call mtemporal(ind)           ! Multiple temporal driver
        elseif (itype .eq. 8) then
          itype = 2
          call mspatial(ind)            ! Multiple spatial driver
        elseif (itype .eq. 9) then
          call mbump(ind)               ! Multiple bump driver
        end if
        
        stop
        end
