!==============================================================================
        subroutine mbump(ind)
!==============================================================================
!
!       Driver for multiple bump receptivity solver
!
!       S. Scott Collis
!
!       Revised: 6-14-97
!==============================================================================
        use stuff
        implicit none

        real :: tmp
        integer :: i, nbody, ind, n
        real, allocatable :: xb(:), delta(:), theta(:)
        character*80 :: base="bump ", name

        integer :: wallbc

        integer :: iparm=11, idelta=10
!==============================================================================
        write(*,"('Enter dalpha_r, dalpha_i ==> ',$)") 
        read(*,*) dalpha_r, dalpha_i
        dalpha = cmplx(dalpha_r,dalpha_i)

        write(*,"('Enter wallbc (0,1) ==> ',$)") 
        read(*,*) wallbc
        if (wallbc .lt. 0 .or. wallbc .gt. 1) then
          write(*,"('Illegal value for wallbc: ',i4)") wallbc
          call exit(1)
        end if

!.... lambda should be the local sweep angle, but I'm currently not using it

!       lambda = zero
!       lambda = lambda * pi / 180.0

!.... read the station information
!.... WARNING:  I used to use delta.dat but now use stat.dat to get local sweep

        open(iparm,file='parm.dat',status='old',err=200)
        read(iparm,*,err=200) n

        open(idelta,file='stat.dat',status='old',err=100)
        nbody = 0
 20     continue
        read(idelta,*,end=30) tmp
        nbody = nbody + 1
        goto 20
 30     continue
        rewind(idelta)
        allocate( xb(nbody), delta(nbody), theta(nbody) )
        do i = 1, nbody
          read(idelta,*) xb(i), delta(i), theta(i)
        end do
        close(idelta)

!.... open the bump output file

        if (wallbc.eq.0) then
          open(50,file='bump.out',status='unknown')
          base = "bump "
        else
          open(50,file='blow.out',status='unknown')
          base = "blow "
        end if
        
        write(*,5)
        do i = 1, n
          read(iparm,*,err=200) ind, x, alphar, alphai

!         read(iparm,*,err=200) ind, betar, alphar, alphai
!         x = xb(ind)
!         beta = cmplx(betar,zero)

          if (ind.lt.1 .or. ind.gt.nbody) then
            write(*,"('ERROR: illegal index...check delta.dat')")
            call exit(1)
          end if
          alpha = cmplx(alphar, alphai)
          yi = two * delta(ind)
          lambda = theta(ind)
          call makename(base, ind, name)
          write(*,10) ind, x, yi, lambda*180/pi, &
                      real(omega), aimag(omega), &
                      real(beta), aimag(beta)
          call bump(name, ind, wallbc)
        end do

        close(50)
        close(iparm)

        deallocate( xb, delta )
        return

5       format(/,110('='),/, &
          'index', '      s       ', '      Yi      ', &
          '    theta_e   ', '    omega_r   ', '    omega_i   ', &
          '    beta_r    ', '    beta_i    ', /,110('='))
10      format(1x,i3,1x,7(1pe13.6,1x))

100     write(*,"('ERROR: stat.dat not found')")
        call exit(1)
200     write(*,"('ERROR: parm.dat not found')")
        call exit(1)

        end
