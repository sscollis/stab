!=============================================================================!
        subroutine mgengrid(y, eta, deta, d2eta, vm)
!=============================================================================!
        use stuff
        implicit none
        
        real vm(ny,ndof,nx), y(ny), eta(ny), deta(ny), d2eta(ny)
        integer j
        real scale
        
        real :: rmax, drmin, b, rc, cm, rr(ny), drr(ny), d2rr(ny), dr
!=============================================================================!

!.... initialize the mean flow

        vm = zero

!.... read in Ted's mean mixing layer profile

!       j = 1
!       open(10,file='profile.dat')
!20     continue
!         read(10,*,end=30) y(j), vm(j,1,1), vm(j,2,1), vm(j,5,1), &
!                            deta(j), d2eta(j)
!         j = j + 1
!         goto 20
!30     continue
        
!       if ( j-1 .ne. ny ) then
!         write(*,*) 'ERROR:  Ny is wrong in mgengrid'
!         call exit(1)
!       end if

!.... renormalize Ted's temperature to t_\infty

!       vm(:,5,1) = vm(:,5,1) * (gamma1 * Ma**2)

!.... make the eta grid

!       dy = one / float(ny-1)      ! quick hack really deta

!.... the scale factor is just the number of points in the truncated
!.... domain minus one.

!       write(*,"('Enter scale ==> ',$)") 
!       read(*,*) scale
        
!       do j = 1, ny
!         eta(j)   = (j-1) * dy
          
!         deta(j) = deta(j) * 60.0 * (scale/250.0)
!         d2eta(j) = d2eta(j) * 60.0 * (scale/250.0)**2
          
!         d2eta(j) = -d2eta(j) * deta(j)**(-3)
!         deta(j)  = one / deta(j)
!         write(30,10) y(j), eta(j), deta(j), d2eta(j)
!       end do

!.... new mapping

        dy = one / real(ny-1)
        write(*,"('Enter rmax, drmin, b, rc ==> ',$)") 
        read(*,*) rmax, drmin, b, rc
        cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * &
             log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) / &
             ( one - (ny-1)*drmin/rmax )
        do j = 1, ny
          eta(j)  = dble(j-1) * dy
          y(j)   = rmax*( cm * (eta(j)-pt5)/pt5 + &
                    log( cosh(b*((eta(j)-pt5)/pt5-rc)) / &
                    cosh(b*((eta(j)-pt5)/pt5+rc)) ) ) / &
                    (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          deta(j)  = rmax*(2*cm - 2*b*(Cosh(b*(-1 + 2*eta(j) - rc)))**(-1)* &
                    (Cosh(b*(-1 + 2*eta(j) + rc)))**(-1)*Sinh(2*b*rc)) / &
                    (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          d2eta(j) = -4*b**2*rmax*Cosh(b*(-1 + 2*eta(j) - rc))**(-2)* &
                    Cosh(b*(-1 + 2*eta(j) + rc))**(-2)*Sinh(2*b - 4*b*eta(j))* &
                    Sinh(2*b*rc) / &
                    (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          d2eta(j) = -d2eta(j) * deta(j)**(-3)
          deta(j)  = one / deta(j)
          write(20,10) y(j), eta(j), deta(j), d2eta(j)
        end do

        return
10      format(4(1pe13.6,1x))
        end
