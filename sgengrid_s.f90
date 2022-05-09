!==============================================================================
        subroutine sgengrid_s(y,eta,deta,d2eta,y_s,eta_s,deta_s,d2eta_s)
!==============================================================================
        use stuff
        implicit none
        
        real y(ny), eta(ny), deta(ny), d2eta(ny), th, dth
        real :: eta_s(ny-1), y_s(ny-1), deta_s(ny-1), d2eta_s(ny-1)

        real :: Lmap
        integer i
        
        real, parameter :: three = 3.0

        real rd1, rd2, dd, rmax, xi, dxi, cm, b, rc, drmin
!==============================================================================
        dth = pi/float(ny-1)

!.... Make mesh in transformed, eta, and Chebyshev space, th

        do i = 1, ny
          th = float(i-1)*dth
          eta(i) = cos(th)
        end do

!.... Also make the staggered mesh

        do i = 1, ny-1
          th = (float(i-1)+pt5)*dth
          eta_s(i) = cos(th)
        end do

!.... make metric transformations

        if (Yi .ne. zero) then   ! use algebraic mapping

!.... infinite domain with algebraic mapping

          if (ymax .eq. 0) then  ! infinite domain with algebraic mapping
            Lmap = Yi            ! for consistency with input
            do i = 1, ny
              deta(i)  = (eta(i)-one)**2/(two*Lmap)
              d2eta(i) = (eta(i)-one)**3/(two*Lmap**2)
              if (eta(i) .ne. one) then
                y(i) = Lmap*(one+eta(i))/(one-eta(i))
              else
                y(i) = 1.0e99    ! approximate infinity
              end if
              if (i.lt.ny) then  ! staggered grid
                deta_s(i)  = (eta_s(i)-one)**2/(two*Lmap)
                d2eta_s(i) = (eta_s(i)-one)**3/(two*Lmap**2)
                y_s(i) = Lmap*(one+eta_s(i))/(one-eta_s(i))
                write(98,10) y_s(i), eta_s(i), deta_s(i), d2eta_s(i)
              end if
              write(99,10) y(i), eta(i), deta(i), d2eta(i)
            end do
          else                   ! Craig Streett's mapping
            do i = 1, ny
              y(i)     = ymax*yi*(one+eta(i)) / (one+two*yi-eta(i)) 
              deta(i)  = (two*yi+one-eta(i))**2 / (two*ymax*yi*(yi+one))
              d2eta(i) = -pt5*(two*yi+one-eta(i))**3 / (ymax*yi*(yi+one))**2
              if (i.lt.ny) then  ! staggered grid
                y_s(i)     = ymax*yi*(one+eta(i)) / (one+two*yi-eta(i)) 
                deta_s(i)  = (two*yi+one-eta(i))**2 / (two*ymax*yi*(yi+one))
                d2eta_s(i) =-pt5*(two*yi+one-eta(i))**3 / (ymax*yi*(yi+one))**2
              end if
            end do
          end if

        else                     ! Hyperbolic tangent mapping

          write(*,"('Enter rmax, drmin, b, rc ==> ',$)") 
          read(*,*) rmax, drmin, b, rc
          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * &
     &           log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) / &
     &         ( one - (ny-1)*drmin/rmax )
          dxi = one / real(ny-1)
          do i = 1, ny
            xi      = one - real(i-1) * dxi
            y(i)    = rmax*( cm * xi + log( cosh(b*(xi-rc)) / &
     &                  cosh(b*(xi+rc)) ) ) / &
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            deta(i) = one/(rmax*(cm + b*tanh(b*(xi-rc)) - b*tanh(b*(xi+rc))) /&
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) )
            d2eta(i) = -(rmax*(-b**2*(tanh(b*(xi-rc)))**2 + &
     &                 b**2*tanh(b*(xi+rc)))/ &
     &                 (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) ) * &
     &                 deta(i)**3
          end do


!.... diagnostic

!         do i = 1, ny
!           xi = real(i-1) * dxi
!           if (i.eq.1) then
!             write(69,10) y(i), xi, deta(i), &
!                          (real(i)-real(i-1))*dxi/(y(i+1)-y(i))
!           else if (i.eq.ny) then
!             write(69,10) y(i), xi, deta(i), &
!                          (real(i-1)-real(i-2))*dxi/(y(i)-y(i-1))
!           else
!             write(69,10) y(i), xi, deta(i), &
!                          (real(i)-real(i-2))*dxi/(y(i+1)-y(i-1))
!           end if
!         end do

!.... correct for the xi -> eta mapping

          do i = 1, ny
            xi = one - real(i-1) * dxi
            d2eta(i) =  pi**2 * cos(pi * xi) * (deta(i))**2 + &
                        pi * sin(pi * xi) * d2eta(i)
            deta(i) =   pi * sin(pi * xi) * deta(i)
          end do

!         open(10,file='metric.out')
!         do i = 1, ny
!           if (i.eq.1) then
!             write(10,10) y(i), eta(i), deta(i), &
!                           (eta(i+1)-eta(i))/(y(i+1)-y(i)), &
!                           d2eta(i), &
!                           (deta(i+1)-deta(i))/(y(i+1)-y(i))
!           else if (i.eq.ny) then
!             write(10,10) y(i), eta(i), deta(i), &
!                           (eta(i)-eta(i-1))/(y(i)-y(i-1)), &
!                           d2eta(i), &
!                           (deta(i)-deta(i-1))/(y(i)-y(i-1))
!           else
!             write(10,10) y(i), eta(i), deta(i), &
!                           (eta(i+1)-eta(i-1))/(y(i+1)-y(i-1)), &
!                           d2eta(i), &
!                           (deta(i+1)-deta(i-1))/(y(i+1)-y(i-1))
!           end if
!         end do
!         close(10)

!.... now do the staggered mesh

          do i = 1, ny-1
            xi    =  one - (real(i-1)+pt5) * dxi
            y_s(i)    = rmax*( cm * xi + log( cosh(b*(xi-rc)) / &
     &                  cosh(b*(xi+rc)) ) ) / &
     &                  (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            deta_s(i) = one/(rmax*(cm + b*tanh(b*(xi-rc)) - &
     &                  b*tanh(b*(xi+rc))) /&
     &                  (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) )
            d2eta_s(i) = -(rmax*(-b**2*(tanh(b*(xi-rc)))**2 + &
     &                 b**2*tanh(b*(xi+rc)))/ &
     &                 (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) ) * &
     &                 deta(i)**3
          end do

!.... correct for the xi -> eta mapping

          do i = 1, ny
            xi = one - (real(i-1)+pt5) * dxi
            d2eta_s(i) =  pi**2 * cos(pi * xi) * (deta_s(i))**2 + &
                        pi * sin(pi * xi) * d2eta_s(i)
            deta_s(i) =   pi * sin(pi * xi) * deta_s(i)
          end do

        end if

        return
 10     format(8(1x,1pe13.6))
        end

