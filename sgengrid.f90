!==============================================================================
        subroutine sgengrid(y,eta,deta,d2eta)
!==============================================================================
        use stuff
        implicit none
        
        real :: y(ny), eta(ny), deta(ny), d2eta(ny), th, dth
        real :: Lmap
        integer :: i
        
        real, parameter :: three = 3.0

        real :: xi, dxi, drmin, b, rc, rmax, cm, rr
!==============================================================================
        dth = pi/float(ny-1)

!.... Make mesh in transformed, eta, and Chebyshev space, th

        do i = 1, ny
          th = float(i-1)*dth
          eta(i) = cos(th)
        end do

!.... make metric transformations

        if (Yi .ne. zero) then   ! use algebraic mapping

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
            end do
          else                   ! Craig Streett's mapping
            do i = 1, ny
              y(i)     = ymax*yi*(one+eta(i)) / (one+two*yi-eta(i)) 
              deta(i)  = (two*yi+one-eta(i))**2 / (two*ymax*yi*(yi+one))
              d2eta(i) = -pt5*(two*yi+one-eta(i))**3 / (ymax*yi*(yi+one))**2
            end do
          end if

        else                     ! Hyperbolic tangent mapping

          write(*,"('Enter ymax ==> ',$)") 
          read(*,*) rmax
          write(*,"('Enter dymin, b, yc ==> ',$)")
          read(*,*) drmin, b, rc

          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * &
                 log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) / &
               ( one - (ny-1)*drmin/rmax )
          
          dxi = one / real(ny-1)
  
          do i = 1, ny
            xi       = one - real(i-1) * dxi
            rr       = real(i-1) * dxi
            y(i)     = rmax*( cm * rr + log( cosh(b*(rr-rc)) / &
                       cosh(b*(rr+rc)) ) ) /  &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            deta(i)  = one/(rmax*(cm + b*tanh(b*(rr-rc)) - &
                       b*tanh(b*(rr+rc))) / &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ))
            d2eta(i) = -(rmax*(-b**2*(tanh(b*(rr-rc)))**2 + &
                       b**2*(tanh(b*(rr+rc)))**2)/ &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )) * &
                       deta(i)**3
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
            xi       = real(i-1) * dxi
            d2eta(i) = -pi**2 * cos(pi * xi) * (deta(i))**2 - &
                        pi * sin(pi * xi) * d2eta(i)
            deta(i)  = -pi * sin(pi * xi) * deta(i)
          end do

          open(10,file='metric.out')
          do i = 1, ny
            if (i.eq.1) then
              write(10,10) y(i), eta(i), deta(i), &
                            (eta(i+1)-eta(i))/(y(i+1)-y(i)), &
                            d2eta(i), &
                            (deta(i+1)-deta(i))/(y(i+1)-y(i))
            else if (i.eq.ny) then
              write(10,10) y(i), eta(i), deta(i), &
                            (eta(i)-eta(i-1))/(y(i)-y(i-1)), &
                            d2eta(i), &
                            (deta(i)-deta(i-1))/(y(i)-y(i-1))
            else
              write(10,10) y(i), eta(i), deta(i), &
                            (eta(i+1)-eta(i-1))/(y(i+1)-y(i-1)), &
                            d2eta(i), &
                            (deta(i+1)-deta(i-1))/(y(i+1)-y(i-1))
            end if
          end do
          close(10)
          
        end if      ! mapping type

        return
 10     format(8(1x,1pe13.6))
        end
