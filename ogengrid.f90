!==============================================================================
        subroutine ogengrid(y,eta,deta,d2eta)
!==============================================================================
        use stuff
        implicit none
        
        real y(ny), eta(ny), deta(ny), d2eta(ny)
        real ym1, etam1, detam1, d2etam1
        integer j

        real :: drmin, b, rc, rmax, cm, rr
!==============================================================================

        write(*,*) "ogengrid: with yi = ", yi

!.... Note that y here is really the "Blasius" variable eta

        dy = (ymax-ymin)/(ny-1)

!.... make the eta grid

        dy = one / float(ny-1)          ! quick hack really deta
        
        if (yi.ne.zero) then

          aa = ymax * yi / ( ymax - two * yi )
          bb = one + aa / ymax
        
          etam1 = -dy
          ym1 = aa * etam1 / (bb - etam1)
          detam1 = (bb - etam1)**2 / (aa * bb)
          d2etam1 = -two * (bb - etam1)**3 / (aa * bb)**2

          write(30,10) ym1, etam1, detam1, d2etam1

          do j = 1, ny
            eta(j)   = (j-1) * dy
            y(j)     = aa * eta(j) / (bb - eta(j))
            deta(j)  = (bb - eta(j))**2 / (aa * bb)
            d2eta(j) = -two * (bb - eta(j))**3 / (aa * bb)**2
            write(30,10) y(j), eta(j), deta(j), d2eta(j)
          end do

        else                           ! Mahesh map

          rmax = ymax
          write(*,"('Enter dymin, b, yc ==> ',$)")
          read(*,*) drmin, b, rc

          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * &
                 log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) / &
               ( one - (ny-1)*drmin/rmax )
          
          do j = 1, ny
            eta(j)   = real(j-1) * dy
            rr       = eta(j)
            y(j)     = rmax*( cm * rr + log( cosh(b*(rr-rc)) / &
                       cosh(b*(rr+rc)) ) ) /  &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            deta(j)  = one/(rmax*(cm + b*tanh(b*(rr-rc)) - &
                       b*tanh(b*(rr+rc))) / &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ))
            d2eta(j) = -(rmax*(-b**2*(tanh(b*(rr-rc)))**2 + &
                       b**2*(tanh(b*(rr+rc)))**2)/ &
                       (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )) * &
                       deta(j)**3
          end do

!.... Uniform mesh

!         detam1 = one / (ymax-ymin)
!         do j = 1, ny
!           eta(j)   = (j-1) * dy
!           y(j)     = (j-1) * (ymax-ymin)/float(ny-1)
!           deta(j)  = one / (ymax-ymin)
!           d2eta(j) = zero
!         end do

        end if

        return
10      format(4(1pe13.6,1x))
        end
