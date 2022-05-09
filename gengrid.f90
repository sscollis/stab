        subroutine gengrid(y,eta,deta,d2eta)
!==============================================================================
        use stuff
        implicit none
        
        real y(ny), eta(ny), deta(ny), d2eta(ny)
        integer i
        
        real, parameter :: three = 3.0

        real rd1, rd2, dd, rmax
!==============================================================================
        
        dy   = one / float(ny-1)        ! quick hack really deta

!.... finite domain with Hyperbolic tangent mapping

        rd1 = 0.0005d0          ! set for R=2400 MSE, r=50
        rd2 = 5.0d0
        dd  = 5.36966703089523d0
        rmax = 40.0d0
          
        do i = 1, ny
          eta(i) =  real(i-1) * dy
          y(i)  =  rmax*(pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 + &
            &        eta(i)))*pt5)/(Sqrt(rd2/rd1) + &
            &        (one - Sqrt(rd2/rd1))* &
            &        (pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 + &
            &         eta(i)))*pt5))
          deta(i) = one/(dd*Sqrt(rd2/rd1)*rmax*Cosh(dd*pt5)* &
            &              one/Cosh(dd*(-pt5 + eta(i)))* &
            &              (Sinh(dd*(one - eta(i))) + Sinh(dd*eta(i)))/ &
            &      (Sqrt(rd2/rd1)*Sinh(dd*(one - eta(i))) + &
            &              Sinh(dd*eta(i)))**2)
          d2eta(i) = -(dd**2*Sqrt(rd2/rd1)*rmax* &
            &               (-Cosh(dd*(one - three*eta(i))) + &
            &               Sqrt(rd2/rd1)*Cosh(dd*(two - three*eta(i))) - &
            &               Cosh(dd*(one - eta(i))) + &
            &               two*Sqrt(rd2/rd1)*Cosh(dd*(one - eta(i))) - &
            &               two*Cosh(dd*eta(i)) + &
            &               Sqrt(rd2/rd1)*Cosh(dd*eta(i)))* &
            &               (one/Cosh(dd*(-pt5 + eta(i))))**2*(Sinh(dd))/ &
            &       (two*(Sqrt(rd2/rd1)*(Sinh(dd*(one - eta(i)))) + & 
            &           (Sinh(dd*eta(i))))**3)) * deta(i)**3
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
          
        return
 10     format(8(1x,1pe13.6))
        end
