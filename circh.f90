!==============================================================================
        subroutine circh( s, n, r, h, dhds, dhdr, dhdsr, dhdrr )
!==============================================================================
!
!       Compute the metrics for a curved wall.  
!       Set the functions xcloc, ycloc and carc for the 
!       geometry under consideration.
!
!       S. Scott Collis
!
!       Revised: 9-23-96
!
!==============================================================================
        implicit none

        integer :: n, i
        real    :: s, r(n), h(n), dhds(n), dhdr(n), dhdsr(n), dhdrr(n)

        real, external :: xcloc, ycloc

        real :: xl, yl, th, bn1, bn2, dxdy, dydx, dxbds, dybds, dx, dy, &
                ddxdx, ddydx, ddxdy, ddydy, d2dxdx2, d2dydx2, d2dxdy2,  &
                d2dydy2, dbn1, dbn2, d2xdy2, d2ydx2, d2xbds2, d2ybds2,  &
                d2bn1, d2bn2

        real :: a(n), b(n), dads(n), dbds(n)
        
        real, parameter :: zero = 0.0, pt5 = 0.5, one = 1.0, onept5 = 1.5, &
                           two = 2.0,  twopt5 = 2.5, three = 3.0, &
                           infty = 1.0e30

        real :: radius
        common /circstuff/ radius
!==============================================================================
        if (s.eq.-one) then
          write(*,*) 'WARNING:  Curvature is turned off!'
          h     = one
          dhdr  = zero
          dhds  = zero
          dhdrr = zero
          dhdsr = zero
          return
        end if

!       write(*,"('Enter the radius ==> ',$)")
!       read(*,*) radius
        radius = s
        
        s = zero

        xl  = xcloc( zero, s )
        yl  = ycloc( xl )
        th  = atan2( -xl, sqrt( radius**2 - xl**2 ) ) 
        bn1 = -sin(th)
        bn2 =  cos(th)

!       write(*,"(8(e13.6,1x))") s, xl, yl, th, bn1, bn2
        
        if (xl .eq. zero) then
          dydx = -infty
        else
          dydx = (-xl)/yl
        end if
        dxbds  = one / sqrt(one + dydx**2)

        if (xl .eq. zero) then
          dxdy = -infty
        else
          dxdy = yl/(-xl)
        end if
        if (xl .le. zero) then
          dybds = one / sqrt( dxdy**2 + one )
        else
          dybds = -one / sqrt( dxdy**2 + one )
        end if

        if (yl .eq. zero) then
          dx = zero
          dy = -xl
          ddxdx = -infty
          ddydx = -one
          ddxdy = one
          ddydy = zero
          d2dxdx2 = -infty
          d2dydx2 = zero
          d2dxdy2 = zero
          d2dydy2 = (xl**2 + yl**2)/xl**3
        else if (xl .eq. zero) then
          dx = yl
          dy = zero
          ddxdx = zero
          ddydx = -one
          ddxdy = one
          ddydy = infty
          d2dxdx2 = -(yl**2 + xl**2)/yl**3
          d2dydx2 = zero
          d2dxdy2 = zero
          d2dydy2 = infty
        else
          dx = yl
          dy = -xl
          ddxdx = -xl / yl
          ddydx = -one
          ddxdy = one
          ddydy = yl / xl
          d2dxdx2 = -(yl**2 + xl**2)/yl**3
          d2dydx2 = zero
          d2dxdy2 = zero
          d2dydy2 = (xl**2 + yl**2)/xl**3
        end if
        
        if ( abs(bn1) .gt. abs(bn2) ) then
          dbn1 = ( -ddydy/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdy + &
                    two*dy*ddydy)/(dx**2 + dy**2)**onept5 ) * dybds
          dbn2 = ( ddxdy/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdy + &
                    two*dy*ddydy)/(dx**2 + dy**2)**onept5 ) * dybds
        else
          dbn1  = ( -ddydx/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdx + &
                    two*dy*ddydx)/(dx**2 + dy**2)**onept5 ) * dxbds
          dbn2  = ( ddxdx/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdx + &
                    two*dy*ddydx)/(dx**2 + dy**2)**onept5 ) * dxbds
        end if

        if ( abs(bn1) .gt. abs(bn2) ) then
          d2xdy2  = -(xl**2+yl**2)/xl**3
          if (xl.le.0) then
            d2ybds2 = -(one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
          else
            d2ybds2 = (one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
          end if
          d2xbds2 = d2xdy2*(dybds)**2 + dxdy*d2ybds2
        else
          d2ydx2  = -(xl**2+yl**2)/yl**3
          d2xbds2 = -(one + dydx**2)**(-onept5) * dydx * d2ydx2 * dxbds
          d2ybds2 = d2ydx2*(dxbds)**2 + dydx*d2xbds2
        end if    

        if ( abs(bn1) .gt. abs(bn2) ) then
          d2bn1 = ((ddxdy*(dy*ddxdy-dx*ddydy)+dx*(dy*d2dxdy2-dx*d2dydy2))/ &
                  (dx**2+dy**2)**(onept5) - &
                  (three*dx*(dy*ddxdy-dx*ddydy)*(dx*ddxdy+dy*ddydy))/ &
                  (dx**2+dy**2)**(twopt5))*(dybds)**2 + &
                  (-ddydy/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdy +  &
                  two*dy*ddydy)/(dx**2 + dy**2)**onept5) * d2ybds2
          d2bn2 = ((ddydy*(dy*ddxdy-dx*ddydy)+dy*(dy*d2dxdy2-dx*d2dydy2))/ &
                  (dx**2+dy**2)**(onept5) - &
                  (three*dy*(dy*ddxdy-dx*ddydy)*(dx*ddxdy+dy*ddydy))/ &
                  (dx**2+dy**2)**(twopt5))*(dybds)**2 + &
                  (ddxdy/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdy +  &
                  two*dy*ddydy)/(dx**2 + dy**2)**onept5) * d2ybds2
        else
          d2bn1 = ((ddxdx*(dy*ddxdx-dx*ddydx)+dx*(dy*d2dxdx2-dx*d2dydx2))/ &
                  (dx**2+dy**2)**(onept5) - &
                  (three*dx*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/ &
                  (dx**2+dy**2)**(twopt5))*(dxbds)**2 + &
                  (-ddydx/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdx +  &
                  two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
          d2bn2 = ((ddydx*(dy*ddxdx-dx*ddydx)+dy*(dy*d2dxdx2-dx*d2dydx2))/ &
                  (dx**2+dy**2)**(onept5) - &
                  (three*dy*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/ &
                  (dx**2+dy**2)**(twopt5))*(dxbds)**2 + &
                  (ddxdx/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdx +  &
                  two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
        end if

        if (xl .eq. zero) d2bn1 = zero

!       write (*,"(6(1pe13.6,1x))") bn1, bn2, dbn1, dbn2, d2bn1, d2bn2
!       write (*,"(6(1pe13.6,1x))") dxbds, dybds, d2xbds2, d2ybds2

!.... now form the actual metric and metric derivatives

        a = dxbds + r * dbn1
        b = dybds + r * dbn2
        
        h = sqrt( a**2 + b**2 )
        
        dads = d2xbds2 + r * d2bn1
        dbds = d2ybds2 + r * d2bn2

        dhds = ( a * dads + b * dbds ) / h

        dhdr = ( a * dbn1 + b * dbn2 ) / h
        
        dhdrr = ( -dhdr**2 + dbn1**2 + dbn2**2 ) / h
        
        dhdsr = -dhds / h**2 * ( a * dbn1 + b * dbn2 ) + &
                ( dads * dbn1 + a * d2bn1 + dbds * dbn2 + b * d2bn2 ) / h

!       do i = 1, n
!         write(70,10) r(i), h(i), dhds(i), dhdr(i), dhdsr(i), dhdrr(i)
!         write(71,10) r(i), h(i), a(i), b(i), dads(i), dbds(i), dhds(i)
!       end do
        
 10     format( 8(1pe13.6,1x) )
 
        return
        end

!=============================================================================c
        function ycloc(x)
!=============================================================================c
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /circstuff/ radius
!=============================================================================c
        ycloc = sqrt( radius**2 - x**2 )

        return
        end

!=============================================================================c
        function xcloc(x,ds)
!=============================================================================c
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /circstuff/ radius
        common /distance/ darc, x1
!=============================================================================c
        darc = ds
        x1   = x

        th1 = atan2( sqrt( radius**2 - x1**2 ), x1 )
        th  = th1 - ds
        xcloc = radius * cos(th)

        return
        end 

!=============================================================================c
        function carc(x1,x2)
!=============================================================================c
        common /circstuff/ radius
!=============================================================================c
        th1 = atan2( sqrt( radius**2 - x1**2 ), x1 )
        th2 = atan2( sqrt( radius**2 - x2**2 ), x2 )

        carc = -(th2 - th1)

        return
        end
