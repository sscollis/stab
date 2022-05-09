!==============================================================================
        subroutine calch( s, n, r, h, dhds, dhdr, dhdsr, dhdrr )
!==============================================================================
!
!       Compute the metrics for a curved wall.  
!       Set the functions xloc, yloc and arc for the 
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

        real, external :: xloc, yloc

        real :: xl, yl, th, bn1, bn2, dxdy, dydx, dxbds, dybds, dx, dy, &
                ddxdx, ddydx, ddxdy, ddydy, d2dxdx2, d2dydx2, d2dxdy2,  &
                d2dydy2, dbn1, dbn2, d2xdy2, d2ydx2, d2xbds2, d2ybds2,  &
                d2bn1, d2bn2

        real :: a(n), b(n), dads(n), dbds(n)
        
        real, parameter :: zero = 0.0, pt5 = 0.5, one = 1.0, onept5 = 1.5, &
                           two = 2.0,  twopt5 = 2.5, three = 3.0, &
                           infty = 1.0e30
!==============================================================================
        if (s.eq.-one) then
!         write(*,*) 'WARNING:  Curvature is turned off!'
          h     = one
          dhdr  = zero
          dhds  = zero
          dhdrr = zero
          dhdsr = zero
          return
        end if

        xl  = xloc( zero, s )
        yl  = yloc( xl )
        th  = atan2( one, yl )
        bn1 = -sin(th)
        bn2 =  cos(th)

!       write(*,"(8(e13.6,1x))") s, xl, yl, th, bn1, bn2
        
        if (xl .eq. zero) then
          dydx = infty
        else
          dydx = one / sqrt( two * xl )
        end if
        dxbds  = one / sqrt(one + dydx**2)

        if (xl .eq. zero) then
          dxdy = zero
        else
          dxdy = sqrt( two * xl )
        end if
        dybds = one / sqrt( dxdy**2 + one )

        if (xl .eq. zero) then
          dx = zero
          dy = one
          ddxdx = infty
          ddydx = zero
          ddxdy = one
          ddydy = zero
          d2dxdx2 = infty
          d2dydx2 = zero
          d2dxdy2 = zero
          d2dydy2 = zero
        else
          dx = yl
          dy = one
          ddxdx = one / yl
          ddydx = zero
          ddxdy = one
          ddydy = zero
          d2dxdx2 = -sqrt(two)/(two**2 * xl**onept5)
          d2dydx2 = zero
          d2dxdy2 = zero
          d2dydy2 = zero
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
          d2xdy2  = one
          d2ybds2 = -(one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
          d2xbds2 = d2xdy2*(dybds)**2 + dxdy*d2ybds2
        else
          d2ydx2  = -one / ( two * sqrt(two) * xl**onept5 )
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
!       end do
        
 10     format( 8(1pe13.6,1x) )
 
        return
        end

!==============================================================================
        function xloc(x, ds)
!==============================================================================
        implicit none
        
        real :: xloc, x, ds
        real, external :: rtflsp, func
        
        real :: darc, x1
        common /distance/ darc, x1
!==============================================================================
        darc = ds
        x1   = x
        xloc = rtflsp(func, x1, x1 + 2.0 * ds, 1.0e-12)
        return
        end 

!==============================================================================
        function yloc(x)
!==============================================================================
!       For a parabolic cylinder
!==============================================================================
        implicit none
        
        real :: yloc, x
!==============================================================================
        yloc = sqrt( 2.0 * x )
        return
        end

!==============================================================================
        function func(x)
!==============================================================================
        implicit none
        
        real :: func, x
        real, external :: arc
        
        real :: darc, x1
        common /distance/ darc, x1
!==============================================================================
        func = darc - arc(x1,x)
        return
        end 

!==============================================================================
        function arc(x1,x2)
!==============================================================================
!       For a parabolic cylinder
!==============================================================================
        implicit none
        
        real :: arc, x1, x2, xi1, xi2
!==============================================================================
        xi1 = sqrt(x1)
        xi2 = sqrt(x2)

        arc = sqrt(2.0)*0.5*xi2*sqrt(1.0+2.0*xi2**2) + &
              0.5*log(sqrt(2.0)*xi2 + sqrt(1.0+2.0*xi2**2))

        arc = arc - ( sqrt(2.0)*0.5*xi1*sqrt(1.0+2.0*xi1**2) + &
              0.5*log(sqrt(2.0)*xi1 + sqrt(1.0+2.0*xi1**2)) )

        return
        end

!==============================================================================
        FUNCTION RTFLSP(FUNC,X1,X2,XACC)
!==============================================================================
!
!       From Numerical Recipes
!
!==============================================================================
        PARAMETER (MAXIT=100)
        EXTERNAL FUNC
!==============================================================================
        FL=FUNC(X1)
        IF(FL.EQ.0.0) THEN
          RTFLSP=X1
          RETURN
        ENDIF
        FH=FUNC(X2)
        IF(FH.EQ.0.0) THEN
          RTFLSP=X2
          RETURN
        ENDIF
        IF(FL*FH.GT.0.) then
          write(*,*) x1, x2, fl, fh
          write(*,*) 'Root must be bracketed for false position.'
          stop
        end if
        IF(FL.LT.0.)THEN
          XL=X1
          XH=X2
        ELSE
          XL=X2
          XH=X1
          SWAP=FL
          FL=FH
          FH=SWAP
        ENDIF
        DX=XH-XL
        DO 11 J=1,MAXIT
          RTFLSP=XL+DX*FL/(FL-FH)
          F=FUNC(RTFLSP)
          IF(F.LT.0.) THEN
            DEL=XL-RTFLSP
            XL=RTFLSP
            FL=F
          ELSE
            DEL=XH-RTFLSP
            XH=RTFLSP
            FH=F
          ENDIF
          DX=XH-XL
          IF(ABS(DEL).LT.XACC.OR.F.EQ.0.)RETURN
11      CONTINUE
        write(*,*) 'RTFLSP exceed maximum iterations'
        stop
        END
