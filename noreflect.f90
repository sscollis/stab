!===============================================================================
        subroutine noreflect( bn1, bn2, rho, u1, u2, u3, t, p, &
                              mu, lm, con, dmu, dlm, dcon, &
                              d2mu, d2lm, d2con, grho, gu, gt, gp, &
                              g11v, g12v, g13v, g22v, g23v, g33v, &
                              A, B, D, G)
!===============================================================================
        use stuff
        implicit none
        
        real :: bn1, bn2, bs1, bs2
        real :: rho, u1, u2, u3, t, p
        real :: grho(3), gu(3,3), gt(3), gp(3)
        real :: mu, lm, con, dmu, dlm, dcon, d2mu, d2lm, d2con
        real :: g11v(5), g12v(5), g13v(5), g22v(5), g23v(5), g33v(5)
        real :: A(5,5), B(5,5), D(5,5), G(5,5)

!.... Local variables

        real :: us, un
        real :: gpn, gps, grhon, grhos, gtn, gts
        real :: gu1n, gu1s, gu2n, gu2s, gu3n, gu3s
        real :: gunn, gusn, guns, guss
        real :: g1mu, g2mu, g3mu, g1lm, g2lm, g3lm, g1con, g2con, g3con
        real :: g1dmu, g2dmu, g3dmu, g1dlm, g2dlm, g3dlm, g1dcon, g2dcon, g3dcon
        
        real :: gmsinv, c, dc, fact1, fact2
        real :: divu, g1divu, g2divu, g3divu, rhoinv
        real :: pinf, cinf, Ly, rk, sigma=0.0
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5
        
!===============================================================================

        return

!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * gt(1)
        g2mu = dmu * gt(2)
        g3mu = dmu * gt(3)

        g1dmu = d2mu * gt(1)
        g2dmu = d2mu * gt(2)
        g3dmu = d2mu * gt(3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gt(1)
        g2con = dcon * gt(2)
        g3con = dcon * gt(3)

        g1dcon = d2con * gt(1)
        g2dcon = d2con * gt(2)
        g3dcon = d2con * gt(3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gt(1)
        g2lm = dlm * gt(2)
        g3lm = dlm * gt(3)

        g1dlm = d2lm * gt(1)
        g2dlm = d2lm * gt(2)
        g3dlm = d2lm * gt(3)
        
!===============================================================================

        pinf = one / (gamma * Ma**2)    ! pressure at infinity
        cinf = one / Ma                 ! speed of sound at infinity
        Ly   = 40.0                     ! length of domain in x
        
        rk = sigma * ( one - Ma**2 ) * cinf / Ly

        bs1 = -bn2
        bs2 =  bn1

!.... compute Us and Un
        
        un = bn1 * u1 + bn2 * u2
        us = bs1 * u1 + bs2 * u2

        c = sqrt( t ) / Ma
        dc = pt5 / ( Ma**2 * c )
        rhoinv = one / rho

!.... compute the divergence of velocity

        divu = gu(1,1) + gu(2,2) + gu(3,3)

!.... compute the gradient of the divergence of um

        g1divu = g11v(2) + g12v(3)
        g2divu = g12v(2) + g22v(3)
        g3divu = zero

!.... compute some derivatives in the boundary normal coordinate system

        gpn = bn1 * gp(1) + bn2 * gp(2)
        gps = bs1 * gp(1) + bs2 * gp(2)
    
        grhon = bn1 * grho(1) + bn2 * grho(2)
        grhos = bs1 * grho(1) + bs2 * grho(2)

        gtn = bn1 * gt(1) + bn2 * gt(2)
        gts = bs1 * gt(1) + bs2 * gt(2)

        gu1n = bn1 * gu(1,1) + bn2 * gu(1,2)
        gu1s = bs1 * gu(1,1) + bs2 * gu(1,2)

        gu2n = bn1 * gu(2,1) + bn2 * gu(2,2)
        gu2s = bs1 * gu(2,1) + bs2 * gu(2,2)

        gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
        gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

        gunn = bn1*bn1 * gu(1,1) + bn2*bn1 * gu(2,1) + &
               bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

        gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
               bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

        guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
               bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

        guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
               bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)
            
!.... compute the Characteristic amplitudes
        
        L1 = ( un - c ) * ( gpn - rho * c * gunn )
        L2 = un * ( c**2 * grhon - gpn )
        L3 = un * gusn
        L4 = un * gu3n
        L5 = ( un + c ) * ( gpn + rho * c * gunn )
        
!.... compute the streamwise derivative terms

        d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
        d2 = pt5 * ( L5 + L1 )
        d3 = pt5 / ( rho * c ) * ( L5 - L1 )
        d4 = L3
        d5 = L4

!.... nonreflecting boundary condition
        
        L1 = rk * ( p - pinf )

!==============================================================================!

!.... Continuity equation

        gmsinv = one / ( gamma * Ma**2 )
        
        if (.false.) then       ! don't do for staggered mesh
        
        A(1,1) = one/c**2*( un*(c**2*bn1-bn1*t*gmsinv) + &
                 pt5*(un+c)*bn1*gmsinv*t ) + us*bs1
        A(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn1 + &
                 rho*bs1*bs1
        A(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn1 + &
                 rho*bs2*bs1
        A(1,4) = zero
        A(1,5) = one/c**2*( un*(-bn1*rho*gmsinv) + &
                 pt5*(un+c)*bn1*rho*gmsinv )

        B(1,1) = one/c**2*( un*(c**2*bn2-bn2*t*gmsinv) + &
                 pt5*(un+c)*bn2*gmsinv*t ) + us*bs2
        B(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
        B(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
        B(1,4) = zero
        B(1,5) = one/c**2*( un*(-bn2*rho*gmsinv) + &
                 pt5*(un+c)*bn2*rho*gmsinv )

        D(1,1) = one/c**2*( un*(-gtn*gmsinv) + &
                 pt5*( (un+c)*(gtn*gmsinv+c*gunn) + rk*t*gmsinv ) ) + &
                 guss + gu(3,3)
        D(1,2) = one/c**2*( bn1*(c**2*grhon-gpn) + &
                 pt5*bn1*(gpn+rho*c*gunn) ) + grhos*bs1
        D(1,3) = one/c**2*( bn2*(c**2*grhon-gpn) + &
                 pt5*bn2*(gpn+rho*c*gunn) ) + grhos*bs2
        D(1,4) = grho(3)
        D(1,5) = -two/c**3*dc*( L2 + pt5*( L5 + L1 ) ) + &
                  one/c**2*( un*(two*c*dc*grhon - grhon*gmsinv) + &
                  pt5*( dc*(gpn+rho*c*gunn) + (un+c)*(grhon*gmsinv + &
                  rho*dc*gunn) + rk*rho*gmsinv ) )

        end if

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        A(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn1*t*gmsinv + &
                 bs1*rhoinv*gmsinv*t*bs1
        A(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn1 + &
                 bs1*un*bs1*bn1 + us*bs1 - &
                 fact1 * g1lm - fact1 * two * g1mu
        A(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn1 + &
                 bs1*un*bs2*bn1 - fact1 * g2mu
        A(2,4) = -fact1 * g3mu
        A(2,5) = bn1/(two*rho*c)*(un+c)*bn1*rho*gmsinv + &
                 bs1*gmsinv*bs1 - &
                 fact1 * dlm * ( gu(1,1) + gu(2,2) + gu(3,3) ) - &
                 fact1 * dmu * ( gu(1,1) + gu(1,1) )

        B(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn2*t*gmsinv + &
                 bs1*rhoinv*gmsinv*t*bs2
        B(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn2 + &
                 bs1*un*bs1*bn2 + us*bs2 - fact1 * g2mu
        B(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn2 + &
                 bs1*un*bs2*bn2 - fact1 * g1lm
        B(2,4) = zero
        B(2,5) = bn1/(two*rho*c)*(un+c)*bn2*rho*gmsinv + &
                 bs1*gmsinv*bs2 - fact1 * dmu * ( gu(1,2) + gu(2,1) )

        D(2,1) = -bn1/(two*rho**2*c)*( L5 - L1 ) + &
                  bn1/(two*rho*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - &
                  rk*t*gmsinv ) + bs1*rhoinv*gmsinv*gts + &
                  rhoinv * ( bn1*d3 + bs1*d4 + us*gu1s + u3*gu(1,3) )
        D(2,2) =  bn1/(two*rho*c)*bn1*(gpn + rho*c*gunn) + &
                  bs1*bn1*gusn + bs1*gu1s
        D(2,3) =  bn1/(two*rho*c)*bn2*(gpn + rho*c*gunn) + &
                  bs1*bn2*gusn + bs2*gu1s
        D(2,4) =  gu(1,3)
        D(2,5) = -bn1*dc/(two*rho*c**2)*( L5 - L1 ) + &
                  bn1/(two*rho*c)*( dc*(gpn + rho*c*gunn) + &
                  (un+c)*(grhon*gmsinv + rho*dc*gunn) - rk*rho*gmsinv ) + &
                  bs1*rhoinv*gmsinv*grhos - &
                  fact1 * ( g1dlm * divu + dlm * g1divu ) -     &
                  fact1 * ( g1dmu * ( gu(1,1) + gu(1,1) ) +     &
                            g2dmu * ( gu(1,2) + gu(2,1) ) +     &
                            g3dmu * ( gu(1,3) + gu(3,1) ) +     &
                              dmu * ( g11v(2) + g11v(2) + &
                                      g22v(2) + g12v(3) ) )

!.... Momentum equation -- x_2

        A(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn1*t*gmsinv + &
                 bs2*rhoinv*gmsinv*t*bs1
        A(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn1 + &
                 bs2*un*bs1*bn1 - fact1 * g2lm
        A(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn1 + &
                    bs2*un*bs2*bn1 + us*bs1 - fact1 * g1mu
        A(3,4) = zero
        A(3,5) = bn2/(two*rho*c)*(un+c)*bn1*rho*gmsinv + &
                 bs2*gmsinv*bs1 - fact1*dmu*( gu(2,1) + gu(1,2) )

        B(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn2*t*gmsinv + &
                 bs2 * rhoinv * gmsinv * t * bs2
        B(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn2 + &
                 bs2*un*bs1*bn2 - fact1 * g1mu
        B(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn2 + &
                 us*bs2 + bs2*un*bs2*bn2 - &
                 fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = bn2/(two*rho*c)*(un+c)*bn2*rho*gmsinv + &
                 bs2*gmsinv*bs2 - fact1 * dlm *                 &
                 (gu(1,1) + gu(2,2) + gu(3,3)) - &
                 fact1 * dmu * ( gu(2,2) + gu(2,2) )

        D(3,1) = -bn2/(two*rho**2*c)*( L5 - L1 ) + &
                  bn2/(two*rho*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - &
                  rk*t*gmsinv ) + bs2*rhoinv*gmsinv*gts + &
                  rhoinv * ( bn2*d3 + bs2*d4 + us*gu2s + u3*gu(2,3) )
        D(3,2) =  bn2/(two*rho*c)*bn1*(gpn + rho*c*gunn) + &
                  bs2*bn1*gusn + bs1*gu2s 
        D(3,3) =  bn2/(two*rho*c)*bn2*(gpn + rho*c*gunn) + &
                  bs2*bn2*gusn + bs2*gu2s
        D(3,4) =  gu(2,3)
        D(3,5) = -bn2*dc/(two*rho*c**2)*( L5 - L1 ) + &
                  bn2/(two*rho*c)*( dc*(gpn + rho*c*gunn) + &
                  (un+c)*(grhon*gmsinv + rho*dc*gunn) - rk*rho*gmsinv ) + &
                  bs2*rhoinv*gmsinv*grhos - &
                  fact1 * ( g2dlm * divu + dlm * g2divu ) - &
                   fact1 * ( g1dmu * ( gu(2,1) + gu(1,2) ) + &
                             g2dmu * ( gu(2,2) + gu(2,2) ) + &
                             g3dmu * ( gu(2,3) + gu(3,2) ) + &
                               dmu * ( g11v(3) + g12v(2) + &
                                       g22v(3) + g22v(3) ) )

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) = -Ma**2*rhoinv*( un*(c**2*bn1-bn1*t*gmsinv) - &
                  pt5*gamma1*(un+c)*bn1*t*gmsinv )
        A(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn1 + &
                  gamma1*t*bs1*bs1 - &
                  fact2 * two * lm * divu - &
                  two * fact2 * mu * ( gu(1,1) + gu(1,1) )
        A(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn1 + &
                  gamma1*t*bs2*bs1 - &
                  two * fact2 * mu * ( gu(2,1) + gu(1,2) )
        A(5,4) = -two * fact2 * mu * ( gu(3,1) + gu(1,3) )
        A(5,5) = -Ma**2*rhoinv*( un*(-bn1*rho*gmsinv) - &
                  pt5*gamma1*(un+c)*bn1*rho*gmsinv ) + us*bs1 - &
                  fact1 * (g1con + dcon * gt(1))

        B(5,1) = -Ma**2*rhoinv*( un*(c**2*bn2-bn2*t*gmsinv) - &
                  pt5*gamma1*(un+c)*bn2*t*gmsinv ) 
        B(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn2 + &
                  gamma1*t*bs1*bs2 - &
                  two*fact2*mu*( gu(1,2) + gu(2,1) )
        B(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn2 + &
                  gamma1*t*bs2*bs2 - fact2 * two * lm * &
                  (gu(1,1) + gu(2,2) + gu(3,3)) - &
                  two * fact2 * mu * ( gu(2,2) + gu(2,2) )
        B(5,4) = -two * fact2 * mu * ( gu(3,2) + gu(2,3) )
        B(5,5) = -Ma**2*rhoinv*( un*(-bn2*rho*gmsinv) - &
                  pt5*gamma1*(un+c)*bn2*rho*gmsinv ) + &
                  us*bs2 - fact1 * (g2con + dcon * gt(2))

        D(5,1) =  Ma**2*rhoinv**2*( L2 - pt5*gamma1*(L5 + L1) ) - &
                  Ma**2*rhoinv*( un*(-gtn*gmsinv) - pt5*gamma1* &
                  ((un+c)*(gtn*gmsinv+c*gunn) + rk*t*gmsinv) ) + &
                  rhoinv * ( gamma*Ma**2*rhoinv*d2 - t*rhoinv*d1 + &
                  us*gts + u3 * gt(3) + gamma1*t*(guss + gu(3,3)) )
        D(5,2) = -Ma**2*rhoinv*( bn1*(c**2*grhon-gpn) - &
                  pt5*gamma1*bn1*(gpn+rho*c*gunn) ) + bs1*gts
        D(5,3) = -Ma**2*rhoinv*( bn2*(c**2*grhon-gpn) - &
                  pt5*gamma1*bn2*(gpn+rho*c*gunn) ) + bs2*gts
        D(5,4) =  gt(3)
        D(5,5) =  -gamma1*rhoinv * ( d1 + us*grhos + u3*grho(3) ) - &
                  Ma**2*rhoinv*( un*(two*c*dc*grhon-grhon*gmsinv) - &
                  pt5*gamma1*( dc*(gpn+rho*c*gunn) + &
                  (un+c)*(grhon*gmsinv + rho*dc*gunn) + rk*rho*gmsinv ) ) - &
                  fact1 * (g1dcon * gt(1) + g2dcon * gt(2) + g3dcon * gt(3) + &
                           dcon * (g11v(5) + g22v(5) + g33v(5)) ) - &
                  fact2 * dlm * ( gu(1,1) + gu(2,2) + gu(3,3) )**2 - &
                  two * fact2 * dmu * ( &
                    ( pt5 * ( gu(1,1) + gu(1,1) ) )**2 + &
                    ( pt5 * ( gu(1,2) + gu(2,1) ) )**2 + &
                    ( pt5 * ( gu(1,3) + gu(3,1) ) )**2 + &
                    ( pt5 * ( gu(2,1) + gu(1,2) ) )**2 + &
                    ( pt5 * ( gu(2,2) + gu(2,2) ) )**2 + &
                    ( pt5 * ( gu(2,3) + gu(2,3) ) )**2 + &
                    ( pt5 * ( gu(3,1) + gu(1,3) ) )**2 + &
                    ( pt5 * ( gu(3,2) + gu(2,3) ) )**2 + &
                    ( pt5 * ( gu(3,3) + gu(3,3) ) )**2 )

!.... convert the matrices back to that used in the stability solver

        A = matmul( G, A)
        B = matmul( G, B)
        D = matmul( G, D)

        return
        end
