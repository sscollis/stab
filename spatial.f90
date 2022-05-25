!==============================================================================
        subroutine spatial(name, ind)
!==============================================================================
!
!       Driver for spatial parallel stability solver with curvature
!       This is the Chebyshev Collocation version.
!
!       I have both IMSL and LAPACK routines.  Currently using LAPACK
!
!       Revised: 10-3-96
!==============================================================================
        use stuff
        use material
        implicit none

        integer i, j, k, ix, ier, ind
        integer i0, idof, j0, jdof
        
        character(80) name
        
        real vm(ny,ndof,nx)
        real eta(ny), y(ny), deta(ny), d2eta(ny)
        
        real D1(ny,ny), D2(ny,ny), Dt1(ny,ny), Dt2(ny,ny)
        
        real A(ny,ndof,ndof), B(ny,ndof,ndof),  C(ny,ndof,ndof)
        real D(ny,ndof,ndof), G(ny,ndof,ndof)
        
        real Vxx(ny,ndof,ndof), Vxy(ny,ndof,ndof), Vyy(ny,ndof,ndof)
        real Vxz(ny,ndof,ndof), Vyz(ny,ndof,ndof), Vzz(ny,ndof,ndof)
        
        real g1vm(ny,ndof),  g2vm(ny,ndof),  g3vm(ny,ndof)
        real g11vm(ny,ndof), g12vm(ny,ndof), g13vm(ny,ndof)
        real g22vm(ny,ndof), g23vm(ny,ndof), g33vm(ny,ndof)
        
        real gum(ny,nsd,nsd), grhom(ny,nsd), gtm(ny,nsd), gpm(ny,nsd)
        real divum(ny), g1divum(ny), g2divum(ny), g3divum(ny)
        
        real rhom(ny),   u1m(ny),    u2m(ny),  u3m(ny),  tm(ny),  pm(ny)

        real rmu(ny),    dmu(ny),    d2mu(ny)
        real rlm(ny),    dlm(ny),    d2lm(ny)
        real con(ny),    dcon(ny),   d2con(ny)
        real g1mu(ny),   g2mu(ny),   g3mu(ny)
        real g1lm(ny),   g2lm(ny),   g3lm(ny)
        real g1con(ny),  g2con(ny),  g3con(ny)
        real g1dmu(ny),  g2dmu(ny),  g3dmu(ny)
        real g1dlm(ny),  g2dlm(ny),  g3dlm(ny)
        real g1dcon(ny), g2dcon(ny), g3dcon(ny)
        
        real S1jj(ny), S2jj(ny), S3jj(ny), S(ny,nsd,nsd), Lapt(ny)
        
        real fact
                
        complex Dh(ny,ndof,ndof), Ah(ny,ndof,ndof), Bh(ny,ndof,ndof)
        complex :: scale
        
        complex A0(2*ndof*ny,2*ndof*ny), B0(2*ndof*ny,2*ndof*ny)
        complex C0(ndof*ny,ndof*ny), C1(ndof*ny,ndof*ny), C2(ndof*ny,ndof*ny)

        complex evec(2*ndof*ny,2*ndof*ny), alp(2*ndof*ny), bet(2*ndof*ny)
        complex omg(2*ndof*ny), cs(2*ndof*ny)
        real    temp1(2*ndof*ny), temp2(2*ndof*ny)
        integer index(2*ndof*ny)

!.... stuff for LAPACK eigensolver

        integer info, lwork
        complex, allocatable :: work(:)
        real, allocatable    :: rwork(:)

!.... stuff for IMSL LU factorization

        integer, allocatable :: ipvt(:)
        real :: rcond

!.... some parameters
        
        real, parameter :: big = 1.0e98

!.... metrics for the curved wall

        real :: h(ny), dhds(ny), dhdr(ny), dhdsr(ny), dhdrr(ny)
!==============================================================================

       write(*,*) "Starting Spatial eigenproblem"

!.... determine if inviscid

        if (Re.ge.big .or. Re.eq.zero) then
          write(*,*) 'A T T E N T I O N:   Inviscid flow'
          Navier = .false.
        end if
        
!.... make the grid

        call sgengrid(y, eta, deta, d2eta)
        
!.... read the mean field

        if (ider) then
          call getmean(vm, y, eta, ny, ind)
        else
          call getmean2(vm, y, eta, g2vm, g22vm, ny, ind)
        end if

        write(*,*) "Finished mean flow"
                    
!.... Initialize the metric terms for curvature

        if (curve.eq.1) then
          call calch( x, ny, y, h, dhds, dhdr, dhdsr, dhdrr )
        else if (curve.eq.2) then
!         h     = x + y         ! note that x is the radius here
!         dhds  = zero
!         dhdr  = one
!         dhdsr = zero
!         dhdrr = zero
          call circh( x, ny, y, h, dhds, dhdr, dhdsr, dhdrr )
        else
          h = one
          dhds = zero
          dhdr = zero
          dhdsr = zero
          dhdrr = zero
        end if
        
!.... loop over the streamwise stations (obsolete)

        do ix = 1, nx
        
!.... initialize

        G   = zero
        A   = zero
        B   = zero
        C   = zero
        D   = zero
        Vxx = zero
        Vxy = zero
        Vyy = zero
        Vxz = zero
        Vyz = zero
        Vzz = zero

!.... setup parallel flow

        rhom = vm(:,1,ix)
        u1m  = vm(:,2,ix)
        u2m  = zero
        u3m  = vm(:,4,ix)
        tm   = vm(:,5,ix)
        pm   = one / (gamma * Ma**2) * rhom * tm

!.... Compute the Chebyshev derivative matrices

        call chebyd(D1, ny-1)   ! use ny-1 since this routine uses 0:ny
        D2  = matmul(D1, D1)    ! compute the second derivative 
        if (wallt.eq.2) then
          Dt1 = D1              ! Derivative operator for temperature
          Dt1(ny,:) = zero      !   adiabatic wall
          Dt2 = matmul(D1, Dt1) ! Second Derivative operator for temperature
        else
          Dt1 = D1
          Dt2 = D2
        end if

!.... Compute derivatives of mean field using the appropriate difference scheme

        g1vm  = zero
        if (ider) then
          do idof = 1, ndof
            g2vm(:,idof) = matmul(D1,vm(:,idof,ix))
          end do
        end if
        g3vm  = zero

        g11vm = zero
        g12vm = zero
        g13vm = zero
        if (ider) then
          do idof = 1, ndof
            g22vm(:,idof) = matmul(D2,vm(:,idof,ix))
          end do
        end if
        g23vm = zero
        g33vm = zero

!.... transform the gradients to physical space

        if (ider) then
          do k = 1, ndof
            g22vm(:,k) = g22vm(:,k) * deta**2 + g2vm(:,k) * d2eta
            g2vm(:,k)  = g2vm(:,k) * deta
          end do
        end if

!.... write out the mean field and its gradients

        if (.true.) then
          open(10,file='rho.out')
          do j = 1, ny
            write (10,13) y(j), vm(j,1,ix), g2vm(j,1), g22vm(j,1)
          end do
          close(10)
          open(10,file='u.out')
          do j = ny, 1, -1
            write (10,13) y(j), vm(j,2,ix), g2vm(j,2), g22vm(j,2)
          end do
          close(10)
          open(10,file='w.out')
          do j = 1, ny
            write (10,13) y(j), vm(j,4,ix), g2vm(j,4), g22vm(j,4)
          end do
          close(10)
          open(10,file='t.out')
          do j = 1, ny
            write (10,13) y(j), vm(j,5,ix), g2vm(j,5), g22vm(j,5)
          end do
          close(10)
13        format(4(1pe20.13,1x))
        end if

!.... initialize gradient of mean velocity

        gum(:,1,1) = g1vm(:,2)
        gum(:,1,2) = g2vm(:,2)
        gum(:,1,3) = g3vm(:,2)
        
        gum(:,2,1) = g1vm(:,3)
        gum(:,2,2) = g2vm(:,3)
        gum(:,2,3) = g3vm(:,3)

        gum(:,3,1) = g1vm(:,4)
        gum(:,3,2) = g2vm(:,4)
        gum(:,3,3) = g3vm(:,4)

!.... compute the divergence (curve*)

        divum = ( gum(:,1,1) + u2m * dhdr ) / h + gum(:,2,2) + gum(:,3,3)
        
!.... initialize gradient of rho and T in the mean

        grhom(:,1) = g1vm(:,1)
        grhom(:,2) = g2vm(:,1)
        grhom(:,3) = g3vm(:,1)

        gtm(:,1)   = g1vm(:,5)
        gtm(:,2)   = g2vm(:,5)
        gtm(:,3)   = g3vm(:,5)
        
        fact = one / (gamma * Ma**2)
        gpm(:,1) = fact * ( grhom(:,1) * tm + rhom * gtm(:,1) )
        gpm(:,2) = fact * ( grhom(:,2) * tm + rhom * gtm(:,2) )
        gpm(:,3) = fact * ( grhom(:,3) * tm + rhom * gtm(:,3) )

!.... compute the gradient of the divergence of um (curve*)

        g1divum = -dhds / h**3 * ( gum(:,1,1) + u2m * dhdr ) + &
                  one / h**2 * ( g11vm(:,2) + gum(:,2,1) * dhdr + &
                  u2m * dhdsr ) + one / h * ( g12vm(:,3) + g13vm(:,4) )

        g2divum = -dhdr / h**2 * ( gum(:,1,1) + u2m * dhdr ) + &
                  one / h * ( g12vm(:,2) + gum(:,2,2) * dhdr + &
                  u2m * dhdrr ) + ( g22vm(:,3) + g23vm(:,4) )

        g3divum = one / h * ( g13vm(:,2) + gum(:,2,3) * dhdr ) + &
                  g23vm(:,3) + g33vm(:,4)
        
!.... compute some stuff that is useful for the viscous terms (curve*)

        S(:,1,1) = ( gum(:,1,1) + u2m * dhdr ) / h
        S(:,1,2) = pt5 * ( ( gum(:,2,1) - u1m * dhdr ) / h + gum(:,1,2) )
        S(:,1,3) = pt5 * ( gum(:,3,1) / h + gum(:,1,3) )
        S(:,2,1) = S(:,1,2)
        S(:,2,2) = gum(:,2,2)
        S(:,2,3) = pt5 * ( gum(:,3,2) + gum(:,2,3) )
        S(:,3,1) = S(:,1,3)
        S(:,3,2) = S(:,2,3)
        S(:,3,3) = gum(:,3,3)

        S1jj = -pt5 * ( dhdr**2 + dhdrr * h ) / h**2 * u1m + &
               pt5 * dhdr * gum(:,1,2) / h + pt5 * g22vm(:,2) - &
               dhds * gum(:,1,1) / h**3 + g11vm(:,2) / h**2 + &
               pt5 * g33vm(:,2) + ( h * dhdsr - dhdr * dhds ) / h**3 * u2m + &
               3.0 * dhdr * gum(:,2,1) / (two * h**2) + &
               pt5 * g12vm(:,3) / h + pt5 * g13vm(:,4) / h 
        
        S2jj = pt5 * ( dhdr * dhds - h * dhdsr ) / h**3 * u1m - &
               3.0 * dhdr * gum(:,1,1) / (two * h**2) + &
               pt5 * g12vm(:,2) / h - dhdr**2 * u2m / h**2 + &
               dhdr * gum(:,2,2) / h + g22vm(:,3) - &
               pt5 * dhds * gum(:,2,1) / h**3 + pt5 * g11vm(:,3) / h**2 + &
               pt5 * g33vm(:,3) + pt5 * g23vm(:,4)
        
        S3jj = pt5 * g13vm(:,2) / h + pt5 * g23vm(:,3) + &
               pt5 * dhdr * gum(:,2,3) / h + pt5 * dhdr * gum(:,3,2) / h + &
               pt5 * g22vm(:,4) + pt5 * g11vm(:,4) / h**2 + g33vm(:,4) - &
               pt5 * dhds * gum(:,3,1) / h**3

!.... Laplacian of Tm (curve*)

        LapT = one/h * ( -dhds/h**2 * gtm(:,1) + one/h * g11vm(:,5) + &
               h * g22vm(:,5) + gtm(:,2) * dhdr + h * g33vm(:,5) )

!.... compute mean material properties
        
        call getmat(tm*te, rmu, rlm, con, dmu, d2mu,    &
                    dlm,   d2lm,   dcon,   d2con)

 !.... nondimensionalize
 
        rmu   = rmu / rmue
        dmu   = dmu * Te / rmue
        d2mu  = d2mu * Te**2 / rmue

        con   = con / cone
        dcon  = dcon * Te / cone
        d2con = d2con * Te**2 / cone
        
        rlm   = rlm / rlme
        dlm   = dlm * Te / rlme
        d2lm  = d2lm * Te**2 / rlme
                
!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * gtm(:,1)
        g2mu = dmu * gtm(:,2)
        g3mu = dmu * gtm(:,3)

        g1dmu = d2mu * gtm(:,1)
        g2dmu = d2mu * gtm(:,2)
        g3dmu = d2mu * gtm(:,3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gtm(:,1)
        g2con = dcon * gtm(:,2)
        g3con = dcon * gtm(:,3)

        g1dcon = d2con * gtm(:,1)
        g2dcon = d2con * gtm(:,2)
        g3dcon = d2con * gtm(:,3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gtm(:,1)
        g2lm = dlm * gtm(:,2)
        g3lm = dlm * gtm(:,3)

        g1dlm = d2lm * gtm(:,1)
        g2dlm = d2lm * gtm(:,2)
        g3dlm = d2lm * gtm(:,3)

!==============================================================================
        
!.... Continuity equation (curve*)

        G(:,1,1) = one
        
        A(:,1,1) = u1m / h
        A(:,1,2) = rhom / h

        B(:,1,1) = u2m
        B(:,1,3) = rhom
        
        C(:,1,1) = u3m
        C(:,1,4) = rhom
        
        D(:,1,1) = divum
        D(:,1,2) = grhom(:,1) / h
        D(:,1,3) = grhom(:,2) + rhom * dhdr / h
        D(:,1,4) = grhom(:,3)
                
!.... Momentum equation -- x_1 (convective + pressure) (curve*)

        G(:,2,2) = rhom
        
        A(:,2,1) = tm/(h * gamma * Ma**2)
        A(:,2,2) = rhom * u1m / h
        A(:,2,5) = rhom/(h * gamma * Ma**2)
        
        B(:,2,2) = rhom * u2m
        
        C(:,2,2) = rhom * u3m
        
        D(:,2,1) = u1m / h * ( gum(:,1,1) + u2m * dhdr ) + &
                   u2m * gum(:,1,2) + u3m * gum(:,1,3)  &
                   + gtm(:,1) / (h * gamma * Ma**2)
        D(:,2,2) = rhom * ( gum(:,1,1) + u2m * dhdr ) / h
        D(:,2,3) = rhom * ( gum(:,1,2) + u1m * dhdr / h )
        D(:,2,4) = rhom * gum(:,1,3)
        D(:,2,5) = grhom(:,1) / (h * gamma * Ma**2)
                
!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,2,2) = A(:,2,2) - fact * ( g1lm / h**2 - rlm / h**3 * dhds )
        A(:,2,3) = A(:,2,3) - fact * rlm / h**2 * dhdr
        A(:,2,5) = A(:,2,5) - fact * dlm * divum / h
        
        B(:,2,3) = B(:,2,3) - fact * ( g1lm / h )
        
        C(:,2,4) = C(:,2,4) - fact * ( g1lm / h )
        
        D(:,2,3) = D(:,2,3) - fact * ( g1lm * dhdr / h**2 - &
                                       rlm / h**3 * dhds * dhdr + &
                                       rlm / h**2 * dhdsr )
        D(:,2,5) = D(:,2,5) - fact * ( g1dlm * divum / h + dlm * g1divum )
        
        Vxx(:,2,2) = fact * rlm / h**2
        
        Vxy(:,2,3) = fact * rlm / h
        
        Vxz(:,2,4) = fact * rlm / h
        
!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,2,2) = A(:,2,2) - fact * ( two * g1mu / h**2 - &
                              two * rmu * dhds / h**3 )
        A(:,2,3) = A(:,2,3) - fact * ( g2mu / h + &
                              rmu * 3.0 * dhdr / h**2 )
        A(:,2,4) = A(:,2,4) - fact * g3mu / h
        A(:,2,5) = A(:,2,5) - fact * dmu * two * S(:,1,1) / h
        
        B(:,2,2) = B(:,2,2) - fact * ( g2mu + rmu * dhdr / h )
        B(:,2,5) = B(:,2,5) - fact * dmu * two * S(:,1,2)
        
        C(:,2,2) = C(:,2,2) - fact * g3mu
        C(:,2,5) = C(:,2,5) - fact * dmu * two * S(:,1,3)
        
        D(:,2,2) = D(:,2,2) - fact * ( g2mu / h * (-dhdr) - &
                              rmu * ( dhdr**2 + dhdrr * h ) / h**2 )
        D(:,2,3) = D(:,2,3) - fact * ( two * g1mu / h**2 * dhdr + &
                              two * rmu * ( dhdsr * h - dhds * dhdr ) / h**3 )
        D(:,2,5) = D(:,2,5) - fact * two * ( g1dmu / h * S(:,1,1) + &
                              g2dmu * S(:,1,2) + g3dmu * S(:,1,3) + &
                              dmu * S1jj )

        Vxx(:,2,2) = Vxx(:,2,2) + fact * two * rmu / h**2
        Vxy(:,2,3) = Vxy(:,2,3) + fact * rmu / h
        Vyy(:,2,2) = Vyy(:,2,2) + fact * rmu
        Vxz(:,2,4) = Vxz(:,2,4) + fact * rmu / h
        Vzz(:,2,2) = Vzz(:,2,2) + fact * rmu

        end if
                                             
!.... Momentum equation -- x_2 (convective + pressure) (curve*)

        G(:,3,3) = rhom
        
        A(:,3,3) = rhom * u1m / h
        
        B(:,3,1) = tm/(gamma * Ma**2)
        B(:,3,3) = rhom * u2m
        B(:,3,5) = rhom/(gamma * Ma**2)
        
        C(:,3,3) = rhom * u3m
        
        D(:,3,1) = u1m / h * ( gum(:,2,1) - u1m * dhdr ) + &
                   u2m * gum(:,2,2) + u3m * gum(:,2,3) + &
                   gtm(:,2) / (gamma * Ma**2)
        D(:,3,2) = rhom * ( gum(:,2,1) - two * u1m * dhdr ) / h
        D(:,3,3) = rhom * gum(:,2,2)
        D(:,3,4) = rhom * gum(:,2,3)
        D(:,3,5) = grhom(:,2) / (gamma * Ma**2)

!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,3,2) = A(:,3,2) - fact * ( g2lm / h - rlm * dhdr / h**2 )
        
        B(:,3,3) = B(:,3,3) - fact * ( g2lm + rlm * dhdr / h )
        B(:,3,5) = B(:,3,5) - fact * dlm * divum
        
        C(:,3,4) = C(:,3,4) - fact * ( g2lm )
        
        D(:,3,3) = D(:,3,3) - fact * ( g2lm / h * dhdr - &
                              rlm * dhdr / h**2 * dhdr + &
                              rlm / h * dhdrr )
        D(:,3,5) = D(:,3,5) - fact * ( g2dlm * divum + dlm * g2divum )

        Vxy(:,3,2) = fact * rlm / h
        
        Vyy(:,3,3) = fact * rlm

        Vyz(:,3,4) = fact * rlm
        
!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,3,2) = A(:,3,2) + fact * rmu * 3.0 * dhdr / h**2
        A(:,3,3) = A(:,3,3) - fact * ( g1mu / h**2 - &
                              rmu * dhds / h**3 )
        A(:,3,5) = A(:,3,5) - fact * dmu * two * S(:,2,1) / h
        
        B(:,3,2) = B(:,3,2) - fact * g1mu / h
        B(:,3,3) = B(:,3,3) - fact * ( two * g2mu + two * rmu * dhdr / h )
        B(:,3,4) = B(:,3,4) - fact * g3mu
        B(:,3,5) = B(:,3,5) - fact * dmu * two * S(:,2,2)
        
        C(:,3,3) = C(:,3,3) - fact * g3mu
        C(:,3,5) = C(:,3,5) - fact * dmu * two * S(:,2,3)
        
        D(:,3,2) = D(:,3,2) - fact * ( g1mu / h**2 * (-dhdr) + &
                              rmu * (dhds * dhdr - h * dhdsr) / h**3 )
        D(:,3,3) = D(:,3,3) + fact * two * rmu * dhdr**2 / h**2
        D(:,3,5) = D(:,3,5) - fact * two * ( g1dmu / h * S(:,2,1) + &
                              g2dmu * S(:,2,2) + g3dmu * S(:,2,3) + &
                              dmu * S2jj )

        Vxx(:,3,3) = Vxx(:,3,3) + fact * rmu / h**2
        Vxy(:,3,2) = Vxy(:,3,2) + fact * rmu / h
        Vyy(:,3,3) = Vyy(:,3,3) + fact * two * rmu
        Vyz(:,3,4) = Vyz(:,3,4) + fact * rmu
        Vzz(:,3,3) = Vzz(:,3,3) + fact * rmu

        end if

!.... Momentum equation -- x_3 (convective + pressure) (curve*)

        G(:,4,4) = rhom
        
        A(:,4,4) = rhom * u1m / h
        
        B(:,4,4) = rhom * u2m
        
        C(:,4,1) = tm/(gamma * Ma**2)
        C(:,4,4) = rhom * u3m
        C(:,4,5) = rhom/(gamma * Ma**2)
        
        D(:,4,1) = u1m * gum(:,3,1) / h + u2m * gum(:,3,2) + &
                   u3m * gum(:,3,3) + gtm(:,3) / (gamma * Ma**2)
        D(:,4,2) = rhom * gum(:,3,1) / h
        D(:,4,3) = rhom * gum(:,3,2)
        D(:,4,4) = rhom * gum(:,3,3)
        D(:,4,5) = grhom(:,3) / (gamma * Ma**2)
        
!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,4,2) = A(:,4,2) - fact * g3lm / h
        
        B(:,4,3) = B(:,4,3) - fact * g3lm
        
        C(:,4,4) = C(:,4,4) - fact * g3lm
        C(:,4,3) = C(:,4,3) - fact * rlm / h * dhdr
        C(:,4,5) = C(:,4,5) - fact * dlm * divum
        
        D(:,4,3) = D(:,4,3) - fact * ( g2lm / h * dhdr )
        D(:,4,5) = D(:,4,5) - fact * ( g3dlm * divum + dlm * g3divum )

        Vxz(:,4,2) = fact * rlm / h
        
        Vyz(:,4,3) = fact * rlm

        Vzz(:,4,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,4,4) = A(:,4,4) - fact * ( g1mu / h**2 - rmu * dhds / h**3 ) 
        A(:,4,5) = A(:,4,5) - fact * dmu * two * S(:,3,1) / h
        
        B(:,4,4) = B(:,4,4) - fact * ( g2mu + rmu * dhdr / h )
        B(:,4,5) = B(:,4,5) - fact * dmu * two * S(:,3,2)
        
        C(:,4,2) = C(:,4,2) - fact * g1mu / h
        C(:,4,3) = C(:,4,3) - fact * ( g2mu + rmu * dhdr / h )
        C(:,4,4) = C(:,4,4) - fact * two * g3mu
        C(:,4,5) = C(:,4,5) - fact * dmu * two * S(:,3,3)
        
        D(:,4,5) = D(:,4,5) - fact * two * ( g1dmu / h * S(:,3,1) + &
                              g2dmu * S(:,3,2) + g3dmu * S(:,3,3) + &
                              dmu * S3jj )

        Vxx(:,4,4) = Vxx(:,4,4) + fact * rmu / h**2
        Vyy(:,4,4) = Vyy(:,4,4) + fact * rmu
        Vxz(:,4,2) = Vxz(:,4,2) + fact * rmu / h
        Vyz(:,4,3) = Vyz(:,4,3) + fact * rmu
        Vzz(:,4,4) = Vzz(:,4,4) + fact * two * rmu

        end if

!.... Energy equation (Advection + pressure) (curve*)

        G(:,5,5) = rhom
        
        A(:,5,2) = rhom * gamma1 * tm / h
        A(:,5,5) = rhom * u1m / h
        
        B(:,5,3) = rhom * gamma1 * tm
        B(:,5,5) = rhom * u2m

        C(:,5,4) = rhom * gamma1 * tm
        C(:,5,5) = rhom * u3m

        D(:,5,1) = u1m / h * gtm(:,1) + u2m * gtm(:,2) + u3m * gtm(:,3) + &
                   gamma1 * tm * divum
        D(:,5,2) = rhom * gtm(:,1) / h
        D(:,5,3) = rhom * gtm(:,2) + rhom * gamma1 * tm * dhdr / h
        D(:,5,4) = rhom * gtm(:,3)

!.... It looks like I'm using the NONLINEAR D(5,5) ??  
!.... It's okay, since d\bar\rho/dt = 0

        D(:,5,5) = rhom * gamma1 * divum

        if (Navier) then
        
!.... diffusion (curve*)

        fact = gamma / (Pr * Re)
        
        A(:,5,5) = A(:,5,5) - fact * (g1con / h**2 + dcon * gtm(:,1) / h**2 - &
                              con * dhds / h**3 )
        B(:,5,5) = B(:,5,5) - fact * (g2con + dcon * gtm(:,2) + &
                              con * dhdr / h )
        C(:,5,5) = C(:,5,5) - fact * (g3con + dcon * gtm(:,3))
        D(:,5,5) = D(:,5,5) - fact * (g1dcon * gtm(:,1) / h**2 + &
                              g2dcon * gtm(:,2) + g3dcon * gtm(:,3) + &
                              dcon * Lapt )

        Vxx(:,5,5) = fact * con / h**2
        Vyy(:,5,5) = fact * con
        Vzz(:,5,5) = fact * con
        
!.... dissipation (lambda) (curve*)

        fact = gamma * gamma1 * Ma**2 * rlme / (Re * rmue)
        
        A(:,5,2) = A(:,5,2) - fact * two * rlm * divum / h
        B(:,5,3) = B(:,5,3) - fact * two * rlm * divum
        C(:,5,4) = C(:,5,4) - fact * two * rlm * divum
        D(:,5,3) = D(:,5,3) - fact * two * rlm * divum * dhdr / h
        D(:,5,5) = D(:,5,5) - fact * dlm * divum * divum
        
!.... dissipation (mu) (curve*)

        fact = gamma * gamma1 * Ma**2 / Re
        
        A(:,5,2) = A(:,5,2) - fact * four * rmu * S(:,1,1) / h
        A(:,5,3) = A(:,5,3) - fact * four * rmu * S(:,2,1) / h
        A(:,5,4) = A(:,5,4) - fact * four * rmu * S(:,3,1) / h

        B(:,5,2) = B(:,5,2) - fact * four * rmu * S(:,1,2)
        B(:,5,3) = B(:,5,3) - fact * four * rmu * S(:,2,2)
        B(:,5,4) = B(:,5,4) - fact * four * rmu * S(:,3,2)

        C(:,5,2) = C(:,5,2) - fact * four * rmu * S(:,1,3)
        C(:,5,3) = C(:,5,3) - fact * four * rmu * S(:,2,3)
        C(:,5,4) = C(:,5,4) - fact * four * rmu * S(:,3,3)

        D(:,5,2) = D(:,5,2) + fact * four * rmu * S(:,2,1) * dhdr / h
        D(:,5,3) = D(:,5,3) - fact * four * rmu * S(:,1,1) * dhdr / h
        D(:,5,5) = D(:,5,5) - fact * two * dmu * ( S(:,1,1)**2 + &
                       S(:,1,2)**2 + S(:,1,3)**2 + S(:,2,1)**2 + &
                       S(:,2,2)**2 + S(:,2,3)**2 + S(:,3,1)**2 + &
                       S(:,3,2)**2 + S(:,3,3)**2)

        end if
        
!==============================================================================

        write(*,*) "Finished building matrices"

!.... initialize

        C0   = zero
        C1   = zero
        C2   = zero

!.... form the equations

        Dh = D + im * beta * C + beta**2 * Vzz
        
        Bh = B - im * beta * Vyz
        
!.... account for the mapping

        do idof = 1, ndof
          do jdof = 1, ndof
            Bh(:,idof,jdof) = Bh(:,idof,jdof) * deta - Vyy(:,idof,jdof) * d2eta
            Vyy(:,idof,jdof) = Vyy(:,idof,jdof) * deta**2
          end do
        end do
        
!.... C0
        
        if (top.eq.1) then
          i = 1
          i0 = (i-1)*ndof
          idof = 1              ! continuity at infinity
            
          do jdof = 1, ndof
            
            do j = 1, ny
              j0 = (j-1)*ndof
              C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + &
                                    Bh(i,idof,jdof) * D1(i,j)
            end do
  
            j0 = (i-1)*ndof
            C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
  
          end do
        end if

        do idof = 1, ndof
          do jdof = 1, ndof

            do i = 2, ny-1
              i0 = (i-1)*ndof

              do j = 1, ny
                j0 = (j-1)*ndof

                C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + &
                                      Bh(i,idof,jdof) * D1(i,j) - &
                                      Vyy(i,idof,jdof) * D2(i,j)
              end do

              j0 = (i-1)*ndof
              C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + Dh(i,idof,jdof)

            end do

          end do
        end do

        i = ny
        i0 = (i-1)*ndof
        idof = 1                ! continuity equation at the wall
          
        do jdof = 1, ndof
          
          do j = 1, ny
            j0 = (j-1)*ndof
            C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * D1(i,j)
          end do

          j0 = (i-1)*ndof
          C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + Dh(i,idof,jdof)

        end do

        if (wallt.eq.2) then
          idof = ndof           ! energy equation at the wall
            
          do jdof = 1, ndof-1
            
            do j = 1, ny
              j0 = (j-1)*ndof
              C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + &
                                    Bh(i,idof,jdof) * D1(i,j) - &
                                    Vyy(i,idof,jdof) * D2(i,j)
            end do
  
            j0 = (i-1)*ndof
            C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
  
          end do
  
          jdof = ndof
            
          do j = 1, ny
            j0 = (j-1)*ndof
            C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * Dt1(i,j) - &
                                  Vyy(i,idof,jdof) * Dt2(i,j)
          end do
  
          j0 = (i-1)*ndof
          C0(i0+idof,j0+jdof) = C0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
        end if

!.... Add the time term into C0 accounting for BC's

        i = 1
        
        if (top.eq.1) then
          C0(1,1) = C0(1,1) - im * omega        ! continuity equation
        else
          C0(1,1) = -one                        ! zero density at infinity
        end if
        C0(2,2) = -one
        C0(3,3) = -one
        C0(4,4) = -one
        C0(5,5) = -one
        
        do idof = 1, ndof
          do jdof = 1, ndof
            do i = 2, ny-1
              i0 = (i-1)*ndof
              C0(i0+idof,i0+jdof) = C0(i0+idof,i0+jdof) - &
                                    im * omega * G(i,idof,jdof)
            end do
          end do
        end do

        i = ny
        i0 = (i-1)*ndof
        idof = 1
        C0(i0+idof,i0+idof) = C0(i0+idof,i0+idof) - im * omega
        idof = 2
        C0(i0+idof,i0+idof) = C0(i0+idof,i0+idof) - one
        idof = 3
        C0(i0+idof,i0+idof) = C0(i0+idof,i0+idof) - one
        idof = 4
        C0(i0+idof,i0+idof) = C0(i0+idof,i0+idof) - one
        idof = ndof
        if (wallt.eq.0) then
          C0(i0+idof,i0+idof) = C0(i0+idof,i0+idof) - one
        else if (wallt.eq.2) then
          do jdof = 1, ndof
            C0(i0+idof,i0+jdof) = C0(i0+idof,i0+jdof) - im * omega * &
                                  G(i,idof,jdof)
          end do
        else
          write(*,"('Illegal value of wallt:  ',i4)") wallt
          call exit(1)
        end if

!.... form the equations

        Dh = im * A + beta * Vxz
        
        Bh = -im * Vxy
        
!.... account for the mapping

        do idof = 1, ndof
          do jdof = 1, ndof
            Bh(:,idof,jdof) = Bh(:,idof,jdof) * deta
          end do
        end do

!.... C1
        
        if (top.eq.1) then
          i = 1
          i0 = (i-1)*ndof
          idof = 1              ! continuity equation at infinity
            
          do jdof = 1, ndof
            
            do j = 1, ny
              j0 = (j-1)*ndof
              C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + &
                                    Bh(i,idof,jdof) * D1(i,j)
            end do
  
            j0 = (i-1)*ndof
            C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + Dh(i,idof,jdof)
  
          end do
        end if

        do idof = 1, ndof
          do jdof = 1, ndof

            do i = 2, ny-1
              i0 = (i-1)*ndof

              do j = 1, ny
                j0 = (j-1)*ndof

                C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + &
                                      Bh(i,idof,jdof) * D1(i,j)
              end do

              j0 = (i-1)*ndof
              C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + Dh(i,idof,jdof)
              
            end do

          end do
        end do
        
        i = ny
        i0 = (i-1)*ndof
        idof = 1                ! density equation at the wall
          
        do jdof = 1, ndof
          
          do j = 1, ny
            j0 = (j-1)*ndof
            C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * D1(i,j)
          end do

          j0 = (i-1)*ndof
          C1(i0+idof,j0+jdof) =  C1(i0+idof,j0+jdof) + Dh(i,idof,jdof)

        end do

        if (wallt.eq.2) then
          i = ny
          i0 = (i-1)*ndof
          idof = ndof           ! energy equation at the wall
            
          do jdof = 1, ndof-1
            
            do j = 1, ny
              j0 = (j-1)*ndof
              C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + &
                                    Bh(i,idof,jdof) * D1(i,j)
            end do
  
            j0 = (i-1)*ndof
            C1(i0+idof,j0+jdof) =  C1(i0+idof,j0+jdof) + Dh(i,idof,jdof)
  
          end do
  
          jdof = ndof
            
          do j = 1, ny
            j0 = (j-1)*ndof
            C1(i0+idof,j0+jdof) = C1(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * Dt1(i,j)
          end do
  
          j0 = (i-1)*ndof
          C1(i0+idof,j0+jdof) =  C1(i0+idof,j0+jdof) + Dh(i,idof,jdof)
        end if

!.... C2 (only viscous contribution)

        do idof = 1, ndof
          do jdof = 1, ndof
            do i = 2, ny-1
              i0 = (i-1)*ndof
              C2(i0+idof,i0+jdof) = Vxx(i,idof,jdof)
            end do
          end do
        end do

        if (wallt.eq.2) then
          i = ny
          i0 = (i-1)*ndof
          idof = ndof           ! energy equation at the wall
          
          do jdof = 1, ndof
            C2(i0+idof,i0+jdof) = Vxx(i,idof,jdof)
          end do
        end if

!.... form the extended system

        write(*,*) "Forming the extended system"

        A0   = zero
        B0   = zero
        evec = zero
        alp  = zero
        bet  = zero
        omg  = zero

        allocate (ipvt(ndof*ny))

#ifdef CRAY
        call CGETRF(ndof*ny, ndof*ny, C0, ndof*ny, ipvt, info)
        if (info.ne.zero) write(*,*) 'CGETRF: ',info
#else
        call ZGETRF(ndof*ny, ndof*ny, C0, ndof*ny, ipvt, info)
        if (info.ne.zero) write(*,*) 'ZGETRF: ',info
#endif
        
        C1 = -C1
        C2 = -C2

#ifdef CRAY
        call CGETRS('N', ndof*ny, ndof*ny, C0, ndof*ny, &
                    ipvt, C1, ndof*ny, info)
        if (info.ne.0) then
          write(*,*) 'CGETRS: ',info
          call exit(1)
        end if
        call CGETRS('N', ndof*ny, ndof*ny, C0, ndof*ny, &
                    ipvt, C2, ndof*ny, info)
        if (info.ne.0) write(*,*) 'CGETRS: ',info
#else
        call ZGETRS('N', ndof*ny, ndof*ny, C0, ndof*ny, &
                    ipvt, C1, ndof*ny, info)
        if (info.ne.0) then
          write(*,*) 'ZGETRS: ',info
          call exit(1)
        end if
        call ZGETRS('N', ndof*ny, ndof*ny, C0, ndof*ny, &
                    ipvt, C2, ndof*ny, info)
        if (info.ne.0) write(*,*) 'ZGETRS: ',info
#endif

        B0(1:ndof*ny,1:ndof*ny)           = C1
        B0(1:ndof*ny,ndof*ny+1:2*ndof*ny) = C2
        
        deallocate(ipvt)

!.... put in the identity matrices

        do i = ndof*ny+1, 2*ndof*ny
          B0(i,i-ndof*ny) = one
        end do
        
        write(*,*) "Solving eigenvalue problem..."

!.... IMSL regular complex eigensolver

!       call EVCCG (2*ndof*ny, B0, 2*ndof*ny, alp, evec, 2*ndof*ny)

!.... LAPACK regular complex eigensolver

        if (.true.) then
        lwork = 2*2*ndof*ny
        allocate (work(lwork), rwork(2*2*ndof*ny), STAT=ier)
        if (ier .ne. 0) then
          write(*,*) 'Error allocating work space'
          call exit(1)
        end if
#ifdef CRAY
        if (ievec.eq.1) then
          call CGEEV('N', 'V', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                      2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
        else
          call CGEEV('N', 'N', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                      2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
        end if
#else
        if (ievec.eq.1) then
          call ZGEEV('N', 'V', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                      2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
        else
          call ZGEEV('N', 'N', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                      2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
        end if
#endif
        if (info.ne.0) then
          if (info.le.0) then
            write(*,*) 'Error in aurgument: ',abs(info)
          else
            write(*,*) 'Warning: Error computing eigenvalues: Error # ',info
          end if
        end if
        deallocate (work, rwork)
        end if

        write(*,*) "Output results"

!.... Note that I have to invert the eigenvalue since I defined the system
!.... Backwards from Bridges and Morris!

        where (alp .ne. 0) 
          alp = 1.0 / alp
        elsewhere
          alp = zero
        end where

!.... sort the eigenvalues by the imaginary part

        do j = 1, 2*ndof*ny
          temp2(j) = AIMAG(alp(j))
          index(j) = j
        end do
        call PIKSR2(2*ndof*ny, temp2, index)
        do j = 1, 2*ndof*ny
          temp1(j) = REAL(alp(index(j)))
          A0(:,j) = evec(:,index(j))
        end do
        
        alp(1:2*ndof*ny) = cmplx(temp1(1:2*ndof*ny),temp2(1:2*ndof*ny))
        evec(:,1:2*ndof*ny) = A0(:,1:2*ndof*ny)
        
!.... end loop on x

        end do
        
!.... compute the phase speed

        where (alp .ne. 0)
          cs = omega / alp
        elsewhere
          cs = zero
        end where

!.... Scale the eigenvectors in a reasonable way

!       if (ievec.eq.1) then
!         do j = 1, 2*ndof*ny
!           scale = zero
!!          do i = ny*ndof+1, 2*ny*ndof
!           do i = 1, ny*ndof
!             if ( abs(evec(i,j)) .gt. abs(scale) ) then
!               scale = evec(i,j)
!             end if
!           end do
!           if (scale .ne. zero) then
!!            do i = ny*ndof+1, 2*ny*ndof
!             do i = 1, ny*ndof
!               evec(i,j) = evec(i,j) / scale
!             end do
!           end if
!         end do
!       end if

!.... write out the eigensystem in an unformatted output file

        open(unit=77,file=name,form='unformatted',status='unknown')
        write (77) ind, ny, ndof, itype, ievec, curve, top, wall, wallt, ider
        write (77) omega, alpha, beta, Re, Ma, Pr
        write (77) x, y, eta, deta, d2eta, yi, ymax
        write (77) alp
        if (ievec.eq.1) write (77) evec
        close (77)

!       do j = 1, ny
!         write(21,*) (vm(j,idof,i), idof=1,ndof)
!       end do

!.... output the eigenvalues and eigenfunctions to the terminal

        if (.false.) then
        
          do j = 1, 2*ndof*ny
            if ( aimag(cs(j)) .gt. zero) then
              write (*,25) j, abs(alp(j)), real(alp(j)), aimag(alp(j)), &
                              real(cs(j)), aimag(cs(j))
            else
              write (*,20) j, abs(alp(j)), real(alp(j)), aimag(alp(j)), &
                              real(cs(j)), aimag(cs(j))
            end if
          end do
        
 100      continue
          write (*,"(/,'Which eigenfunction ==> ',$)")
          read (*,*) j
        
          if ( j .lt. 0 .or. j .gt. 2*ndof*ny ) goto 100

          if (j .ne. 0) then
            open (unit=20,file='space.dat',form='formatted',status='unknown')
            do i = 1, ny
!             i0 = (i+ny-1)*ndof
              i0 = (i-1)*ndof
              write (20,50)eta(i), &
                            real(evec(i0+1,j)), &
                            aimag(evec(i0+1,j)), &
                            real(evec(i0+2,j)), &
                            aimag(evec(i0+2,j)), &
                            real(evec(i0+3,j)), &
                            aimag(evec(i0+3,j)), &
                            real(evec(i0+4,j)), &
                            aimag(evec(i0+4,j)), &
                            real(evec(i0+5,j)), &
                            aimag(evec(i0+5,j))
            end do
            close (20)
            goto 100
          end if
        
        end if

        return
        
 10     format(1p,7(e20.13,1x))
 20     format(1p,i5,5x,e13.6,5x,2(e13.6,1x),5x,2(e13.6,1x))
 25     format(1p,i5,5x,e13.6,5x,2(e13.6,1x),5x,2(e13.6,1x),' <==')
 50     format(1p,11(e20.13,1x))

        end
