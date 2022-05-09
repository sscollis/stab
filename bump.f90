!==============================================================================
        subroutine bump(name, ind, wallbc)
!==============================================================================
!
!       Driver for receptivity solver with curvature.  
!       This routine uses the extended system approach with Chebyshev 
!       Collocation on a nonstaggered mesh.
!
!       It could be extensively cleaned up as could all of the STAB routines.
!
!       Uses LAPACK routines
!
!       S. Scott Collis
!
!       Revised: 10-2-96
!==============================================================================
        use stuff
        use material
        implicit none

        integer i, j, k, ix, ier, ind, wallbc
        integer i0, idof, j0, jdof, neq
        
        character*80 name
        
        real vm(ny,ndof,nx)

        real eta(ny), y(ny), deta(ny), d2eta(ny)
        
        real D1(ny,ny), D2(ny,ny), Dt1(ny,ny), Dt2(ny,ny)
        
        real    :: A(ny,ndof,ndof), C(ny,ndof,ndof)
        real    :: D(ny,ndof,ndof), G(ny,ndof,ndof)
        complex :: B(ny,ndof,ndof)
        
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
        
        real fact, fact1(ny), dke

        complex, allocatable :: Eh(:,:,:), Fh(:,:,:)
        complex, allocatable :: A0(:,:,:,:), rhs(:,:), duda(:,:)

        real :: yeff, yl(ny-1)
        complex :: alphap, eff, pro(ny-1)
        integer :: ialpha, info

!.... some parameters
        
        real, parameter :: big = 1.0e98

!.... metrics for the curved wall

        real :: h(ny), dhds(ny), dhdr(ny), dhdsr(ny), dhdrr(ny)
        real :: hinv(ny), hinv2(ny), hinv3(ny)

!.... stuff for LU factorization

        integer, allocatable :: ipiv(:)
        real :: rcond, cpu

        real, external :: second
        complex, external :: findmax
!==============================================================================

!       cpu = second()

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
        hinv = one / h
        hinv2 = hinv**2
        hinv3 = hinv**3
        
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
        
        ix = 1          ! obsolete

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

!       if (.false.) then
!         open(10,file='rho.out')
!         do j = 1, ny
!           write (10,13) y(j), vm(j,1,ix), g2vm(j,1), g22vm(j,1)
!         end do
!         close(10)
!         open(10,file='u.out')
!         do j = 1, ny
!           write (10,13) y(j), vm(j,2,ix), g2vm(j,2), g22vm(j,2)
!         end do
!         close(10)
!         open(10,file='w.out')
!         do j = 1, ny
!           write (10,13) y(j), vm(j,4,ix), g2vm(j,4), g22vm(j,4)
!         end do
!         close(10)
!         open(10,file='t.out')
!         do j = 1, ny
!           write (10,13) y(j), vm(j,5,ix), g2vm(j,5), g22vm(j,5)
!         end do
!         close(10)
!       end if

!.... initialize gradient of mean velocity (not the velocity gradient tensor)

        gum(:,1,1) = g1vm(:,2)
        gum(:,1,2) = g2vm(:,2)
        gum(:,1,3) = g3vm(:,2)
        
        gum(:,2,1) = g1vm(:,3)
        gum(:,2,2) = g2vm(:,3)
        gum(:,2,3) = g3vm(:,3)

        gum(:,3,1) = g1vm(:,4)
        gum(:,3,2) = g2vm(:,4)
        gum(:,3,3) = g3vm(:,4)

!.... compute the mean divergence (curve*)

        divum = ( gum(:,1,1) + u2m * dhdr ) * hinv + gum(:,2,2) + gum(:,3,3)
        
!.... initialize gradient of rho and T in the mean
!.... not the gradient for curvilinear coordinates

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

        g1divum = -dhds * hinv3 * ( gum(:,1,1) + u2m * dhdr ) + &
                  hinv2 * ( g11vm(:,2) + gum(:,2,1) * dhdr + &
                  u2m * dhdsr ) + hinv * ( g12vm(:,3) + g13vm(:,4) )

        g2divum = -dhdr * hinv2 * ( gum(:,1,1) + u2m * dhdr ) + &
                  hinv * ( g12vm(:,2) + gum(:,2,2) * dhdr + &
                  u2m * dhdrr ) + ( g22vm(:,3) + g23vm(:,4) )

        g3divum = hinv * ( g13vm(:,2) + gum(:,2,3) * dhdr ) + &
                  g23vm(:,3) + g33vm(:,4)
        
!.... compute some stuff that is useful for the viscous terms (curve*)

        S(:,1,1) = ( gum(:,1,1) + u2m * dhdr ) * hinv
        S(:,1,2) = pt5 * ( ( gum(:,2,1) - u1m * dhdr ) * hinv + gum(:,1,2) )
        S(:,1,3) = pt5 * ( gum(:,3,1) * hinv + gum(:,1,3) )
        S(:,2,1) = S(:,1,2)
        S(:,2,2) = gum(:,2,2)
        S(:,2,3) = pt5 * ( gum(:,3,2) + gum(:,2,3) )
        S(:,3,1) = S(:,1,3)
        S(:,3,2) = S(:,2,3)
        S(:,3,3) = gum(:,3,3)

        S1jj = -pt5 * ( dhdr**2 + dhdrr * h ) * hinv2 * u1m + &
               pt5 * dhdr * gum(:,1,2) * hinv + pt5 * g22vm(:,2) - &
               dhds * gum(:,1,1) * hinv3 + g11vm(:,2) * hinv2 + &
               pt5 * g33vm(:,2) + ( h * dhdsr - dhdr * dhds ) * hinv3 * u2m + &
               3.0 * dhdr * gum(:,2,1) / (two * h**2) + &
               pt5 * g12vm(:,3) * hinv + pt5 * g13vm(:,4) * hinv 
        
        S2jj = pt5 * ( dhdr * dhds - h * dhdsr ) * hinv3 * u1m - &
               3.0 * dhdr * gum(:,1,1) / (two * h**2) + &
               pt5 * g12vm(:,2) * hinv - dhdr**2 * u2m * hinv2 + &
               dhdr * gum(:,2,2) * hinv + g22vm(:,3) - &
               pt5 * dhds * gum(:,2,1) * hinv3 + pt5 * g11vm(:,3) * hinv2 + &
               pt5 * g33vm(:,3) + pt5 * g23vm(:,4)
        
        S3jj = pt5 * g13vm(:,2) * hinv + pt5 * g23vm(:,3) + &
               pt5 * dhdr * gum(:,2,3) * hinv + &
               pt5 * dhdr * gum(:,3,2) * hinv + &
               pt5 * g22vm(:,4) + pt5 * g11vm(:,4) * hinv2 + g33vm(:,4) - &
               pt5 * dhds * gum(:,3,1) * hinv3

!.... Laplacian of Tm (curve*)

        LapT = hinv * ( -dhds/h**2 * gtm(:,1) + hinv * g11vm(:,5) + &
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
        
        A(:,1,1) = u1m * hinv
        A(:,1,2) = rhom * hinv

        B(:,1,1) = u2m
        B(:,1,3) = rhom
        
        C(:,1,1) = u3m
        C(:,1,4) = rhom
        
        D(:,1,1) = divum
        D(:,1,2) = grhom(:,1) * hinv
        D(:,1,3) = grhom(:,2) + rhom * dhdr * hinv
        D(:,1,4) = grhom(:,3)
                
!.... Momentum equation -- x_1 (convective + pressure) (curve*)

        G(:,2,2) = rhom
        
        A(:,2,1) = tm/(h * gamma * Ma**2)
        A(:,2,2) = rhom * u1m * hinv
        A(:,2,5) = rhom/(h * gamma * Ma**2)
        
        B(:,2,2) = rhom * u2m
        
        C(:,2,2) = rhom * u3m
        
        D(:,2,1) = u1m * hinv * ( gum(:,1,1) + u2m * dhdr ) + &
                   u2m * gum(:,1,2) + u3m * gum(:,1,3)  &
                   + gtm(:,1) / (h * gamma * Ma**2)
        D(:,2,2) = rhom * ( gum(:,1,1) + u2m * dhdr ) * hinv
        D(:,2,3) = rhom * ( gum(:,1,2) + u1m * dhdr * hinv )
        D(:,2,4) = rhom * gum(:,1,3)
        D(:,2,5) = grhom(:,1) / (h * gamma * Ma**2)
                
!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,2,2) = A(:,2,2) - fact * ( g1lm * hinv2 - rlm * hinv3 * dhds )
        A(:,2,3) = A(:,2,3) - fact * rlm * hinv2 * dhdr
        A(:,2,5) = A(:,2,5) - fact * dlm * divum * hinv
        
        B(:,2,3) = B(:,2,3) - fact * ( g1lm * hinv )
        
        C(:,2,4) = C(:,2,4) - fact * ( g1lm * hinv )
        
        D(:,2,3) = D(:,2,3) - fact * ( g1lm * dhdr * hinv2 - &
                                       rlm * hinv3 * dhds * dhdr + &
                                       rlm * hinv2 * dhdsr )
        D(:,2,5) = D(:,2,5) - fact * ( g1dlm * divum * hinv + dlm * g1divum )
        
        Vxx(:,2,2) = fact * rlm * hinv2
        
        Vxy(:,2,3) = fact * rlm * hinv
        
        Vxz(:,2,4) = fact * rlm * hinv
        
!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,2,2) = A(:,2,2) - fact * ( two * g1mu * hinv2 - &
                              two * rmu * dhds * hinv3 )
        A(:,2,3) = A(:,2,3) - fact * ( g2mu * hinv + &
                              rmu * 3.0 * dhdr * hinv2 )
        A(:,2,4) = A(:,2,4) - fact * g3mu * hinv
        A(:,2,5) = A(:,2,5) - fact * dmu * two * S(:,1,1) * hinv
        
        B(:,2,2) = B(:,2,2) - fact * ( g2mu + rmu * dhdr * hinv )
        B(:,2,5) = B(:,2,5) - fact * dmu * two * S(:,1,2)
        
        C(:,2,2) = C(:,2,2) - fact * g3mu
        C(:,2,5) = C(:,2,5) - fact * dmu * two * S(:,1,3)
        
        D(:,2,2) = D(:,2,2) - fact * ( g2mu * hinv * (-dhdr) - &
                              rmu * ( dhdr**2 + dhdrr * h ) * hinv2 )
        D(:,2,3) = D(:,2,3) - fact * ( two * g1mu * hinv2 * dhdr + &
                              two * rmu * ( dhdsr * h - dhds * dhdr ) * hinv3 )
        D(:,2,5) = D(:,2,5) - fact * two * ( g1dmu * hinv * S(:,1,1) + &
                              g2dmu * S(:,1,2) + g3dmu * S(:,1,3) + &
                              dmu * S1jj )

        Vxx(:,2,2) = Vxx(:,2,2) + fact * two * rmu * hinv2
        Vxy(:,2,3) = Vxy(:,2,3) + fact * rmu * hinv
        Vyy(:,2,2) = Vyy(:,2,2) + fact * rmu
        Vxz(:,2,4) = Vxz(:,2,4) + fact * rmu * hinv
        Vzz(:,2,2) = Vzz(:,2,2) + fact * rmu

        end if
                                             
!.... Momentum equation -- x_2 (convective + pressure) (curve*)

        G(:,3,3) = rhom
        
        A(:,3,3) = rhom * u1m * hinv
        
        B(:,3,1) = tm/(gamma * Ma**2)
        B(:,3,3) = rhom * u2m
        B(:,3,5) = rhom/(gamma * Ma**2)
        
        C(:,3,3) = rhom * u3m
        
        D(:,3,1) = u1m * hinv * ( gum(:,2,1) - u1m * dhdr ) + &
                   u2m * gum(:,2,2) + u3m * gum(:,2,3) + &
                   gtm(:,2) / (gamma * Ma**2)
        D(:,3,2) = rhom * ( gum(:,2,1) - two * u1m * dhdr ) * hinv
        D(:,3,3) = rhom * gum(:,2,2)
        D(:,3,4) = rhom * gum(:,2,3)
        D(:,3,5) = grhom(:,2) / (gamma * Ma**2)

!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,3,2) = A(:,3,2) - fact * ( g2lm * hinv - rlm * dhdr * hinv2 )
        
        B(:,3,3) = B(:,3,3) - fact * ( g2lm + rlm * dhdr * hinv )
        B(:,3,5) = B(:,3,5) - fact * dlm * divum
        
        C(:,3,4) = C(:,3,4) - fact * ( g2lm )
        
        D(:,3,3) = D(:,3,3) - fact * ( g2lm * hinv * dhdr - &
                              rlm * dhdr * hinv2 * dhdr + &
                              rlm * hinv * dhdrr )
        D(:,3,5) = D(:,3,5) - fact * ( g2dlm * divum + dlm * g2divum )

        Vxy(:,3,2) = fact * rlm * hinv
        
        Vyy(:,3,3) = fact * rlm

        Vyz(:,3,4) = fact * rlm
        
!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,3,2) = A(:,3,2) + fact * rmu * 3.0 * dhdr * hinv2
        A(:,3,3) = A(:,3,3) - fact * ( g1mu * hinv2 - &
                              rmu * dhds * hinv3 )
        A(:,3,5) = A(:,3,5) - fact * dmu * two * S(:,2,1) * hinv
        
        B(:,3,2) = B(:,3,2) - fact * g1mu * hinv
        B(:,3,3) = B(:,3,3) - fact * ( two * g2mu + two * rmu * dhdr * hinv )
        B(:,3,4) = B(:,3,4) - fact * g3mu
        B(:,3,5) = B(:,3,5) - fact * dmu * two * S(:,2,2)
        
        C(:,3,3) = C(:,3,3) - fact * g3mu
        C(:,3,5) = C(:,3,5) - fact * dmu * two * S(:,2,3)
        
        D(:,3,2) = D(:,3,2) - fact * ( g1mu * hinv2 * (-dhdr) + &
                              rmu * (dhds * dhdr - h * dhdsr) * hinv3 )
        D(:,3,3) = D(:,3,3) + fact * two * rmu * dhdr**2 * hinv2
        D(:,3,5) = D(:,3,5) - fact * two * ( g1dmu * hinv * S(:,2,1) + &
                              g2dmu * S(:,2,2) + g3dmu * S(:,2,3) + &
                              dmu * S2jj )

        Vxx(:,3,3) = Vxx(:,3,3) + fact * rmu * hinv2
        Vxy(:,3,2) = Vxy(:,3,2) + fact * rmu * hinv
        Vyy(:,3,3) = Vyy(:,3,3) + fact * two * rmu
        Vyz(:,3,4) = Vyz(:,3,4) + fact * rmu
        Vzz(:,3,3) = Vzz(:,3,3) + fact * rmu

!.... now use continuity to remove Vyy(:,3,3) term

        fact1 = Vyy(:,3,3) / rhom
        
        A(:,3,1) = A(:,3,1) + fact1 * ( gum(:,1,2) * hinv - &
                                        u1m * dhdr * hinv2 )
        A(:,3,2) = A(:,3,2) + fact1 * ( grhom(:,2) * hinv - &
                                        rhom * dhdr * hinv2 )

        B(:,3,1) = B(:,3,1) + fact1 * ( divum + gum(:,2,2) - im * omega )
        B(:,3,2) = B(:,3,2) + fact1 * ( grhom(:,1) * hinv )
        B(:,3,3) = B(:,3,3) + fact1 * ( two * grhom(:,2) + rhom * dhdr * hinv )
        B(:,3,4) = B(:,3,4) + fact1 * ( grhom(:,3) )
        
        C(:,3,1) = C(:,3,1) + fact1 * ( gum(:,3,2) )
        C(:,3,4) = C(:,3,4) + fact1 * ( grhom(:,2) )
        
        D(:,3,1) = D(:,3,1) + fact1 * ( g2divum )
        D(:,3,2) = D(:,3,2) + fact1 * ( g12vm(:,1) * hinv - &
                                        grhom(:,1) * dhdr * hinv2 )
        D(:,3,3) = D(:,3,3) + fact1 * ( grhom(:,2) * dhdr * hinv - &
                              rhom * dhdr**2 * hinv2 + rhom * dhdrr * hinv + &
                              g22vm(:,1) )
        D(:,3,4) = D(:,3,4) +  fact1 * ( g23vm(:,1) )
        
        Vxy(:,3,1) = Vxy(:,3,1) - fact1 * u1m * hinv
        Vxy(:,3,2) = Vxy(:,3,2) - fact1 * rhom * hinv
        Vyy(:,3,1) = Vyy(:,3,1) - fact1 * u2m
        Vyz(:,3,1) = Vyz(:,3,1) - fact1 * u3m
        Vyz(:,3,4) = Vyz(:,3,4) - fact1 * rhom

!.... zero the offending term

        Vyy(:,3,3) = zero

        end if

!.... Momentum equation -- x_3 (convective + pressure) (curve*)

        G(:,4,4) = rhom
        
        A(:,4,4) = rhom * u1m * hinv
        
        B(:,4,4) = rhom * u2m
        
        C(:,4,1) = tm/(gamma * Ma**2)
        C(:,4,4) = rhom * u3m
        C(:,4,5) = rhom/(gamma * Ma**2)
        
        D(:,4,1) = u1m * gum(:,3,1) * hinv + u2m * gum(:,3,2) + &
                   u3m * gum(:,3,3) + gtm(:,3) / (gamma * Ma**2)
        D(:,4,2) = rhom * gum(:,3,1) * hinv
        D(:,4,3) = rhom * gum(:,3,2)
        D(:,4,4) = rhom * gum(:,3,3)
        D(:,4,5) = grhom(:,3) / (gamma * Ma**2)
        
!.... (viscous lambda) (curve*)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,4,2) = A(:,4,2) - fact * g3lm * hinv
        
        B(:,4,3) = B(:,4,3) - fact * g3lm
        
        C(:,4,4) = C(:,4,4) - fact * g3lm
        C(:,4,3) = C(:,4,3) - fact * rlm * hinv * dhdr
        C(:,4,5) = C(:,4,5) - fact * dlm * divum
        
        D(:,4,3) = D(:,4,3) - fact * ( g2lm * hinv * dhdr )
        D(:,4,5) = D(:,4,5) - fact * ( g3dlm * divum + dlm * g3divum )

        Vxz(:,4,2) = fact * rlm * hinv
        
        Vyz(:,4,3) = fact * rlm

        Vzz(:,4,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re
        
        A(:,4,4) = A(:,4,4) - fact * ( g1mu * hinv2 - rmu * dhds * hinv3 ) 
        A(:,4,5) = A(:,4,5) - fact * dmu * two * S(:,3,1) * hinv
        
        B(:,4,4) = B(:,4,4) - fact * ( g2mu + rmu * dhdr * hinv )
        B(:,4,5) = B(:,4,5) - fact * dmu * two * S(:,3,2)
        
        C(:,4,2) = C(:,4,2) - fact * g1mu * hinv
        C(:,4,3) = C(:,4,3) - fact * ( g2mu + rmu * dhdr * hinv )
        C(:,4,4) = C(:,4,4) - fact * two * g3mu
        C(:,4,5) = C(:,4,5) - fact * dmu * two * S(:,3,3)
        
        D(:,4,5) = D(:,4,5) - fact * two * ( g1dmu * hinv * S(:,3,1) + &
                              g2dmu * S(:,3,2) + g3dmu * S(:,3,3) + &
                              dmu * S3jj )

        Vxx(:,4,4) = Vxx(:,4,4) + fact * rmu * hinv2
        Vyy(:,4,4) = Vyy(:,4,4) + fact * rmu
        Vxz(:,4,2) = Vxz(:,4,2) + fact * rmu * hinv
        Vyz(:,4,3) = Vyz(:,4,3) + fact * rmu
        Vzz(:,4,4) = Vzz(:,4,4) + fact * two * rmu

        end if

!.... Energy equation (Advection + pressure) (curve*)

        G(:,5,5) = rhom
        
        A(:,5,2) = rhom * gamma1 * tm * hinv
        A(:,5,5) = rhom * u1m * hinv
        
        B(:,5,3) = rhom * gamma1 * tm
        B(:,5,5) = rhom * u2m

        C(:,5,4) = rhom * gamma1 * tm
        C(:,5,5) = rhom * u3m

        D(:,5,1) = u1m * hinv * gtm(:,1) + u2m * gtm(:,2) + u3m * gtm(:,3) + &
                   gamma1 * tm * divum
        D(:,5,2) = rhom * gtm(:,1) * hinv
        D(:,5,3) = rhom * gtm(:,2) + rhom * gamma1 * tm * dhdr * hinv
        D(:,5,4) = rhom * gtm(:,3)
        D(:,5,5) = rhom * gamma1 * divum

        if (Navier) then
        
!.... diffusion (curve*)

        fact = gamma / (Pr * Re)
        
        A(:,5,5) = A(:,5,5) - fact * (g1con * hinv2 + &
                              dcon * gtm(:,1) * hinv2 - &
                              con * dhds * hinv3 )
        B(:,5,5) = B(:,5,5) - fact * (g2con + dcon * gtm(:,2) + &
                              con * dhdr * hinv )
        C(:,5,5) = C(:,5,5) - fact * (g3con + dcon * gtm(:,3))
        D(:,5,5) = D(:,5,5) - fact * (g1dcon * gtm(:,1) * hinv2 + &
                              g2dcon * gtm(:,2) + g3dcon * gtm(:,3) + &
                              dcon * Lapt )

        Vxx(:,5,5) = fact * con * hinv2
        Vyy(:,5,5) = fact * con
        Vzz(:,5,5) = fact * con
        
!.... dissipation (lambda) (curve*)

        fact = gamma * gamma1 * Ma**2 * rlme / (Re * rmue)
        
        A(:,5,2) = A(:,5,2) - fact * two * rlm * divum * hinv
        B(:,5,3) = B(:,5,3) - fact * two * rlm * divum
        C(:,5,4) = C(:,5,4) - fact * two * rlm * divum
        D(:,5,3) = D(:,5,3) - fact * two * rlm * divum * dhdr * hinv
        D(:,5,5) = D(:,5,5) - fact * dlm * divum * divum
        
!.... dissipation (mu) (curve*)

        fact = gamma * gamma1 * Ma**2 / Re
        
        A(:,5,2) = A(:,5,2) - fact * four * rmu * S(:,1,1) * hinv
        A(:,5,3) = A(:,5,3) - fact * four * rmu * S(:,2,1) * hinv
        A(:,5,4) = A(:,5,4) - fact * four * rmu * S(:,3,1) * hinv

        B(:,5,2) = B(:,5,2) - fact * four * rmu * S(:,1,2)
        B(:,5,3) = B(:,5,3) - fact * four * rmu * S(:,2,2)
        B(:,5,4) = B(:,5,4) - fact * four * rmu * S(:,3,2)

        C(:,5,2) = C(:,5,2) - fact * four * rmu * S(:,1,3)
        C(:,5,3) = C(:,5,3) - fact * four * rmu * S(:,2,3)
        C(:,5,4) = C(:,5,4) - fact * four * rmu * S(:,3,3)

        D(:,5,2) = D(:,5,2) + fact * four * rmu * S(:,2,1) * dhdr * hinv
        D(:,5,3) = D(:,5,3) - fact * four * rmu * S(:,1,1) * dhdr * hinv
        D(:,5,5) = D(:,5,5) - fact * two * dmu * ( S(:,1,1)**2 + &
                       S(:,1,2)**2 + S(:,1,3)**2 + S(:,2,1)**2 + &
                       S(:,2,2)**2 + S(:,2,3)**2 + S(:,3,1)**2 + &
                       S(:,3,2)**2 + S(:,3,3)**2)
        end if

!       write(*,*) 'Form system ',second()-cpu
!       cpu = second()
!==============================================================================
!
!.... Setup the extended system
!
!       I did it this way because I originally wanted to solve the adjoint
!       problem using this code.  However, it soon became clear that 
!       enforcing boundary conditions on the global adjoint problem was
!       quite difficult.  So I converted this code to solve the receptivity
!       problem with the curvature terms included.
!
!==============================================================================

        neq = 8         ! now BIGGER for the extended system
        
        allocate( Eh(ny,neq,neq), Fh(ny,neq,neq) )
        allocate( A0(ny,neq,ny,neq), rhs(ny,neq), duda(ny,neq) )
        allocate( ipiv(ny*neq) )

        do ialpha = 1, 2        ! do two alpha's

        if (ialpha.eq.1) then
          alphap = alpha + dalpha
        else
          alphap = alpha - dalpha
        end if

        Eh = zero
        Fh = zero
        
        Eh(:,1:5,1:5) = B - im * alphap * Vxy - im * beta * Vyz
        Eh(:,2,6)     = -Vyy(:,2,2)
        Eh(:,4,7)     = -Vyy(:,4,4)
        Eh(:,5,8)     = -Vyy(:,5,5)
        Eh(:,6,2)     = one
        Eh(:,7,4)     = one
        Eh(:,8,5)     = one
        
        Fh(:,1:5,1:5) = D + im * alphap * A + im * beta * C + &
                        alphap**2 * Vxx + alphap * beta * Vxz + &
                        beta**2 * Vzz - im * omega * G
        Fh(:,6,6)     = -one
        Fh(:,7,7)     = -one
        Fh(:,8,8)     = -one
        
!==============================================================================
!       L H S
!==============================================================================

        A0 = zero
        
!.... interior (note the mapping metric)

        do idof = 1, neq
          do jdof = 1, neq
            do i = 1, ny
              do j = 1, ny
                A0(i,idof,j,jdof) =  A0(i,idof,j,jdof) + &
                                     Eh(i,idof,jdof) * deta(i) * D1(i,j)
              end do
            end do
          end do
        end do

        do idof = 1, neq
          do jdof = 1, neq
            do i = 1, ny
              A0(i,idof,i,jdof) =  A0(i,idof,i,jdof) + Fh(i,idof,jdof)
            end do
          end do
        end do

!.... freestream boundary conditions

        if (top .eq. 0) then
          i = 1
          A0(i,2,:,:) = zero
          A0(i,2,i,2) = one
          A0(i,3,:,:) = zero
          A0(i,3,i,3) = one
          A0(i,4,:,:) = zero
          A0(i,4,i,4) = one
          A0(i,5,:,:) = zero
          A0(i,5,i,5) = one
        else
          write(*,*) 'Illegal value of top = ',top
          call exit(1)
        end if
        
!.... wall boundary conditions (on the regular solution)

        i = ny

        A0(i,3,:,:) = zero
        A0(i,3,i,3) = one
        if (Navier) then
          A0(i,2,:,:) = zero
          A0(i,2,i,2) = one
          A0(i,4,:,:) = zero
          A0(i,4,i,4) = one
          if (wallt.eq.0) then
            A0(i,5,:,:) = zero
            A0(i,5,i,5) = one
          else
            write(*,"('Illegal value for wallt: ',i4)") wallt
            call exit(1)
          end if
        end if

!==============================================================================
!       R H S
!==============================================================================

        rhs = zero

!... Freestream boundary conditions

        if (top.eq.0) then
          i = 1
          rhs(i,2) = zero
          rhs(i,3) = zero
          rhs(i,4) = zero
          rhs(i,5) = zero
        else
          write(*,*) 'Illegal value of top = ',top
          call exit(1)
        end if
        
!.... Wall boundary conditions (inhomogeneous)

        i = ny

        if (wallbc .eq. 0) then        ! bump

          rhs(i,3) = -gum(i,2,2)

          if (Navier) then
            rhs(i,2) = -gum(i,1,2)
            rhs(i,4) = -gum(i,3,2)
            if (wallt.eq.0) then
              rhs(i,5) = -gtm(i,2)
            else
              write(*,"('Illegal value for wallt: ',i4)") wallt
              call exit(1)
            end if
          end if

        elseif (wallbc .eq. 1) then    ! blow

          rhs(i,2) = zero
          rhs(i,3) = one
          rhs(i,4) = zero
          
          if (wallt.eq.0) then
            rhs(i,5) = zero
          else
            write(*,"('Illegal value for wallt: ',i4)") wallt
            call exit(1)
          end if
          
        else

          write(*,"('Illegal value for wallbc: ',i4)") wallbc
          call exit(1)

        end if


!==============================================================================
!       S O L V E
!==============================================================================
#ifdef CRAY
        call CGESV( ny*neq, 1, A0, ny*neq, ipiv, rhs, ny*neq, info )
#else
        call ZGESV( ny*neq, 1, A0, ny*neq, ipiv, rhs, ny*neq, info )
#endif
        if (info .ne. 0) then
          if (info .lt. 0) then
            write(*,*) 'Argument ',abs(info),' had an illegal value'
            call exit(1)
          else
            write(*,*) 'Factor ', info , &
                       ' is exactly zero so the matrix is singular!'
            call exit(1)
          end if
        end if

        if (ialpha.eq.1) then   ! compute derivative
          duda = zero
          do idof = 1, neq
            do j = 2, ny-1
              duda(j,idof) = one / rhs(j,idof)
            end do
            if (idof.eq.1) then
              j = ny
              duda(j,idof) = one / rhs(j,idof)
            end if
          end do
        else                    ! finish derivative
          do idof = 1, neq
            do j = 2, ny-1
              duda(j,idof) = im * sqrt(two * pi) * two * dalpha / &
                              ( duda(j,idof) - one / rhs(j,idof) )
            end do
            if (idof.eq.1) then
              j = ny
              duda(j,idof) = im * sqrt(two * pi) * two * dalpha / &
                             ( duda(j,idof) - one / rhs(j,idof) )
            end if
          end do
        end if

        end do          ! ialpha

!       write(*,*) 'Solve system ',second()-cpu

        open(30,file=name,form='formatted')
        do j = 1, ny
          write (30,50) y(j), &
                        real (duda(j,1)), &
                        aimag(duda(j,1)), &
                        real (cos(lambda) * duda(j,2) + &
                              sin(lambda) * duda(j,4)), &
                        aimag(cos(lambda) * duda(j,2) + &
                              sin(lambda) * duda(j,4)), &
                        real (duda(j,3)), &
                        aimag(duda(j,3)), &
                        real (-sin(lambda) * duda(j,2) + &
                               cos(lambda) * duda(j,4)), &
                        aimag(-sin(lambda) * duda(j,2) + &
                               cos(lambda) * duda(j,4)), &
                        real (duda(j,5)), &
                        aimag(duda(j,5)), &
                        abs( cos(lambda) * duda(j,2) + &
                             sin(lambda) * duda(j,4))
        end do
        close(30)

!.... compute the disturbance kinetic energy integral
!.... Recall that the data is stored backwards

        dke = zero
        do j = ny, 2, -1
          if ( y(j-1) .le. ymax ) then
          dke = dke + pt5 * ( abs(duda(j,2))**2   + abs(duda(j,3))**2 + &
                              abs(duda(j,4))**2   + abs(duda(j-1,2))**2 + &
                              abs(duda(j-1,3))**2 + abs(duda(j-1,4))**2 ) * &
                              abs( y(j-1) - y(j) )
!         write(*,*) j, y(j), dke
          end if
        end do

!.... Since the velocities are already in the body-fitted coordinate system
!.... I will also compute receptivity efficiency based on the v_s velocity.  
!.... NOTE:  this is NOT the local streamline velocity, u_s, unless lambda is
!.... set to the local sweep angle (which it now is in mbump)

        do j = 1, ny-1
          pro(j) = cos(lambda) * duda(ny-j+1,2) + sin(lambda) * duda(ny-j+1,4)
          yl(j)  = y(ny-j+1)
        end do
#ifdef CRAY
        eff = findmax( ny-1, yl, pro, yeff )
#endif

!.... Write to file:
!.... x (or s) location, disturbance kinetic energy, sqrt of disturbance KE
!.... |Lambda| based on |v_s|_max, and arg(Lambda) based on |v_s|_max

        write(50,"(6(1pe20.13,1x))") x, dke, sqrt(dke), yeff, abs(eff), &
                                     atan2( aimag(eff), real(eff) )

!.... print out the solution

!       if (.false.) then
!         do i = 1, ny
!           write (21,50) y(i), real(rhs(i,1)), aimag(rhs(i,1))
!           write (22,50) y(i), real(rhs(i,2)), aimag(rhs(i,2))
!           write (23,50) y(i), real(rhs(i,3)), aimag(rhs(i,3))
!           write (24,50) y(i), real(rhs(i,4)), aimag(rhs(i,4))
!           write (25,50) y(i), real(rhs(i,5)), aimag(rhs(i,5))
!           write (26,50) y(i), real(rhs(i,6)), aimag(rhs(i,6))
!           write (27,50) y(i), real(rhs(i,7)), aimag(rhs(i,7))
!           write (28,50) y(i), real(rhs(i,8)), aimag(rhs(i,8))
!         end do
!         close(21)
!         close(22)
!         close(23)
!         close(24)
!         close(25)
!         close(26)
!         close(27)
!         close(28)
!       end if

!.... clean up memory

        deallocate( Eh, Fh )
        deallocate( A0, rhs, duda )
        deallocate( ipiv )

        return

!.... Some handy format statements

10      format(1p,7(e20.13,1x))
13      format(4(1pe20.13,1x))
20      format(1p,i5,5x,e13.6,5x,2(e13.6,1x),5x,2(e13.6,1x))
25      format(1p,i5,5x,e13.6,5x,2(e13.6,1x),5x,2(e13.6,1x),' <==')
50      format(1p,12(e20.13,1x))

        end

!.... I only have access to the imsl libraries on the Cray right now

#ifdef CRAY
!
! I should really use Chebyschev interpolation here, check out the nconvert
!
!=============================================================================!
!           F I N D M A X   O F   A   C O M P L E X   F U N C T I O N
!=============================================================================!
        module bspline

        integer           :: n, kord = 5
        real, allocatable :: knot(:), bs(:)

        real, external :: bsder

        end module bspline
!=============================================================================!
        function getval( nl, x, f, xmax)
        
!=============================================================================!
        use bspline
        implicit none
        
        integer :: nl, i, u, imax
        real :: getval, x(nl), f(nl), xmax
!=============================================================================!
        n = nl
        
        allocate( knot(n+kord), bs(n) )
        call BSNAK( n, x, kord, knot )
        call BSINT( n, x, f, kord, knot, bs )
        
        getval = BSDER( 0, xmax, kord, knot, n, bs )

        deallocate( knot, bs )

        return
        end     
!=============================================================================!
        function findmax( nl, x, cf, xmax  )
!
!       Find the first maximum mgnitude of a complex function
!
!=============================================================================!
        use bspline
        implicit none
        
        integer :: nl, i, u, imax
        real :: x(nl), xmax, f(nl), fmaxr, fmaxi
        complex :: findmax, cf(nl)

        real, external :: RTSAFE
        
        external fmax
!=============================================================================!
        f = abs(cf)
        n = nl
        
        allocate( knot(n+kord), bs(n) )
        call BSNAK( n, x, kord, knot )
        call BSINT( n, x, f, kord, knot, bs )
        
        do i = 1, n-1
          if ( f(i+1) .lt. f(i) ) goto 10
        end do
        write(*,*) 'Error in findmax'
        findmax = 1.0
        goto 20
        
  10    continue
        imax = i
        xmax = RTSAFE( fmax, x(i-2), x(i+2), 1.0e-14 )

        f = real( cf )
        call BSINT( n, x, f, kord, knot, bs )
        fmaxr = BSDER( 0, xmax, kord, knot, n, bs )
        f = aimag( cf )
        call BSINT( n, x, f, kord, knot, bs )
        fmaxi = BSDER( 0, xmax, kord, knot, n, bs )
        findmax = cmplx( fmaxr, fmaxi )

  20    deallocate( knot, bs )

        return
        end     
!=============================================================================!
        subroutine fmax( x, g, d )
!=============================================================================!
        use bspline
        implicit none

        real :: x, f, g, d
!=============================================================================!

!.... maximum value of function f

        f = BSDER( 0, x, kord, knot, n, bs )
        g = BSDER( 1, x, kord, knot, n, bs )
        d = BSDER( 2, x, kord, knot, n, bs )

        return
        end 
#endif
