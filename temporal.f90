!==============================================================================
        subroutine temporal(name, ind)
!==============================================================================
!
!       Driver for temporal parallel stability solver
!       This is the Chebyshev Collocation version.
!
!       I have both IMSL and LAPACK routines.  Currently using LAPACK
!
!       Revised: 7-25-96
!       Revised: 9-2-96
!==============================================================================
        use stuff
        use material
        implicit none

        integer i, j, k, ix, ier, ind
        integer i0, idof, j0, jdof
        
        real vm(ny,ndof,nx)
        real dvm(ny,ndof,nx)
        real d2vm(ny,ndof,nx)

        real eta(ny), y(ny), deta(ny), d2eta(ny)
        
        real D1(ny,ny), D2(ny,ny), Dt1(ny,ny), Dt2(ny,ny)
        
        real A(ny,ndof,ndof), B(ny,ndof,ndof),  C(ny,ndof,ndof)
        real D(ny,ndof,ndof), G(ny,ndof,ndof)
        
        real Vxx(ny,ndof,ndof), Vxy(ny,ndof,ndof), Vyy(ny,ndof,ndof)
        real Vxz(ny,ndof,ndof), Vyz(ny,ndof,ndof), Vzz(ny,ndof,ndof)
        
        real g1vm(ny,ndof),  g2vm(ny,ndof),  g3vm(ny,ndof)
        real g11vm(ny,ndof), g12vm(ny,ndof), g13vm(ny,ndof)
        real g21vm(ny,ndof), g22vm(ny,ndof), g23vm(ny,ndof)
        real g31vm(ny,ndof), g32vm(ny,ndof), g33vm(ny,ndof)
        
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
        
        real S1jj(ny),  S2jj(ny),  S3jj(ny),  S(ny,nsd,nsd)
        
        real fact
                
        complex Dh(ny,ndof,ndof), Ah(ny,ndof,ndof), Bh(ny,ndof,ndof)
        complex :: scale
        
        complex A0(ndof*ny,ndof*ny), B0(ndof*ny,ndof*ny), &
                C0(ndof*ny,ndof*ny)
        
        complex evec(ndof*ny,ndof*ny), alp(ndof*ny), bet(ndof*ny)
        complex omg(ndof*ny), cs(ndof*ny)
        real    temp1(ndof*ny), temp2(ndof*ny)
        integer index(ndof*ny)
        
        character(80) name

        complex up(ny), vp(ny)

        real, parameter :: big = 1.0e98
!       logical, parameter :: ider = .true.     ! internal derivatives
!       logical :: Navier = .true.
!       integer :: wallt  = 0                   ! wall temperature BC

!.... stuff for LAPACK eigensolver

        integer info, lwork
        complex, allocatable :: work(:)
        real, allocatable    :: rwork(:)

!.... stuff for LAPACK linear solver

        integer :: ipiv(ndof*ny)
!==============================================================================
        if (Re.ge.big .or. Re.eq.zero) then
          write(*,*) 'A T T E N T I O N:   Inviscid flow'
          Navier = .false.
        end if

!.... grid

        call sgengrid(y, eta, deta, d2eta)
        
!.... read the mean field

        if (ider) then
          call getmean(vm, y, eta, ny, ind)
        else
          call getmean2(vm, y, eta, g2vm, g22vm, ny, ind)
        end if
                    
!.... loop over the streamwise stations

        do ix = 1, nx
        
!       write (*,*) 
!       write (*,*) 'F O R M I N G   E I G E N S Y S T E M'

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
        
        rhom = vm(:,1,ix)
        u1m  = vm(:,2,ix)
        u2m  = vm(:,3,ix)
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
        
        g21vm = zero
        if (ider) then
          do idof = 1, ndof
            g22vm(:,idof) = matmul(D2,vm(:,idof,ix))
          end do
        end if
        g23vm = zero
        
        g31vm = zero
        g32vm = zero
        g33vm = zero

!.... transform the gradients to physical space

        if (ider) then
          do k = 1, ndof
            g22vm(:,k) = g22vm(:,k) * deta**2 + g2vm(:,k) * d2eta
            g2vm(:,k)  = g2vm(:,k) * deta
          end do
        end if
        
!.... write out the mean field and its gradients

        open(10,file='rho.out')
        do j = 1, ny
          write (10,13) y(j), vm(j,1,ix), g2vm(j,1), g22vm(j,1)
        end do
        close(10)
        open(10,file='u.out')
        do j = 1, ny
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
                        
 13     format(4(1pe20.13,1x))

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
        
        divum = gum(:,1,1) + gum(:,2,2) + gum(:,3,3)
        
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

!.... compute the gradient of the divergence of um

        g1divum = g11vm(:,2) + g12vm(:,3)
        g2divum = g12vm(:,2) + g22vm(:,3)
        g3divum = zero
        
!.... compute some stuff that is useful for the viscous terms

        do i = 1, nsd
          do j = 1, nsd
            S(:,i,j) = pt5 * ( gum(:,i,j) + gum(:,j,i) )
          end do
        end do

        S1jj = pt5 * ( g11vm(:,2) + g11vm(:,2) + g22vm(:,2) + &
                       g12vm(:,3) + g33vm(:,2) + g13vm(:,4) )
               
        S2jj = pt5 * ( g11vm(:,3) + g21vm(:,2) + g22vm(:,3) + &
                       g22vm(:,3) + g33vm(:,3) + g23vm(:,4) )
        
        S3jj = pt5 * ( g11vm(:,4) + g31vm(:,2) + g22vm(:,4) + &
                       g32vm(:,3) + g33vm(:,4) + g33vm(:,4) )

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
        
!.... Continuity equation

        G(:,1,1) = one
        
        A(:,1,1) = u1m
        A(:,1,2) = rhom

        B(:,1,1) = u2m
        B(:,1,3) = rhom
        
        C(:,1,1) = u3m
        C(:,1,4) = rhom
        
        D(:,1,1) = divum
        D(:,1,2) = grhom(:,1)
        D(:,1,3) = grhom(:,2)
        D(:,1,4) = grhom(:,3)
                
!.... Momentum equation -- x_1 (convective + pressure)

        G(:,2,2) = rhom
        
        A(:,2,1) = tm/(gamma * Ma**2)
        A(:,2,2) = rhom * u1m
        A(:,2,5) = rhom/(gamma * Ma**2)
        
        B(:,2,2) = rhom * u2m
        
        C(:,2,2) = rhom * u3m
        
        D(:,2,1) = u1m * gum(:,1,1) + u2m * gum(:,1,2) + u3m * gum(:,1,3) &
                   + gtm(:,1) / (gamma * Ma**2)
        D(:,2,2) = rhom * gum(:,1,1)
        D(:,2,3) = rhom * gum(:,1,2)
        D(:,2,4) = rhom * gum(:,1,3)
        D(:,2,5) = grhom(:,1) / (gamma * Ma**2)
                
!.... (viscous lambda)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,2,2) = A(:,2,2) - fact * g1lm
        A(:,2,5) = A(:,2,5) - fact * dlm * divum
        
        B(:,2,3) = B(:,2,3) - fact * g1lm
        
        C(:,2,4) = C(:,2,4) - fact * g1lm
        
        D(:,2,5) = D(:,2,5) - fact * ( g1dlm * divum + dlm * g1divum )
        
        Vxx(:,2,2) = fact * rlm
        
        Vxy(:,2,3) = fact * rlm
        
        Vxz(:,2,4) = fact * rlm
                
!.... (viscous mu)

        fact = one / Re
        
        A(:,2,2) = A(:,2,2) - fact * two * g1mu
        A(:,2,3) = A(:,2,3) - fact * g2mu
        A(:,2,4) = A(:,2,4) - fact * g3mu
        A(:,2,5) = A(:,2,5) - fact * dmu * two * S(:,1,1)
        
        B(:,2,2) = B(:,2,2) - fact * g2mu
        B(:,2,5) = B(:,2,5) - fact * dmu * two * S(:,1,2)
        
        C(:,2,2) = C(:,2,2) - fact * g3mu
        C(:,2,5) = C(:,2,5) - fact * dmu * two * S(:,1,3)
        
        D(:,2,5) = D(:,2,5) - fact * two * ( g1dmu * S(:,1,1) + &
                              g2dmu * S(:,1,2) + g3dmu * S(:,1,3) + &
                              dmu * S1jj )

        Vxx(:,2,2) = Vxx(:,2,2) + fact * two * rmu
        Vxy(:,2,3) = Vxy(:,2,3) + fact * rmu
        Vyy(:,2,2) = Vyy(:,2,2) + fact * rmu
        Vxz(:,2,4) = Vxz(:,2,4) + fact * rmu
        Vzz(:,2,2) = Vzz(:,2,2) + fact * rmu

        end if
                                             
!.... Momentum equation -- x_2 (convective + pressure)

        G(:,3,3) = rhom
        
        A(:,3,3) = rhom * u1m
        
        B(:,3,1) = tm/(gamma * Ma**2)
        B(:,3,3) = rhom * u2m
        B(:,3,5) = rhom/(gamma * Ma**2)
        
        C(:,3,3) = rhom * u3m
        
        D(:,3,1) = u1m * gum(:,2,1) + u2m * gum(:,2,2) + u3m * gum(:,2,3) &
                   + gtm(:,2) / (gamma * Ma**2)
        D(:,3,2) = rhom * gum(:,2,1)
        D(:,3,3) = rhom * gum(:,2,2)
        D(:,3,4) = rhom * gum(:,2,3)
        D(:,3,5) = grhom(:,2) / (gamma * Ma**2)

!.... (viscous lambda)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,3,2) = A(:,3,2) - fact * g2lm
        
        B(:,3,3) = B(:,3,3) - fact * g2lm
        B(:,3,5) = B(:,3,5) - fact * dlm * divum
        
        C(:,3,4) = C(:,3,4) - fact * g2lm
        
        D(:,3,5) = D(:,3,5) - fact * ( g2dlm * divum + dlm * g2divum )

        Vxy(:,3,2) = fact * rlm
        
        Vyy(:,3,3) = fact * rlm

        Vyz(:,3,4) = fact * rlm
        
!.... (viscous mu)

        fact = one / Re
        
        A(:,3,3) = A(:,3,3) - fact * g1mu
        A(:,3,5) = A(:,3,5) - fact * dmu * two * S(:,2,1)

        B(:,3,2) = B(:,3,2) - fact * g1mu
        B(:,3,3) = B(:,3,3) - fact * two * g2mu
        B(:,3,4) = B(:,3,4) - fact * g3mu
        B(:,3,5) = B(:,3,5) - fact * dmu * two * S(:,2,2)
        
        C(:,3,3) = C(:,3,3) - fact * g3mu
        C(:,3,5) = C(:,3,5) - fact * dmu * two * S(:,2,3)
        
        D(:,3,5) = D(:,3,5) - fact * two * ( g1dmu * S(:,2,1) + &
                              g2dmu * S(:,2,2) + g3dmu * S(:,2,3) + &
                              dmu * S2jj )
        Vxx(:,3,3) = Vxx(:,3,3) + fact * rmu
        Vxy(:,3,2) = Vxy(:,3,2) + fact * rmu
        Vyy(:,3,3) = Vyy(:,3,3) + fact * two * rmu
        Vyz(:,3,4) = Vyz(:,3,4) + fact * rmu
        Vzz(:,3,3) = Vzz(:,3,3) + fact * rmu

        end if

!.... Momentum equation -- x_3 (convective + pressure)

        G(:,4,4) = rhom
        
        A(:,4,4) = rhom * u1m
        
        B(:,4,4) = rhom * u2m
        
        C(:,4,1) = tm/(gamma * Ma**2)
        C(:,4,4) = rhom * u3m
        C(:,4,5) = rhom/(gamma * Ma**2)
        
        D(:,4,1) = u1m * gum(:,3,1) + u2m * gum(:,3,2) + u3m * gum(:,3,3)  &
                   + gtm(:,3) / (gamma * Ma**2)
        D(:,4,2) = rhom * gum(:,3,1)
        D(:,4,3) = rhom * gum(:,3,2)
        D(:,4,4) = rhom * gum(:,3,3)
        D(:,4,5) = grhom(:,3) / (gamma * Ma**2)
        
!.... (viscous lambda)

        if (Navier) then

        fact = rlme / (rmue * Re)
        
        A(:,4,2) = A(:,4,2) - fact * g3lm
        
        B(:,4,3) = B(:,4,3) - fact * g3lm
        
        C(:,4,4) = C(:,4,4) - fact * g3lm
        C(:,4,5) = C(:,4,5) - fact * dlm * divum
        
        D(:,4,5) = D(:,4,5) - fact * ( g3dlm * divum + dlm * g3divum )

        Vxz(:,4,2) = fact * rlm
        
        Vyz(:,4,3) = fact * rlm

        Vzz(:,4,4) = fact * rlm

!.... (viscous mu)

        fact = one / Re
        
        A(:,4,4) = A(:,4,4) - fact * g1mu
        A(:,4,5) = A(:,4,5) - fact * dmu * two * S(:,3,1)
        
        B(:,4,4) = B(:,4,4) - fact * g2mu
        B(:,4,5) = B(:,4,5) - fact * dmu * two * S(:,3,2)
        
        C(:,4,2) = C(:,4,2) - fact * g1mu
        C(:,4,3) = C(:,4,3) - fact * g2mu
        C(:,4,4) = C(:,4,4) - fact * two * g3mu
        C(:,4,5) = C(:,4,5) - fact * dmu * two * S(:,3,3)
        
        D(:,4,5) = D(:,4,5) - fact * two * ( g1dmu * S(:,3,1) + &
                              g2dmu * S(:,3,2) + g3dmu * S(:,3,3) + &
                              dmu * S3jj )
        Vxx(:,4,4) = Vxx(:,4,4) + fact * rmu
        Vyy(:,4,4) = Vyy(:,4,4) + fact * rmu
        Vxz(:,4,2) = Vxz(:,4,2) + fact * rmu
        Vyz(:,4,3) = Vyz(:,4,3) + fact * rmu
        Vzz(:,4,4) = Vzz(:,4,4) + fact * two * rmu

        end if
                              
!.... Energy equation (Advection + pressure)

        G(:,5,1) = -gamma1 * tm / gamma
        G(:,5,5) = rhom / gamma
        
        A(:,5,1) = -gamma1 * u1m * tm / gamma
        A(:,5,5) = rhom * u1m / gamma
        
        B(:,5,1) = -gamma1 * u2m * tm / gamma
        B(:,5,5) = rhom * u2m / gamma

        C(:,5,1) = -gamma1 * u3m * tm / gamma
        C(:,5,5) = rhom * u3m / gamma

        D(:,5,1) = one/gamma * ( u1m * gtm(:,1) + u2m * gtm(:,2) + &
                                 u3m * gtm(:,3) )
        D(:,5,2) = rhom * gtm(:,1) - gamma1 * Ma**2 * gpm(:,1)
        D(:,5,3) = rhom * gtm(:,2) - gamma1 * Ma**2 * gpm(:,2)
        D(:,5,4) = rhom * gtm(:,3) - gamma1 * Ma**2 * gpm(:,3)
        D(:,5,5) = -gamma1/gamma * ( u1m * grhom(:,1) + u2m * grhom(:,2) + &
                                     u3m * grhom(:,3) )

        if (Navier) then

!.... diffusion

        fact = one / (Pr * Re)
        
        A(:,5,5) = A(:,5,5) - fact * (g1con + dcon * gtm(:,1))
        B(:,5,5) = B(:,5,5) - fact * (g2con + dcon * gtm(:,2))
        C(:,5,5) = C(:,5,5) - fact * (g3con + dcon * gtm(:,3))
        D(:,5,5) = D(:,5,5) - fact * (g1dcon * gtm(:,1) + &
                              g2dcon * gtm(:,2) + g3dcon * gtm(:,3) + &
                              dcon * ( g11vm(:,5) + g22vm(:,5) + &
                                       g33vm(:,5) ) )

        Vxx(:,5,5) = fact * con
        Vyy(:,5,5) = fact * con
        Vzz(:,5,5) = fact * con
        
!.... dissipation (lambda)

        fact = gamma1 * Ma**2 * rlme / (Re * rmue)
        
        A(:,5,2) = A(:,5,2) - fact * two * rlm * divum
        B(:,5,3) = B(:,5,3) - fact * two * rlm * divum
        C(:,5,4) = C(:,5,4) - fact * two * rlm * divum
        D(:,5,5) = D(:,5,5) - fact * dlm * divum * divum
        
!.... dissipation (mu)

        fact = two * gamma1 * Ma**2 / Re
        
        A(:,5,2) = A(:,5,2) - fact * two * rmu * S(:,1,1)
        A(:,5,3) = A(:,5,3) - fact * two * rmu * S(:,2,1)
        A(:,5,4) = A(:,5,4) - fact * two * rmu * S(:,3,1)

        B(:,5,2) = B(:,5,2) - fact * two * rmu * S(:,1,2)
        B(:,5,3) = B(:,5,3) - fact * two * rmu * S(:,2,2)
        B(:,5,4) = B(:,5,4) - fact * two * rmu * S(:,3,2)

        C(:,5,2) = C(:,5,2) - fact * two * rmu * S(:,1,3)
        C(:,5,3) = C(:,5,3) - fact * two * rmu * S(:,2,3)
        C(:,5,4) = C(:,5,4) - fact * two * rmu * S(:,3,3)

        D(:,5,5) = D(:,5,5) - fact * dmu * (   S(:,1,1)**2 + &
                   S(:,1,2)**2 + S(:,1,3)**2 + S(:,2,1)**2 + &
                   S(:,2,2)**2 + S(:,2,3)**2 + S(:,3,1)**2 + &
                   S(:,3,2)**2 + S(:,3,3)**2)

        end if
        
!==============================================================================

!.... form the equations

        Dh = D + im * alpha * A + im * beta * C &
             + alpha**2 * Vxx + alpha * beta * Vxz + beta**2 * Vzz
        
!       Ah = A - two * im * alpha * Vxx - im * beta * Vxz
        
        Bh = B - im * alpha * Vxy - im * beta * Vyz
        
!.... account for the mapping

        do idof = 1, ndof
          do jdof = 1, ndof
            Bh(:,idof,jdof) = Bh(:,idof,jdof) * deta - Vyy(:,idof,jdof) * d2eta
            Vyy(:,idof,jdof) = Vyy(:,idof,jdof) * deta**2
          end do
        end do
        
!.... initialize

        A0   = zero
        B0   = zero
        evec = zero
        alp  = zero
        bet  = zero

!.... B0

        do idof = 1, ndof
          B0(idof,idof) = im * one
        end do
        
        do idof = 1, ndof
          do jdof = 1, ndof
            do i = 2, ny-1
              i0 = (i-1)*ndof
              B0(i0+idof,i0+jdof) = im * G(i,idof,jdof)
            end do
          end do
        end do

!.... wall boundary

        i = ny
        i0 = (i-1)*ndof
        do idof = 1, ndof-1
          B0(i0+idof,i0+idof) = im * one
        end do

        idof = ndof
        i0 = (i-1)*ndof
        if (wallt.eq.0) then
          B0(i0+idof,i0+idof) = im * one
        else if (wallt.eq.2) then
          do jdof = 1, ndof
            B0(i0+idof,i0+jdof) = im * G(i,idof,jdof)
          end do
        else
          write(*,"('Illegal value of wallt:  ',i4)") wallt
          stop
        end if  
        
!.... A0  (modified for Chebyshev)
        
        idof = 1
        i = 1
        i0 = (i-1)*ndof
          
        do jdof = 1, ndof
          
          do j = 1, ny
            j0 = (j-1)*ndof         
            A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * D1(i,j)
          end do
          
          j0 = (i-1)*ndof 
          A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
        
        end do
                        
        do idof = 1, ndof
          do jdof = 1, ndof

            do i = 2, ny-1
              i0 = (i-1)*ndof

              do j = 1, ny
                j0 = (j-1)*ndof
                
                A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                       Bh(i,idof,jdof) * D1(i,j) - &
                                       Vyy(i,idof,jdof) * D2(i,j)
              end do
              
              j0 = (i-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
              
            end do
            
          end do
        end do

!.... enforce continuity

        i = ny
        idof = 1                ! density equation
        i0 = (i-1)*ndof
          
        do jdof = 1, ndof
          
          do j = 1, ny
            j0 = (j-1)*ndof         
            A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + &
                                  Bh(i,idof,jdof) * D1(i,j)
          end do
          
          j0 = (i-1)*ndof 
          A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + Dh(i,idof,jdof)

        end do

        if (wallt.eq.2) then
          i = ny
          idof = ndof           ! energy equation
          i0 = (i-1)*ndof
            
          do jdof = 1, ndof-1
            
            do j = 1, ny
              j0 = (j-1)*ndof       
                  A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                          Bh(i,idof,jdof) * D1(i,j) - &
                                          Vyy(i,idof,jdof) * D1(i,j)
            end do
            
            j0 = (i-1)*ndof 
            A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
  
          end do
  
          jdof = ndof
          do j = 1, ny
            j0 = (j-1)*ndof         
                A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                        Bh(i,idof,jdof) * Dt1(i,j) - &
                                        Vyy(i,idof,jdof) * Dt2(i,j)
          end do
          j0 = (i-1)*ndof 
          A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + Dh(i,idof,jdof)
        end if

!       write (*,*)
!       write (*,*) 'S O L V I N G   E I G E N S Y S T E M'

!       call WRCRN('A0',ndof,ndof,A0,ndof*ny,0)         ! print out

!.... IMSL regular complex eigensolver (method 1)

!       call LINCG (ndof*ny, B0, ndof*ny, B0, ndof*ny)
        
!       call MCRCR (ndof*ny, ndof*ny, B0, ndof*ny, &
!                   ndof*ny, ndof*ny, A0, ndof*ny, &
!                   ndof*ny, ndof*ny, C0, ndof*ny)
        
!       call EVCCG (ndof*ny, C0, ndof*ny, omg, evec, ndof*ny)

!.... Regular complex eigensolver (method 2)

#ifdef CRAY
        call CGESV( ndof*ny, ndof*ny, B0, ndof*ny, ipiv, A0, ndof*ny, info )
#else
        call ZGESV( ndof*ny, ndof*ny, B0, ndof*ny, ipiv, A0, ndof*ny, info )
#endif
        if (info .ne. 0) then
          if (info .lt. 0) then
            write(*,*) 'Argument ',abs(info),' had an illegal value'
            stop
          else
            write(*,*) 'Factor ',info, &
                       ' is exactly zero so the matrix is singular'
            stop
          end if
        end if

!.... IMSL regular complex eigensolver

!       call EVCCG (ndof*ny, A0, ndof*ny, omg, evec, ndof*ny)

!.... LAPACK regular complex eigensolver

        lwork = 2*ndof*ny
        allocate (work(lwork), rwork(2*ndof*ny), STAT=ier)
        if (ier .ne. 0) then
          write(*,*) 'Error allocating work space'
          stop
        end if
#ifdef CRAY
        call CGEEV('N', 'V', ndof*ny, A0, ndof*ny, omg, evec, ndof*ny, &
                   evec, ndof*ny, work, lwork, rwork, info)
#else
        call ZGEEV('N', 'V', ndof*ny, A0, ndof*ny, omg, evec, ndof*ny, &
                   evec, ndof*ny, work, lwork, rwork, info)
#endif
        if (info.ne.0) then
          write(*,*) 'Error computing eigenvalues: Error # ',info
          stop
        end if
        deallocate (work, rwork)
                
!.... IMSL generalized complex eigensolver

!       CALL GVCCG (ndof*ny, A0, ndof*ny, B0, ndof*ny, &
!                   alp, bet, evec, ndof*ny)

!.... Lapack generalized complex eigensolver
!
!       lwork = 8*ndof*ny
!       allocate (work(lwork), rwork(lwork), STAT=ier)
!       if (ier .ne. 0) then
!         write(*,*) 'Error allocating work space'
!         call exit(1)
!       end if
!       
!       call CGEGV ( 'N', 'V', ndof*ny, A0, ndof*ny, B0, ndof*ny,       &
!                     alp, bet, evec, ndof*ny, evec, ndof*ny,           &
!                    work, lwork, rwork, info )
!                    
!       write (*,*) 'Info = ', info
!       
!       deallocate (work, rwork)
!
!.... compute the frequency (temporal)

!       where (bet .ne. 0) 
!         omg = alp / bet
!       elsewhere
!         omg = zero
!       end where

!.... sort the eigenvalues by the imaginary part

        do j = 1, ndof*ny
          temp2(j) = AIMAG(omg(j))
          index(j) = j
        end do
        call PIKSR2(ndof*ny, temp2, index)
        do j = 1, ndof*ny
          temp1(j) = REAL(omg(index(j)))
          A0(:,j) = evec(:,index(j))
        end do
        
        omg = cmplx(temp1,temp2)
        evec = A0

!.... end loop on x
        
        end do

!.... compute the phase speed

        cs = omg / sqrt(alpha**2 + beta**2)

!.... Scale the eigenvectors in a reasonable way

        do j = 1, ndof*ny
          scale = zero
          do i = 1, ny*ndof
            if ( abs(evec(i,j)) .gt. abs(scale) ) then
              scale = evec(i,j)
            end if
          end do
          if (scale .ne. zero) then
            do i = 1, ny*ndof
              evec(i,j) = evec(i,j) / scale
            end do
          end if
        end do
        
!.... write out the eigensystem in an unformatted output file

        open(unit=77,file=name,form='unformatted',status='unknown')
        write (77) ind, ny, ndof, itype, ievec, curve, &
                   top, wall, wallt, ider
        write (77) omega, alpha, beta, Re, Ma, Pr
        write (77) x, y, eta, deta, d2eta, yi, ymax
        write (77) omg
        write (77) evec
        close (77)

!.... output the eigenvalues and eigenfunctions

        if (.false.) then
        
        do j = 1, ndof*ny
          if ( aimag(omg(j)) .gt. zero) then
            write (*,25) j, real(omg(j)), aimag(omg(j)), &
                            real(cs(j)),  aimag(cs(j))
          else
            write (*,20) j, real(omg(j)), aimag(omg(j)), &
                            real(cs(j)),  aimag(cs(j))
          end if
        end do

 100    continue
        write (*,"(/,'Which eigenfunction ==> ',$)")
        read (*,*) j

        if ( j .lt. 0 .or. j .gt. ndof*ny ) goto 100

        if (j .ne. 0) then
          open (unit=20, file='eig.dat', form='formatted', &
                status='unknown')
          do i = 1, ny
            i0 = (i-1)*ndof
            write (20,50) eta(i), &
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

!.... compute the divergence of the eigenfunction

          do i = 1, ny
            i0 = (i-1)*ndof + 2
            up(i) = (0.0,1.0) * alpha * evec(i0,j)
            i0 = (i-1)*ndof + 3
            vp(i) = evec(i0,j) 
          end do

          vp = matmul(D1, vp)

          vp = up + vp * deta

          do i = 1, ny
            write(55,50) eta(i), real(vp(i)), aimag(vp(i))
          end do

!.... do it again

          goto 100
        end if

        end if
        
        return
        
 10     format(1p,7(e20.13,1x))
 20     format(1p,i5,1x,2(e20.13,1x),5x,2(e20.13,1x))
 25     format(1p,i5,1x,2(e20.13,1x),5x,2(e20.13,1x),' <==')
 50     format(1p,11(e20.13,1x))

        end
