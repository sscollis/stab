!==============================================================================
        subroutine stokes(ind)
!==============================================================================
!
!       Driver for Stokes flow solver.  
!       This is the staggered grid version.
!
!       Revised: 2-19-96
!==============================================================================
        use stuff
        use material
        implicit none

        integer i, j, k, ix, ier, ind
        integer i0, i1, i2, idof, j0, jdof
        
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
        
        complex A0(ndof*ny,ndof*ny), rhs(ndof*ny), sol(ndof*ny)
        
        real, parameter :: big = 1.0e98
!       logical :: Navier = .true.
!       integer :: wallt = 0

!.... stuff for LAPACK linear solver

        integer :: info, ipiv(ndof*ny)

!.... staggered grid variables

        real :: eta_s(ny-1), y_s(ny-1), deta_s(ny-1), d2eta_s(ny-1)
        real :: D1G(ny-1,ny-1), D2G(ny-1,ny-1)
        real :: GLtoG(ny-1,ny), GtoGL(ny,ny-1)
        real :: vm_s(ny-1,ndof,nx), g2vm_s(ny-1,ndof), g22vm_s(ny-1,ndof)
        
        real :: D1_GtoGL(ny,ny-1), D1G_GLtoG(ny-1,ny)   
        real :: D2_GtoGL(ny,ny-1), D2G_GLtoG(ny-1,ny)

!.... sponge variables
        
        real :: spg(ny), spg_s(ny-1), ys = 10.0, ye = 40.0, As = 0.1
        integer :: Ns = 3, ispg = 1

        character*1 :: ans
        
        complex :: u(ny), v(ny), div(ny), vor(ny), t(ny), rho(ny), p(ny)
        complex :: div_s(ny-1), vor_s(ny-1), p_s(ny-1), rho_s(ny-1)
!==============================================================================
        write (*,*) 'S T A G G E R E D   M E S H   S O L V E R'

        write(*,"('Enter ispg ==> ',$)")
        read(*,*) ispg
        if (ispg.ne.0) then
          write(*,"('Enter As ==> ',$)")
          read(*,*) As
          write(*,"('Enter ys ==> ',$)")
          read(*,*) ys
        end if
          
        if (Re.ge.big .or. Re.eq.zero) then
          write(*,*) 'A T T E N T I O N:   Inviscid flow'
          Navier = .false.
        end if

!.... generate the normal and staggered meshes

        write (*,*) 
        write (*,*) 'G E N E R A T I N G   G R I D ( S )'
        call sgengrid_s(y, eta, deta, d2eta, y_s, eta_s, deta_s, d2eta_s)
        
!.... make the sponges

        ye = y(1)

        spg = zero
        spg_s = zero
        if (ispg.eq.1) then
          do j = 1, ny
            if ( y(j) .ge. ys ) then
              spg(j) = As * ( ( y(j) - ys ) / ( ye - ys ) )**Ns
            else
              spg(j) = zero
            end if
            write (30,10) y(j), eta(j), spg(j)
          end do
          do j = 1, ny-1
            if ( y_s(j) .ge. ys ) then
              spg_s(j) = As * ( ( y_s(j) - ys ) / ( ye - ys ) )**Ns
            else
              spg_s(j) = zero
            end if
            write (31,10) y_s(j), eta_s(j), spg_s(j)
          end do
        end if
 
!.... generate the mean field

        write(*,"('(P)rofile, (U)niform, (Z)ero mean flow ==> ',$)")
        read(*,"(a1)") ans

        if (ans.eq.'P' .or. ans.eq.'p') then
          if (ider) then
            write(*,*) 'Using internal derivatives'
            call getmean(vm, y, eta, ny, ind)
            call getmean(vm_s, y_s, eta_s, ny-1, ind)
          else
            call getmean2(vm, y, eta, g2vm, g22vm, ny, ind)
            call getmean2(vm_s, y_s, eta_s, g2vm_s, g22vm_s, ny-1, ind)
          end if
        else if (ans.eq.'U' .or. ans.eq.'u') then
          write(*,*) 'W A R N I N G:   Uniform mean flow'
          vm(:,1,:) = one
          vm(:,2,:) = one
          vm(:,3,:) = zero
          vm(:,4,:) = zero
          vm(:,5,:) = one
          g2vm = zero
          g22vm = zero
          vm_s(:,1,:) = one
          vm_s(:,2,:) = one
          vm_s(:,3,:) = zero
          vm_s(:,4,:) = zero
          vm_s(:,5,:) = one
          g2vm_s = zero
          g22vm_s = zero
        else
          write(*,*) 'W A R N I N G:   Zero mean flow'
          vm(:,1,:) = one
          vm(:,2,:) = zero
          vm(:,3,:) = zero
          vm(:,4,:) = zero
          vm(:,5,:) = one
          g2vm = zero
          g22vm = zero
          vm_s(:,1,:) = one
          vm_s(:,2,:) = zero
          vm_s(:,3,:) = zero
          vm_s(:,4,:) = zero
          vm_s(:,5,:) = one
          g2vm_s = zero
          g22vm_s = zero
        end if
 
!.... loop over the streamwise stations

        do ix = 1, nx
        
        write (*,*) 
        write (*,*) 'F O R M I N G   E Q U A T I O N S'

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

!.... Compute the Chebyshev-Gauss-Lobatto derivative matrices

        call chebyd(D1, ny-1)   ! use ny-1 since this routine uses 0:ny
        D2  = matmul(D1, D1)    ! compute the second derivative
        Dt1 = D1                ! Derivative operator for temperature
        Dt1(ny,:) = zero        !   adiabatic wall
        Dt2 = matmul(D1, Dt1)   ! Second Derivative operator for temperature

!.... Compute the Chebyshev-Gauss matrices

        call chebydG(D1G, ny-2) ! derivative at the Gauss points
        D2G = matmul(D1G, D1G)  ! second derivative at Gauss points

        call chebyint( GLtoG, GtoGL, ny-1 )
        
        D1_GtoGL  = matmul( D1, GtoGL )
        D1G_GLtoG = matmul( D1G, GLtoG )
        D2_GtoGL  = matmul( D2, GtoGL )
        D2G_GLtoG = matmul( D2G, GLtoG )

!.... Compute derivatives of mean field using the appropriate difference scheme

        g1vm  = zero
        if (ider) then
          do idof = 1, ndof
            g2vm(:,idof)   = matmul(D1,vm(:,idof,ix))
            g2vm_s(:,idof) = matmul(D1G,vm_s(:,idof,ix))
          end do
        end if
        g3vm  = zero
        
        g11vm = zero
        g12vm = zero
        g13vm = zero
        
        g21vm = zero
        if (ider) then
          do idof = 1, ndof
            g22vm(:,idof)   = matmul(D2,vm(:,idof,ix))
            g22vm_s(:,idof) = matmul(D2G,vm_s(:,idof,ix))
          end do
        end if
        g23vm = zero
        
        g31vm = zero
        g32vm = zero
        g33vm = zero

!.... transform the gradients to physical space

        if (ider) then
          do k = 1, ndof
            g22vm(:,k)   = g22vm(:,k)   * deta**2   + g2vm(:,k)   * d2eta
            g22vm_s(:,k) = g22vm_s(:,k) * deta_s**2 + g2vm_s(:,k) * d2eta_s
            g2vm(:,k)    = g2vm(:,k)   * deta
            g2vm_s(:,k)  = g2vm_s(:,k) * deta_s
          end do
        end if
        
!.... write out the mean field and its gradients

        open(10,file='rho.out')
        do j = 1, ny-1
          write (10,13) y_s(j), vm_s(j,1,ix), g2vm_s(j,1), g22vm_s(j,1)
        end do
        close(10)
        open(10,file='u.out')
        do j = 1, ny
          write (10,13) y(j), vm(j,2,ix), g2vm(j,2), g22vm(j,2)
        end do
        close(10)
        open(10,file='v.out')
        do j = 1, ny
          write (10,13) y(j), vm(j,3,ix), g2vm(j,3), g22vm(j,3)
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
        
!.... Continuity equation (now on the staggered mesh)

        do j = 1, ny-1
          G(j,1,1) = one
          
          A(j,1,1) = vm_s(j,2,ix)       !       u1m
          A(j,1,2) = vm_s(j,1,ix)       !       rhom
  
          B(j,1,1) = vm_s(j,3,ix)       !       u2m
          B(j,1,3) = vm_s(j,1,ix)       !       rhom
          
          C(j,1,1) = vm_s(j,4,ix)       !       u3m
          C(j,1,4) = vm_s(j,1,ix)       !       rhom
          
          D(j,1,1) = g2vm_s(j,3)        !       divum
          D(j,1,2) = zero               !       grhom(j,1)
          D(j,1,3) = g2vm_s(j,1)        !       grhom(j,2)
          D(j,1,4) = zero               !       grhom(j,3)
        end do
                
!.... Momentum equation -- x_1 (convective + pressure)

        G(:,2,2) = rhom
        
        A(:,2,1) = tm/(gamma * Ma**2)
        A(:,2,2) = rhom * u1m
        A(:,2,5) = rhom/(gamma * Ma**2)
        
        B(:,2,2) = rhom * u2m
        
        C(:,2,2) = rhom * u3m
        
        D(:,2,1) = u1m * gum(:,1,1) + u2m * gum(:,1,2) + u3m * gum(:,1,3)  &
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
        
        D(:,3,1) = u1m * gum(:,2,1) + u2m * gum(:,2,2) + u3m * gum(:,2,3)  &
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
        
        Ah = A - two * im * alpha * Vxx - im * beta * Vxz
        
        Bh = B - im * alpha * Vxy - im * beta * Vyz
        
!.... account for the mapping

        idof = 1
        do jdof = 1, ndof
          do j = 1, ny-1
            Bh(j,idof,jdof)  = Bh(j,idof,jdof)  * deta_s(j) - &
                               Vyy(j,idof,jdof) * d2eta_s(j)
            Vyy(j,idof,jdof) = Vyy(j,idof,jdof) * deta_s(j)**2
          end do
        end do

        do idof = 2, ndof
          do jdof = 1, ndof
            Bh(:,idof,jdof) = Bh(:,idof,jdof) * deta - Vyy(:,idof,jdof) * d2eta
            Vyy(:,idof,jdof) = Vyy(:,idof,jdof) * deta**2
          end do
        end do

!==============================================================================
!       L H S
!==============================================================================

        A0 = zero

!.... interior

        idof = 1
        jdof = 1
        do i = 1, ny-1
          i0 = (i-1)*ndof
          do j = 1, ny-1
            j0 = (j-1)*ndof
            A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                   Bh(i,idof,jdof) * D1G(i,j) - &
                                   Vyy(i,idof,jdof) * D2G(i,j)
          end do
          j0 = (i-1)*ndof
          A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + Dh(i,idof,jdof) - &
                                 im * omega * G(i,idof,jdof) + spg_s(i)
        end do

        idof = 1
        do jdof = 2, ndof
          do i = 1, ny-1
            i0 = (i-1)*ndof
            do j = 1, ny
              j0 = (j-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Bh(i,idof,jdof) * D1G_GLtoG(i,j) - &
                                     Vyy(i,idof,jdof) * D2G_GLtoG(i,j)

              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Dh(i,idof,jdof) * GLtoG(i,j) - &
                                     im * omega * G(i,idof,jdof) * GLtoG(i,j) 
            end do
          end do
        end do

        jdof = 1
        do idof = 2, ndof
          do i = 1, ny
            i0 = (i-1)*ndof
            do j = 1, ny-1
              j0 = (j-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Bh(i,idof,jdof) * D1_GtoGL(i,j) - &
                                     Vyy(i,idof,jdof) * D2_GtoGL(i,j)

              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Dh(i,idof,jdof) * GtoGL(i,j) - &
                                     im * omega * G(i,idof,jdof) * GtoGL(i,j)
            end do
          end do
        end do

        do idof = 2, ndof
          do jdof = 2, ndof
            do i = 1, ny
              i0 = (i-1)*ndof
              do j = 1, ny
                j0 = (j-1)*ndof
                A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                       Bh(i,idof,jdof) * D1(i,j) - &
                                       Vyy(i,idof,jdof) * D2(i,j)
              end do
              j0 = (i-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + Dh(i,idof,jdof) - &
                                     im * omega * G(i,idof,jdof)
              if (jdof.eq.idof) then
                A0(i0+idof,j0+jdof) = A0(i0+idof,j0+jdof) + spg(i)
              end if
            end do
          end do
        end do

!.... Account for the extra pressure node

        i = ny
        i0 = (i-1)*ndof
        
        A0(i0+1,i0+1) = one

!.... Free-stream boundary conditions

        i = 1
        i0 = (i-1)*ndof

        do idof = 2, ndof
          A0(i0+idof,:)       = zero
          A0(i0+idof,i0+idof) = one
        end do
        
!.... First-order extrapolation of v-velocity in the far-field

        if (.true.) then
          i = 1
          idof = 3
          
          i0 = (i-1)*ndof
          A0(i0+idof,:) = zero
          A0(i0+idof,i0+idof) = one
          
          i1 = i*ndof
          A0(i0+idof,i1+idof) = -two
          
          i2 = (i+1)*ndof
          A0(i0+idof,i2+idof) = one
        end if

!.... Wall boundary conditions

        i = ny
        i0 = (i-1)*ndof

        A0(i0+3,:)    = zero
        A0(i0+3,i0+3) = one

        if (Navier) then
          A0(i0+2,:)    = zero
          A0(i0+2,i0+2) = one
  
          A0(i0+4,:)    = zero
          A0(i0+4,i0+4) = one

          idof = ndof                   ! energy equation
  
          write(*,"('Enter wallt ==> ',$)") 
          read(*,*) wallt
          
          if (wallt.eq.0) then          ! direct enforcement of t = zero

            A0(i0+idof,:)       = zero
            A0(i0+idof,i0+idof) = one

          elseif (wallt.eq.1) then      ! modified energy equation

            A0(i0+idof,:) = zero

            jdof = 1
            do j = 1, ny-1
              j0 = (j-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Bh(i,idof,jdof) * D1_GtoGL(i,j) - &
                                     Vyy(i,idof,jdof) * D2_GtoGL(i,j)

              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Dh(i,idof,jdof) * GtoGL(i,j) - &
                                     im * omega * G(i,idof,jdof) * &
                                     GtoGL(i,j) + spg(i) * GtoGL(i,j)
            end do

            do jdof = 2, ndof-1
              do j = 1, ny
                j0 = (j-1)*ndof
                A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                       Bh(i,idof,jdof) * D1(i,j) - &
                                       Vyy(i,idof,jdof) * D2(i,j)
              end do
              j0 = (i-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + Dh(i,idof,jdof) - &
                                     im * omega * G(i,idof,jdof) + spg(i)
            end do

            jdof = ndof
            do j = 1, ny
              j0 = (j-1)*ndof
              A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + &
                                     Bh(i,idof,jdof) * Dt1(i,j) - &
                                     Vyy(i,idof,jdof) * Dt2(i,j)
            end do
            j0 = (i-1)*ndof
            A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + Dh(i,idof,jdof) - &
                                   im * omega * G(i,idof,jdof) + spg(i)

          elseif (wallt.eq.2) then      ! direct enforcement of dt/dy = zero

            A0(i0+idof,:) = zero
            jdof = ndof
            do j = 1, ny
              j0 = (j-1)*ndof       
                  A0(i0+idof,j0+jdof) =  A0(i0+idof,j0+jdof) + D1(i,j)
            end do
        
          end if        ! wallt
          
        end if          ! Navier

!==============================================================================
!       R H S
!==============================================================================

        rhs = zero

        if (ispg.eq.1) then
          do i = 1, ny-1
            i0 = (i-1)*ndof
            rhs(i0+1) = spg_s(i) * Ma**2
          end do
          do i = 1, ny
            i0 = (i-1)*ndof
            rhs(i0+2) = spg(i) * Ma
            rhs(i0+3) = spg(i) * zero
            rhs(i0+4) = spg(i) * zero
            rhs(i0+5) = spg(i) * gamma1 * Ma**2
          end do
        end if

!... Freestream boundary conditions

        i = 1
        i0 = (i-1)*ndof

        rhs(i0+3) = zero

!.... Ymax = 40, zero flow, M=0.1, R=2200, t=0

!       rhs(i0+3) = (-7.821309E-05,7.821309E-05)

!.... Ymax=40, zero flow, M=0.1, R=2200, dtdy=0

!       rhs(i0+3) = (-5.2991378691201E-05,5.2936951723434E-05)

!.... Ymax=40, BL profile, M=0.1, R=2200, dtdy=0

!       rhs(i0+3) = (-4.0342853313555E-05,2.4520137783262E-04)

        rhs(i0+2) = Ma
        rhs(i0+4) = zero
        rhs(i0+5) = gamma1 * Ma**2      

!.... Wall boundary conditions

        i = ny
        i0 = (i-1)*ndof

        rhs(i0+3) = zero
        
!       rhs(i0+3) = -(-1.9711074057057E-03,1.9553900205955E-03) - &
!                     ( 8.4792985940422E-06,3.6110911422545E-08) - &
!                     (-1.8228275321512E-08,-1.8393903194749E-08)

        if (Navier) then
          rhs(i0+2) = zero
          rhs(i0+4) = zero
          rhs(i0+5) = zero
        end if

!==============================================================================
!       S O L V E
!==============================================================================

        write (*,*)
        write (*,*) 'S O L V I N G   E Q U A T I O N S'

        if (.false.) then
        
!.... IMSL complex linear solver

        call LSACG (ndof*ny, A0, ndof*ny, rhs, 1, sol)
        do i = 1, ndof*ny
          rhs(i) = sol(i)
        end do
        
        else
        
!.... LAPACK complex linear solver

#ifdef CRAY
        call CGESV( ndof*ny, 1, A0, ndof*ny, ipiv, rhs, ndof*ny, info )
#else
        call ZGESV( ndof*ny, 1, A0, ndof*ny, ipiv, rhs, ndof*ny, info )
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

        end if

!.... print out the solution

        do i = 1, ny-1
          i0 = (i-1)*ndof
          write (21,50) y_s(i), &
                        real(rhs(i0+1)), aimag(rhs(i0+1))
        end do
        do i = 1, ny
          i0 = (i-1)*ndof
          write (22,50) y(i), &
                        real(rhs(i0+2)), aimag(rhs(i0+2))
        end do
        do i = 1, ny
          i0 = (i-1)*ndof
          write (23,50) y(i), &
                        real(rhs(i0+3)), aimag(rhs(i0+3))
        end do
        do i = 1, ny
          i0 = (i-1)*ndof
          write (24,50) y(i), &
                        real(rhs(i0+4)), aimag(rhs(i0+4))
        end do
        do i = 1, ny
          i0 = (i-1)*ndof
          write (25,50) y(i), &
                        real(rhs(i0+5)), aimag(rhs(i0+5))
        end do

!.... compute divergence and vorticity
        
        do i = 1, ny
          i0 = (i-1)*ndof
          u(i) = rhs(i0+2)
          v(i) = rhs(i0+3)
          t(i) = rhs(i0+5)
        end do
        
        do i = 1, ny-1
          i0 = (i-1)*ndof
          rho_s(i) = rhs(i0+1)
        end do
        rho(:) = matmul( GtoGL, rho_s )

        div(:) = im * alpha * u(:) + deta(:) * matmul( D1, v(:) )
        vor(:) = deta(:) * matmul( D1, u(:) ) - im * alpha * v(:)
        p(:)   = one/(gamma*Ma**2) * ( tm(:) * rho(:) + rhom(:) * t(:) )
        
        do i = 1, ny
          write (26,50) y(i), &
                        real(div(i)), aimag(div(i)), &
                        real(vor(i)), aimag(vor(i)), &
                        real(p(i)),   aimag(p(i))
        end do

        div_s(:) = im * alpha * matmul( GLtoG, u ) + &
                   deta_s(:) * matmul( D1G_GLtoG, v )
        vor_s(:) = deta_s(:) * matmul( D1G_GLtoG, u ) - im * alpha * &
                   matmul( GLtoG, v )
        p_s(:) = one/(gamma*Ma**2) * (  matmul( GLtoG, tm ) * rho_s(:) + &
                 matmul( GLtoG, rhom * t ) )

        do i = 1, ny-1
          write(27,50) y_s(i), &
                       real(div_s(i)), aimag(div_s(i)), &
                       real(vor_s(i)), aimag(vor_s(i)), &
                       real(p_s(i)), aimag(p_s(i))
        end do
        
        end do          ! loop on ix
        
        return
        
 10     format(1p,7(e14.6,1x))
 20     format(1p,i5,1x,2(e20.13,1x),5x,2(e20.13,1x))
 25     format(1p,i5,1x,2(e20.13,1x),5x,2(e20.13,1x),' <==')
 50     format(1p,11(e20.13,1x))

        end
