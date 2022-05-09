!==============================================================================
        subroutine mspatial(ind)
!==============================================================================
!
!       Driver for multiple spatial parallel stability solver
!
!       Data is stored in unformated eigenmode files.  Use the postprocessor
!       getab or getax to extract data.
!==============================================================================
        use stuff
        implicit none
        real :: omin, omax, oinc
        real :: bmin, bmax, binc
        real :: tmp
        integer :: i, itmp, nbody
        real,allocatable :: xb(:), delta(:)
        integer :: io, ib, no, nb, iver, ind, ind1, ind2, ind_inc, dtype
        character*80 :: base="eig ", name
!==============================================================================
        write (*,"('Enter omega_min, omega_max, omega_inc ==> ',$)")
        read (*,*) omin, omax, oinc
        write (*,"('Enter beta_min, beta_max, beta_inc ==> ',$)")
        read (*,*) bmin, bmax, binc
        write (*,"('Enter ind1, ind2 , ind_inc ==> ',$)")
        read (*,*) ind1, ind2, ind_inc
        write (*,"('Enter dtype ==> ',$)")      ! 0 = Cheby or 1 = FD
        read (*,*) dtype

!.... Use the boundary layer thickness data to scale the mesh spacing

        open(10,file='delta.dat',status='old')
        nbody = 0
 20     continue
        read(10,*,end=30) tmp
        nbody = nbody + 1
        goto 20
 30     continue
        rewind(10)
        if (ind1.lt.1 .or. ind2.gt.nbody) then
          write(*,"('ERROR: illegal index...check body.dat')")
          call exit(1)
        end if
        allocate( xb(nbody), delta(nbody) )
        do i = 1, nbody
          read(10,*) xb(i), delta(i)
        end do
        close(10)

        if (oinc .eq. zero) oinc = one
        if (binc .eq. zero) binc = one
        if (ind_inc .eq. 0) ind_inc = 1
        
        no = nint((omax-omin)/oinc)+1
        nb = nint((bmax-bmin)/binc)+1
        iver = zero

        write(*,5)
        do ind = ind1, ind2, ind_inc
          x = xb(ind)
          yi = two * delta(ind)
!         write(*,"('xb = ',1pe13.6,'   Yi = ',1pe13.6)") x, yi
          do io = 1, no
            do ib = 1, nb             
              iver = iver + 1
              omega = omin + float(io-1) * oinc
              beta  = bmin + float(ib-1) * binc
              call makename(base,iver,name)
              write(*,10) ind, iver, x, yi, real(omega), aimag(omega), &
                          real(beta), aimag(beta)
              if (dtype.eq.0) then
                call spatial(name, ind)
              else
                call fd_spatial(name, ind)
              end if
            end do
          end do
        end do

        return
  5     format(/,110('='),/, &
          'Index',' iver', &
          '      s       ', &
          '     Y_i      ', &
          '    omega_r   ', &
          '    omega_i   ', &
          '    beta_r    ', &
          '    beta_i    ', &
          /,110('='))
 10     format(2(1x,i3,1x),6(1pe13.6,1x))
        end
