!=============================================================================!
        program getevec
!
!  This program reads in a eigensystem data file and allows the user
!  to select the eigenvectors to output in ASCII format.
!
!  Usage:      getevec filename
!
!  Author:     S. Scott Collis
!
!  Revised:    6-13-97
!
!  Notes:      I have switched the order of printing for the eigenvectors
!
!  Copyright:  S. Scott Collis
!              Department of Mechanical Engineering and Materials Science
!              Rice University, MS 321
!              Houston, TX 77005-1892
!              (713) 527-8101 x3617
!=============================================================================!
        integer :: nx, ny, ndof, itype, ievec, icurve, top, wall, wallt
        real, allocatable :: y(:), eta(:), deta(:), d2eta(:)
        complex, allocatable :: eval(:), evec(:,:)
        complex :: omega, alpha, beta, k
        real :: x, Re, Ma, Pr, yi, ymax, scale
        character(80) :: base, fname
        integer :: iver, iloc, nmax
#if 0 
        integer, external :: iargc
#endif
        integer :: narg, iarg
!=============================================================================!
        narg = iargc()
        if (narg.eq.0) then
          write(*,"('Enter the Eigensystem filename ==> ',$)")
          read(*,"(a)") fname
        else
          call getarg(1,fname)
        end if

!.... open the eigenvector file

        open(unit=10,file=fname,form='unformatted',status='old',err=1000)
        read(10) nx, ny, ndof, itype, ievec, icurve, top, wall, wallt
        allocate( y(ny), eta(ny), deta(ny), d2eta(ny) )
        write(*,*) "itype = ", itype
        if (itype.eq.1) then
          nmax = ndof*ny
        else
          nmax = 2*ndof*ny
        end if
        allocate( eval(nmax), evec(nmax,nmax) )
        read(10) omega, alpha, beta, Re, Ma, Pr
        read(10) x, y, eta, deta, d2eta, yi, ymax
        read(10) eval
        if (ievec.eq.1) read(10) evec
        close(10)

!.... compute the total wave-number

        k = sqrt( alpha**2 + beta**2 )

!.... write out the eigenvalues

20      continue
        if (itype.eq.1) then
          base = 'time'
          write(*,"(/,'Temporal Eigensystem:',/)")
          write(*,"('s = ',1pe13.6,' Yi = ',1pe20.13,', Ymax = ',1pe13.6)") &
            x, yi, ymax
          write(*,"('Alpha = (',1pe13.6,',',1pe13.6,')')") alpha
          write(*,"('Beta  = (',1pe13.6,',',1pe13.6,')')") beta
          write(*,"('Re = ',1pe13.6,', Ma = ',1pe13.6,', Pr = ',1pe13.6)") &
            Re, Ma, Pr
          write(*,"(100('='))")
          write(*,"('  Index',9x,'omega_r',16x,'omega_i',18x,'c_r',20x,'c_i')")
          write(*,"(100('='))")
          do j = 1, nmax
            write (*,"(i5,4x,4(1pe21.13E3,2x))") j, real(eval(j)), &
            aimag(eval(j)), real(eval(j)/k), aimag(eval(j)/k)
          end do
        else if (itype.eq.2) then
          base = 'space'
          write(*,"(/,'Spatial Eigensystem:',/)")
          write(*,"('s = ',1pe13.6,' Yi = ',1pe20.13,', Ymax = ',1pe13.6)") &
            x, yi, ymax
          write(*,"('Omega = (',1pe13.6,',',1pe13.6,')')") omega
          write(*,"('Beta  = (',1pe13.6,',',1pe13.6,')')") beta
          write(*,"('Re = ',1pe13.6,', Ma = ',1pe13.6,', Pr = ',1pe13.6)") &
            Re, Ma, Pr
          write(*,"(100('='))")
          write(*,"('  Index',9x,'alpha_r',16x,'alpha_i',18x,'c_r',20x,'c_i')")
          write(*,"(100('='))")
          do j = 1, nmax
            if ( eval(j) .ne. 0 ) then
              write (*,"(i5,4x,4(1pe21.13E3,2x))") j, real(eval(j)), &
                aimag(eval(j)), real(omega/eval(j)), aimag(omega/eval(j))
            else
              write (*,"(i5,4x,4(1pe21.13E3,2x))") j, real(eval(j)), &
                aimag(eval(j)), 0.0, 0.0
            end if
          end do
        else
          write(*,*) 'Itype = ',itype,' is not supported'
          stop
        end if
        write(*,"(100('='),/)")
        
!.... output selected eigenvectors
          
        iver = 0
 10     continue
        write (*,"('Which eigenfunction ==> ',$)")
        read (*,*) j
        
        if ( j.eq.0 ) stop

        if ( j.eq.-1 ) goto 20  ! reprint the eigenvalues

        if ( j.lt.-1 .or. j.gt.nmax ) goto 10   ! illegal input

!.... Scale the eigenvectors in a reasonable way

        scale = 0.0
        do i = 1, ny*ndof
          if ( abs(evec(i,j)) .gt. abs(scale) ) then
            scale = evec(i,j)
          end if
        end do
        if (scale .ne. 0.0) then
          do i = 1, ny*ndof
            evec(i,j) = evec(i,j) / scale
          end do
        end if

!.... output the eigenfunction

        iver = iver + 1
        call makename(base,iver,fname)
        open (unit=20, file=fname, form='formatted', status='unknown')
        write(20,"('# Re = ',1pe13.6,', Ma = ',1pe13.6,', Pr = ',1pe13.6)") &
          Re, Ma, Pr
        if (itype.eq.1 .or. itype.eq.3) then
          write(20,"('# Omega = (',1pe21.13E3,',',1pe21.13E3,')')") eval(j)
          write(20,"('# Alpha = (',1pe21.13E3,',',1pe21.13E3,')')") alpha
          write(20,"('# Beta  = (',1pe21.13E3,',',1pe21.13E3,')')") beta
        elseif (itype.eq.2 .or. itype.eq.4) then
          write(20,"('# Omega = (',1pe21.13E3,',',1pe21.13E3,')')") omega 
          write(20,"('# Alpha = (',1pe21.13E3,',',1pe21.13E3,')')") eval(j) 
          write(20,"('# Beta  = (',1pe21.13E3,',',1pe21.13E3,')')") beta
        endif 
        do i = ny, 1, -1
!         i0 = (i+ny-1)*ndof
          i0 = (i-1)*ndof
          write (20,50) y(i), &
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
        iloc = index(fname,' ')
        write(*,"('  Eigenfunction saved in: ',a)") fname(1:iloc)
        goto 10
        
 50     format(1p,11(1pe21.13E3,1x))

        stop

1000    write(*,"('Error opening file...')")
        stop
        end

!=============================================================================!
    subroutine makename(base,iver,fname)
!
!.... put a version number on the filename
!
!=============================================================================!
      character(80) base, fname

      length = index(base,' ')
      fname = base
      if (iver .lt. 10) then
        write(fname(length:80),"('.',i1)") iver
      else if (iver .lt. 100) then
        write(fname(length:80),"('.',i2)") iver
      else
        write(fname(length:80),"('.',i3)") iver
      end if

      return
      end

