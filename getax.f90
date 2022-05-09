!=============================================================================!
        program getalpha_x
!
!  This program reads in a eigensystem data file and allows the user
!  to select the eigenvectors to output in ASCII format.
!
!  This version is for streamwise sweeps and uses linear extrapolation to 
!  find alpha.
!
!  Author:     S. Scott Collis
!
!  Revised:    6-13-97
!
!  Copyright:  S. Scott Collis
!              Department of Mechanical Engineering and Materials Science
!              Rice University, MS 321
!              Houston, TX 77005-1892
!              (713) 527-8101 x3617
!=============================================================================!
        implicit none
        
        integer :: nx, ny, ndof, itype, ievec, icurve, top, wall, wallt, ider
        real, allocatable :: y(:), eta(:), deta(:), d2eta(:)
        complex, allocatable :: eval(:), evec(:,:)
        complex :: omega, alpha, beta, k
        real :: x, Re, Ma, Pr, Yi, Ymax
        character*80 :: base, fname
        integer :: j, iloc, nmax, icount
        integer :: iver, iver1, iver2, iver_inc
        integer :: ind, ind1, ind2, ind_inc
        integer :: i, nbody, itmp, jloc
        real :: error
        real, allocatable :: xb(:)
        
        complex :: eval1, eval2, evalg
        complex :: parm1, parm2
!=============================================================================!
        
        write (*,"('Enter iver1, iver2, iver_inc ==> ',$)")
        read (*,*) iver1, iver2, iver_inc

!.... open the body.dat file to get the streamwise locations

        open(10,file='body.dat',status='old')
        nbody = 0
 20     continue
        read(10,*,end=30) itmp
        nbody = nbody + 1
        goto 20
 30     continue
        rewind(10)
        allocate( xb(nbody) )
        do i = 1, nbody
          read(10,*) itmp, xb(i)
        end do
        close(10)

        icount = 0
        base = 'eig '
        do iver = iver1, iver2, iver_inc

        icount = icount + 1
        
!.... open the eigensystem file

        call makename(base,iver,fname)
        open(unit=10,file=fname,form='unformatted',status='unknown')
        read(10) nx, ny, ndof, itype, ievec, icurve, top, wall, wallt, ider
        if (.not. allocated( y ) ) then
          allocate( y(ny), eta(ny), deta(ny), d2eta(ny) )
        end if
        if (itype.eq.1) then
          nmax = ndof*ny
        else
          nmax = 2*ndof*ny
        end if
        if (.not. allocated( eval ) ) allocate( eval(nmax) )
        read(10) omega, alpha, beta, Re, Ma, Pr
        read(10) x, y, eta, deta, d2eta, yi, ymax
        read(10) eval
        close(10)
          
!.... compute the total wave-number

        k = sqrt( alpha**2 + beta**2 )

!.... write out the eigenvalues

        if (icount .eq. 1 .or. icount .eq. 2) then
        
        if (itype.eq.1) then
          write(*,"(/,'Temporal Eigensystem:',/)")
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
          write(*,"(/,'Spatial Eigensystem:',/)")
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
        
        write (*,"('Which eigenvalue ==> ',$)")
        read (*,*) jloc

        if (icount.eq.2) then
          eval2 = eval1
          parm2 = parm1
        end if
        
        eval1 = eval(jloc)
        parm1 = x

        open(20,file='stab.out')
        write(20,50) nx, x, real(eval(jloc)), aimag(eval(jloc))
        
        else    ! extrapolate for the eigenvalue
        
          evalg = eval1 + (eval2-eval1)/(parm2-parm1)*(x-parm1)

          error = 1.0e99
          do j = 1, nmax
            if ( abs(eval(j)-evalg) .lt. error ) then
              jloc = j
              error = abs(eval(j)-evalg)
            end if
          end do
          
          write(20,50) nx, x, real(eval(jloc)), aimag(eval(jloc))

          eval2 = eval1
          parm2 = parm1
          eval1 = eval(jloc)
          parm1 = x
          
        end if          ! iver
        
        end do          ! ind
        
        close(20)

        stop
 50     format(i4,3(1x,1pe20.13))
 60     format(3(1x,1pe20.13))
        end

!=============================================================================!
    subroutine makename(base,iver,fname)
!
!.... put a version number on the filename
!
!=============================================================================!
      character*80 base, fname

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

