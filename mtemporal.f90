!==============================================================================
        subroutine mtemporal(ind)
!==============================================================================
!
!       Driver for multiple temporal parallel stability solver
!
!       Can only do alpha and beta sweeps right now.
!
!       S. Scott Collis
!
!       Revised: 6-14-97
!==============================================================================
        use stuff
        implicit none
        real :: amin, amax, ainc
        real :: bmin, bmax, binc
        integer :: ia, ib, na, nb, iver, ind
        character*80 :: base="eig ", name
!==============================================================================
        write (*,"('Enter alpha_min, alpha_max, alpha_inc ==> ',$)")
        read (*,*) amin, amax, ainc
        write (*,"('Enter beta_min, beta_max, beta_inc ==> ',$)")
        read (*,*) bmin, bmax, binc

        na = max( nint((amax-amin)/ainc), 1 )
        nb = max( nint((bmax-bmin)/binc), 1 )
        iver = 0
        
        do ia = 1, na
          do ib = 1, nb
            iver = iver + 1
            alpha = amin + float(ia-1) * ainc
            beta  = bmin + float(ib-1) * binc
            call makename(base,iver,name)
            write(*,10) ia, ib, real(alpha), aimag(alpha), &
                        real(beta), aimag(beta)
            call temporal(name, ind)
          end do
        end do
        
        return
  5     format(/,70('='),/, &
          '  ia ','  ib ', &
          '    alpha_r   ', &
          '    alpha_i   ', &
          '    beta_r    ', &
          '    beta_i    ', &
          /,70('='))
 10     format(2(1x,i3,1x),4(1pe13.6,1x))
        end

!=============================================================================!
        subroutine makename(base,iver,fname)
!
!.... put a version number on a filename
!
!=============================================================================!
        character*80 base, fname
  
        length = index(base,' ')
        fname = base
        if (iver .lt. 10) then
          write(fname(length:80),"('.',i1)") iver
        else if (iver .lt. 100) then
          write(fname(length:80),"('.',i2)") iver
        else if (iver .lt. 1000) then
          write(fname(length:80),"('.',i3)") iver
        else if (iver .lt. 10000) then
          write(fname(length:80),"('.',i4)") iver
        else
          write(*,*) 'Error in MakeName:  iver too large'
          call exit(1)
        end if
  
        return
        end
