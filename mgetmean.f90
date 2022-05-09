!=============================================================================!
        subroutine mgetmean( vm, y, eta)
!
!  Mixing layer version for Ted
!
!  Input the mean flow profile assuming parallel flow 
!  (i.e. the v-velocity is forced to zero)
!
!=============================================================================!
        use stuff
        implicit none

        real vm(ny,ndof,nx), eta(ny), y(ny), ymaxm
        real, allocatable :: ym(:,:), vt(:,:,:), vs(:,:,:)
        
        integer i, j, k, nxm, nym, ndofm, ier
!=============================================================================!

!.... read the mean field and spline to the new grid

        open (unit=10, file='profile.dat', form='formatted', status='unknown')
        
        read (10,*) nxm, nym, ndofm, ymaxm
        
        if (ndofm .ne. ndof) then
          write (*,*) 'NDOF in mean field is incorrect ',ndofm
          call exit(1)
        end if
        
        allocate( ym(nym,nxm), vt(nym,nxm,ndof), vs(nym,nxm,ndof), STAT=ier )
        if (ier .ne. 0) then
          write(*,*) 'Error allocating mean field'
          call exit(1)
        end if
        
        do i = 1, nxm
          do j = 1, nym
            read (10,*) ym(j,i), vt(j,i,1), vt(j,i,2), vt(j,i,5)
            vt(j,i,3) = zero
            vt(j,i,4) = zero
            vt(j,i,5) = vt(j,i,5) * (gamma1 * Ma**2)
          end do
          ym(j,i) = ym(j,i)
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        
        close (10)
        
!.... Evaluate the mean field on the disturbance grid

        do i = 1, nxm
          do j = 1, ny
            if ( abs(y(j)) .le. ymaxm) then
              call SPEVAL(nym, ym(1,i), vt(1,i,1), vs(1,i,1), y(j), vm(j,1,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,2), vs(1,i,2), y(j), vm(j,2,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,3), vs(1,i,3), y(j), vm(j,3,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,4), vs(1,i,4), y(j), vm(j,4,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,5), vs(1,i,5), y(j), vm(j,5,i))
            else if (y(j) .lt. -ymaxm) then
              vm(j,1,i) = vt(nym,i,1)
              vm(j,2,i) = vt(nym,i,2)
              vm(j,3,i) = vt(nym,i,3)
              vm(j,4,i) = vt(nym,i,4)
              vm(j,5,i) = vt(nym,i,5)
            else if (y(j) .gt. ymaxm) then
              vm(j,1,i) = vt(1,i,1)
              vm(j,2,i) = vt(1,i,2)
              vm(j,3,i) = vt(1,i,3)
              vm(j,4,i) = vt(1,i,4)
              vm(j,5,i) = vt(1,i,5)
            end if
!           write (69,10) eta(j), (vm(j,k,i), k = 1, ndof)
          end do
        end do
        
        deallocate( ym, vt, vs )
        
        return
10      format(8(1pe13.6,1x))
        end

