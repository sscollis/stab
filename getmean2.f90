!=============================================================================!
        subroutine getmean2( vm, y, eta, g2vm, g22vm, lny, iver)
!
!  Input the mean flow profile assuming parallel flow.
!  (i.e. v-velocity is forced to zero).
!
!  This version of getmean also reads in the first and second derivatives
!  of the boundary-layer profiles.
!=============================================================================!
        use stuff
        implicit none

        integer :: lny
        real :: vm(lny,ndof,nx), eta(lny), y(lny), ymaxm
        real :: g2vm(lny,ndof,nx), g22vm(lny,ndof,nx)
        real, allocatable :: ym(:,:), vt(:,:,:), vs(:,:,:)
        
        integer :: i, j, k, nxm, nym, ndofm, ier, iver
        real :: tmp

        character(80) :: base, fname
!=============================================================================!
        
!.... read the mean field and spline to the new grid

        base = 'profile'
        call makename(base,iver,fname)
!       write(*,"(/,'Reading mean flow from:  ',a)") fname
        open (unit=10, file=fname, form='formatted', status='old',err=1000)
        
!       read (10,*) nxm, nym, ndofm, ymaxm
!       if (ndofm .ne. ndof) then
!         write (*,*) 'NDOF in mean field is incorrect ',ndofm
!         call exit(1)
!       end if
        
        nxm = 1
        nym = 0
 20     continue
        read(10,*,end=30) tmp
        nym = nym + 1
        goto 20
 30     continue
        rewind(10)

        allocate( ym(nym,nxm), vt(nym,nxm,ndof), vs(nym,nxm,ndof), &
                  STAT=ier )
        if (ier .ne. 0) then
          write(*,*) 'Error allocating mean field'
          call exit(1)
        end if
        
        do i = 1, nxm
          do j = 1, nym
            read (10,*) ym(j,i), (vt(j,i,k),k=1,ndof)
            vt(j,i,3) = zero ! parallel flow assumption
          end do
          ym(j,i) = ym(j,i)
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        
        close (10)

        ymaxm = ym(nym,1)
        
!.... Evaluate the mean field on the disturbance grid

        do i = 1, nxm
          do j = 1, lny
            if (y(j) .le. ymaxm) then
              call SPEVAL(nym, ym(1,i), vt(1,i,1), vs(1,i,1), y(j), vm(j,1,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,2), vs(1,i,2), y(j), vm(j,2,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,3), vs(1,i,3), y(j), vm(j,3,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,4), vs(1,i,4), y(j), vm(j,4,i))
              call SPEVAL(nym, ym(1,i), vt(1,i,5), vs(1,i,5), y(j), vm(j,5,i))
            else
              vm(j,1,i) = vt(nym,i,1)
              vm(j,2,i) = vt(nym,i,2)
              vm(j,3,i) = vt(nym,i,3)
              vm(j,4,i) = vt(nym,i,4)
              vm(j,5,i) = vt(nym,i,5)
            end if
!           write (69,10) eta(j), (vm(j,k,i), k = 1, ndof)
          end do
        end do

!.... read the first derivative and spline to the new grid

        base = 'first'
        call makename(base,iver,fname)
        open (unit=10, file=fname, form='formatted', status='old',err=1000)

        do i = 1, nxm
          do j = 1, nym
            read (10,*) ym(j,i), (vt(j,i,k),k=1,ndof)
            vt(j,i,3) = zero ! parallel flow assumption
          end do
          ym(j,i) = ym(j,i)
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        close (10)
        
!.... Evaluate the mean derivative field on the disturbance grid

        do i = 1, nxm
          do j = 1, lny
            if (y(j) .le. ymaxm) then
              call SPEVAL(nym,ym(1,i), vt(1,i,1), vs(1,i,1), y(j), g2vm(j,1,i))
              call SPEVAL(nym,ym(1,i), vt(1,i,2), vs(1,i,2), y(j), g2vm(j,2,i))
              call SPEVAL(nym,ym(1,i), vt(1,i,3), vs(1,i,3), y(j), g2vm(j,3,i))
              call SPEVAL(nym,ym(1,i), vt(1,i,4), vs(1,i,4), y(j), g2vm(j,4,i))
              call SPEVAL(nym,ym(1,i), vt(1,i,5), vs(1,i,5), y(j), g2vm(j,5,i))
            else
              g2vm(j,1,i) = vt(nym,i,1)
              g2vm(j,2,i) = vt(nym,i,2)
              g2vm(j,3,i) = vt(nym,i,3)
              g2vm(j,4,i) = vt(nym,i,4)
              g2vm(j,5,i) = vt(nym,i,5)
            end if
          end do
        end do

!.... read the second derivative and spline to the new grid

        base = 'second'
        call makename(base,iver,fname)
        open (unit=10, file=fname, form='formatted', status='old',err=1000)

        do i = 1, nxm
          do j = 1, nym
            read (10,*) ym(j,i), (vt(j,i,k),k=1,ndof)
            vt(j,i,3) = zero ! parallel flow assumption
          end do
          ym(j,i) = ym(j,i)
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        close (10)
        
!.... Evaluate the mean derivative field on the disturbance grid

        do i = 1, nxm
          do j = 1, lny
            if (y(j) .le. ymaxm) then
              call SPEVAL(nym,ym(1,i),vt(1,i,1), vs(1,i,1), y(j), g22vm(j,1,i))
              call SPEVAL(nym,ym(1,i),vt(1,i,2), vs(1,i,2), y(j), g22vm(j,2,i))
              call SPEVAL(nym,ym(1,i),vt(1,i,3), vs(1,i,3), y(j), g22vm(j,3,i))
              call SPEVAL(nym,ym(1,i),vt(1,i,4), vs(1,i,4), y(j), g22vm(j,4,i))
              call SPEVAL(nym,ym(1,i),vt(1,i,5), vs(1,i,5), y(j), g22vm(j,5,i))
            else
              g22vm(j,1,i) = vt(nym,i,1)
              g22vm(j,2,i) = vt(nym,i,2)
              g22vm(j,3,i) = vt(nym,i,3)
              g22vm(j,4,i) = vt(nym,i,4)
              g22vm(j,5,i) = vt(nym,i,5)
            end if
          end do
        end do

        deallocate( ym, vt, vs )
        
        return
10      format(8(1pe13.6,1x))

1000    write(*,*) 'Error reading profile...'
        call exit(1)

        end
