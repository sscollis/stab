!=============================================================================!
        subroutine getmean( vm, y, eta, lny, iver)
!
!  Input the mean flow profile assuming parallel flow 
!  (i.e. the v-velocity is forced to zero)
!
!=============================================================================!
        use stuff
        implicit none

        integer :: lny
        real :: vm(lny,ndof,nx), eta(lny), y(lny), ymaxm
        real, allocatable :: ym(:,:), vt(:,:,:), vs(:,:,:)
        
        integer :: i, j, k, nxm, nym, ndofm, ier, iver
#if 0
        real :: tmp
#else
        character(256) :: tmp
#endif
        
        character(80) :: base, fname
!=============================================================================!

!.... read the mean field and spline to the new grid
        
        base = 'profile'
        call makename(base,iver,fname)
        if (verbose) then
          write(*,"(/,'Reading mean flow from:  ',a)") fname
        endif
        open (unit=10, file=fname, form='formatted', status='old',err=1000)
        
!       read (10,*) nxm, nym, ndofm, ymaxm
!       if (ndofm .ne. ndof) then
!         write (*,*) 'NDOF in mean field is incorrect ',ndofm
!         call exit(1)
!       end if
        
        nxm = 1
        nym = 0
 20     continue
#if 0
          read(10,*,end=30) tmp
          nym = nym + 1
#else
          read(10,'(a)',end=30) tmp
          if (tmp(1:1).ne.'#') then
            nym = nym + 1
          endif
#endif
        goto 20
 30     continue
        rewind(10)

        allocate( ym(nym,nxm), vt(nym,nxm,ndof), vs(nym,nxm,ndof), &
                  STAT=ier)
        if (ier .ne. 0) then
          write(*,*) 'Error allocating mean field'
          stop
        end if
        
        do i = 1, nxm
#if 0
          do j = 1, nym
            read (10,*) ym(j,i), (vt(j,i,k),k=1,ndof)
            vt(j,i,3) = zero ! parallel flow assumption
          end do
#else
          j = 1
 40       continue
            read (10,'(a)',end=50) tmp
            if (tmp(1:1).ne.'#') then
              read (tmp,*) ym(j,i), (vt(j,i,k),k=1,ndof)
              vt(j,i,3) = zero ! parallel flow assumption
              j = j + 1
            endif
            goto 40
 50       continue
#endif
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
        
        deallocate( ym, vt, vs )
        
        return
10      format(8(1pe13.6,1x))

1000    write(*,*) 'Error reading profile...'
        stop

        end subroutine getmean
