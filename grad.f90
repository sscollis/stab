        subroutine g2(v, g2v, ny, dy)
        
!.... updated to fourth order accurate differencing

        implicit none
        integer ny
        
        real v(ny), g2v(ny), dy
        
!.... compute the gradient in y

!       g2v(1)      = (-3.0*v(1) + 4.0*v(2) - v(3))     &
!                     / (2.0 * dy)

        g2v(1)      = (-25.0*v(1) + 48.0*v(2) - 36.0*v(3) &
                       +16.0*v(4) -  3.0*v(5))/(12.0*dy)

!       g2v(2)      = (v(3) - v(1))/(2.0 * dy)

        g2v(2)      = ( -3.0*v(1) - 10.0*v(2) + &
                        18.0*v(3) -  6.0*v(4) + &
                         1.0*v(5) ) / (12.0*dy)

        g2v(3:ny-2) = (-v(5:ny) + 8.0 * v(4:ny-1)        &
                      -8.0 * v(2:ny-3) + v(1:ny-4) ) /   &
                      (12.0 * dy)

!       g2v(ny-1)   = (v(ny) - v(ny-2))/(2.0 * dy)

        g2v(ny-1)   = -(-3.0*v(ny)   - 10.0*v(ny-1) + &
                        18.0*v(ny-2) -  6.0*v(ny-3) + &
                         1.0*v(ny-4) ) / (12.0*dy)

!       g2v(ny)     = (3.0*v(ny) - 4.0*v(ny-1) +         &
!                     v(ny-2)) / (2.0 * dy)

        g2v(ny)      = -(-25.0*v(ny)   + 48.0*v(ny-1) - &
                          36.0*v(ny-2) + 16.0*v(ny-3) - & 
                           3.0*v(ny-4)) /(12.0*dy)
        
        return
        end

!-----------------------------------------------------------------------------
        subroutine grad( ndof, nx, ny, nz, v, g1v, g2v, g3v, ix, dx, dy, dz)

!.... updated to fourth order accurate differencing

        integer nx, ny, nz
        
        real v(ny,ndof,nx), g1v(ny,ndof), g2v(ny,ndof), g3v(ny,ndof)
        
!.... compute the gradient in x

        if ( nx .eq. 1) g1v = 0.0
        
!.... compute the gradient in y

!       g2v(1,:)      = (-3.0*v(1,:,ix) + 4.0*v(2,:,ix) - v(3,:,ix)) &
!                       / (2.0 * dy)

        g2v(1,:)      = (-25.0*v(1,:,ix) + 48.0*v(2,:,ix) - 36.0*v(3,:,ix) &
                         +16.0*v(4,:,ix) -  3.0*v(5,:,ix))/(12.0*dy)

!       g2v(2,:)      = (v(3,:,ix) - v(1,:,ix))/(2.0 * dy)

        g2v(2,:)      = ( -3.0*v(1,:,ix) - 10.0*v(2,:,ix) + &
                          18.0*v(3,:,ix) -  6.0*v(4,:,ix) + &
                           1.0*v(5,:,ix) ) / (12.0*dy)

        g2v(3:ny-2,:) = (-v(5:ny,:,ix) + 8.0 * v(4:ny-1,:,ix)        &
                        -8.0 * v(2:ny-3,:,ix) + v(1:ny-4,:,ix) ) /   &
                        (12.0 * dy)

!       g2v(ny-1,:)   = (v(ny,:,ix) - v(ny-2,:,ix))/(2.0 * dy)

        g2v(ny-1,:)   = -(-3.0*v(ny,:,ix)   - 10.0*v(ny-1,:,ix) + &
                          18.0*v(ny-2,:,ix) -  6.0*v(ny-3,:,ix) + &
                           1.0*v(ny-4,:,ix) ) / (12.0*dy)

!       g2v(ny,:)     = (3.0*v(ny,:,ix) - 4.0*v(ny-1,:,ix) +         &
!                       v(ny-2,:,ix)) / (2.0 * dy)

        g2v(ny,:)      = -(-25.0*v(ny,:,ix)   + 48.0*v(ny-1,:,ix) - &
                            36.0*v(ny-2,:,ix) + 16.0*v(ny-3,:,ix) - & 
                             3.0*v(ny-4,:,ix)) /(12.0*dy)

!.... compute the gradient in z

        if ( nz .eq. 1) g3v = 0.0
        
        return
        end
