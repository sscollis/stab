!-----------------------------------------------------------------------------
        subroutine getmat(t, rmu, rlm, con, &
                          dmu, d2mu, dlm, d2lm, dcon, d2con)
                           
!.... Calculate the Navier-Stokes material properties

        use stuff
        implicit none
        real t(:), rmu(:), rlm(:), con(:), dmu(:), d2mu(:)
        real dlm(:), d2lm(:), dcon(:), d2con(:)
        
        if (mattyp .eq. 0) then         ! constant visc
          rmu = datmat(1)
          dmu = 0.0
          d2mu = 0.0
          dlm = 0.0
          d2lm = 0.0
        else                            ! Sutherland's law
          rmu  = datmat(1) * (t)/datmat(2) * sqrt((t)/datmat(2)) *  &
                 (datmat(2) + datmat(3)) / (t + datmat(3))
          dmu  = (datmat(1) * (3.0*datmat(3)+t)*(datmat(3)+datmat(2)) * &
                 sqrt(t/datmat(2)))/(2.0 * (datmat(3)+t)**2 * datmat(2))
          d2mu = (datmat(1) * (3.0*datmat(3)**2 - 6.0*datmat(3)*t - t**2) * &
                 (datmat(3)+datmat(2)))/(4.0*(datmat(3)+t)**3 * &
                 sqrt(t/datmat(2)) * datmat(2)**2)
        end if

        con = rmu * cp / Pr
        dcon = dmu * cp / Pr
        d2con = d2mu * cp / Pr
        
        rlm = -pt66 * rmu
        dlm = -pt66 * dmu
        d2lm = -pt66 * d2mu
        
        return
        end

!-----------------------------------------------------------------------------
        subroutine sgetmat(t, rmu, rlm, con, &
                          dmu, d2mu, dlm, d2lm)
                           
!.... Calculate the Navier-Stokes material properties

        use stuff
        implicit none
        real t, rmu, rlm, con, dmu, d2mu
        real dlm, d2lm, dcon, d2con
        
        if (mattyp .eq. 0) then         ! constant visc
          rmu = datmat(1)
          dmu = 0.0
          d2mu = 0.0
          dlm = 0.0
          d2lm = 0.0
        else                            ! Sutherland's law
          rmu  = datmat(1) * (t)/datmat(2) * sqrt((t)/datmat(2)) *  &
                 (datmat(2) + datmat(3)) / (t + datmat(3))
          dmu  = (datmat(1) * (3.0*datmat(3)+t)*(datmat(3)+datmat(2)) * &
                 sqrt(t/datmat(2)))/(2.0 * (datmat(3)+t)**2 * datmat(2))
          d2mu = (datmat(1) * (3.0*datmat(3)**2 - 6.0*datmat(3)*t - t**2) * &
                 (datmat(3)+datmat(2)))/(4.0*(datmat(3)+t)**3 * &
                 sqrt(t/datmat(2)) * datmat(2)**2)
        end if

        con = rmu * cp / Pr
        dcon = dmu * cp / Pr
        d2con = d2mu * cp / Pr
        
        rlm = -pt66 * rmu
        dlm = -pt66 * dmu
        d2lm = -pt66 * d2mu
        
        return
        end

