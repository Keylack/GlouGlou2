module algo
    use model
    implicit none
    contains

        subroutine SIMPLE 
            integer :: i,j,k,l
            i = 0
            j = 0
            k = 0
            l = 0
            residualUx = 10
            residualUy = 10
            residualp = 10
            rescontinuity = 10
        
            do while (residualUx > convergence)
                i = i+1
                call XMomentum
                call TDMA(2, imax, 1, jmax, ux, aeu, awu, asu, anu, scu, a0u, apu, ux_old)
                residualUx = Residuals(ux0, ux)
            enddo
            do while (residualUy > convergence)
                j = j + 1
                call YMomentum
                call TDMA(1, imax, 2, jmax, uy, aev, awv, asv, anv, scv, a0v, apv, uy_old)
                residualUy = Residuals(uy0, uy)
            enddo
            pcorr = 0
            pcorr0 = 0
            do while (residualp > convergence)
                k = k + 1
                call Continuity
                call GaussSeidel(1, imax, 1, jmax, pcorr, aep, awp, asp, anp, scp, a0p, app, zero)
                residualp = Residuals(pcorr, pcorr0)
                if( mod(k, 5000) == 0) then
                    !print*,'Iteration pression : ', k, '    Résidu :', residualp
                endif
            enddo
            print*, 'Equation P converge en ', k, 'itérations'
            !print*,'app',app
            !print*,'aep',aep
            !print*,'awp',awp
            !print*,'asp',asp
            !print*,'anp',anp
            !print*,'scp',scp
            !print*,'ux',ux 
            !print*,'uy',uy
            !print*,pcorr
            !print*,'ux avant correction', ux
            !print*,'uy avant correction', uy
            call correctionChamps
            call verifContinuite
            !print*,rescontinuity
            !print*,'ux après correction', ux
            !print*,'uy après correction', uy

        return
        end subroutine

        subroutine GaussSeidel(istart, iend, jstart, jend, Phi, ae, aw, as, an, sc, a0, ap, Phi0)
            integer :: i,j, istart, iend, jstart, jend
            real, dimension(imax,jmax) :: Phi, ae, aw, an, as, sc, a0, ap, Phi0
        
            do i = istart, iend 
                
                do j = jstart, jend 
                    Phi(i,j) = a0(i,j)*Phi0(i,j) + sc(i,j)
                    if(i /= 1) then
                        Phi(i,j) = Phi(i,j) + aw(i,j)*Phi(i - 1, j)
                    endif
        
                    if(i /= imax) then
                        Phi(i,j) = Phi(i,j) + ae(i,j)*Phi(i + 1,j)
                    endif
        
                    if(j /= 1) then
                        Phi(i,j) = Phi(i,j) + as(i,j)*Phi(i, j-1)
                    endif
        
                    if(j /= jmax) then
                        Phi(i,j) = Phi(i,j) + an(i,j)*Phi(i, j+1)
                    endif
        
                    Phi(i,j) = Phi(i,j) / ap(i,j)
                enddo
            enddo
        return
        end subroutine

        subroutine TDMA(istart, iend, jstart, jend, Phi, ae, aw, as, an, sc, a0, ap, Phi0)
        
            integer :: i,j, istart, iend, jstart, jend
            real, dimension(imax,jmax) :: Phi, ae, aw, an, as, sc, a0, ap, Phi0, Phidt_2
            real,dimension(imax,jmax) :: a, P, Q, b, c, d

            real :: dt_2

            dt_2 = dt/2


            ! Résolution ligne par ligne alternée, première itération

            do i =istart, iend
                do j = jstart, jend
                    c(istart,j) = 0
                    b(iend,j) = 0
                    a(i,j) = ap(i,j) * dt / dt_2
                    if(i /= iend) b(i,j) = ae(i,j)
                    if(i /= istart) c(i,j) = aw(i,j)
                    d(i,j) = a0(i,j) * dt / dt_2 * Phi0(i,j) + sc(i,j)
                    if(j /= jstart)d(i,j) = d(i,j) + as(i,j)*Phi0(i,j - 1)
                    if(j /= jend)d(i,j) = d(i,j) + an(i,j)*Phi0(i,j + 1)
                enddo
            enddo

            !1re étape init descente
            do j = jstart, jend
                P(istart,j) = b(istart,j)/a(istart,j)
                Q(istart,j) = d(istart,j)/a(istart,j)
            enddo
            !2e étape : descente
            do i = istart + 1, iend 
                do j = jstart, jend
                    P(i,j) = b(i,j) / (a(i,j) - c(i,j) * P(i - 1,j))
                    Q(i,j) = (d(i,j) + c(i,j)*Q(i - 1,j))/(a(i,j) - c(i,j)*P(i - 1,j))
                enddo
            enddo
            
            !3ème étape : init remontée
            do j = jstart, jend
                Phidt_2(iend, j) = Q(iend, j)
            enddo

            !Remontée
            do i = iend - 1, istart, - 1
                do j = jstart, jend 
                    Phidt_2(i,j) = P(i,j) * Phidt_2(i + 1, j) + Q(i,j)
                enddo
            enddo

            !Deuxième itération

            do i = istart, iend
                do j = jstart, jend
                    c(i,jstart) = 0
                    b(i,jend) = 0
                    if(i /= iend) b(i,j) = an(i,j)
                    if(i /= istart) c(i,j) = as(i,j)
                    d(i,j) = a0(i,j) * dt / dt_2 * Phidt_2(i,j) + sc(i,j)
                    if(i /= istart)d(i,j) = d(i,j) + aw(i,j)*Phidt_2(i - 1,j)
                    if(i /= iend)d(i,j) = d(i,j) + ae(i,j)*Phidt_2(i + 1,j)
                enddo
            enddo

            !init descente
            do i = istart, iend
                P(i,jstart) = b(i,jstart)/a(i,jstart)
                Q(i,jstart) = d(i,jstart)/a(i,jstart)
            enddo

            !2e étape : descente
            do i = istart, iend 
                do j = jstart + 1, jend
                    P(i,j) = b(i,j) / (a(i,j) - c(i,j) * P(i,j - 1))
                    Q(i,j) = (d(i,j) + c(i,j)*Q(i,j - 1))/(a(i,j) - c(i,j)*P(i,j - 1))
                enddo
            enddo

            !3ème étape : init remontée
            do i = istart, iend
                Phi(i, jend) = Q(i, jend)
            enddo

            !Remontée
            do i = istart, iend
                do j = jend - 1, jstart, -1
                    Phi(i,j) = P(i,j) * Phi(i, j + 1) + Q(i,j)
                enddo
            enddo
        return
        end subroutine

        subroutine GradConj(istart, iend, jstart, jend, Phi, ae, aw, as, an, sc, a0, ap, Phi0)
            integer :: istart, iend, jstart, jend
            real, dimension(imax,jmax) :: Phi, ae, aw, an, as, sc, a0, ap, Phi0
            integer :: i,j,k 
            real, dimension(imax, jmax) :: rk, rk2, pk, ApK
            real :: rksum, alphak, betak, pkApk
            k = 0
            rksum = 0
            alphak = 0
            betak = 0
            pkApk = 0
            do i = istart,iend 
                do j = jstart, jend 
                    rk(i,j) = sc(i,j) + Phi0(i,j)*a0(i,j) - ap(i,j)*Phi(i,j)
                    if(i /= istart) rk(i,j) = rk(i,j) + aw(i,j)*Phi(i - 1,j)
                    if(i /= iend) rk(i,j) = rk(i,j) + ae(i,j)*Phi(i + 1,j)
                    if(j /= jstart) rk(i,j) = rk(i,j) + as(i,j) * Phi(i, j - 1)
                    if(j /= jend) rk(i,j) = rk(i,j) + an(i,j) * Phi(i, j + 1)
                    rksum = rksum + abs(rk(i,j))
                enddo
            enddo
            pk = rk
            do while (rksum > 1e-10)

                do i = istart, iend 
                    do j = jstart, jend 
                        alphak = alphak + rk(i,j) * rk(i,j)
                        Apk(i,j) = ap(i,j) * pk(i,j)
                        if(i /= istart) Apk(i,j) = Apk(i,j) - aw(i,j)*pk(i - 1,j)
                        if(i /= iend) Apk(i,j) = Apk(i,j) - ae(i,j) * pk(i + 1,j)
                        if(j /= jstart) Apk(i,j) = Apk(i,j) - as(i,j) * pk(i, j - 1)
                        if(j /= jend) Apk(i,j) = Apk(i,j) - an(i,j) * pk(i, j + 1)
                        pkApk = pkApk + pk(i,j) * Apk(i,j)
                    enddo 
                enddo 

                alphak = alphak / pkApk
                k = k + 1
                do i = istart, iend 
                    do j = jstart, jend 
                        Phi(i,j) = Phi(i,j) + alphak * pk(i,j)
                        rk2(i,j) = rk(i,j )- alphak * Apk(i,j)
                    enddo 
                enddo
                rksum = 0
                betak = 0
                do i = istart, iend 
                    do j = jstart, jend 
                        rksum = rksum + abs(rk2(i,j))
                    enddo
                enddo

                do i = istart, iend 
                    do j = jstart, jend 
                        if(rk(i,j) /= 0) then
                            betak = betak + rk2(i,j) * rk2(i,j) / (rk(i,j) * rk(i,j))
                        endif
                    enddo 
                enddo 
                print*, 'rk2', rk2
                print*, 'rk', rk
                print*,'betak', betak
                pk = rk2 + betak * pk
                rk = rk2
                rk2 = 0
            enddo
        
        end subroutine
end module algo