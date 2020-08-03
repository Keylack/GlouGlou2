module model
    implicit none

    include 'var'

    contains
    subroutine maillage

        integer :: i,j

        dx = height/imax
        dy = width/jmax 

        do i = 1, imax 
            do j = 1, jmax
                x(i,j) = dx/2 + (i-1)*dx
                y(i,j) = dy/2 + (j-1)*dy
            enddo
        enddo 
    return
    end subroutine

    subroutine XMomentum
        integer :: i, j
        real :: uw, ue, vn, vs

        a0u = 0
        aeu = 0
        awu = 0
        asu = 0
        scu = 0
        apu = 1

        do i = 2, imax 
            do j = 1, jmax 
                a0u(i,j) = dx*dy/dt 
                apu(i,j) = dx*dy/dt 
                if(i /= 1) then
                    uw = (ux(i, j) + ux(i - 1, j))/2
                    awu(i,j) = max(uw, 0.) * dy + nu*dy/dx
                    apu(i,j) = apu(i,j) + max(-uw, 0.) * dy + nu*dy/dx
                    scu(i,j) = (p(i - 1,j) - p(i,j))*dy
                endif

                if(i /= imax) then
                    ue = (ux(i, j) + ux(i + 1, j))/2
                    aeu(i, j) = max(-ue, 0.) * dy + nu *dy/dx 
                    apu(i,j) = apu(i,j) + max(ue, 0.) * dy + nu *dy/dx 
                endif

                if(j /= 1) then
                    vs = (uy(i ,j) + uy(i - 1, j))/2
                    asu(i,j) = max(vs, 0.) * dx + nu*dx/dy
                    apu(i,j) = apu(i,j) + max(-vs, 0.) * dx + nu*dx/dy
                endif

                if(j /= jmax) then
                    vn = (uy(i, j + 1) + uy(i - 1, j + 1)) / 2
                    anu(i,j) = max(-vn, 0.) * dx + nu*dx/dy
                    apu(i,j) = apu(i,j) + max(vn, 0.) * dx + nu*dx/dy
                endif


                scu(i,jmax) = scu(i, jmax) + 2 * uwall * nu * dx/dy
                apu(i, jmax) = apu(i, jmax) + 2*nu*dx/dy
                
                apu(i, 1) = apu(i, 1) + 2*nu*dx/dy
                 
                apu(imax, j) = apu(imax, j) + nu*dy/dx
            enddo
        enddo
        do i = 2, imax 
            do j = 1, jmax
                ap_x(i,j) = apu(i,j) - (awu(i,j) + aeu(i,j) + asu(i,j) + anu(i,j))
            enddo
        enddo

        ux0 = ux

    return
    end subroutine


    subroutine YMomentum
        integer :: i, j
        real :: uw, ue, vn, vs

        a0v = 0
        aev = 0
        awv = 0
        asv = 0
        scv = 0
        apv = 1

        do i = 1, imax 
            do j = 2, jmax 
                a0v(i,j) = dx*dy/dt 
                apv(i,j) = apv(i,j) + a0v(i,j)

                if(i /= 1) then
                    uw = (ux(i, j) + ux(i, j - 1))/2
                    awv(i,j) = max(uw, 0.) * dy + nu*dy/dx
                    apv(i,j) = apv(i,j) + max(-uw, 0.) * dy + nu*dy/dx
                endif

                if(i /= imax) then
                    ue = (ux(i + 1, j) + ux(i + 1, j - 1))/2
                    aev(i, j) = max(-ue, 0.) * dy + nu *dy/dx 
                    apv(i,j) = apv(i,j) + max(ue, 0.) * dy + nu *dy/dx 
                endif

                if(j /= 1) then
                    vs = (uy(i ,j) + uy(i, j - 1))/2
                    asv(i,j) = max(vs, 0.) * dx + nu*dx/dy
                    scv(i,j) = (p(i,j - 1) - p(i,j))*dx
                    apv(i,j) = apv(i,j) + max(-vs, 0.) * dx + nu*dx/dy
                endif

                if(j /= jmax) then
                    vn = (uy(i, j) + uy(i, j + 1)) / 2
                    anv(i,j) = max(-vn, 0.) * dx + nu*dx/dy
                    apv(i,j) = apv(i,j) + max(vn, 0.) * dx + nu*dx/dy
                endif

                apv(i, jmax) = apv(i, jmax) + nu*dx/dy

                apv(1, j) = apv(1, j) + 2*nu*dy/dx
                apv(imax, j) = apv(imax, j) + 2*nu*dy/dx
            enddo
        enddo
        do i = 1, imax 
            do j = 2, jmax
                ap_y(i,j) = apv(i,j) - (awv(i,j) + aev(i,j) + asv(i,j) + anv(i,j))
            enddo
        enddo

        uy0 = uy

    return
    end subroutine

    subroutine Continuity

        integer :: i, j, k 
        real, dimension(imax, jmax) :: dnb_x, dnb_y
        a0p = 0
        aep = 0
        awp = 0
        asp = 0
        scp = 0
        app = 1

        do i = 1, imax 
            do j = 1, jmax 
                dnb_x(i,j) = dy/ap_x(i,j)
                dnb_y(i,j) = dx/ap_y(i,j)
            enddo
        enddo

        do i = 1, imax 
            do j = 1, jmax 

                if(i /= 1) then
                    awp(i,j) = dnb_x(i,j)*dy
                    scp(i,j) = scp(i,j) + ux(i,j) * dy
                endif

                if(i /= imax) then
                    aep(i,j) = dnb_x(i + 1, j)*dy
                    scp(i,j) = scp(i,j) - ux(i + 1,j) * dy
                endif

                if(j /= 1) then
                    asp(i,j) = dnb_y(i,j)*dx
                    scp(i,j) = scp(i,j) + uy(i, j) * dx
                endif

                if(j /= jmax) then
                    anp(i,j) = dnb_y(i,j + 1)*dx
                    scp(i,j) = scp(i,j) - uy(i, j + 1) * dx
                endif
                app(i,j) = aep(i,j) + awp(i,j) + anp(i,j) + asp(i,j)
            enddo
        enddo

        pcorr0 = pcorr


    return
    end subroutine

    real function Residuals(Phi, Phi0)

        integer :: i, j
        real, dimension(imax, jmax) :: Phi, Phi0

        Residuals = 0
        do i = 1, imax 
            do j = 1, jmax
                Residuals = Residuals + (abs(Phi(i,j) - Phi0(i,j)))/(max(abs(Phi0(i,j)), 1e-10))
            enddo
        enddo

    return 
    end function


    subroutine correctionChamps

        integer :: i, j 
        real, dimension(imax, jmax) :: dnb_x, dnb_y
        real :: alpha

        do i = 1, imax 
            do j = 1, jmax 
                dnb_x(i,j) = dy/ap_x(i,j)
                dnb_y(i,j) = dx/ap_y(i,j)
            enddo
        enddo
        alpha = 0.6
        p = p + pcorr

        do i = 2, imax
            do j = 1, jmax
                ux(i,j) = ux(i,j) + alpha * dnb_x(i,j) * (pcorr(i-1,j) - pcorr(i,j))
            enddo
        enddo

        do i = 1, imax
            do j = 2, jmax
                uy(i,j) = uy(i,j) + alpha * dnb_y(i,j) * (pcorr(i,j - 1) - pcorr(i,j))
            enddo
        enddo

    return
    end subroutine

    subroutine verifContinuite
        integer :: i,j
        real, dimension(imax,jmax) :: resCont, sumU
        real :: umax

        rescontinuity = 0.
        sumU = 0
        resCont = 0.
        do i = 1, imax
            do j = 1, jmax
                if(i /= 1) then
                    resCont(i,j) = resCont(i,j) + ux(i,j)
                    sumU(i,j) = sumU(i,j) + abs(ux(i,j))
                endif

                if(i /= imax) then
                    resCont(i,j) = resCont(i,j) - ux(i+1, j)
                    sumU(i,j) = sumU(i,j) + abs(ux(i+1, j))
                endif

                if(j /= 1) then
                    resCont(i,j) = resCont(i,j) + uy(i,j)
                    sumU(i,j) = sumU(i,j) + abs(uy(i,j))
                endif

                if(j /= jmax)then
                    resCont(i,j) = resCont(i,j) - uy(i, j + 1)
                    sumU(i,j) = sumU(i,j) + abs(uy(i, j + 1))
                endif
            enddo
        enddo

        do i = 1, imax 
            do j = 1,jmax 
                rescontinuity = rescontinuity + abs(resCont(i,j))
            enddo
        enddo
    return
    end subroutine

end module model