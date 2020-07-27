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

subroutine init_var
    integer :: i, j

    do i = 1, imax 
        do j = 1, jmax 
            p(i,j) = 0.
            ux(i,j) = 0.
            uy(i,j) = 0.
            sc(i,j) = 0.
            ae(i,j) = 0.
            aw(i,j) = 0.
            as(i,j) = 0.
            an(i,j) = 0.
            pcorr(i,j) = 0.
        enddo
    enddo
return
end subroutine

subroutine XMomentum
    integer :: i, j
    real :: uw, ue, vn, vs
    real, dimension(imax, jmax) :: ux0

    a0 = 0
    ae = 0
    aw = 0
    as = 0
    sc = 0
    ap = 1

    do i = 1, imax 
        do j = 1, jmax 
            a0(i,j) = dx*dy/dt 

            if(i /= 1) then
                uw = (ux(i, j) + ux(i - 1, j))/2
                aw(i,j) = max(uw, 0.) * dy + nu*dy/dx
                sc(i,j) = (p(i - 1,j) - p(i,j))*dy
            endif

            if(i /= imax) then
                ue = (ux(i, j) + ux(i + 1, j))/2
                ae(i, j) = -max(-ue, 0.) * dy + nu *dy/dx 
            endif

            if(j /= 1) then
                vs = (uy(i ,j) + uy(i - 1, j))/2
                as(i,j) = max(vs, 0.) * dx + nu*dx/dy
            endif

            if(j /= jmax) then
                vn = (uy(i, j + 1) + uy(i - 1, j + 1)) / 2
                an(i,j) = -max(-vn, 0.) * dx + nu*dx/dy
            endif

            ap(i,j) = a0(i,j) + aw(i,j) + ae(i,j) + an(i,j) + as(i,j)

            if(j == jmax) then
                sc(i,jmax) = sc(i, jmax) + 2 * uwall * nu * dx/dy
                ap(i, jmax) = ap(i, jmax) + 2*nu*dx/dy
            endif

            if(j == 1) then
                ap(i, 1) = ap(i, jmax) + 2*nu*dx/dy
            endif
            if(i == imax) then
                ap(i, 1) = ap(i, jmax) + 2*nu*dy/dx
            endif
        enddo
    enddo
    do i = 1, imax 
        do j = 1, jmax
            ap_x(i,j) = ap(i,j)
        enddo
    enddo
    ux0 = ux
    call GaussSeidel(2, imax, 1, jmax, ux, ae, aw, as, an, sc, a0, ap, ux_old)
    residualUx = Residuals(ux0, ux)
!    print*,'ae', ae
!    print*,'aw',aw
!    print*,'ap',ap
!    print*,'as',as
!    print*,'an',an
!    print*,'sc',sc
return
end subroutine

subroutine YMomentum
    integer :: i, j
    real :: uw, ue, vn, vs
    real, dimension(imax, jmax) :: uy0

    a0 = 0
    ae = 0
    aw = 0
    as = 0
    sc = 0
    ap = 1

    do i = 1, imax 
        do j = 1, jmax 
            a0(i,j) = dx*dy/dt 

            if(i /= 1) then
                uw = (ux(i, j) + ux(i, j - 1))/2
                aw(i,j) = max(uw, 0.) * dy + nu*dy/dx
            endif

            if(i /= imax) then
                ue = (ux(i + 1, j) + ux(i + 1, j - 1))/2
                ae(i, j) = -max(-ue, 0.) * dy + nu *dy/dx 
            endif

            if(j /= 1) then
                vs = (uy(i ,j) + uy(i, j - 1))/2
                as(i,j) = max(vs, 0.) * dx + nu*dx/dy
                sc(i,j) = (p(i,j - 1) - p(i,j))*dy
            endif

            if(j /= jmax) then
                vn = (uy(i, j) + uy(i, j + 1)) / 2
                an(i,j) = -max(-vn, 0.) * dx + nu*dx/dy
            endif

            ap(i,j) = a0(i,j) + aw(i,j) + ae(i,j) + an(i,j) + as(i,j)

            if(j == jmax) then
                ap(i, jmax) = ap(i, jmax) + 2*nu*dx/dy
            endif

            if(i == 1) then
                ap(1, j) = ap(1, j) + 2*nu*dy/dx
            endif
            if(i == imax) then
                ap(1, j) = ap(1, j) + 2*nu*dy/dx
            endif
        enddo
    enddo
    do i = 1, imax 
        do j = 1, jmax
            ap_y(i,j) = ap(i,j)
        enddo
    enddo
   ! print*,'ap',ap
   ! print*,'ae', ae
   ! print*,'aw',aw
   ! print*,'as',as
   ! print*,'an',an
   ! print*,'sc',sc        
   ! print*,'a0',a0    
   ! print*,ux
    uy0 = uy
    call GaussSeidel(1, imax, 2, jmax, uy, ae, aw, as, an, sc, a0, ap, uy_old)
    residualUy = Residuals(uy0, uy)

return
end subroutine

subroutine Continuity

    integer :: i, j, k 
    real, dimension(imax, jmax) :: dnb_x, dnb_y
    a0 = 0
    ae = 0
    aw = 0
    as = 0
    sc = 0
    ap = 1
    

    do i = 1, imax 
        do j = 1, jmax 
            dnb_x(i,j) = dy/ap_x(i,j)
            dnb_y(i,j) = dx/ap_y(i,j)
        enddo
    enddo

    do i = 1, imax 
        do j = 1, jmax 

            if(i /= 1) then
                aw(i,j) = dnb_x(i,j)*dy
                sc(i,j) = sc(i,j) + ux(i,j) * dy
            endif

            if(i /= imax) then
                ae(i,j) = dnb_x(i + 1, j)*dy
                sc(i,j) = sc(i,j) - ux(i + 1,j) * dy
            endif

            if(j /= 1) then
                as(i,j) = dnb_y(i,j)*dx
                sc(i,j) = sc(i,j) + uy(i, j) * dx
            endif

            if(j /= jmax) then
                an(i,j) = dnb_y(i,j + 1)*dx
                sc(i,j) = sc(i,j) - uy(i + 1, j) * dx
            endif
            ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j)
        enddo
    enddo
!    print*,'ap',ap
!    print*,'ae', ae
!    print*,'aw',aw
!    print*,'as',as
!    print*,'an',an
!    print*,'sc',sc        
!    print*,'a0',a0    
!    print*,uy
    pcorr0 = pcorr
    call GaussSeidel(1, imax, 1, jmax, pcorr, ae, aw, as, an, sc, a0, ap, pcorr0)
    residualp = Residuals(pcorr, pcorr0)

return
end subroutine

real function Residuals(Phi, Phi0)

    integer :: i, j
    real, dimension(imax, jmax) :: Phi, Phi0

    Residuals = 0
    do i = 1, imax 
        do j = 1, jmax
            Residuals = Residuals + (abs(Phi(i,j) - Phi0(i,j)))/(max(abs(Phi(i,j)), 1e-10))
        enddo
    enddo

return 
end function


subroutine GaussSeidel(istart, iend, jstart, jend, Phi, ae, aw, as, an, sc, a0, ap, Phi0)
    integer :: i,j, istart, iend, jstart, jend
    real, dimension(imax,jmax) :: Phi, ae, aw, an, as, sc, a0, ap, Phi0

    do i = istart, iend lboundl
        
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
    p = p + alpha * pcorr

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
    real, dimension(imax,jmax) :: resCont

    rescontinuity = 0.
    resCont = 0.
    do i = 1, imax
        do j = 1, jmax
            if(i /= 1) then
                resCont(i,j) = resCont(i,j) + ux(i,j)
            endif

            if(i /= imax) then
                resCont(i,j) = resCont(i,j) - ux(i+1, j)
            endif

            if(j /= 1) then
                resCont(i,j) = resCont(i,j) + uy(i,j)
            endif

            if(j /= jmax)then
                resCont(i,j) = resCont(i,j) - uy(i, j + 1)
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