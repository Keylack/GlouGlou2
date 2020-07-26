program main
    use model
    !use results

    implicit none
    integer :: i
    real :: residualU, residualV
    i = 0

    print*,'Entrer dimensions du maillage (longueur et largeur)'
    read*, height, width
    print*, 'Entrer taille du maillage (X et Y)'
    read*,imax,jmax
    print*,'Pas de temps'
    read*,dt
    print*,'Temps final'
    read*,tmax

    allocate(p(1:imax,1:jmax), ux(1:imax,1:jmax), uy(1:imax,1:jmax), ae(1:imax,1:jmax), aw(1:imax,1:jmax), &
&        an(1:imax,1:jmax), as(1:imax,1:jmax), x(1:imax,1:jmax), y(1:imax,1:jmax), sc(1:imax,1:jmax), &
&        a0(1:imax,1:jmax), ap(1:imax,1:jmax), ap_x(1:imax, 1:jmax), ap_y(1:imax,1:jmax), ux_old(1:imax,1:jmax), &
&        uy_old(1:imax,1:jmax), pcorr(1:imax, 1:jmax), ux00(1:imax, 1:jmax), uy00(1:imax, 1:jmax))

call maillage
!call init_var

t = 0
do while (t < tmax)

    print*, 'Temps :', t

    ux_old = ux
    uy_old = uy
    residual = 10
    i = 0
    do while (residual > convergence)
        ux00 = ux
        uy00 = uy
        i = i + 1
        call PISO
        residualU = Residuals(ux, ux00)
        residualV = Residuals(uy, uy00)
        residual = max(residualU, residualV)
    enddo
    print*,'convergence de PISO en ', i, 'iterations'

t = t + dt

enddo

end program main


