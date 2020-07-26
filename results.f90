module results
    use model
    implicit none

    
    contains

    subroutine resultat
        character*3 :: cartime
        character*1 :: ZERO
        integer :: tRes, k, i, j
        print*,'Writing results t = :', nint(t)
10      format('i = ', I3, '  j = ', I3)
20      format(5(1X, E12.5))
        ZERO = '0'
        tRes = aint(t)
        write(carTime, '(I3)')tRes 

        k = 1
        do while(carTime(k:k) == ' ' .and. k < 3)
            k = k+1
        enddo
        carTime = carTime(k:3)
        if(k == 2) carTime = ZERO//carTime
        if (k == 3) carTime = ZERO//ZERO//carTime

        open(unit = 1, file='res.'//cartime)
        write(1, '(A)')'Variables= "x(m)" "y(m)" "Ux(m/s)" "Uy(m/s)" "p(m2/s2)"'
        write(1, 10)imax,jmax

        do i = 1, imax
            do j = 1, jmax
                write(1, 20)x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j)
            enddo
        enddo

        close(unit=1)

    end subroutine
end module results