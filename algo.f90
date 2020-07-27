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
            enddo
        
        !    print*,ux
        
            !print*,'Equation Ux converge en : ', i, 'iterations'
        
            do while (residualUy > convergence)
                j = j + 1
                call YMomentum
            enddo
        
            !print*,'Equation Uy converge en : ', j, 'iterations'

            pcorr = 0
            pcorr0 = 0
            do while (residualp > convergence)
                k = k + 1
                call Continuity
                if( mod(k, 5000) == 0) then
                    !print*,'Iteration pression : ', k, '    Résidu :', residualp
                endif
            enddo
            !print*, 'Equation P converge en ', k, 'itérations'
            call correctionChamps
            !print*,ux
            
        
        return
        end subroutine


end module algo