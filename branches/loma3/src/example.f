        program main
        implicit none

        x(1)=0.0  
        Dx=0.01
        N=1000;
        do j=2,N
            x(j) = x(j) + Dx
