!----------------------------------------------------------------------!
!      Giving a Runge-Kutta timestep                                   !
!      -----------------------------                                   !
!      Solving non linear system for Low Mach number aproximation      !
!      using Helmholtz decomposition.                                  !
!                                                                      !
!      Explanation and description included on ../doc folder.          !
!                                                                      !
!   variables:                                                         !
!     vor: Omega_y (vertical component of momentum vector curl)        !
!     phi: phi,vertical component of lapl. of div.free momentum comp.  !
!     nxwk,nywk,nzwk,mxwk,mywk,mzwk: working buffers                   !
!     chwk,rhst,drho: more working buffers                             !
!     u00,w00:  00 modes                                               !
!     rf0u,rf0w,u00wk,w00wk: wk buffers for 00modes                    !
!                                                                      !
!    new cross.f for LoMa --AAF (26/02/2014)                           !
!                                                                      !
!   Warning : hvect subroutine in real*8                               !
!             00 modes in real*8                                       !
!                                                                      !
!**********************************************************************!

      subroutine interpol(vor,
     .                  phi,
     .                  psi,
     .                  scal,
     .                  u00,w00,v00,
     .                  u00wk,w00wk,v00wk,
     .                  vor2,
     .                  phi2,
     .                  psi2,
     .                  scal2,
     .                  myid)
      use tem
      use fis
      use wave
      use point
      use ficheros
      use timers
      use MPI_GROUPS
      use input
     

      implicit none

      include "mpif.h"
      include "ctes3D"

      !------------- Parameters for SMR R-K ------------------------!
!!!RK Spalart      
!      real(8)  :: gama(3), alpha(4), beta(3), ibeta(3), xi(3)
!      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
!      parameter (alpha=(/ 29d0/96d0, -3d0/40d0, 1d0/6d0, 29d0/96d0/))
!      parameter (beta =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /))
!      parameter (ibeta=(/ 160d0/37d0, 24d0/5d0,   6d0     /))
!      parameter (xi   =(/        0d0,-17d0/60d0, -5d0/12d0 /)) !aaf 

!RK3LS A. Almagro
c ------------------- Variables ----------------------------------------

      complex(realprec) :: phi(0:my-1,0:mz1,pb:pe),
     .          vor(0:my-1,0:mz1,pb:pe),
     .          psi(0:my-1,0:mz1,pb:pe),scal(0:my-1,0:mz1,pb:pe)
      complex(realprec) :: phi2(0:mye-1,0:mz1,pb:pe),
     .          vor2(0:mye-1,0:mz1,pb:pe),
     .          psi2(0:mye-1,0:mz1,pb:pe),scal2(0:mye-1,0:mz1,pb:pe)

      real(realprec) :: u00(0:my-1),w00(0:my-1),
     .       u00wk(0:mye-1),w00wk(0:mye-1),
     .       v00(0:my-1),v00wk(0:mye-1)
     
      real(realprec)  :: trapz 
      
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,kk,
     .        rksteptest
      
      
      real(realprec) :: iter_time,write_time
      complex(realprec) :: temp1,temp2
      real(4) ::  temporal(3),comun(3)
      
      real(4)  enerdis,H00u,H00w,dumu,dumw,cteu,ctew
      real(4)  massu,massw,massu1,massu2,massw1,massw2,dum
      real(4)  uner(9)
      real(4)  hbot,htop,zbot,ztop,rhobot,rhotop
      integer :: j1,j2
          

      character(100) fname,fdebug
      complex(realprec),allocatable :: tempa(:,:),tempb(:,:)

!debugging variables
      character(3) ext1,ext2 
      integer :: ideb
      
c    ---------------------- Program ------------------------------------

      allocate(tempa(0:mz1,pb:pe),tempb(0:mz1,pb:pe))
      !DEBUG
      !rksteptest=1

c    Defining combinations of wavenumbers useful for calculations.      
      if (myid.eq.0) then

         do k=1,mz1
            tempa(k,1) = 1d0/(bet2(k))
         enddo
         do i=pb+1,pe
            do k=0,mz1
               tempa(k,i) = 1d0/(alp2(i-1) + bet2(k))
            enddo
         enddo
         tempa(0,1)=0d0 !00 mode
      else
         do i=pb,pe
            do k=0,mz1
               tempa(k,i) = 1d0/(alp2(i-1) + bet2(k))
            enddo
         enddo
      endif
      do i=pb,pe
         do k=0,mz1
            tempb(k,i) = tempa(k,i)*xbet(k)
            tempa(k,i) = tempa(k,i)*xalp(i-1)
         enddo
      enddo


c ========================================================================
c                    THIS IS THE PROGRAM
c ========================================================================
      !call interp_linear(1,my,y,u00, mye,ye,u00wk)
      !call interp_linear(1,my,y,v00, mye,ye,v00wk)
      !call interp_linear(1,my,y,w00, mye,ye,w00wk)

      !call interp_nearest(1,my,y,u00, mye,ye,u00wk)
      !call interp_nearest(1,my,y,v00, mye,ye,v00wk)
      !call interp_nearest(1,my,y,w00, mye,ye,w00wk)

      if (myid.lt.numerop) then     ! Comp. Proc.
        call interp_linear(1,my,y,u00, mye,ye,u00wk)
        call interp_linear(1,my,y,v00, mye,ye,v00wk)
        call interp_linear(1,my,y,w00, mye,ye,w00wk)
      write(*,*) "comp proc making interpolation..."
      !-------------------------------------------
      !Before interpolating clean boundaries
      !Mode 00 of scal is special
!!!      !Scalar doesnt need cleaning
!      j1 = 0
!      !j1 = 30
!      !Clean bottom
!      do i=pb,pe
!         do k=0,mz1
!            do j=0,j1
!               vor(j,k,i)= 0.0
!               phi(j,k,i)= 0.0
!               psi(j,k,i)= 0.0
!             enddo
!         enddo
!      enddo
!      !Clean TOP   
!      do i=pb,pe
!         do k=0,mz1
!            do j=my1-j1,my1
!               vor(j,k,i)= 0.0
!               phi(j,k,i)= 0.0
!               psi(j,k,i)= 0.0
!             enddo
!         enddo
!      enddo
!      !-------------------------------------------


      do i=pb,pe
        do k=0,mz1
          call interp_complex(vor(0,k,i),vor2(0,k,i))
        enddo
      enddo
      do i=pb,pe
        do k=0,mz1
          call interp_complex(phi(0,k,i),phi2(0,k,i))
        enddo
      enddo
      do i=pb,pe
        do k=0,mz1
          call interp_complex(psi(0,k,i),psi2(0,k,i))
        enddo
      enddo
      do i=pb,pe
        do k=0,mz1
          call interp_complex(scal(0,k,i),scal2(0,k,i))
        enddo
      enddo
      endif

      !
      call escru(vor2,phi2,psi2,scal2,u00wk,w00wk,v00wk,myid)
      !call escru(vor,phi,psi,scal,mfz,u00,w00,v00,myid)
! ------------------- finished writing image --------------------------!
      !call MPI_BARRIER(MPI_COMM_CALC,ierr)

      end


      subroutine interp_complex(vec_inp,vec_out)
      use input
      use fis
      use MPI_GROUPS
      implicit none
      include "mpif.h"
      include "ctes3D"
     
      integer i,j,k 
      real(realprec),dimension(2,my):: vec_inp
      real(realprec),dimension(2,mye):: vec_out
        
      !call interp_nearest(2,my,y,vec_inp, mye,ye, vec_out)
      call interp_linear(2,my,y,vec_inp, mye,ye, vec_out)
      !call interp_nearest(1,my,y,vec_inp(2,1), mye,ye, vec_out(2,1))
      !call interp_linear(1,my,y,vec_inp(1,1), mye,ye, vec_out(1,1))
      !call interp_linear(1,my,y,vec_inp(2,1), mye,ye, vec_out(2,1))
      

      end subroutine interp_complex
       






