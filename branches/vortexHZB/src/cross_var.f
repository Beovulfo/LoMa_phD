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

      subroutine cross1(vor,
     .                  phi,
     .                  psi,
     .                  scal,
     .                  u00,w00,v00,
     .                  rf0u,rf0w,rf0v,
     .                  u00wk,w00wk,v00wk,
     .                  mxwk,
     .                  mywk,
     .                  mzwk,
     .                  nxwk,
     .                  nywk,
     .                  nzwk,
     .                  ten12,
     .                  ten13,
     .                  ten23,
     .                  rhst,
     .                  drho,
     .                  chwk2,
     .                  chwk,
     .                  myid)
      use diag
      use tem
      use timacc
      use fis
      use point
      use wave
      use ficheros
      use timers
      use cnan
      use matrices
      use MPI_GROUPS
     

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

!!RK Spalart-Lucky Toni (alpha=xi, beta= gamma)
      real(8)  :: gama(3), alpha(3), beta(3), ibeta(3), xi(3)
      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
      parameter (alpha=(/      0d0, -17d0/60d0, -5d0/12d0 /))!aaf
      parameter (beta =(/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))!aaf
      parameter (ibeta=(/ 15d0/8d0,   12d0/5d0,   4d0/3d0 /))!aaf
      parameter (xi   =(/      0d0, -17d0/60d0, -5d0/12d0 /)) !aaf 



!RK appropiate to use with PSI updating
!      real*8     gama(3), alpha(3), beta(3), ibeta(3), xi(3)
!      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
!      parameter (xi   =(/        0d0,-17d0/60d0, -5d0/12d0 /))  
!      parameter (alpha=(/ 4d0/15d0, 1d0/15d0, 1d0/6d0/))
!      parameter (beta =(/ 4d0/15d0, 1d0/15d0,  1d0/6d0 /))
!     parameter (ibeta=(/ 15d0/4d0, 15d0,   6d0     /))
!this RKstep works worse than spalart!
   
c ------------------- Variables ----------------------------------------

      complex*8 :: phi(0:my-1,0:mz1,pb:pe),vor(0:my-1,0:mz1,pb:pe),
     .          psi(0:my-1,0:mz1,pb:pe),scal(0:my-1,0:mz1,pb:pe)

      complex*8 ::  nxwk(0:my-1,0:mgalz-1,pb:pe),
     .          nzwk(0:my-1,0:mgalz-1,pb:pe),
     .          nywk(0:my-1,0:mgalz-1,pb:pe),
     .          mxwk(0:my-1,0:mgalz-1,pb:pe),
     .          mywk(0:my-1,0:mgalz-1,pb:pe),
     .          mzwk(0:my-1,0:mgalz-1,pb:pe),
     .          rhst(0:my-1,0:mgalz-1,pb:pe),
     .          drho(0:my-1,0:mgalz-1,pb:pe),
     .          chwk(0:my-1,0:mgalz-1,pb:pe),
     .          chwk2(0:my-1,0:mgalz-1,pb:pe),
     .          ten12(0:my-1,0:mgalz-1,pb:pe),
     .          ten13(0:my-1,0:mgalz-1,pb:pe),
     .          ten23(0:my-1,0:mgalz-1,pb:pe)

      real(8) :: u00(0:my-1),w00(0:my-1),rf0u(0:my-1),
     .       rf0w(0:my-1),u00wk(0:my-1),w00wk(0:my-1),
     .       v00(0:my-1),v00wk(0:my-1),rf0v(0:my-1)
     
      real(8)  :: trapz 
      real(8)  :: drho00(0:my1) 
      
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,kk
      
      real(8) :: rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire,
     .        dtbeta
      
      real(8)  iter_time,write_time
      complex*8 temp1,temp2
      real(4) ::  temporal(3),comun(3)
      
      real(4)  reynota,H00u,H00w,dumu,dumw,cteu,ctew
      real(4)  massu,massw,massu1,massu2,massw1,massw2,dum
      real(4)  bcpsi00

      character(100) fname
      complex*8,allocatable :: tempa(:,:),tempb(:,:)

!debugging variables
      character(3) ext1,ext2 
      integer :: ideb
      
c    ---------------------- Program ------------------------------------
c debug
      !ideb = 100     
      ideb = 99 

      allocate(tempa(0:mz1,pb:pe),tempb(0:mz1,pb:pe))

      
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
      
      commtimer =0.0D0
      transtimer=0.0D0
      totaltimer=0.0D0

      ire      = 1./re
      ihist    = 0
      irun     = 0   ! first time step is special in tim3rkp
      icfl     = 1   ! first time step always needs a step size
      nanerror = 0   ! Check NAN in the 50 firsts steps
      nanproblem = 0 ! initializes NAN controller

      if(myid.eq.0) then
         fname=filstt(1:index(filstt,' ')-1)//'.cf'
         write (*,*) fname
         open(29,file=fname,status='unknown')
      endif


!Initialization: 
      if (myid.eq.0) then
          do j=0,my1
             rf0u(j)=0d0
             rf0w(j)=0d0
             rf0v(j)=0d0
          enddo
      endif
c ========================================================================
c                 THIS IS THE TIME LOOP
c ========================================================================

      do 30 istep=1,nstep
      !save istep into ext1
       write(ext1,'(i3.3)') istep

         if (myid.eq.0) then
           totaltimer = totaltimer-MPI_WTIME()
           iter_time=-MPI_WTIME()
         endif

         if (mod(istep-1,nhist).eq.0) ihist=1
         if (mod(istep-1,ncfl).eq.0) icfl = 1
         
! ------------------- write image to a file ---------------------------!

         IF (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then

            if (myid.eq.0) then
               write_time = -MPI_WTIME()
            endif
            do j=0,my1
               u00(j) = u00wk(j) 
               w00(j) = w00wk(j)
               v00(j) = v00wk(j)
            enddo
            !call escru(vor,phi,u00,w00,sp,spwk,myid)
            call escru(vor,phi,u00,w00,v00,myid)
!write scalar field
            call escruscal(scal,1,myid)
!aaf Here we need to call a escru for scalar field

            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
            endif
         ENDIF

! ------------------- finished writing image --------------------------!

      if (myid.gt.numerop-1) goto 30    ! only for save procs

      if (irun.eq.0) then ! this is done only for the first step
           irun = 1
         if (myid.eq.0) then
            do j=0,my1
              u00wk(j)=u00(j)
              w00wk(j)=w00(j)
              v00wk(j)=v00(j)
            enddo
!--------------------------------------------------------------!
!!debug print to file
            open(unit=40,file="u00init.txt",status="unknown")
            open(unit=38,file="T00init.txt",status="unknown")
                  do j=0,my-1
                      write(40,*) y(j+1), u00wk(j)
                      write(38,*) y(j+1), scal(j,0,pb)
                  enddo
             close(40)
             close(38)
!--------------------------------------------------------------!
         endif
      endif ! end special first step

               !  Runge-Kutta third order  !
!=======================================================================!                
!----------- RUNGE-KUTTA SUBSTEPS START---------------------------------!
!=======================================================================!                

      do 10 rkstep=1,3
      !save rkstep to char ext2
          write(ext2,'(i3.3)') rkstep
!at rkstep=1 dtxi=0, for rkstep>1 we have the proper DT for the istep
!and it shouldn't change through the istep.
         dtxi=Deltat*xi(rkstep)

c-------------save vor and psi in nxwk/nywk to work---------------
         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  nxwk(j,k,i)=vor(j,k,i)   !Omega_y
                  nywk(j,k,i)=psi(j,k,i)   !PSI
               enddo
            enddo
         enddo

!-----------------(1) Laplacian phi solver--------------------------!
!  Solve laplacian to find my and dmy/dy
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               !call Lapvdv(phi(0,k,i),mywk(0,k,i),nzwk(0,k,i),rK)
               call Lapvdvhom(phi(0,k,i),mywk(0,k,i),nzwk(0,k,i),rK)
            enddo
         enddo
            
 
!*******************************************************************!
!After State(1):
!   vor  : Omega_y
!   phi  : phi
!   psi  : psi
!   scal : T
!   nxwk : Omega_y
!   nzwk : dmy/dy
!   nywk : psi
!   mxwk : -
!   mywk : my
!   mzwk : -
!   ten12: RHS(Omega_y)(n-1)
!   ten13: RHS(phi)(n-1) 
!   ten23: d(psi)/dy (from last substep)
!   rhst : RHS(T)(n-1)
!   drho : -
!*******************************************************************!

!-----------------(2) Solver for mx, mz ----------------------------!
! Solving eq system 
         do i=pb,pe
            do k=0,mz1
               temp1 = tempa(k,i)
               temp2 = tempb(k,i)
               do j=0,my-1
                  ! Solving mx and mz
                  ! mx=tempa*dmy/dy-tempb*Omega_y
                  mxwk(j,k,i)  = nzwk(j,k,i)*temp1-nxwk(j,k,i)*temp2
                  ! mz=tempb*dmy/dy+tempa*Omega_y
                  mzwk(j,k,i)  = nzwk(j,k,i)*temp2+nxwk(j,k,i)*temp1
               enddo
            enddo
         enddo
!
!Save scal into drho
         do i=pb,pe
           do k=0,mz1
             do j=0,my-1
                drho(j,k,i)=scal(j,k,i)  !save T in drho too
             enddo
           enddo
         enddo
!        Now build u_i+RHS(1-1)*Dt*eps_i
!        -------------------------------
         do i=pb,pe
           do k=0,mz1
             do j=0,my-1
               vor(j,k,i)=vor(j,k,i)+ten12(j,k,i)*dtxi!Omega_y
               phi(j,k,i)=phi(j,k,i)+ten13(j,k,i)*dtxi!PSI
               scal(j,k,i)=scal(j,k,i)+rhst(j,k,i)*dtxi!scal 
!save buffer with RHS of Scal
               chwk2(j,k,i)= rhst(j,k,i)*dtxi
             enddo
           enddo
         enddo
!Save scal into drho
         do i=pb,pe
           do k=0,mz1
             do j=0,my-1
               rhst(j,k,i)=drho(j,k,i)  !save T in rhst too
             enddo
           enddo
         enddo
!
!.............................................................!
!        Prepare d(psi)/dy 
          do i=pb,pe
             call deryr2(nywk(0,0,i),ten23(0,0,i),mz)
             !!ten23= d(psi)/dy= d(nywk)/dy -- F-F
          enddo
!
!*******************************************************************!
!After State(2):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T**
!   nxwk : Omega_y
!   nzwk : dmy/dy
!   nywk : psi
!   mxwk : mx
!   mywk : my
!   mzwk : mz
!   chwk:  
!   chwk2: RHS(RKstep of T)  
!   ten12: RHS(Omega_y)(n-1)
!   ten13: RHS(phi)(n-1)
!   ten23: d(psi)/dy
!   rhst : T
!   drho : T
!*******************************************************************!

!-----------------(3) Obtain rho*u ----------------------------!
!
!Warning:
!Noted that its really important to use good value for d(psi)/dy
         do i=pb,pe
             do k=0,mz1
                do j=0,my1
                   !rhou=mx+ikx*psi
                   mxwk(j,k,i)=mxwk(j,k,i)+xalp(i-1)*nywk(j,k,i)    
                   !rhow=mz+ikz*psi
                   mzwk(j,k,i)=mzwk(j,k,i)+xbet(k)*nywk(j,k,i)    
                   !rhov=my+d(psi)/dy
                   mywk(j,k,i)=mywk(j,k,i)+ten23(j,k,i)    
               enddo
            enddo
         enddo
!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=16, file='u00T00.'//ext1//ext2)
                do j=0,my1
                   write(16,'(3e20.8)') y(j+1), u00wk(j),
     .                                  real(rhst(j,0,1))
                enddo
                close(unit=16)
             endif
!!------------------------------------------------------------------! 




!!---------- Adding modes 00 to rhou and rhow------------------------!
!---------- Only master computes 00  modes
        if (myid.eq.0) then
           do j=0,my1  
              mxwk(j,0,1) = u00wk(j)
              mywk(j,0,1) = v00wk(j) !rhov00 
              mzwk(j,0,1) = w00wk(j)
           enddo

!       !Now build in u00/w00 the variable u00**/w00**
! Set BC rho u remains unchanged at boundaries
           do j=0,my1
              u00wk(j)=u00wk(j)+dtxi*rf0u(j)
              w00wk(j)=w00wk(j)+dtxi*rf0w(j)
           enddo
        endif
!!*******************************************************************!
!After State(3):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T**
!   nxwk : Omega_y
!   nzwk : dmy/dy
!   nywk : psi
!   mxwk : rhou
!   mywk : rhov
!   mzwk : rhow
!   chwk :  -  
!   ten12: RHS(Omega_y)(n-1)
!   ten13: RHS(phi)(n-1)
!   ten23: d(psi)/dy
!   rhst : T
!   drho : T
!*******************************************************************!
!-----------------(4) Obtain T derivatives---------------------------!
!aaf      calculate dT/dy: ten13=dT/dy -- F-F
         do i=pb,pe
            call deryr2(rhst(0,0,i),ten13(0,0,i),mz)
         enddo
        !calculate laplacian of T and save it in "nywk"
         do i=pb,pe
            call laplacian(rhst(0,0,i),nywk(0,0,i),mz,i)
        enddo 

!---------!Calculate dT/dz,dT/dx------------------------------------!
         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  !dT/dz=i kz T
                  ten23(j,k,i)=xbet(k)*rhst(j,k,i)
                  !dT/dx=i kx T
                  ten12(j,k,i)=xalp(i-1)*rhst(j,k,i)
               enddo
            enddo
        enddo

!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=14, file='foupin.'//ext1//ext2)
                do j=0,my1
                   write(14,'(7e20.8)') y(j+1),real(mxwk(j,0,1)),
     .           real(mywk(j,0,1)),real(mzwk(j,0,1)),
     .           real(mxwk(j,1,1)),real(mywk(j,1,1)),real(mzwk(j,1,1))
                enddo
                close(unit=14)
             endif
!!------------------------------------------------------------------! 




!*******************************************************************!
!After State(4):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T**
!   nxwk : Omega_y
!   nzwk : dmy/dy
!   nywk : lap(T)
!   mxwk : rhou
!   mywk : rhov
!   mzwk : rhow
!   chwk:  d(psi)/dy
!   ten12: dT/dx
!   ten13: dT/dy
!   ten23: dT/dz
!   rhst : T
!   drho : T
!*******************************************************************!
!------------(5) FOU-PHYS(CALCULATIONS)-FOU--------------------------!
!Here we compute the new Dt                                          !
        call foup(mxwk,mywk,mzwk,nxwk,nywk,nzwk,
     .          rhst,ten12,ten13,ten23,chwk,myid,rkstep)
! outputs: 
! mxwk = u
! mywk = v
! mzwk = w
! nxwk = rhouu
! nywk = rhovv
! nzwk = rhoww
! ten12= rhouv
! ten13= rhouw
! ten23= rhovw
! rhst = RHS(T)
!---------------------------------------------------------------------!
!Computes Dt*gamma, Dt*xi before calculating this substep CFL(calculated
!within fouphysfou subroutine
         dtgamma = Deltat*gama(rkstep)
         dtxi    = Deltat*xi(rkstep)
         dtbeta  = Deltat*beta(rkstep)

!-------------(6a) Evolving Temperature ----------------------------!       
!              call rkstepexp(scal(0,0,i),rhst(0,0,i), mz,i,dtgamma)
          do i=pb,pe
             do k=0,mz1
                do j=0,my1
                  scal(j,k,i) = scal(j,k,i)  + dtgamma*rhst(j,k,i)
                  chwk2(j,k,i)= chwk2(j,k,i) + dtgamma*rhst(j,k,i)
                enddo
             enddo
          enddo

!*******************************************************************!
!After State(6a):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T(i+1) <<--<
!   nxwk : rhouu
!   nzwk : rhoww
!   nywk : rhovv
!   mxwk : u
!   mywk : v
!   mzwk : w
!   chwk:  -
!   ten12: rhouv
!   ten13: rhouw
!   ten23: rhovw
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!   drho : T                                                        
!*******************************************************************!
!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=18, file='foupout.'//ext1//ext2)
                do j=0,my1
                   write(18,'(7e20.8)') y(j+1),real(nxwk(j,0,1)),
     .                         real(ten12(j,0,1)),real(ten13(j,0,1)),
     .           real(nxwk(j,1,1)),real(ten12(j,0,1)),real(ten13(j,0,1))
                enddo
                close(unit=18)
             endif
!!------------------------------------------------------------------! 




!-------------(6b) Calculate Nx,Ny,Nz -  -------------------------!       
!        6b.1)Calculate d(rhovv)/dy
         do i=pb,pe
            call deryr2(nywk(0,0,i),nywk(0,0,i),mz)
         enddo
!        6b.2.1)Calculate Ny=-d(rhouv)/dx-d(rhovv)/dy-d(rhovw)/dz
         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                   !Ny
                   nywk(j,k,i)=-xalp(i-1)*ten12(j,k,i)-nywk(j,k,i)-
     .                         xbet(k)*ten23(j,k,i)
               enddo
            enddo
         enddo 
!        6b.2.2) Calculate new derivatives on y
         do i=pb,pe
            call deryr2(ten12(0,0,i),ten12(0,0,i),mz)!d(rhouv)/dy
            call deryr2(ten23(0,0,i),ten23(0,0,i),mz)!d(rhovw)/dy
         enddo
!       6b.3) Calculate Nx,Nz
          do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                   !Nx
                   nxwk(j,k,i)=-xalp(i-1)*nxwk(j,k,i)-ten12(j,k,i)-
     .                         xbet(k)*ten13(j,k,i)
                   !Nz
                   nzwk(j,k,i)=-xalp(i-1)*ten13(j,k,i)-ten23(j,k,i)-
     .                         xbet(k)*nzwk(j,k,i)
               enddo
            enddo
         enddo 

!----------------------00 MODES-------------------------------------!
!Save / compute terms needed for evolution of 00 modes
         if (myid.eq.0) then !only proc 0 does this
            do j=0,my1
               rf0u(j)= real(nxwk(j,0,pb))
               rf0w(j)= real(nzwk(j,0,pb))
            enddo
         endif  

!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=18, file='nxwk00.'//ext1//ext2)
                do j=0,my1
                   write(18,'(2e20.8)') y(j+1), rf0u(j)
                enddo
                close(unit=18)
             endif
!!------------------------------------------------------------------! 


!-------------------------------------------------------------------!
        if (rkstep.eq.3.and.ihist.eq.1.and.myid.eq.0) then
             H00u=0.
             H00w=0.
             do j=2,my-2
                H00u = H00u + rf0u(j-1)*(y(j+1)-y(j-1))
                H00w = H00w + rf0w(j-1)*(y(j+1)-y(j-1))
             enddo
             H00u = 0.25*H00u
             H00w = 0.25*H00w
         endif
!        
!*******************************************************************!
!After State(6b):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T(i+1) <<--< RIP buffer
!   nxwk : Nx
!   nzwk : Nz
!   nywk : Ny
!   mxwk : u
!   mywk : v
!   mzwk : w
!   chwk:  -
!   ten12: -
!   ten13: -
!   ten23: -
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!

!-------------(6c) Calculate Sij-----------------------------------!       
! Calculate derivatives on y
         do i=pb,pe
            call deryr2(mxwk(0,0,i),ten12(0,0,i),mz)!du/dy
            call deryr2(mzwk(0,0,i),ten23(0,0,i),mz)!dw/dy
         enddo
!!!-------------------------------------------------------------------! 
!debug
!Write all before tauij
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=19, file='Sij.'//ext1//ext2)
                do j=0,my1
                   write(19,'(7e20.8)') y(j+1), real(mxwk(j,0,1)),
     .     real(mywk(j,0,1)),real(mzwk(j,0,1)),real(mxwk(j,1,1)),
     .     real(mywk(j,1,1)),real(mzwk(j,1,1))
                enddo
                close(unit=19)
             endif
!
!!-------------------------------------------------------------------!


! Calculate all tensor Sij
          !S12, S13 and S23
          do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  !S12
                  ten12(j,k,i)=(ten12(j,k,i)+
     .                         xalp(i-1)*mywk(j,k,i))/2.0D0
                  !S13
                  ten13(j,k,i)=(xbet(k)*mxwk(j,k,i)+
     .                          xalp(i-1)*mzwk(j,k,i))/2.0D0
                  !S23
                  ten23(j,k,i)=(ten23(j,k,i)+
     .                         xbet(k)*mywk(j,k,i))/2.0D0
               enddo
            enddo
          enddo
          !Now the S11 and S33
          do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  !S11
                  mxwk(j,k,i)=xalp(i-1)*mxwk(j,k,i) 
                  !S33
                  mzwk(j,k,i)=xbet(k)*mzwk(j,k,i) 
               enddo
            enddo
          enddo
          ! and last: S22
          do i=pb,pe
             call deryr2(mywk(0,0,i),mywk(0,0,i),mz)
          enddo

!*******************************************************************!

!*******************************************************************!
!After State(6c):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T(i+1) <<--< RIP buffer
!   nxwk : Nx
!   nzwk : Nz
!   nywk : Ny
!   mxwk : S11
!   mywk : S22
!   mzwk : S33
!   chwk:  -
!   chwk2: T(n+1)-T(n)
!   ten12: S12
!   ten13: S13
!   ten23: S23
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!
!*******************************************************************!
!-------------------------------------------------------------------! 
!debug
!Write all before tauij
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=19, file='tauijin.'//ext1//ext2)
                do j=0,my1
                   write(19,'(7e20.8)') y(j+1), real(mxwk(j,0,1)),
     .     real(mywk(j,0,1)),real(mzwk(j,0,1)),real(ten12(j,0,1)),
     .     real(ten13(j,0,1)),real(ten23(j,0,1))
                enddo
                close(unit=19)
             endif
!
!!-------------------------------------------------------------------!

!!------------(7) Go to Phys and compute TAUij and Drho (and back)----!
        call tauij(mxwk,mywk,mzwk,ten12,ten13,ten23,
     .                  drho,chwk2,chwk,myid,rkstep)

!! outputs: 
!! mxwk = tau11 
!! mywk = tau22
!! mzwk = tau33
!! ten12= tau12
!! ten13= tau13
!! ten23= tau23
!! drho = drho
!!---------------------------------------------------------------------!
!!!-------------------------------------------------------------------! 
!debug
!Write all after tauij
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=19, file='tauijout.'//ext1//ext2)
                do j=0,my1
                   write(19,'(7e20.8)') y(j+1), real(mxwk(j,0,1)),
     .     real(mywk(j,0,1)),real(mzwk(j,0,1)),real(ten12(j,0,1)),
     .     real(ten13(j,0,1)),real(ten23(j,0,1))
                enddo
                close(unit=19)
             endif
!
!!-------------------------------------------------------------------!



!----------(8a) Calculate Mx,My,Mz------------------------------------!
!         tau22 only needed as d(tau22)/dy and goes for My
!         First we can calculate My
       do i=pb,pe
          call deryr2(mywk(0,0,i),mywk(0,0,i),mz) !now mywk=d(tau22)/dy
       enddo
       do i=pb,pe
          do k=0,mz1
             do j=0,my-1
                 !My=ikxTau12+d(tau22)/dy+ikztau23
                  mywk(j,k,i)=xalp(i-1)*ten12(j,k,i)+mywk(j,k,i)+
     .                        xbet(k)*ten23(j,k,i) !My ready
             enddo
          enddo
       enddo 
       !compute derivatives on y needed for Mx and Mz
       do i=pb,pe
          call deryr2(ten12(0,0,i),ten12(0,0,i),mz) !now ten12=d(tau12)/dy
          call deryr2(ten23(0,0,i),ten23(0,0,i),mz) !now ten23=d(tau23)/dy
       enddo
       
       !Compute Mx, Mz
       do i=pb,pe
          do k=0,mz1
             do j=0,my-1
                mxwk(j,k,i)=xalp(i-1)*mxwk(j,k,i)+ten12(j,k,i)+
     .                      xbet(k)*ten13(j,k,i) !Mx
                mzwk(j,k,i)=xbet(k)*mzwk(j,k,i)+ten23(j,k,i)+
     .                      xalp(i-1)*ten13(j,k,i) !Mz
             enddo
          enddo
       enddo 

!!----------------------00 MODES-------------------------------------!
!!!Save / compute terms needed for evolution of 00 modes
         if (myid.eq.0) then !only proc 0 does this
            do j=0,my1
                rf0u(j)=rf0u(j)+ire*real(mxwk(j,0,pb)) !RHS(rhou00)
                rf0w(j)=rf0w(j)+ire*real(mzwk(j,0,pb))!RHS(rhow00)
            enddo
         endif  

!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=19, file='rhsu00.'//ext1//ext2)
                do j=0,my1
                   write(19,'(2e20.8)') y(j+1), rf0u(j)
                enddo
                close(unit=19)
             endif
!
!!-------------------------------------------------------------------!

!*******************************************************************!
!After State(8a):
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi
!   scal : T(i+1) <<--< RIP buffer
!   nxwk : Nx
!   nzwk : Nz
!   nywk : Ny
!   mxwk : Mx
!   mywk : My
!   mzwk : Mz
!   chwk:  -
!   ten12: -
!   ten13: -
!   ten23: -
!   drho : drho                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!



!------------------------------------------------------------------!
!------------(8b) Obtain psi(n+1)-----------------------------------!
!1) First solve 00 mode of psi
!   Only for myid.eq.0
         if (myid.eq.0) then
!        1.1) Integrate drho00
              !F(y) = F(-Ly) + int(f(y),y,[-Ly,y])
              !cumtrapz calculates the integral part of last formula
              do j=0,my1
                  drho00(j) = real(drho(j,0,1))
              enddo 
              !CUMTRAPZ:
              rf0v(0) = 0d0
              do j=0,my1-1
                 rf0v(j+1) = rf0v(j)+0.5d0*(y(j+2)-y(j+1))*
     .                (drho00(j+1)+drho00(j))
              enddo
              trapz = rf0v(my1) !last element is the area of f!
              !now rf0v keeps the primitive from -Ly to y of drho
              !We will force that rhov00(-Ly) = -rhov00(Ly):
              !    + Implying that:
              !         int(f(y),[-Ly,Ly])=F(Ly)-F(-Ly)=-2*F(-Ly)
              !         So that:
              !             F(-Ly) =-0.5d0* trapz(f(y),-Ly,Ly)
              do j=0,my1
                  rf0v(j)=rf0v(j)-0.5d0*trapz
              enddo
!rf0v  keeps the integral of drho00

!    1.2) Update v00
             do j=0,my1 
                v00wk(j)=-alpha(rkstep)/beta(rkstep)*v00wk(j)
     .                                        -rf0v(j)/dtbeta
             enddo
         endif
!----------------------------------------------------------------------
!!!-------------------------------------------------------------------! 
!debug
             if (myid.eq.0.and.istep.le.ideb) then
                open(unit=12, file='rf0v.'//ext1//ext2)
                do j=0,my1
                   write(12,'(4e24.16)') y(j+1), real(drho(j,0,1)),
     .             rf0v(j),real(drho(j,1,1))
                enddo
                close(unit=12)
             endif
!!------------------------------------------------------------------! 



!2) Solve other modes...
!New variable:
!   2.1) Lap(eta)=drho
!   2.2) solve Lapvdv of Eta (BC's??)
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               call Lapvdvhom(drho(0,k,i),drho(0,k,i),ten23(0,k,i),rK)
            enddo
         enddo
!   2.3)Update psi(n+1), using eta instead of drho
          do i=pb,pe
             do k=0,mz1
               do j=0,my1
                  psi(j,k,i)=-drho(j,k,i)/dtbeta -
     .                       psi(j,k,i)*alpha(rkstep)/beta(rkstep)
               enddo
             enddo
          enddo 
! Now:
! ten23 = not usesful...
! psi   = psi(n+1)
! v00   = rhov00(n+1)
! v00wk = rhov00(n+1)

!----------(9) Construct the RHS of Omega_y and phi--------------------!
! 
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               do j=0,my-1
                 mywk(j,k,i)= -rK*(ire*mywk(j,k,i)+nywk(j,k,i)) !aux1
                 nywk(j,k,i)=xalp(i-1)*(nxwk(j,k,i)+ire*mxwk(j,k,i))+
     .                       xbet(k)*(nzwk(j,k,i)+ire*mzwk(j,k,i))
!                  !aux2
               enddo
            enddo
         enddo
!we need the derivative on y of aux2
         do i=pb,pe
              call deryr2(nywk(0,0,i),nywk(0,0,i),mz)
         enddo 
! RHS(phi)=aux1-d(aux2)/dy
! RHS(Omega_y)=ikzNx-ikxNx+ire*(ikzMx-ikxMz) 
         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
!                  !RHS(Omega_y)
                  ten12(j,k,i)=xbet(k)*(nxwk(j,k,i)+ire*mxwk(j,k,i))-
     .                       xalp(i-1)*(nzwk(j,k,i)+ire*mzwk(j,k,i))
!                  !RHS(phi)
                  ten13(j,k,i)=mywk(j,k,i)-nywk(j,k,i)
               enddo
            enddo
         enddo

!*******************************************************************!
!After state(9)
!   vor  : Omega_y**
!   phi  : phi**
!   psi  : psi(n+1)
!   scal : T(n+1) <<--< RIP buffer
!   nxwk : -
!   nzwk : -
!   nywk : -
!   mxwk : -
!   mywk : -
!   mzwk : -
!   chwk:  -
!   ten12: RHS(Omega_y)
!   ten13: RHS(phi)
!   ten23: dpsi(n+1)/dy
!   drho : -                                                       
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!
!
!-------------(10) Evolving Omega_y and phi-------------------------!       
              !call rkstepexp(vor(0,0,i),ten12(0,0,i), mz,i,dtgamma)
          do i=pb,pe
             do k=0,mz1
                do j=0,my1
                   vor(j,k,i) = vor(j,k,i) + dtgamma*ten12(j,k,i)
                enddo
             enddo
          enddo

              !call rkstepexp(phi(0,0,i),ten13(0,0,i), mz,i,dtgamma)
          do i=pb,pe
             do k=0,mz1
                do j=0,my1
                   phi(j,k,i) = phi(j,k,i) + dtgamma*ten13(j,k,i)
                enddo
             enddo
          enddo


!!*******************************************************************!
!After state(10)
!   vor  : Omega_y(n+1)
!   phi  : phi(n+1)
!   psi  : psi(n+1)
!   scal : T(n+1) 
!   nxwk : -
!   nzwk : -
!   nywk : -
!   mxwk : -
!   mywk : -
!   mzwk : -
!   chwk:  -
!   ten12: RHS(Omega_y)
!   ten13: RHS(phi)
!   ten23: dpsi(n+1)/dy
!   drho : -                                                       
!   rhst : RHS(T) 
!*******************************************************************!
!===================================================================!
! ----------------------------------------------------------------------!
!             this is the block for the 00 mode
!                   Only master does it
! ----------------------------------------------------------------------!
        if (myid.eq.0) then
! evolve modes 00 por rhou, rhow
           do j=0,my1
             u00wk(j)=u00wk(j)+dtgamma*rf0u(j)
             w00wk(j)=w00wk(j)+dtgamma*rf0w(j)
           enddo
!compute mass
           massu=0.
           massw=0.
           do j=0,my-1
              massu = massu+u00wk(j)*hy(j+1)
              massw = massw+w00wk(j)*hy(j+1)
           enddo

        endif   ! 00 modes


 10    continue  !!! End of R-K sub-step


         if (myid.eq.0) then
            if (ihist.eq.1) then
               if (istep.eq.1) then
                   ctm2=0.
                   ttm =0.
               endif
324            format(8(d22.14))
               write(29,324) time, Deltat,dm,massu,massw,
     .         ire ,MPI_WTIME()+iter_time, commtimer-ctm2
               call flush(29)
325            format(i5,12(d14.6))
               write(*,325)istep,time,Deltat,dm,massu,massw,
     .                 ire,H00u,H00w,MPI_WTIME()+iter_time,
     .                    transtimer-ttm,commtimer-ctm2
            endif
            ctm2 = commtimer
            ttm  = transtimer 
         endif

c          	! time:
         time=time+Deltat

         if(icfl .eq.1)   icfl   = 0
         if(ihist.eq.1)  ihist  = 0

         if (myid.eq.0) then
           totaltimer=totaltimer+MPI_WTIME()
           print *,MPI_WTIME()+iter_time
         endif

         if (istep.eq.1) ctm=commtimer
         if (istep.gt.1.and.istep.lt.5) then
            temporal(istep-1) = MPI_WTIME()+ iter_time
            comun   (istep-1) = commtimer-ctm
            ctm = commtimer
         endif
         call MPI_BARRIER(MPI_COMM_CALC,ierr)
         
         call MPI_allreduce(nanerror,nanproblem,1,MPI_INTEGER,MPI_MAX,
     .                              MPI_COMM_CALC,ierr)
         if (nanproblem.gt.0) go to 100


 30   continue
 

      if (myid.eq.0) then
         open (11,file='stemp.dat',status='unknown') 

c-----  means

        temporal(1)  = (temporal(2)+temporal(3)+temporal(1))/3
        comun(1)     = (comun(1)+comun(2)+comun(3) )/3 
        write(11,*) temporal(1),comun(1),
     .              commtimer/totaltimer,totaltimer
        print *,"Total time: ",totaltimer
        
        print *,"Trans. time: ",transtimer
        print *,"Comm. time: ",commtimer
        print *,"Comm/Total: ",commtimer/totaltimer
        print *,"Trans/Total: ",transtimer/totaltimer
        close(11)      
      endif
      
      write(*,*) nanproblem

 100   if (nanproblem.gt.0) then
         write(*,*) 'Corrupted Memory'
         write(*,*) 'OH MY GOODNESS...'
      endif


      !deallocate(tempa,tempb)
      

      end

!******************************************************************!
!    FOUPHYSFOU1                                                  ! 
!    ------------                                                  ! 
!    From Inputs->go to Phys-> do operations->go back with outputs !     
!                                                                  !
!       VARIABLE       INPUT          OUTPUT                       !
!       --------------------------------------                     !
!       up1     =      rhou             u                          !
!       up2     =      rhov             v                          !
!       up3     =      rhow             w                          !
!       ten11   =       -             rhouu                        !
!       ten22   =      lap(T)         rhovv                        !
!       ten33   =       -             rhoww                        !
!       rhst    =        T            RHS(T)                       !
!       ten12   =      dT/dx          rhouv                        !
!       ten13   =      dT/dy          rhouw                        !
!       ten23   =      dT/dz          rhovw                        !
!                                                                  !
!******************************************************************!
      subroutine foup(up1,up2,up3,ten11,ten22,ten33,
     .           rhst,ten12,ten13,ten23,chwk,myid,rkstep)
      use MPI_GROUPS
      use point
      use diag
      use tem
      use fis
      use timacc
      use wave
      use timers
      use wkhvect
      use cnan

       
      implicit none

      include "mpif.h"
      include "ctes3D"


c     ---------------------- Variables ----------------------
      
      complex*8, dimension(0:my1,0:mgalz-1,pb:pe):: 
     .          up1,up2,up3,ten11,ten22,ten33,ten12,
     .          ten13,ten23,rhst,chwk

      real(4) uner(9),aa
      integer myid,rkstep,i,j,k,kk,ierr
      
c     ---------------------- Program ----------------------      


!  local transpose yz to zy, and adding zeros between positive
!  and negative wavenumbers.
      call localyz2zy(up1,up1,chwk) !rhou
      call localyz2zy(up2,up2,chwk) !rhov
      call localyz2zy(up3,up3,chwk) !rhow
!      !TEMPERATURE
      call localyz2zy(ten22,ten22,chwk) !Lap(T) 
      call localyz2zy(ten12,ten12,chwk) !dT/dx         
      call localyz2zy(ten13,ten13,chwk) !dT/dy  
      call localyz2zy(ten23,ten23,chwk) !dT/dz        
      call localyz2zy(rhst,rhst,chwk) !T      
       
! inverse transform (fou --> fis @z),does the compacting aswell
      do i = pb,pe
         call fourz(up1(0,0,i),1)    !rhou
         call fourz(up2(0,0,i),1)    !rhov
         call fourz(up3(0,0,i),1)    !rhow
         !TEMPERATURE
         call fourz(ten22(0,0,i),1)  !lap(T)
         call fourz(ten12(0,0,i),1)  ! dT/dx
         call fourz(ten13(0,0,i),1)  ! dT/dy
         call fourz(ten23(0,0,i),1)  ! dT/dz
         call fourz(rhst(0,0,i),1)  ! T
      enddo

!      nanerror = 0
!
!      call check_nan(up1,my*mgalz*mmp*2,nanerror)
!      call check_nan(up2,my*mgalz*mmp*2,nanerror)
!      call check_nan(up3,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten22,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten12,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten13,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten23,my*mgalz*mmp*2,nanerror)
!      call check_nan(rhst,my*mgalz*mmp*2,nanerror)
!
!      if (nanerror.eq.1) then
!         write(*,*) 'NAN Found in ',myid,' before changes'
!         stop
!      endif

!  change plane to line(1:mx,lb:le)
!  and distribute package of lines to processors
      call chpl2ln(up1,up1,chwk,myid)
      call chpl2ln(up2,up2,chwk,myid)
      call chpl2ln(up3,up3,chwk,myid)
!temperature
      call chpl2ln(ten22,ten22,chwk,myid)
      call chpl2ln(ten12,ten12,chwk,myid)
      call chpl2ln(ten13,ten13,chwk,myid)
      call chpl2ln(ten23,ten23,chwk,myid)
      call chpl2ln(rhst,rhst,chwk,myid)

!==========================================================================!
! before this call upi,ten22,ten12,ten13, ten23 and rhst
! are still lines in fourier space x (phys z).
      call phys1(up1,up2,up3,ten11,ten22,ten33,
     .     ten12,ten13,ten23,rhst,myid,rkstep)
!==========================================================================!

c  ---------- back to the yz planes an fft 

      call chln2pl(up1,up1,chwk,myid)     !u
      call chln2pl(up2,up2,chwk,myid)     !v
      call chln2pl(up3,up3,chwk,myid)     !w
      call chln2pl(ten11,ten11,chwk,myid) !rhouu
      call chln2pl(ten22,ten22,chwk,myid) !rhovv  
      call chln2pl(ten33,ten33,chwk,myid) !rhoww
      call chln2pl(ten12,ten12,chwk,myid) !rhouv  
      call chln2pl(ten13,ten13,chwk,myid) !rhouw
      call chln2pl(ten23,ten23,chwk,myid) !rhovw
      call chln2pl(rhst,rhst,chwk,myid)   !RHS(T)
!
!!
!      call check_nan(rhst,my*mgalz*mmp*2,nanerror)
!      call check_nan(up1,my*mgalz*mmp*2,nanerror)
!      call check_nan(up2,my*mgalz*mmp*2,nanerror)
!      call check_nan(up3,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten11,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten22,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten33,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten12,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten13,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten23,my*mgalz*mmp*2,nanerror)
!
!      if (nanerror.eq.1) then
!         write(*,*) 'NAN Found in ',myid,' after changes'
!         stop
!      endif

!  convert to fourier Z now.
      do i=pb,pe
         call fourz(up1(0,0,i),-1)   ! u
         call fourz(up2(0,0,i),-1)   ! v
         call fourz(up3(0,0,i),-1)   ! w
         call fourz(ten11(0,0,i),-1)  !rhouu
         call fourz(ten22(0,0,i),-1)  !rhovv
         call fourz(ten33(0,0,i),-1)  !rhoww
         call fourz(ten12(0,0,i),-1)  !rhouv
         call fourz(ten13(0,0,i),-1)  !rhouw
         call fourz(ten23(0,0,i),-1)  !rhovw
         call fourz(rhst(0,0,i),-1)   !RHS(T)
      enddo
!transpose back
      call localzy2yz(up1,up1,chwk)
      call localzy2yz(up2,up2,chwk)
      call localzy2yz(up3,up3,chwk) 
      call localzy2yz(ten11,ten11,chwk) 
      call localzy2yz(ten22,ten22,chwk) 
      call localzy2yz(ten33,ten33,chwk) 
      call localzy2yz(ten12,ten12,chwk)   
      call localzy2yz(ten13,ten13,chwk)   
      call localzy2yz(ten23,ten23,chwk)   
      call localzy2yz(rhst,rhst,chwk)

!------------------------------------------------------!
!!aaf------- momentum thickness calculation---------!
!! dm will keep the momentum thickness
              !dm =int(1/4-(um/Du)**2)
         if (myid.eq.0.and.rkstep.eq.1) then
            if (ihist.eq.1) then
              dm=0d0
              do j=1,my-2
                  dm=dm+trp(j+1)*((real(up1(j+1,0,1)))**2+
     .                          (real(up1(j-1,0,1)))**2)
              enddo
              dm=dm+trp(1)*((real(up1(1,0,1)))**2+
     .        (real(up1(0,0,1)))**2)+
     .        trp(my)*((real(up1(my1,0,1)))**2+
     .                  (real(up1(my-2,0,1)))**2) 

              dm =1d0/4.0d0*(y(my)-y(1))-
     .           dm/(real(up1(my1,0,1))-real(up1(0,0,1)))**2
             endif
          endif
!------------------------------------------------------!

           
      end


!*********************************************************************!
!                                                                     !
!         transform to physical and do operations                     !
!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!                                                                     !
!                                                                     !
!*********************************************************************!
      subroutine phys1(up1,up2,up3,ten11,ten22,ten33,
     .           ten12,ten13,ten23,rhst,myid,rkstep)
      use point
      use tem
      use fis
      use timacc
      use diag 
      use wkhvect
      use MPI_GROUPS

      implicit none

      include "mpif.h"
      include "ctes3D"

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
      real(8)  cflx,cfly,cflz,hxalp,hzbet,hyy,cfl0,reigmx1,aa
      real(4) ::  up1(0:mx-1,lb:le),up2(0:mx-1,lb:le), 
     .        up3(0:mx-1,lb:le),ten11(0:mx-1,lb:le),
     .        ten22(0:mx-1,lb:le),ten33(0:mx-1,lb:le),
     .        ten12(0:mx-1,lb:le),ten13(0:mx-1,lb:le),
     .        ten23(0:mx-1,lb:le),
     .        rhst (0:mx-1,lb:le)

      real(8) duma,dumb,dumc,duma1,dumb1,dumc1
      real(8) pi
      !aaf adding Peclet
      real(8) peclet,sigma

!      ----------------------- Programa ------------------------      
!initialization
      !aaf Peclet
      peclet=re     !for debuggin consider peclet=reynodls
      !sigma=0.7D0
      sigma=0D0
      cflx = 0.
      cfly = 0.
      cflz = 0.
      cfl0 = 0.
      uv0=0.
      uvL=0.
      pi = 4d0*atan(1d0)
!calculate geometric staff
      !1/Dx and 1/Dz
      hxalp=alp*(mgalx-1)/pi !2/Dx
      hzbet=bet*(mgalz-1)/pi !2/Dz
      hyy=hy(1)
      do j=2,my
         hyy = min(hyy,hy(j))
      enddo
      hyy = 1./hyy
! computing viscous Dt:
! computing viscous Dt inverse:
! Dt_visc=1/(2 nu)*Dy_min^2
      hyy=2/Re*(max(hyy,max(hxalp/2,hzbet/2)))**2

c   -------------------------------------------------

c        Move everything to pPP, line by line 


      duma = 0d0
      dumb = 0d0
      dumc = 0d0

      do 40 j = lb,le
       
c     copy lines...
         do i=0,mx-1
            up1wk(i) = up1(i,j) !rhou
         enddo
         do i=0,mx-1           
            up2wk(i) = up2(i,j) !rhov
         enddo
         do i=0,mx-1
            up3wk(i) = up3(i,j) !rhow
         enddo
         do i=0,mx-1
            ten22wk(i) = ten22(i,j) !lap(T)
         enddo
         do i=0,mx-1
            ten12wk(i) = ten12(i,j) !dT/dx
         enddo
         do i=0,mx-1
            ten13wk(i) = ten13(i,j) !dT/dy
         enddo
         do i=0,mx-1
            ten23wk(i) = ten23(i,j) !dT/dz
         enddo
         do i=0,mx-1
            rhstwk(i) = rhst(i,j)   !T
         enddo
!        !registering the mpi_wtime...
         if (myid.eq.0) then
            duma = duma+MPI_WTIME()
         endif

! convert lines in fourier to lines in physical for x
! upiwk8 means real 8 format
         call fourx(up1wk,up1wk8,1)    !rhou
         call fourx(up2wk,up2wk8,1)    !rhov
         call fourx(up3wk,up3wk8,1)    !rhow
         call fourx(ten22wk,ten22wk8,1) ! lap(T)
         call fourx(ten12wk,ten12wk8,1) ! dT/dx
         call fourx(ten13wk,ten13wk8,1) ! dT/dy
         call fourx(ten23wk,ten23wk8,1) ! dT/dz
         call fourx(rhstwk,rhstwk8,1) ! T
! ----------------------------------------------------------------------!
!                                                                       !
!         computing in physical domain ( line by line)                  !
!                                                                       !
! ----------------------------------------------------------------------!

!         !Calculate u,v,w
         do i=0,mgalx-1
            tmp1wk8(i)= up1wk8(i)*rhstwk8(i) !u for wk 
            up1wk(i)  = tmp1wk8(i)           !u for output

            tmp2wk8(i)= up2wk8(i)*rhstwk8(i) !v
            up2wk(i)  = tmp2wk8(i)           !v for output

            tmp3wk8(i)= up3wk8(i)*rhstwk8(i) !w for wk
            up3wk(i)  = tmp3wk8(i)           !w for output
         enddo

!         !Calculate rhs(T)
         do i=0,mgalx-1
            rhstwk(i)=-tmp1wk8(i)*ten12wk8(i)-tmp2wk8(i)*ten13wk8(i)-
     .      tmp3wk8(i)*ten23wk8(i)+ !-u_i*dT/dxi
     .      sigma/peclet*(rhstwk8(i)**sigma)*(ten12wk8(i)**2+
     .      ten13wk8(i)**2+ten23wk8(i)**2) + !(1/Pe)s*T^s*(dT/dxi dT/dxi)
     .      (1/peclet)*(rhstwk8(i))**(sigma+1)*ten22wk8(i)
!            1/pe*T**(sigma+1)*lap(T)
         enddo         
 
!===================================================================!
!Now we do not need Lap(T),DT/dxi anymore and can replaces variables
!===================================================================!

!         !Calculate rhouu(ten11),rhovv(ten22),rhoww(ten33)
!         !rhou*u_i*u_i=rhou_i*u_i
         do i=0,mgalx-1
            ten11wk(i)=up1wk8(i)*tmp1wk8(i) !rhouu ready for output
            ten22wk(i)=up2wk8(i)*tmp2wk8(i) !rhovv ready for output
            ten33wk(i)=up3wk8(i)*tmp3wk8(i) !rhoww ready for output
         enddo

!         !Calculate the cross products rhou_i*u_j
         do i=0,mgalx-1
            ten12wk(i)=up1wk8(i)*tmp2wk8(i) !rhouv ready for output
            ten13wk(i)=up1wk8(i)*tmp3wk8(i) !rhouw ready for output
            ten23wk(i)=up2wk8(i)*tmp3wk8(i) !rhovw ready for ouput
         enddo
!
! now lines of mgalx, each mgalz one height of y.
!========================================================!
!IMPORTANT POINT: all quantities are now PHYSICAL
!========================================================!
c    ------------------------------------------------
c    ----------- triple products 
! each sub-block of lines are in the same y. So that
! the total number of lines in each processor are
! sections in y. (j is the line index).
         jj = (j-1)/mgalz +1
! obtain the equivalent y index (height of block lines).
         if (jj.lt.0) jj=1
!      estimates maximum time step     !
! still within loop for each line
! now tmpiwk8 keeps the velocities
         if (rkstep.eq.1.and.icfl.eq.1) then
            cfl0=0.
            do i = 0,mgalx-1
               cflx = max(cflx,abs(tmp1wk8(i)))
               cfl0 = max(cfl0,abs(tmp2wk8(i)))
               cflz = max(cflz,abs(tmp3wk8(i)))
!cflx and cflz have same hx,hz for all domain
!but cfly needs to account for hy. 
            enddo
            cfly=max(cfly,cfl0/hy(jj))
         endif

c--------------------- back to F-P-P  

         call fourx(up1wk,up1wk8,-1)  ! u
         call fourx(up2wk,up1wk8,-1)  ! v
         call fourx(up3wk,up1wk8,-1)  ! w
         call fourx(ten11wk,up1wk8,-1)  ! rhouu
         call fourx(ten22wk,up1wk8,-1)  ! rhovv
         call fourx(ten33wk,up1wk8,-1)  ! rhoww
         call fourx(ten12wk,up1wk8,-1)  ! rhouv
         call fourx(ten13wk,up1wk8,-1)  ! rhouw
         call fourx(ten23wk,up1wk8,-1)  ! rhovw
         call fourx(rhstwk,up1wk8,-1)  ! RHS(T)
c --------    back to lines. We throw away the high modes (Dealiasing)

         if (myid.eq.0) then
            dumc = dumc-MPI_WTIME()
         endif
         do i=0,mx-1            
            rhst(i,j) = rhstwk(i) !RHS(T)
         enddo
         do i=0,mx-1
            up1(i,j) = up1wk(i)
         enddo
         do i=0,mx-1
            up2(i,j) = up2wk(i)
         enddo
         do i=0,mx-1            
            up3(i,j) = up3wk(i)
         enddo
         do i=0,mx-1            
            ten11(i,j) = ten11wk(i)
         enddo
         do i=0,mx-1            
            ten22(i,j) = ten22wk(i)
         enddo
         do i=0,mx-1            
            ten33(i,j) = ten33wk(i)
         enddo
         do i=0,mx-1            
            ten12(i,j) = ten12wk(i)
         enddo
         do i=0,mx-1            
            ten13(i,j) = ten13wk(i)
         enddo
         do i=0,mx-1            
            ten23(i,j) = ten23wk(i)
         enddo

         if (myid.eq.0) then
            dumc = dumc+MPI_WTIME()
         endif

 40    continue

      if (myid.eq.-2) then
         write(*,*) 'primer bucle',duma
         write(*,*) 'segun  bucle',dumb
         write(*,*) 'tercer bucle',dumc
      endif


      if (rkstep.eq.1.and.icfl.eq.1) then 
c-------------------- computes Deltat
!hxalp, hzbet are 1/Dx,1/Dz
         cflx = cflx*hxalp
         cflz = cflz*hzbet
!using a 3 factor in order to account for 3D mesh
!ex:Dt<CFL/(u/Dx+v/Dy+w/Dz) (it's a safety factor)
         cfl0 = max(cflx,max(cfly,cflz))

         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL8,MPI_MAX,
     .                     MPI_COMM_CALC,ierr)

!now take max of (Dt_diff^(-1), Dt_visc^(-1))
         reigmx1 = max(reigmx1,hyy)

!check viscous time vs diffusive time
         if (hyy.eq.reigmx1.and.myid.eq.0) then
            write(*,*) '------------Viscous Dt used----------' 
            write(*,*) 'Dt_visc = ',CFL/hyy,'myid= ', myid
         endif
!compute Deltat and dtr
         Deltat=CFL/reigmx1
         dtr=Re/Deltat

         if (reigmx1.lt.1e-1) then
            write(*,*) 'UGGGHH',myid,ihist,lb,le,pb,pe
         endif

      endif

      end
      
      
c----------------------------------------------------------------------
c    Advances the 00 modes (linear term explicit, order 3)
c----------------------------------------------------------------------

      subroutine rk00exp(f,rhs,dtgamma)
      use fis

      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real(8)  f(my),rhs(my)
      real(8)  dtgamma
c     ----------------------- work --------------------
      integer i,j
! ----------------------------------------------------------
          do j=1,my 
            f(j) = f(j)+dtgamma*rhs(j) 
            !new u00/w00:
          enddo
! not evolving the BC
      end



!************************************************************************
!    Subroutine rkstepexp                                               !
!    by aaf                                                             !
!    makes the explicit rkstep                                          !
!                                                                       !
!    size(u)=(2,my,m)                                                   !
!    adds the rhs term multiplied by dtgamma(i)                         ! 
!    Input:                                                             !
!         u : u(i)+Dt*xi(i)*RHS(u,i-1)                                  !
!        rhs : RHS(i)                                                   !
!                                                                       !
!    Output:                                                            !
!         u: u(i+1)                                                     !
!                                                                       !
!************************************************************************

      subroutine rkstepexp(u,rhs,m,xpl,dtgamma)
      use wave

      implicit none
      include "ctes3D"
      integer m,i,j,k,xpl

      real(4) :: u(2,my,m),rhs(2,my,m),tmp
      real(8) dtgamma

!------------------------------------------------------------!
       do k=1,m
          do j=1,my
!         !update the phi/vor
              u(1,j,k)=u(1,j,k)+dtgamma*rhs(1,j,k)
          enddo
!         imaginary part       
          do j=1,my
          !tmp keeps the nl rhs
              u(2,j,k)=u(2,j,k)+dtgamma*rhs(2,j,k)
          enddo
      enddo
 
      end subroutine rkstepexp

!-----------------------------------------------------------------------!


!***********************************************************************!
!   Make laplacian from u and give it through du
!                                                                       !
!************************************************************************

      subroutine laplacian(u,du,m,xpl)
      use matrices
      use wave

      implicit none
      include "ctes3D"

!--------------Variables----------------------!
      integer m,i,j,k,xpl
      real(4) :: u(2,my,m),du(2,my,m),tmp
      real(8) :: wk1(my),wk2(my)
      real(8) coef

!--------------Program----------------------!
      do k=1,m

        wk1(1) = dt22(3,1)*u(1,1,k)+dt22(4,1)*u(1,2,k)+ 
     $    dt22(5,1)*u(1,3,k)
       
        wk2(1) = dt22(3,1)*u(2,1,k)+dt22(4,1)*u(2,2,k)+ 
     $      dt22(5,1)*u(2,3,k)
       
       wk1(2) = dt22(2,2)*u(1,1,k)+dt22(3,2)*u(1,2,k)+ 
     $      dt22(4,2)*u(1,3,k)+dt22(5,2)*u(1,4,k)
       
       wk2(2) = dt22(2,2)*u(2,1,k)+dt22(3,2)*u(2,2,k)+ 
     $      dt22(4,2)*u(2,3,k)+dt22(5,2)*u(2,4,k)

       do j=3,my-2
          wk1(j)= dt22(1,j)*u(1,j-2,k)
          wk2(j)= dt22(1,j)*u(2,j-2,k)
          do i=2,5
             wk1(j) = wk1(j) + dt22(i,j)*u(1,i+j-3,k)
             wk2(j) = wk2(j) + dt22(i,j)*u(2,i+j-3,k)
          enddo
       enddo

       wk1(my-1)=dt22(1,my-1)*u(1,my-3,k)+dt22(2,my-1)*u(1,my-2,k)+ 
     $      dt22(3,my-1)*u(1,my-1,k)+dt22(4,my-1)*u(1,my  ,k)
       
       wk2(my-1)=dt22(1,my-1)*u(2,my-3,k)+dt22(2,my-1)*u(2,my-2,k)+ 
     $      dt22(3,my-1)*u(2,my-1,k)+dt22(4,my-1)*u(2,my  ,k)

       wk1(my)=  dt22(1,my  )*u(1,my-2,k)+dt22(2,my  )*u(1,my-1,k)+ 
     $       dt22(3,my  )*u(1,my,k)

        wk2(my)=  dt22(1,my  )*u(2,my-2,k)+dt22(2,my  )*u(2,my-1,k)+ 
     $       dt22(3,my  )*u(2,my,k)

! second derivative of u
! real part
        call banbks(prem3,my,wk1)
! imaginary part       
        call banbks(prem3,my,wk2)
!------------------------------------------------------------!
        coef=bet2(icx(k-1))+alp2(xpl-1) !kh^2
        !Compute Laplacian and save in "du"
          !real part 
          do j=1,my
           du(1,j,k)=wk1(j)-coef*u(1,j,k)
          enddo
         !imaginary part       
          do j=1,my
            du(2,j,k)=wk2(j)-coef*u(2,j,k)
          enddo
      enddo
 
      end subroutine laplacian

!-----------------------------------------------------------------------!

!!******************************************************************!
!!    FOUPHYSFOU2-tauij                           ! 
!!    ------------                                                  ! 
!!    From Inputs->go to Phys-> do operations->go back with outputs !     
!!                                                                  !
!!       VARIABLE       INPUT          OUTPUT                       !
!!       --------------------------------------                     !
!!       ten11   =       S11           tau11                        !
!!       ten22   =       S22           tau22                        !
!!       ten33   =       S33           tau33                        !
!!       ten12   =       S12           tau12                        !
!!       ten13   =       S13           tau13                        !
!!       ten23   =       S23           tau23                        !
!!      drho    =        T             Drho                         !
!!       tnext   =   T(n+1)-T(n)        -                           !
!!******************************************************************!
      subroutine tauij(ten11,ten22,ten33,ten12,ten13,ten23,
     .           drho,tnext,chwk,myid,rkstep)

      use MPI_GROUPS
      use point
      use tem
      use fis
      use timacc
      use wave
      use timers
      use wkhvect
      use cnan
   
      implicit none

      include "mpif.h"
      include "ctes3D"


c     ---------------------- Variables ----------------------
      
      complex*8 ::  ten11 (0:my-1,0:mgalz-1,pb:pe),
     .          ten22 (0:my-1,0:mgalz-1,pb:pe),
     .          ten33 (0:my-1,0:mgalz-1,pb:pe),
     .          ten12 (0:my-1,0:mgalz-1,pb:pe),
     .          ten13 (0:my-1,0:mgalz-1,pb:pe),
     .          ten23 (0:my-1,0:mgalz-1,pb:pe),
     .          tnext (0:my-1,0:mgalz-1,pb:pe),
     .          drho  (0:my-1,0:mgalz-1,pb:pe)
               
!      real*8 rf0u(0:my-1),rf0w(0:my-1),dum

      real(4) ::  chwk(*),uner(9),aa
      integer myid,rkstep,i,j,k,kk,ierr

! -------------------Programa-------------------------------!    
      
!  local transpose yz to zy, and adding zeros between positive
!  and negative wavenumbers.
      call localyz2zy(ten11,ten11,chwk) !S11
      call localyz2zy(ten22,ten22,chwk) !S22
      call localyz2zy(ten33,ten33,chwk) !S33
      call localyz2zy(ten12,ten12,chwk) !S12         
      call localyz2zy(ten13,ten13,chwk) !S13  
      call localyz2zy(ten23,ten23,chwk) !S23        
      call localyz2zy(drho ,drho ,chwk) !T      
      call localyz2zy(tnext,tnext,chwk) !T(n+1)-T(n)    
       
! inverse transform (fou --> fis @z),does the compacting aswell
      do i = pb,pe
         call fourz(ten11(0,0,i),1)    
         call fourz(ten22(0,0,i),1)    
         call fourz(ten33(0,0,i),1)    
         call fourz(ten12(0,0,i),1)  
         call fourz(ten13(0,0,i),1)  
         call fourz(ten23(0,0,i),1)  
         call fourz(drho (0,0,i),1) 
         call fourz(tnext(0,0,i),1) 
      enddo

!      nanerror = 0
!
!      call check_nan(ten11,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten22,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten33,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten12,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten13,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten23,my*mgalz*mmp*2,nanerror)
!      call check_nan(drho,my*mgalz*mmp*2,nanerror)
!      call check_nan(tnext,my*mgalz*mmp*2,nanerror)
!
!      if (nanerror.eq.1) then
!         write(*,*) 'NAN Found in ',myid,' before tauij changes'
!         stop
!      endif


!  change plane to line(1:mx,lb:le)
!  and distribute package of lines to processors
      call chpl2ln(ten11,ten11,chwk,myid)
      call chpl2ln(ten22,ten22,chwk,myid)
      call chpl2ln(ten33,ten33,chwk,myid)
      call chpl2ln(ten12,ten12,chwk,myid)
      call chpl2ln(ten13,ten13,chwk,myid)
      call chpl2ln(ten23,ten23,chwk,myid)
      call chpl2ln(drho,drho,chwk,myid)
      call chpl2ln(tnext,tnext,chwk,myid)


!==========================================================================!
! before this call ten11,ten22,ten33,ten12,ten13, ten23,drho and tnext
! are still lines in fourier space x (phys z).
      call tensor(ten11,ten22,ten33,ten12,ten13,ten23,
     .           drho,tnext,myid,rkstep)
!==========================================================================!
! they are still lines in fourier x after hvect

c  ---------- back to the yz planes and fft 
      call chln2pl(drho,drho,chwk,myid)   !drho 
      call chln2pl(ten11,ten11,chwk,myid) !tau11
      call chln2pl(ten22,ten22,chwk,myid) !tau22
      call chln2pl(ten33,ten33,chwk,myid) !tau33
      call chln2pl(ten12,ten12,chwk,myid) !tau12
      call chln2pl(ten13,ten13,chwk,myid) !tau13
      call chln2pl(ten23,ten23,chwk,myid) !tau23
!
!
!!debug CHECK NAN
!      call check_nan(ten11,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten22,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten33,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten12,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten13,my*mgalz*mmp*2,nanerror)
!      call check_nan(ten23,my*mgalz*mmp*2,nanerror)
!      call check_nan(drho,my*mgalz*mmp*2,nanerror)
!
!      if (nanerror.eq.1) then
!         write(*,*) 'NAN Found in ',myid,' after tauij changes'
!         stop
!      endif


!  convert to fourier Z now.
      do i=pb,pe
         call fourz(drho (0,0,i),-1)  !drho
         call fourz(ten11(0,0,i),-1)  !tau11
         call fourz(ten22(0,0,i),-1)  !tau22
         call fourz(ten33(0,0,i),-1)  !tau33
         call fourz(ten12(0,0,i),-1)  !tau12
         call fourz(ten13(0,0,i),-1)  !tau13
         call fourz(ten23(0,0,i),-1)  !tau23
      enddo

!transpose back
      call localzy2yz(drho , drho,chwk) 
      call localzy2yz(ten11,ten11,chwk) 
      call localzy2yz(ten22,ten22,chwk) 
      call localzy2yz(ten33,ten33,chwk) 
      call localzy2yz(ten12,ten12,chwk)   
      call localzy2yz(ten13,ten13,chwk)   
      call localzy2yz(ten23,ten23,chwk)   
           

      end



!!*********************************************************************!
!!                                                                     !
!!         transform to physical and do operations                     !
!!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!!                                                                     !
!!                                                                     !
!!*********************************************************************!
      subroutine tensor(ten11,ten22,ten33,ten12,ten13,ten23,
     .                 drho,tnext,myid,rkstep)

      use point
      use tem
      use fis
      use timacc 
      use wkhvect
      use cnan
      use MPI_GROUPS

      implicit none

      include "mpif.h"
      include "ctes3D"

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
      real(8)  cflx,cfly,cflz,hxalp,hzbet,hyy,cfl0,reigmx1,aa
      real(4) :: tnext(0:mx-1,lb:le),
     .        ten11(0:mx-1,lb:le),
     .        ten22(0:mx-1,lb:le),ten33(0:mx-1,lb:le),
     .        ten12(0:mx-1,lb:le),ten13(0:mx-1,lb:le),
     .        ten23(0:mx-1,lb:le),
     .        drho (0:mx-1,lb:le)
      real(8)  sigma,const1
      real(8) duma,dumb,dumc,duma1,dumb1,dumc1
      
!      ----------------------- Programa ------------------------      
      !sigma=0.7D0
      sigma=0D0
      const1=2/3

c        Move everything to pPP, line by line 
      duma= 0d0
      dumb= 0d0
      dumc= 0d0
      do 50 j = lb,le
         do i=0,mx-1
            up1wk(i) = drho(i,j)!T(n)
         enddo
         do i=0,mx-1           
            up2wk(i) = tnext(i,j) !T(n+1)-T(n)
         enddo
         do i=0,mx-1
            ten11wk(i) = ten11(i,j) !S11
         enddo
         do i=0,mx-1
            ten22wk(i) = ten22(i,j) !S22
         enddo
         do i=0,mx-1
            ten33wk(i) = ten33(i,j) !S33
         enddo
         do i=0,mx-1
            ten12wk(i) = ten12(i,j) !S12
         enddo
         do i=0,mx-1
            ten13wk(i) = ten13(i,j) !S13
         enddo
         do i=0,mx-1
            ten23wk(i) = ten23(i,j) !S23
         enddo
        
         if (myid.eq.0) then
            duma =duma+MPI_WTIME()
         endif

! convert lines in fourier to lines in physical for x
! upiwk8 means real 8 format
         call fourx(up1wk,up1wk8,1)    !T
         call fourx(up2wk,up2wk8,1)    !T(n+1)-T(n)
         call fourx(ten11wk,ten11wk8,1) ! S11
         call fourx(ten22wk,ten22wk8,1) ! S22
         call fourx(ten33wk,ten33wk8,1) ! S33
         call fourx(ten12wk,ten12wk8,1) ! S12
         call fourx(ten13wk,ten13wk8,1) ! S13
         call fourx(ten23wk,ten23wk8,1) ! S23
! ----------------------------------------------------------------------!
!                                                                       !
!         computing in physical domain ( line by line)                  !
!                                                                       !
! ----------------------------------------------------------------------!
         !Calculate drho=1/T(n+1)-1/T
         do i=0,mgalx-1
             up1wk(i)  = -up2wk8(i)/(up1wk8(i)*(up1wk8(i)+up2wk8(i)))
             !drho ready
         enddo
         do i=0,mgalx-1
             tmp1wk8(i) =const1*(ten11wk8(i)+ten22wk8(i)+ten33wk8(i))
             !2/3Skk
         enddo
         !S11,S22,S33 like: tau11=mu*(2S11-2/3Skk)
         !Sij like        : tauij=2*mu*Sij
         do i=0,mgalx-1
!Full formulation:
            ten11wk(i)=up1wk8(i)**sigma*(2*ten11wk8(i)
     .                                   -tmp1wk8(i))

            ten22wk(i)=up1wk8(i)**sigma*(2*ten22wk8(i)
     .                                   -tmp1wk8(i))

            ten33wk(i)=up1wk8(i)**sigma*(2*ten33wk8(i)
     .                                   -tmp1wk8(i))

            ten12wk(i)=2*up1wk8(i)**sigma*ten12wk8(i)

            ten13wk(i)=2*up1wk8(i)**sigma*ten13wk8(i)

            ten23wk(i)=2*up1wk8(i)**sigma*ten23wk8(i)
        enddo
!!Incompressible formulation:
!            ten11wk(i)=2*ten11wk8(i)-up2wk8(i)
!            ten22wk(i)=2*ten22wk8(i)-up2wk8(i)
!            ten33wk(i)=2*ten33wk8(i)-up2wk8(i)
!            ten12wk(i)=2*ten12wk8(i)
!            ten13wk(i)=2*ten13wk8(i)
!            ten23wk(i)=2*ten23wk8(i)
 

!========================================================!
!IMPORTANT POINT: all quantities are now PHYSICAL
!========================================================!
c--------------------- back to F-P-P  

         call fourx(ten11wk,up1wk8,-1)  ! tau11
         call fourx(ten22wk,up1wk8,-1)  ! tau22
         call fourx(ten33wk,up1wk8,-1)  ! tau33
         call fourx(ten12wk,up1wk8,-1)  ! tau12
         call fourx(ten13wk,up1wk8,-1)  ! tau13
         call fourx(ten23wk,up1wk8,-1)  ! tau23
         call fourx(up1wk  ,up1wk8,-1)  ! drho

c --------    back to lines. We throw away the high modes (Dealiasing)

         if (myid.eq.0) then
            dumc = dumc-MPI_WTIME()
         endif

         do i=0,mx-1
            drho(i,j) = up1wk(i) !drho
         enddo
         do i=0,mx-1            
            ten11(i,j) = ten11wk(i)
         enddo
         do i=0,mx-1            
            ten22(i,j) = ten22wk(i)
         enddo
         do i=0,mx-1            
            ten33(i,j) = ten33wk(i)
         enddo
         do i=0,mx-1            
            ten12(i,j) = ten12wk(i)
         enddo
         do i=0,mx-1            
            ten13(i,j) = ten13wk(i)
         enddo
         do i=0,mx-1            
            ten23(i,j) = ten23wk(i)
         enddo

         if (myid.eq.0) then
            dumc = dumc+MPI_WTIME()
         endif

 50    continue

      call check_nan(drho,mml*mx,nanerror)
      
      if (nanerror.eq.1) then
         write(*,*) 'NAN found in',myid,'on drho lines'
         stop
      endif

      end
     

