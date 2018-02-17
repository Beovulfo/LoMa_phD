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
     .                  spwk,
     .                  sp,myid)
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
      use statis
      use MPI_GROUPS
      use spectra
      use rkcoef
     

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
!spectra
      real*4 sp  (0:nz1,1:  nspec+1,12,pb:pe),
     .       spwk(0:nz1,1:  nspec+1,12,pb:pe)
      
      real(8) :: rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire
      
      real(8)  iter_time,write_time
      complex*8 temp1,temp2
      real(4) ::  temporal(3),comun(3)
      
      real(4)  enerdis,H00u,H00w,dumu,dumw,cteu,ctew
      real(4)  massu,massw,massu1,massu2,massw1,massw2,dum
      real(4)  uner(9) 

      character(100) fname,fdebug
      complex*8,allocatable :: tempa(:,:),tempb(:,:)

!debugging variables
      character(3) ext1,ext2 
      integer :: ideb
      
c    ---------------------- Program ------------------------------------

      allocate(tempa(0:mz1,pb:pe),tempb(0:mz1,pb:pe))

c    Defining combinations of wavenumbers useful for calculations.      
      if (myid.eq.0) then
         !write RK coefs:
          write(*,*) "RK3 coefficients (alpha,beta,gama,xi)"
          write(*,*) (alpha(j),beta(j),gama(j),xi(j),j=1,3)

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
!     Timers and other stuff initialized
      commtimer =0.0D0
      transtimer=0.0D0
      totaltimer=0.0D0
      writetimer=0.0D0

      ire      = 1.0d0/re
      ihist    = 1   !write first step
      irun     = 0   ! first time step is special in tim3rkp
      icfl     = 1   ! first time step always needs a step size
      nanerror = 0   ! Check NAN in the 50 firsts steps
      nanproblem = 0 ! initializes NAN controller

      if(myid.eq.0) then
         fname=filstt(1:index(filstt,' ')-1)//'.cf'
         write (*,*) fname
         open(29,file=fname,status='unknown')
      endif

!!!!    For testing s=1.0 with temperature constant
!!!       !*********
!         do i=pb,pe
!            do k=0,mz1
!               do j=0,my-1
!                  scal(j,k,i)=0.0  !Temperature set to 1
!!                  !WARNING!!!!!!!!!!!!!!!!!!!!!TAKE THIS OUT!
!               enddo
!            enddo
!         enddo
!!!         !*********
!!!         !WARNING!!!!!!!!!!!
!         if (myid.eq.0) then
!         do j=0,my1
!            scal(j,0,1) = 1.0
!         enddo
!         endif
!!!         !*********
!!!




!Initialization: 
       nacum = 0
       nacumsp = 0
      if (myid.eq.0) then
      !Save densities at boundaries
         rhobot =1d0/real(scal(0,0,1)) !mode 00,bot
         rhotop = 1d0/real(scal(my-1,0,1)) !mode 00,top
         rhocrit = min(rhobot,rhotop) !minimum density,for Dt
          do j=0,my1
             rf0u(j)=0d0
             rf0w(j)=0d0
             rf0v(j)=0d0
             drho00(j)=0d0
          enddo
      endif

      call MPI_BCAST(rhobot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhotop,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhocrit,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      !Clean all stats
      call resetstats(sp,myid) !reset sp and stats    



c ========================================================================
c                 THIS IS THE TIME LOOP
c ========================================================================
      do 30 istep=1,nstep
         if (myid.eq.0) then
           totaltimer = totaltimer-MPI_WTIME()
           iter_time=-MPI_WTIME()
         endif
         !-----------------------------------------------------------! 
         !ISTEP CHECKINGS/FLAGS UPDATINGS
         !-----------------------------------------------------------! 
         !stats accumulations
         !accumulating "ntimes" steps continously just before 
         !writing stats time
         if (mod(istep-1,nhist).eq.(nhist-ntimes).and.nstart.ne.0) then
                istati = 1
         endif
         if (mod(istep,nhist).eq.0)  ihist = 1 
         if (mod(istep-1,ncfl).eq.0) icfl = 1 !update Deltat

         ! ------------------- write image to a file ---------------------------!

         IF (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then
            if (myid.eq.0) then
               write_time = -MPI_WTIME()
               writetimer = writetimer - MPI_WTIME()
            endif
            do j=0,my1
               u00(j) = u00wk(j) 
               w00(j) = w00wk(j)
               v00(j) = v00wk(j)
            enddo
            !time need to be broadcasted to save procs
            !call MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
            call escru(vor,phi,psi,scal,u00,w00,v00,myid)

            if (myid.eq.0) then
               writetimer = writetimer + MPI_WTIME()
            endif
         ENDIF
! ------------------- finished writing image --------------------------!

      if (myid.gt.numerop-1) goto 30    ! only for save procs


      if (irun.eq.0) then ! this is done only for the first step
           irun = 1
         if (myid.eq.0) then
            do j=0,my1
              u00wk(j)=u00(j)
              v00wk(j)=v00(j)
              w00wk(j)=w00(j)
            enddo
         endif
      endif ! end special first step

               !  Runge-Kutta third order  !
!=======================================================================!                
!----------- RUNGE-KUTTA SUBSTEPS START---------------------------------!
!=======================================================================!                

      do 10 rkstep=1,3
!at rkstep=1 dtxi=0, for rkstep>1 we have the proper Dt for the istep
!and it shouldn't change through the istep.
         dtxi=Deltat*xi(rkstep)
c-------------save vor and psi in nxwk/nywk to work---------------
!Clean useless modes 00 of vor and phi
         if (myid.eq.0) then
            do j=0,my-1
             vor(j,0,1)=0.0   !Omega_y
             phi(j,0,1)=0.0   !PHI
             psi(j,0,1)=0.0   !PSI
            enddo
         endif

         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  nxwk(j,k,i)=vor(j,k,i)   !Omega_y
                  nywk(j,k,i)=psi(j,k,i)   !PSI
               enddo
            enddo
         enddo

!Setting BC'S
         !BC on Omegay
         do i=pb,pe
            do k=0,mz1
                !These BCs are not well imposed for CFD scheme
                !d(Omega_y)/dy @ Bound's = 0
               !testing adding this
                !nxwk(1  ,k,i)   = nxwk(2,k,i)
                !nxwk(my1-1,k,i) = nxwk(my1-2,k,i)
                !one point
                nxwk(0  ,k,i) = nxwk(1,k,i)   
                nxwk(my1,k,i) = nxwk(my1-1,k,i)   
            enddo
         enddo
         !BC on phi
         do i=pb,pe
            do k=0,mz1
                 !phi = 0 @ Bound's
                 phi(0  ,k,i) = cmplx(0.0,0.0)
                 phi(my1,k,i) = cmplx(0.0,0.0)
            enddo
         enddo
         !BC on T
         if (myid.eq.0) then
            do i=pb,pe
               do k=0,mz1
                  if (pb.eq.1.and.k.eq.0) then
                  !00 mode
                     scal(0,0,1) = cmplx(1.0/rhobot,0.0) 
                     scal(my1,0,1) = cmplx(1.0/rhotop,0.0)
                  else 
                    !phi = 0 @ Bound's
                     scal(0  ,k,i) = cmplx(0.0,0.0)
                     scal(my1,k,i) = cmplx(0.0,0.0)
                  endif
               enddo
            enddo
         else 
           do i=pb,pe
               do k=0,mz1
                 !phi = 0 @ Bound's
                 scal(0  ,k,i) = cmplx(0.0,0.0)
                 scal(my1,k,i) = cmplx(0.0,0.0)
               enddo
            enddo
         endif
 
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
!        Save scal into drho
         do i=pb,pe
           do k=0,mz1
             do j=0,my-1
                drho(j,k,i)=scal(j,k,i) 
             enddo
           enddo
         enddo

!        Now build u_i+RHS(i-1)*Dt*eps_i
!        -------------------------------
         do i=pb,pe
           do k=0,mz1
             do j=0,my1
               vor(j,k,i) = vor(j,k,i) +ten12(j,k,i)*dtxi!Omega_y
               phi(j,k,i) = phi(j,k,i) +ten13(j,k,i)*dtxi!PSI
               scal(j,k,i)= scal(j,k,i)+ rhst(j,k,i)*dtxi!scal 
               !save buffer with RHS of Scal
               chwk2(j,k,i)= rhst(j,k,i)*dtxi
             enddo
           enddo
         enddo
!        Save scal into drho, once it has been used.
         do i=pb,pe
           do k=0,mz1
             do j=0,my1
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

!!---------- Adding modes 00 to rhou and rhow------------------------!
!---------- Only master computes 00  modes
        if (myid.eq.0) then
           do j=0,my1  
              mxwk(j,0,1) = u00wk(j)
              mywk(j,0,1) = v00wk(j) !rhov00 
              mzwk(j,0,1) = w00wk(j)
           enddo

!       !Now build in u00/w00 the variable u00**/w00**
!        Set BC rho u remains unchanged at boundaries
           do j=1,my1-1
              u00wk(j)=u00wk(j)+dtxi*rf0u(j)
              w00wk(j)=w00wk(j)+dtxi*rf0w(j)
           enddo
        endif


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
!         do i=pb,pe
!            call laplacian(rhst(0,0,i),nywk(0,0,i),mz,i)
!        enddo 

         do i=pb,pe
            call deryyr2(rhst(0,0,i),nywk(0,0,i),mz)
         enddo
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               do j=0,my1
                 nywk(j,k,i) = nywk(j,k,i) - rK*rhst(j,k,i)
               enddo
            enddo
         enddo


!---------!Calculate dT/dz,dT/dx------------------------------------!
         do i=pb,pe
            do k=0,mz1
               do j=0,my1
                  !dT/dz=i kz T
                  ten23(j,k,i)=xbet(k)*rhst(j,k,i)
                  !dT/dx=i kx T
                  ten12(j,k,i)=xalp(i-1)*rhst(j,k,i)
               enddo
            enddo
        enddo


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
!======================================================!
      ! STATIS
      !if istati==1 means we need to stack stats
      if (istati.eq.1.and.rkstep.eq.1) then
          !Accumulating stats while istati==1
          !Computing stats for RHOUi means and T
          do i=pb,pe
             call addstats1(mxwk(0,0,i),mywk(0,0,i),mzwk(0,0,i),
     .                     rhst(0,0,i),sp(0,1,1,i),i)
          enddo
      endif
!======================================================!

!------------(5) FOU-PHYS(CALCULATIONS)-FOU--------------------------!
!Here we compute the new Dt                                          !
        call foup2(mxwk,mywk,mzwk,nxwk,nywk,nzwk,
     .          rhst,ten12,ten13,ten23,chwk,sp,myid,rkstep)
!
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
!   chwk2: DT = T(i+1) - T(i)
!   ten12: rhouv
!   ten13: rhouw
!   ten23: rhovw
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!   drho : T                                                        
!*******************************************************************!


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
!   chwk2: DT = T(i+1) - T(i)
!   ten12: -
!   ten13: -
!   ten23: -
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!
!======================================================!
      ! STATIS
      !adding velocity and vorticity stats
      !if istati==1 means we need to stack stats
      if (istati.eq.1.and.rkstep.eq.1) then
!=================================================================
!                   Calculating VORITICITIES
!                 Only for STATISTICS proposes
!==================================================================
! Calculate derivatives on y
         do i=pb,pe
            call deryr2(mxwk(0,0,i),ten12(0,0,i),mz)!du/dy
            call deryr2(mzwk(0,0,i),ten23(0,0,i),mz)!dw/dy
         enddo
!We need to build vor_y=ikzu - ikxw
          do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  !vorx
                  ten23(j,k,i)= ten23(j,k,i)
     .                        - xbet(k)*mywk(j,k,i)
                  !vory
                  ten13(j,k,i)= xbet(k)*mxwk(j,k,i)
     .                        - xalp(i-1)*mzwk(j,k,i)
                  !vorz
                  ten12(j,k,i)= xalp(i-1)*mywk(j,k,i)
     .                        - ten12(j,k,i)
               enddo
            enddo
          enddo

!NOW WE HAVE VORTICITES!!
!!*******************************************************************!
!After State(6c):
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
!   chwk2: DT = T(i+1) - T(i)
!   ten12: vorx
!   ten13: vory
!   ten23: vorz
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!


          !Accumulating stats while istati==1
          nacum = nacum +1
          nacumsp = nacumsp +1
          !STATS FOR U,V,W, VORx VORy VORz
          !reset uner 
          do kk = 1,9
              uner(kk) = 0.
          enddo
          !Calling for accumulation and dissipation
          do i=pb,pe
             call addstats2(mxwk(0,0,i),mywk(0,0,i),mzwk(0,0,i),
     .       ten23(0,0,i),ten13(0,0,i),ten12(0,0,i),chwk(0,0,i),
     .       uner,sp(0,1,1,i),i)
          enddo
          !If time to write hist REDUCE uner
          if (ihist.eq.1) then
               call MPI_ALLREDUCE(uner,energy,9,MPI_REAL,MPI_SUM,
     .                            MPI_COMM_CALC,ierr)
               do i=1,9
                  energy(i) = sqrt(abs(energy(i)))
               enddo
          endif
      endif


!-------------(7) Calculate Sij-----------------------------------!       
! Calculate derivatives on y
         do i=pb,pe
            call deryr2(mxwk(0,0,i),ten12(0,0,i),mz)!du/dy
            call deryr2(mzwk(0,0,i),ten23(0,0,i),mz)!dw/dy
         enddo
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

!!*******************************************************************!
!After State(7):
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
!   chwk2: DT = T(i+1) - T(i)
!   ten12: S12
!   ten13: S13
!   ten23: S23
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!
        if (sigma.ne.0.0) then
!CALL TAUIJ2 in order to compute DRHO
        call tauij2(drho,chwk2,chwk,
     .       mxwk,mywk,mzwk,ten12,ten13,ten23,
     .       myid,rkstep)
        else
          do i=pb,pe
             do k=0,mz1
                do j=0,my1
                    !TauIJ = 2.0*Sij
                    mxwk(j,k,i) = 2.0* mxwk(j,k,i)
                    mywk(j,k,i) = 2.0* mywk(j,k,i)
                    mzwk(j,k,i) = 2.0* mzwk(j,k,i)
                    ten12(j,k,i) = 2.0* ten12(j,k,i)
                    ten13(j,k,i) = 2.0* ten13(j,k,i)
                    ten23(j,k,i) = 2.0* ten23(j,k,i)
                    drho(j,k,i) = 0d0
               enddo
             enddo
          enddo

        endif
 
!*******************************************************************!
!! outputs: 
!! mxwk = tau11 
!! mywk = tau22
!! mzwk = tau33
!! ten12= tau12
!! ten13= tau13
!! ten23= tau23
!! drho = drho
!
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
!!*******************************************************************!
!After State(8):
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
!   chwk2: DT = T(i+1) - T(i)
!   ten12: -
!   ten13: -
!   ten23: -
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!

!!----------(9) Construct the RHS of Omega_y and phi--------------------!
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

!-------------(9) Evolving Omega_y and phi-------------------------!       
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
!*******************************************************************!
!After State(9):
!   vor  : Omega_y(i+1)
!   phi  : phi(i+1)
!   psi  : psi(i+1)
!   scal : T(i+1) 
!   nxwk : Nx
!   nzwk : Nz
!   nywk : Ny
!   mxwk : -
!   mywk : -
!   mzwk : -
!   chwk:  -
!   chwk2: DT = T(i+1) - T(i)
!   ten12: RHS(omegay)
!   ten13: RHS(phi)
!   ten23: -
!   drho : drho                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!!*******************************************************************!
! ----------------------------------------------------------------------!
!             this is the block for the 00 mode
!                   Only master does it
! ----------------------------------------------------------------------!
        if (myid.eq.0) then
! evolve modes 00 por rhou, rhow
           do j=1,my1-1
             u00wk(j)=u00wk(j)+dtgamma*rf0u(j)
             w00wk(j)=w00wk(j)+dtgamma*rf0w(j)
           enddo
!compute momentum 
           massu=0.
           massw=0.
           do j=0,my-1
              massu = massu+u00wk(j)*hy(j+1)
              massw = massw+w00wk(j)*hy(j+1)
           enddo
        endif   ! 00 modes
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

!------------(8) Obtain psi(n+1)-----------------------------------!
!1) First solve 00 mode of psi
!   Only for myid.eq.0
         if (myid.eq.0) then
!        1.1) Integrate drho00
              !F(y) = F(-Ly) + int(f(y),y,[-Ly,y])
              do j=0,my1
              !    drho00(j) = real(drho(j,0,1))/(dtgamma+dtxi)
                  drho00(j) =-real(drho(j,0,1))
!     .                   -alpha(rkstep)/beta(rkstep)*drho00(j)
              enddo 
              !CUMTRAPZ:
              !rf0v(0) = 0d0
              !do j=0,my1-1
              !   rf0v(j+1) = rf0v(j)+0.5d0*(y(j+2)-y(j+1))*
     .        !        (drho00(j+1)+drho00(j))
              !enddo
              !compact
              call cumtrapz(drho00,rf0v,0d0)
              !call inty8(drho00,rf0v) !first element is zero
              trapz = rf0v(my1) !last element is the area of f!
              !now rf0v keeps the primitive from -Ly to y of drho
              do j=0,my1
                  !A) condition v(-Ly) = -v(Ly)
                  !--------------------------------
                  ! v00wk(j) = -rf0v(j) + rhobot*trapz/(rhotop+rhobot)

                  !B) condition rhov(-Ly) = -rhov(Ly)
                  !..................................
                  !v00wk(j) = -rf0v(j) + trapz/2.0d0
                  !C) Pressure Condition pbot=ptop
                  !             rhov(Ly) = trapz/(1+s^(-0.5))
                  !             rhov(-Ly) = trapz*(1/(1+s^(-0.5))-1)
                  v00wk(j) = rf0v(j) - trapz*
     .              (1d0/(1d0+(rhotop/rhobot)**0.5))
     .                   -alpha(rkstep)/beta(rkstep)*v00wk(j)
              enddo
!rf0v  keeps the integral of drho00
         endif
!----------------------------------------------------------------------
!2) Solve other modes...
!         !first we need to obtain LAP of PSI
!         do i=pb,pe
!            call deryyr2(psi(0,0,i),nywk(0,0,i),mz)
!         enddo
!         do i=pb,pe
!            do k=0,mz1
!               k1 = icx(k)
!               rK = bet2(k1)+alp2(i-1)
!               do j=0,my1
!                 nywk(j,k,i) = nywk(j,k,i) - rK*psi(j,k,i)
!               enddo
!            enddo
! construct -drho/dtbeta in order to solve laplacian
         do i=pb,pe
            do k=0,mz1
               do j=0,my1
                   nywk(j,k,i) = -drho(j,k,i)
               enddo
            enddo
         enddo
!solve laplacian of alpha/beta*psi(n)+psi(n+1)
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               !call Lapvdv(nywk(0,k,i),drho(0,k,i),ten23(0,k,i),rK)
               call Lapvdvhom(nywk(0,k,i),drho(0,k,i),ten23(0,k,i),rK)
            enddo
         enddo

         do i=pb,pe
            do k=0,mz1
               do j=0,my1
                   psi(j,k,i) =drho(j,k,i)-
     .                alpha(rkstep)/beta(rkstep)*psi(j,k,i)
               enddo
            enddo
         enddo


!!   2.1)Update psi(n+1)
!          do i=pb,pe
!             do k=0,mz1
!               do j=0,my1
!                  !psi(j,k,i)=-drho(j,k,i)/(dtgamma+dtxi)!Euler Implicit
!                  nywk(j,k,i)=-drho(j,k,i)/dtbeta-
!     .                        alpha(rkstep)/beta(rkstep)*nywk(j,k,i)
!                  !nywk(j,k,i)=drho(j,k,i)/dtbeta 
!               enddo
!             enddo
!          enddo 
!
!!New variable:
!!   2.2) solve Lapvdv of PSI
!         do i=pb,pe
!            do k=0,mz1
!               k1 = icx(k)
!               rK = bet2(k1)+alp2(i-1)
!               call Lapvdvhom(nywk(0,k,i),psi(0,k,i),ten23(0,k,i),rK)
!               !call Lapvdv(psi(0,k,i),psi(0,k,i),ten23(0,k,i),rK)
!               !call solvepsi(nywk(0,k,i),psi(0,k,i),ten23(0,k,i),rK)
!            enddo
!         enddo
!

! Now:
! ten23 = not useful...
! psi   = psi(n+1)
! v00   = rhov00(n+1)
! v00wk = rhov00(n+1)
!***************************************************************************!


 10    continue  !!! End of R-K sub-step

       if (ihist.eq.1.and.istati.eq.1) then
            if (myid.eq.0) then 
               writetimer = writetimer - MPI_WTIME()
            endif
            call writestats(sp,spwk,myid,istep) !added spectrum
            !writestats also makes ALLREDUCE for ep->wkstats11
            !CALCULATE ENERGY DISSIPATION
            if (myid.eq.0) then
                 writetimer = writetimer + MPI_WTIME()
                 enerdis =0.0
                 !ep reduced is on wkstats(j,11)
                 do j=1,my
                      enerdis = enerdis + wkstats(j,11)/nacum*hy(j)
                 enddo
            endif
            call resetstats(sp,myid) !reset sp and stats    
            istati = 0 !not more accumulation
       endif
!WRITE CF AND OUTPUT FILE------------------------
         if (myid.eq.0) then !ONLY MASTER
             if (ihist.eq.1) then !WHEN WRITING STATS/CF/OUTPUT
                 if (istep.eq.1) then
                     ctm2=0.
                     ttm =0.
                 endif
324              format(17(d22.14))
                 write(29,324)  time,dm,rhocrit,enerdis, energy,
     .            massu,massw,
     .           MPI_WTIME()+iter_time, commtimer-ctm2
                 call flush(29) !for XLF
                 !call flush(29)

325              format(i5,13(d14.6))
                 write(*,325) istep,time,Deltat,dm,rhocrit,massu,massw,
     .                   ire,H00u,H00w,MPI_WTIME()+iter_time,
     .                    transtimer-ttm,commtimer-ctm2
!================================================================!
        !clean stats now that ENER and EP have been used!
            !clean stats after writing history files (out of rkstep
            !loop)
            endif
!==================================================================!
            ctm2 = commtimer
            ttm  = transtimer 
         endif

c          	! time:
         time=time+Deltat

         if(icfl .eq.1)   icfl   = 0
         if(ihist.eq.1)  ihist  = 0

         if (myid.eq.0) then
           totaltimer=totaltimer+MPI_WTIME()
           !print *,MPI_WTIME()+iter_time
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
        !write(11,*) temporal(1),comun(1),
     .  !            commtimer/totaltimer,totaltimer
        print *,"Total time: ",totaltimer
        
        print *,"Trans. time: ",transtimer
        print *,"Comm. time: ",commtimer
        print *,"Writ. time: ",writetimer
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
      subroutine foup2(up1,up2,up3,ten11,ten22,ten33,
     .           rhst,ten12,ten13,ten23,chwk,sp,myid,rkstep)
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
      use statis
      use spectra

       
      implicit none

      include "mpif.h"
      include "ctes3D"


c     ---------------------- Variables ----------------------
      
      complex*8, dimension(0:my1,0:mgalz-1,pb:pe):: 
     .          up1,up2,up3,ten11,ten22,ten33,ten12,
     .          ten13,ten23,rhst,chwk

      real(4) uner(9),aa
      integer myid,rkstep,i,j,k,kk,ierr
c-------spectra----------
      real*4 sp(0:nz1,1:nspec+1,12,pb:pe)

      
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
      call phys2(up1,up2,up3,ten11,ten22,ten33,
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
                  dm=dm+trp(j+1)*((real(up1(j+1,0,1)))**2.0+
     .                          (real(up1(j-1,0,1)))**2.0)
              enddo
              dm=dm+trp(1)*((real(up1(1,0,1)))**2.0+
     .        (real(up1(0,0,1)))**2.0)+
     .        trp(my)*((real(up1(my1,0,1)))**2.0+
     .                  (real(up1(my-2,0,1)))**2.0) 

              dm =1d0/4.0d0*(y(my)-y(1))-
     .           dm/(real(up1(my1,0,1))-real(up1(0,0,1)))**2.0
             endif
          endif
      end subroutine
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!         transform to physical and do operations                     !
!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!                                                                     !
!                                                                     !
!*********************************************************************!
      subroutine phys2(up1,up2,up3,ten11,ten22,ten33,
     .           ten12,ten13,ten23,rhst,myid,rkstep)
      use point
      use tem
      use fis
      use timacc
      use diag 
      use wkhvect
      use MPI_GROUPS
      use statis

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

      real(8) pi
      !aaf adding Peclet
      real(8) peclet
!debugging
      real*8 duma,dumb,dumc,duma1,dumb1,dumc1

c     ------------------------- Indices -----------------------

      integer vec(2,3),vecrc(2,3),vecy(2)         ! i,j
      integer resulm(2,3)

!      ----------------------- Programa ------------------------      
!initialization
      !aaf Peclet
      peclet=0.7d0*re  
      !sigma=0.7D0
      cflx = 0d0
      cfly = 0d0
      cflz = 0d0
      cfl0 = 0d0
      pi = 4.0d0*atan(1.0d0)
!calculate geometric staff
      !1/Dx and 1/Dz
      hxalp=alp*(mgalx-1)/pi !2/Dx
      hzbet=bet*(mgalz-1)/pi !2/Dz
      hyy=hy(1)
      do j=2,my
         hyy = min(hyy,hy(j))
      enddo
      hyy = 1.0d0/hyy
! computing viscous Dt:
! computing viscous Dt inverse:
! Dt_visc=1/(2 nu)*Dy_min^2
! nu = mu/rho = 1/(Re*rho00)
      hyy=2.0d0/(peclet*rhocrit**(sigma+1d0))*
     .   (hyy**2  +(hxalp/2.0d0)**2 + (hzbet/2.0d0)**2)
!      hyy=2.0d0/(Re*rhocrit)*
!     .   (max(hyy,max(hxalp/2.0d0,hzbet/2.0d0)))**2
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
!========================================================!
!IMPORTANT POINT: all quantities are now PHYSICAL
!========================================================!
c    ------------------------------------------------
c        compute statistics in PPP
c    ----------- triple products 
! each sub-block of lines are in the same y. So that
! the total number of lines in each processor are
! sections in y. (j is the line index).
         jj = (j-1)/mgalz +1
! obtain the equivalent y index (height of block lines).
         if (jj.lt.0) jj=1
         if (istati.eq.1 .and. rkstep.eq.1) then
            !obtain cross products rho u_i u_j
            do i=0,mgalx-1
               aa = up1wk8(i)*rhstwk8(i) !u
               ruu(jj) = ruu(jj) +aa*up1wk8(i) !ru*u
               ruv(jj) = ruv(jj) +aa*up2wk8(i) !rv*u
               ruw(jj) = ruw(jj) +aa*up3wk8(i) !rw*u
            enddo
            do i=0,mgalx-1
               aa = up2wk8(i)*rhstwk8(i) !v
               rvv(jj) = rvv(jj) +aa*up2wk8(i) !rv*v
               rvw(jj) = rvw(jj) +aa*up3wk8(i) !rw*v
            enddo
            do i=0,mgalx-1
               aa = up3wk8(i)*rhstwk8(i) !w
               rww(jj) = rww(jj) +aa*up3wk8(i) !rw*w
            enddo
            do i=0,mgalx-1
               rhom(jj) = rhom(jj) +1.0d0/rhstwk8(i) !rho=1/T.
               mum(jj) = mum(jj) +rhstwk8(i)**(sigma) !mu/mu0=T^sigma.
               num(jj) = num(jj) +rhstwk8(i)**(sigma+1.0) !nu/nu0=T^sigma+1.
            enddo
         endif



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
     .      (1.0d0/peclet)*(rhstwk8(i))**(sigma+1d0)*ten22wk8(i)
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
! each sub-block of lines are in the same y. So that
! the total number of lines in each processor are
! sections in y. (j is the line index).
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
!hxalp, hzbet are 2*(1/Dx,1/Dz)
         cflx = cflx*hxalp*0.5
         cflz = cflz*hzbet*0.5
!using a 3 factor in order to account for 3D mesh
!ex:Dt<CFL/(u/Dx+v/Dy+w/Dz) (it's a safety factor)
         cfl0 = max(cflx,max(cfly,cflz))

!!      !debug
!!      !-----------------------------------------------
!!         write(*,*) "myid = ",myid, "cfl0 =",cfl0,
!!     .      "cflx=",cflx,"cflz=",cflz,"hyy=",hyy,
!!     .      "hxalp=",hxalp,"hzbet=",hzbet,"cfly=",cfly
!!      !-----------------------------------------------


         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL8,MPI_MAX,
     .                     MPI_COMM_CALC,ierr)

!now take max of (Dt_conv^(-1), Dt_visc^(-1))
         reigmx1 = max(reigmx1,hyy)

!check viscous time vs diffusive time
!         if (hyy.eq.reigmx1.and.myid.eq.0) then
!            write(*,*) '------------Viscous Dt used----------' 
!            write(*,*) 'Dt_visc = ',CFL/hyy,'myid= ', myid
!         endif
!compute Deltat and dtr
         Deltat=CFL/reigmx1
         dtr=Re/Deltat


         if (reigmx1.lt.1e-1) then
            write(*,*) 'UGGGHH',myid,ihist,lb,le,pb,pe
         endif

!!Debugging
!      !debug
!      !-----------------------------------------------
!         if(myid.ne.0) then
!            call MPI_SEND(cflx,1,MPI_REAL8,0,
!     &                 myid,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(cfly,1,MPI_REAL8,0,
!     &                    myid,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(cflz,1,MPI_REAL8,0,
!     &                 myid,MPI_COMM_WORLD,ierr)
!!            call MPI_SEND(vec,6,MPI_INTEGER,0,
!!     .                   myid,MPI_COMM_WORLD,ierr)
!         else
!            duma1 = cflx
!            dumb1 = cfly
!            dumc1 = cflz
!!            resulm(1,1) = vec(1,1)
!!            resulm(2,1) = vec(2,1)
!!            resulm(1,2) = vec(1,2)
!!            resulm(2,2) = vec(2,2)
!!            resulm(1,3) = vec(1,3)
!!            resulm(2,3) = vec(2,3)                   
!            do iproc=1,numerop-1
!               call MPI_RECV(duma,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call MPI_RECV(dumb,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call MPI_RECV(dumc,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!!               call MPI_RECV(vecrc,6,MPI_INTEGER,iproc,
!!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               
!               if (duma.gt.duma1) then
!                   duma1 = duma
!!                   resulm(1,1) = vecrc(1,1)
!!                   resulm(2,1) = vecrc(2,1)
!               endif
!               
!               if (dumb.gt.dumb1) then
!                   dumb1 = dumb
!!                   resulm(1,2) = vecrc(1,2)
!!                   resulm(2,2) = vecrc(2,2)
!               endif
!               
!               if (dumc.gt.dumc1) then
!                   dumc1 = dumc
!!                   resulm(1,3) = vecrc(1,3)
!!                   resulm(2,3) = vecrc(2,3)
!               endif
!            
!            enddo
!            write(*,*)'max de u',duma1
!            write(*,*)'max de v',dumb1
!            write(*,*)'max de w',dumc1
!         endif
!      !-----------------------------------------------
!
!endif of Deltat updating
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
!!      drho     =       T             Drho                         !
!!       tnext   =   T(n+1)-T(n)        -                           !
!!******************************************************************!
      subroutine tauij2(drho,tnext,chwk,
     .           ten11,ten22,ten33,ten12,ten13,ten23,
     .           myid,rkstep)
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
      
      complex*8 :: chwk (0:my-1,0:mgalz-1,pb:pe),
     .          tnext (0:my-1,0:mgalz-1,pb:pe),
     .          drho  (0:my-1,0:mgalz-1,pb:pe),
     .          ten11 (0:my-1,0:mgalz-1,pb:pe),
     .          ten22 (0:my-1,0:mgalz-1,pb:pe),
     .          ten33 (0:my-1,0:mgalz-1,pb:pe),
     .          ten12 (0:my-1,0:mgalz-1,pb:pe),
     .          ten13 (0:my-1,0:mgalz-1,pb:pe),
     .          ten23 (0:my-1,0:mgalz-1,pb:pe)
               
!      real*8 rf0u(0:my-1),rf0w(0:my-1),dum

      real(4) :: uner(9),aa
      integer myid,rkstep,i,j,k,kk,ierr

! -------------------Programa-------------------------------!    
      
!  local transpose yz to zy, and adding zeros between positive
!  and negative wavenumbers.
      call localyz2zy(drho ,drho ,chwk) !T      
      call localyz2zy(tnext,tnext,chwk) !T(n+1)-T(n)    
      call localyz2zy(ten11,ten11,chwk) !S11
      call localyz2zy(ten22,ten22,chwk) !S22
      call localyz2zy(ten33,ten33,chwk) !S33
      call localyz2zy(ten12,ten12,chwk) !S12         
      call localyz2zy(ten13,ten13,chwk) !S13  
      call localyz2zy(ten23,ten23,chwk) !S23  
! inverse transform (fou --> fis @z),does the compacting aswell
      do i = pb,pe
         call fourz(drho (0,0,i),1) 
         call fourz(tnext(0,0,i),1) 
         call fourz(ten11(0,0,i),1)    
         call fourz(ten22(0,0,i),1)    
         call fourz(ten33(0,0,i),1)    
         call fourz(ten12(0,0,i),1)  
         call fourz(ten13(0,0,i),1)  
         call fourz(ten23(0,0,i),1)
      enddo
!  change plane to line(1:mx,lb:le)
!  and distribute package of lines to processors
      call chpl2ln(drho,drho,chwk,myid)
      call chpl2ln(tnext,tnext,chwk,myid)
      call chpl2ln(ten11,ten11,chwk,myid)
      call chpl2ln(ten22,ten22,chwk,myid)
      call chpl2ln(ten33,ten33,chwk,myid)
      call chpl2ln(ten12,ten12,chwk,myid)
      call chpl2ln(ten13,ten13,chwk,myid)
      call chpl2ln(ten23,ten23,chwk,myid)
!==========================================================================!
! before this call ten11,ten22,ten33,ten12,ten13, ten23,drho and tnext
! are still lines in fourier space x (phys z).
      call tensor2(drho,tnext,
     . ten11,ten22,ten33,ten12,ten13,ten23,
     . myid,rkstep)
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
      subroutine tensor2(drho,tnext,
     .          ten11,ten22,ten33,ten12,ten13,ten23, 
     .           myid,rkstep)

      use point
      use tem
      use fis
      use timacc 
      use wkhvect
      use cnan
      use rkcoef
      use MPI_GROUPS

      implicit none

      include "mpif.h"
      include "ctes3D"

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
      real(4) :: tnext(0:mx-1,lb:le),
     .        drho (0:mx-1,lb:le),
     .        ten11(0:mx-1,lb:le),
     .        ten22(0:mx-1,lb:le),ten33(0:mx-1,lb:le),
     .        ten12(0:mx-1,lb:le),ten13(0:mx-1,lb:le),
     .        ten23(0:mx-1,lb:le)
      real(8)  const1
      real(8) duma,dumb,dumc,duma1,dumb1,dumc1
      
!      ----------------------- Programa ------------------------      
      !sigma=0.7D0
      const1=2.0d0/3.0d0

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
         !DELTARHO
         do i=0,mgalx-1
             up1wk(i)  = -up2wk8(i)/(up1wk8(i)*(up1wk8(i)+up2wk8(i))
     .                   *dtbeta)
             call check_nan(up1wk(i),1,nanerror)
             !drho ready
         enddo
         !We include the normal stress component in the pressure
         !do i=0,mgalx-1
         !    tmp1wk8(i) =const1*(ten11wk8(i)+ten22wk8(i)+ten33wk8(i))
         !    !2/3Skk
         !enddo
         !S11,S22,S33 like: tau11=mu*(2S11-2/3Skk)
         !Sij like        : tauij=2*mu*Sij
         do i=0,mgalx-1
!Full formulation:
            ten11wk(i)=up1wk8(i)**sigma*2d0*ten11wk8(i)
!     .                                   -tmp1wk8(i))

            ten22wk(i)=up1wk8(i)**sigma*2d0*ten22wk8(i)
!     .                                   -tmp1wk8(i))

            ten33wk(i)=up1wk8(i)**sigma*2d0*ten33wk8(i)
!     .                                   -tmp1wk8(i))

            ten12wk(i)=2*up1wk8(i)**sigma*ten12wk8(i)

            ten13wk(i)=2*up1wk8(i)**sigma*ten13wk8(i)

            ten23wk(i)=2*up1wk8(i)**sigma*ten23wk8(i)
        enddo

!========================================================!
!IMPORTANT POINT: all quantities are now PHYSICAL
!========================================================!
c--------------------- back to F-P-P  
         call fourx(up1wk  ,up1wk8,-1)  ! drho
         call fourx(ten11wk,up1wk8,-1)  ! tau11
         call fourx(ten22wk,up1wk8,-1)  ! tau22
         call fourx(ten33wk,up1wk8,-1)  ! tau33
         call fourx(ten12wk,up1wk8,-1)  ! tau12
         call fourx(ten13wk,up1wk8,-1)  ! tau13
         call fourx(ten23wk,up1wk8,-1)  ! tau23

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
      
      if (nanerror.eq.1) then
         write(*,*) 'NAN found in',myid,'on drho lines'
         stop
      endif

      end
     
!--------------------------------------------------------!
!addstats1 subroutine
! AAF June 2014
!--------------------------------------------------------!
      subroutine addstats1(u1c,u2c,u3c,Tc,spaux,plane)
      
      use statis
      use wave
      use point
      use spectra
     
      implicit none
      include "mpif.h"
      include "ctes3D"

      !Modules to include
!.........................................................!
!.........................................................!

      complex*8 u1c(my,0:mz1),u2c(my,0:mz1),u3c(my,0:mz1),
     .           Tc(my,0:mz1)

      real*4 spaux(0:nz1,nspec+1,12)
      integer j,k,iy,kk,plane
      real(8) sigma
 

c     -----  spectra
c     ---- rhou(9),rhov(10),rhow(11),T(12)
      sigma = 0.7d0

      if (plane.eq.1)  then          ! i = 0
         do iy=1,nspec+1
            do kk = 0,mz1
               k = icx(kk)
               spaux(k,iy,9) = spaux(k,iy,9)+u1c(jsptot(iy),kk)*
     &                        conjg(u1c(jsptot(iy),kk))
               spaux(k,iy,10) = spaux(k,iy,10)+u2c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk))
               spaux(k,iy,11) = spaux(k,iy,11)+u3c(jsptot(iy),kk)*
     &                        conjg(u3c(jsptot(iy),kk))
               spaux(k,iy,12) = spaux(k,iy,12)+Tc(jsptot(iy),kk)*
     &                        conjg(Tc(jsptot(iy),kk))
            enddo
         enddo
      else
      do iy=1,nspec+1
            do kk = 0,mz1
               k = icx(kk)
               spaux(k,iy,9) = spaux(k,iy,9)+2.*u1c(jsptot(iy),kk)*
     &                        conjg(u1c(jsptot(iy),kk))
               spaux(k,iy,10) = spaux(k,iy,10)+2.*u2c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk))
               spaux(k,iy,11) = spaux(k,iy,11)+2.*u3c(jsptot(iy),kk)*
     &                        conjg(u3c(jsptot(iy),kk))
               spaux(k,iy,12) = spaux(k,iy,12)+2.*Tc(jsptot(iy),kk)*
     &                        conjg(Tc(jsptot(iy),kk))
            enddo
         enddo
      endif


      if (plane.eq.1) then             ! i=1
         do k=0,mz1
            do j=1,my
               Tp(j) = Tp(j) + Tc(j,k)*conjg(Tc(j,k))
            enddo
         enddo
      else
         do k=0,mz1
            do j=1,my
               Tp(j) = Tp(j) + 2.*Tc(j,k)*conjg(Tc(j,k))
            enddo
         enddo
      endif
       
!c ------------ MEAN only accumulated by MASTER node (plane = 1)
      if (plane.eq.1) then 
         do j=1,my
            rum(j) = rum(j)+u1c(j,0)
            rvm(j) = rvm(j)+u2c(j,0)
            rwm(j) = rwm(j)+u3c(j,0)
            Tm(j) = Tm(j)+Tc(j,0)
         enddo
      endif

      end subroutine

!--------------------------------------------------------------------!



!--------------------------------------------------------!
!addstats2 subroutine
! AAF June 2014
!--------------------------------------------------------!
      subroutine addstats2(u1c,u2c,u3c,o1c,o2c,o3c,
     .                     chwkc,uner,spaux,plane)
      
      !Modules to include
      use point
      use statis
      use diag
      use wave
      use fis
      use spectra
!.........................................................!
      implicit none
      include "mpif.h"
      include "ctes3D"

!.........................................................!

      complex*8 u1c(my,0:mz1),u2c(my,0:mz1),u3c(my,0:mz1),
     .        o1c(my,0:mz1),o2c(my,0:mz1),o3c(my,0:mz1),chwkc(my,0:mz1),
     .        chwkc2(my,0:mz1)

      integer j,k,iy,kk,plane
      real*4 uner(9)
      real*4 hyy
      real*8 aa
      complex*16 cc
      real*4 spaux(0:nz1,nspec+1,12)

!-------------Program-------------------------------------------!
c     -----  spectra
c     --u(1),v(2),w(3),uv_re(4),uv_im(5),o1(6),o2(7),o3(8)

      if (plane.eq.1)  then          ! i = 0
         do iy=1,nspec+1
            do kk = 0,mz1
               k = icx(kk)
               spaux(k,iy,1) = spaux(k,iy,1)+u1c(jsptot(iy),kk)*
     &                        conjg(u1c(jsptot(iy),kk))
               spaux(k,iy,2) = spaux(k,iy,2)+u2c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk))
               spaux(k,iy,3) = spaux(k,iy,3)+u3c(jsptot(iy),kk)*
     &                        conjg(u3c(jsptot(iy),kk))
               spaux(k,iy,4) = spaux(k,iy,4)+real(u1c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk)))
               spaux(k,iy,5) = spaux(k,iy,5)+imag(u1c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk)))
               spaux(k,iy,6) = spaux(k,iy,6)+o1c(jsptot(iy),kk)*
     &                        conjg(o1c(jsptot(iy),kk))
               spaux(k,iy,7) = spaux(k,iy,7)+o2c(jsptot(iy),kk)*
     &                        conjg(o2c(jsptot(iy),kk))
               spaux(k,iy,8) = spaux(k,iy,8)+o3c(jsptot(iy),kk)*
     &                        conjg(o3c(jsptot(iy),kk))
            enddo
         enddo
      else
      do iy=1,nspec+1
            do kk = 0,mz1
               k = icx(kk)
               spaux(k,iy,1) = spaux(k,iy,1)+2.*u1c(jsptot(iy),kk)*
     &                        conjg(u1c(jsptot(iy),kk))
               spaux(k,iy,2) = spaux(k,iy,2)+2.*u2c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk))
               spaux(k,iy,3) = spaux(k,iy,3)+2.*u3c(jsptot(iy),kk)*
     &                        conjg(u3c(jsptot(iy),kk))
               spaux(k,iy,4) = spaux(k,iy,4)+2.*real(u1c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk)))
               spaux(k,iy,5) = spaux(k,iy,5)+2.*imag(u1c(jsptot(iy),kk)*
     &                        conjg(u2c(jsptot(iy),kk)))
               spaux(k,iy,6) = spaux(k,iy,6)+2.*o1c(jsptot(iy),kk)*
     &                        conjg(o1c(jsptot(iy),kk))
               spaux(k,iy,7) = spaux(k,iy,7)+2.*o2c(jsptot(iy),kk)*
     &                        conjg(o2c(jsptot(iy),kk))
               spaux(k,iy,8) = spaux(k,iy,8)+2.*o3c(jsptot(iy),kk)*
     &                        conjg(o3c(jsptot(iy),kk))
            enddo
         enddo
      endif



      !Compute dV/dy
      call deryr2(u2c(1,0),chwkc(1,0),mz) !chwk=dv/dy       
      if (plane.eq.1) then             ! i=1
         do k=0,mz1
            do j=1,my
!              INTENSITIES
               ener(1) = u1c(j,k)*conjg(u1c(j,k))
               up(j) = up(j) + ener(1)
              
               ener(2) = u2c(j,k)*conjg(u2c(j,k))
               vp(j) = vp(j) + ener(2) 

               ener(3) = u3c(j,k)*conjg(u3c(j,k))
               wp(j) = wp(j) + ener(3)

               ener(4) = real(u1c(j,k)*conjg(u2c(j,k)))
               uvr(j)= uvr(j)+ ener(4)

               ener(5) = real(u1c(j,k)*conjg(u3c(j,k)))
               uwr(j)= uwr(j)+ ener(5)

               ener(6) = real(u3c(j,k)*conjg(u2c(j,k)))
               vwr(j)= vwr(j)+ ener(6)

               ener(7) = o1c(j,k)*conjg(o1c(j,k))
               w1p(j) = w1p(j) + ener(7)

               ener(8) = o2c(j,k)*conjg(o2c(j,k))
               w2p(j) = w2p(j) + ener(8)

               ener(9) = o3c(j,k)*conjg(o3c(j,k))
               w3p(j) = w3p(j) + ener(9)

               !compressibility effects
               cc = chwkc(j,k) + xbet(k)*u3c(j,k)
               thep(j) = thep(j) + cc*conjg(cc)
               !theup = theta*v
               theup(j) = theup(j) + real(cc*conjg(u2c(j,k)))


!              DISSIPATION
               aa =  bet2(k) *
     &               ( u1c(j,k)*conjg(u1c(j,k)) +
     &                 u2c(j,k)*conjg(u2c(j,k)) +
     &                 u3c(j,k)*conjg(u3c(j,k)) ) +
     &               chwkc(j,k)*conjg(chwkc(j,k))

               cc = o1c(j,k) + xbet(k)*u2c(j,k)
               aa = aa + cc*conjg(cc)
               cc = o3c(j,k)
               aa = aa + cc*conjg(cc)

               ep(j) = ep(j) + aa 


cc --------------- add this plane energy
               hyy = hy(j)
               do kk = 1,9
                  uner(kk) = uner(kk) + ener(kk)*hyy
               enddo

            enddo !end j loop
         enddo!end k loop

      else
         do k=0,mz1
            do j=1,my

               ener(1) = 2.*u1c(j,k)*conjg(u1c(j,k))
               up(j) = up(j) + ener(1)

               ener(2) = 2.*u2c(j,k)*conjg(u2c(j,k))
               vp(j) = vp(j) + ener(2)

               ener(3) = 2.*u3c(j,k)*conjg(u3c(j,k))
               wp(j) = wp(j) + ener(3)

               ener(4) = 2.*real(u1c(j,k)*conjg(u2c(j,k)))
               uvr(j) = uvr(j) + ener(4)

               ener(5)= 2.*real(u1c(j,k)*conjg(u3c(j,k)))
               uwr(j)= uwr(j) + ener(5)

               ener(6)= 2.*real(u3c(j,k)*conjg(u2c(j,k)))
               vwr(j)= vwr(j) + ener(6)
               
               ener(7)= 2.*o1c(j,k)*conjg(o1c(j,k))
               w1p(j)= w1p(j) + ener(7)

               ener(8) = 2.*o2c(j,k)*conjg(o2c(j,k))
               w2p(j)= w2p(j) + ener(8)

               ener(9) = 2.*o3c(j,k)*conjg(o3c(j,k))
               w3p(j)= w3p(j) + ener(9)

               !compressibility effects
               cc = xalp(plane-1)*u1c(j,k) +
     .           chwkc(j,k) + xbet(k)*u3c(j,k)
               thep(j) = thep(j) + 2.0*cc*conjg(cc)
               !theup = theta*v
               theup(j) = theup(j) + 2.0*real(cc*conjg(u2c(j,k)))


c           dissipation  ----------------

               aa = ( alp2(plane-1) + bet2(k) ) *
     &              ( u1c(j,k)*conjg(u1c(j,k)) +
     &                u2c(j,k)*conjg(u2c(j,k)) +
     &                u3c(j,k)*conjg(u3c(j,k)) ) +
     &              chwkc(j,k)*conjg(chwkc(j,k))

               cc = o1c(j,k) + xbet(k)*u2c(j,k)
               aa = aa + cc*conjg(cc)
               cc = o3c(j,k) - xalp(plane-1)*u2c(j,k)
               aa = aa + cc*conjg(cc)

               ep(j) = ep(j) + 2.*aa

cc ------------ add this plane energy
  
               hyy = hy(j)
               do kk = 1,9
                  uner(kk) = uner(kk) + ener(kk)*hyy
               enddo

            enddo !end j loop
         enddo !end k loop
      endif !end IF plane.eq.0
c ------------ MEAN only accumulated by MASTER node (plane = 1)
      if (plane.eq.1) then 
         do j=1,my
            um(j) = um(j)+u1c(j,0)
            vm(j) = vm(j)+u2c(j,0)
            wm(j) = wm(j)+u3c(j,0)
            w1m(j) = w1m(j)+o1c(j,0)
            w2m(j) = w2m(j)+o2c(j,0)
            w3m(j) = w3m(j)+o3c(j,0)
            them(j) = them(j)+chwkc(j,0) !diverg mean profile
         enddo
      endif




      end subroutine

!--------------------------------------------------------------------!

           
      subroutine cumtrapz(f,ff,init)
      use fis
      implicit none
      include "ctes3D"

      real(8),intent(in)::f(my)
      real(8),intent(in)::init
      real(8),intent(out)::ff(my)
      integer j
      !.------
      ff(1) =init
      do j=1,my-1
         ff(j+1)=ff(j)+0.5d0*(y(j+1)-y(j))*(f(j+1)+f(j))
      enddo
      end subroutine cumtrapz


        
  
 
     

