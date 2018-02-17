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
!     sp: spectral                                                     !
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
     .                  u00,w00,
     .                  rf0u,rf0w,u00wk,w00wk,
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
!     .                  spwk,
     .                  chwk2,
     .                  chwk,
     .                  myid)
!     .                  sp,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

      !------------- Parameters for SMR R-K ------------------------!
      
      real*8     gama(3), alpha(3), beta(3), ibeta(3), xi(3)
      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
!      parameter (alpha=(/ 29d0/96d0, -3d0/40d0, 1d0/6d0/))
!      parameter (beta =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /))
!      parameter (ibeta=(/ 160d0/37d0, 24d0/5d0,   6d0     /))
      parameter (alpha=(/ 4d0/15d0, 1d0/15d0, 1d0/6d0/))
      parameter (beta =(/ 4d0/15d0, 1d0/15d0,  1d0/6d0 /))
      parameter (ibeta=(/ 15d0/4d0, 15d0,   6d0     /))
      parameter (xi   =(/        0d0,-17d0/60d0, -5d0/12d0 /)) !aaf 
   
 
!     ---------------------   commons ------------------------------- 

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .        ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save   /tem/

      integer nimag,nstep,nhist,ihist,icfl,ncfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
      save   /timacc/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(0:my-1),hy(my),fmap(my),
     .            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/
      
      integer nacumsp,jsptot
      common /spectra/   nacumsp,jsptot(2*nspec+1)
      save   /spectra/
      
      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save /ficheros/
      character*84 fname
      
      real*8 commtimer,transtimer,totaltimer,ctm,ctm2,ttm
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

      integer nanerror,nanproblem
      common /cnan/ nanerror
      save   /cnan/


      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/

      
c ------------------- Variables ----------------------------------------

      complex*8 phi(0:my-1,0:mz1,pb:pe),vor(0:my-1,0:mz1,pb:pe),
     .          psi(0:my-1,0:mz1,pb:pe),scal(0:my-1,0:mz1,pb:pe)

      complex*8 nxwk  (0:my-1,0:mgalz-1,pb:pe),
     .          nzwk  (0:my-1,0:mgalz-1,pb:pe),
     .          nywk  (0:my-1,0:mgalz-1,pb:pe),
     .          mxwk  (0:my-1,0:mgalz-1,pb:pe),
     .          mywk  (0:my-1,0:mgalz-1,pb:pe),
     .          mzwk  (0:my-1,0:mgalz-1,pb:pe),
     .          rhst  (0:my-1,0:mgalz-1,pb:pe),
     .          drho  (0:my-1,0:mgalz-1,pb:pe),
     .          chwk  (0:my-1,0:mgalz-1,pb:pe),
     .          chwk2 (0:my-1,0:mgalz-1,pb:pe),
     .          ten12 (0:my-1,0:mgalz-1,pb:pe),
     .          ten13 (0:my-1,0:mgalz-1,pb:pe),
     .          ten23 (0:my-1,0:mgalz-1,pb:pe)

      real*8 u00(0:my-1),w00(0:my-1),rf0u(0:my-1),
     .       rf0w(0:my-1),u00wk(0:my-1),w00wk(0:my-1)
  
!      real*4 sp  (0:nz1,1:2*nspec+1,8,pb:pe),
!     .       spwk(0:nz1,1:  nspec+1,8,pb:pe)
      
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,kk
      
      real*8  rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire,
     .        dtbeta
      
      real*8  iter_time,write_time
      complex*8 temp1,temp2
      real*4  temporal(3),comun(3)
      
      real*4  reynota,H00u,H00w,dumu,dumw,cteu,ctew
      real*4  massu,massw,massu1,massu2,massw1,massw2,dum
      real*4  bcpsi00

      complex*8,allocatable :: tempa(:,:),tempb(:,:)

!debugging variables
      real*8  kkita
      character*6 varname
      complex*8 v1(0:my1),v2(0:my1),v3(0:my1),v4(0:my1) 
      integer box
      
c    ---------------------- Program ------------------------------------
    
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
      istati   = 0
      irun     = 0   ! first time step is special in tim3rkp
      icfl     = 1   ! first time step always needs a step size
      nanerror = 0   ! Check NAN in the 50 firsts steps
      nanproblem = 0 ! initializes NAN controller

      if(myid.eq.0) then
         fname=filstt(1:index(filstt,' ')-1)//'.cf'
         write (*,*) fname
         open(29,file=fname,status='unknown')
      endif
!debug
      !write(*,*) 'WITHIN CROSS',myid,'scal(1)',scal(0,0,pb)
c ========================================================================
c                 THIS IS THE TIME LOOP
c ========================================================================

      do 30 istep=1,nstep
         
         
         if (myid.eq.0) then
           totaltimer = totaltimer-MPI_WTIME()
           iter_time=-MPI_WTIME()
         endif

         if (mod(istep-1,nhist).eq.0) ihist=1
         if (mod(istep-1,ntimes).eq.0.and.nstart.ne.0) istati=1
!        nstart is a flag, if 0 no stats, otherwise make stats each
!        ntimes
         if (mod(istep-1,ncfl).eq.0) icfl = 1
         
! ------------------- write image to a file ---------------------------!

         IF (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then

            if (myid.eq.0) then
               write_time = -MPI_WTIME()
            endif
            do j=0,my1
               u00(j) = u00wk(j) 
               w00(j) = w00wk(j)
            enddo
!aaf disabling spectra 
            !call escru(vor,phi,u00,w00,sp,spwk,myid)
            call escru(vor,phi,u00,w00,myid)
!write scalar field
!!            call escruscal(scal,1,myid)
!aaf Here we need to call a escru for scalar field

            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
            endif
         ENDIF

! ------------------- finished writing image --------------------------!

      if (myid.gt.numerop-1) goto 30    ! only for save procs

      if(irun.eq.0) then ! this is done only for the first step
           irun = 1
         if (myid.eq.0) then
            do j=0,my1
              u00wk(j)=u00(j)
              w00wk(j)=w00(j)
            enddo
         endif
      endif ! end special first step

               !  Runge-Kutta third order  !
!=======================================================================!                
!----------- RUNGE-KUTTA SUBSTEPS START---------------------------------!
!=======================================================================!                

      do 10 rkstep=1,3
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
               call Lapvdvhom(phi(0,k,i),mywk(0,k,i),nzwk(0,k,i),rK)
!               call Lapvdv(phi(0,k,i),mywk(0,k,i),nzwk(0,k,i),rK)
            enddo
         enddo
!=====================DEBUGGING!!!!=================================!
!debugging
         if (rkstep.eq.1.and.istep.eq.500) then   
           do i=pb,pe
             !mode kx=1
             box=2;
             if (i.eq.box) then
                do j=0,my1
                  v1(j)=mywk(j,1,i)
                enddo
                !Write the velocity of this mode
                varname='000v11'
                call write2file(v1,varname,my,1,1,myid)

                !Computate and write dvdy
                call deryr2(v1,v2,1)
                varname='00dv11'
                call write2file(v2,varname,my,1,1,myid)

                !Computate and write dvdvdy
                call deryr2(v2,v3,1)
                varname='dvdv11'
                call write2file(v3,varname,my,1,1,myid)

                !Computate and write dvdvdvdy
                call deryr2(v3,v2,1)
                varname='0d3v11'
                call write2file(v2,varname,my,1,1,myid)
                  
                !Compute directly second derivative
                call deryyr2(v1,v3,1) 
                varname='0d2v11'
                call write2file(v3,varname,my,1,1,myid)
 
                !Compute derivative of second derivative
                call deryr2 (v3,v1,1) 
                varname='d2dv11'
                call write2file(v1,varname,my,1,1,myid)

             write(*,*) 'v11 written by myid=', myid

             endif

             box=4;
             !mode(4,4)
             if (i.eq.box) then
                do j=0,my1
                  v1(j)=mywk(j,3,i)
                enddo
                !Write the velocity of this mode
                varname='000v44'
                call write2file(v1,varname,my,1,1,myid)

                !Computate and write dvdy
                call deryr2(v1,v2,1)
                varname='00dv44'
                call write2file(v2,varname,my,1,1,myid)

                !Computate and write dvdvdy
                call deryr2(v2,v3,1)
                varname='dvdv44'
                call write2file(v3,varname,my,1,1,myid)

                !Computate and write dvdvdvdy
                call deryr2(v3,v2,1)
                varname='0d3v44'
                call write2file(v2,varname,my,1,1,myid)
                  
                !Compute directly second derivative
                call deryyr2(v1,v3,1) 
                varname='0d2v44'
                call write2file(v3,varname,my,1,1,myid)
 
                !Compute derivative of second derivative
                call deryr2 (v3,v1,1) 
                varname='d2dv44'
                call write2file(v1,varname,my,1,1,myid)

             write(*,*) 'v44 written by myid=', myid

             endif
             !mode(15,15)
             box=16;
             if (i.eq.box) then
                do j=0,my1
                  v1(j)=mywk(j,box-1,i)
                enddo
                !Write the velocity of this mode
                varname='000v15'
                call write2file(v1,varname,my,1,1,myid)

                !Computate and write dvdy
                call deryr2(v1,v2,1)
                varname='00dv15'
                call write2file(v2,varname,my,1,1,myid)

                !Computate and write dvdvdy
                call deryr2(v2,v3,1)
                varname='dvdv15'
                call write2file(v3,varname,my,1,1,myid)

                !Computate and write dvdvdvdy
                call deryr2(v3,v2,1)
                varname='0d3v15'
                call write2file(v2,varname,my,1,1,myid)
                  
                !Compute directly second derivative
                call deryyr2(v1,v3,1) 
                varname='0d2v15'
                call write2file(v3,varname,my,1,1,myid)
 
                !Compute derivative of second derivative
                call deryr2 (v3,v1,1) 
                varname='d2dv15'
                call write2file(v1,varname,my,1,1,myid)
             write(*,*) 'v15 written by myid=', myid, pb,pe


             endif

           enddo
         endif



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

!        Now build u_i+RHS(1-1)*Dt*eps_i
!        -------------------------------
         do i=pb,pe
           do k=0,mz1
             do j=0,my-1
               vor(j,k,i)=vor(j,k,i)+ten12(j,k,i)*dtxi!Omega_y
               phi(j,k,i)=phi(j,k,i)+ten13(j,k,i)*dtxi!PSI
!              Similarly for Temperature(scal and rhst)
!              Save T in chwk first then change scal
!commenting all related with Temperature
!               !drho(j,k,i)=scal(j,k,i)  !save T in drho too
!               !scal(j,k,i)=scal(j,k,i)+rhst(j,k,i)*dtxi 
!               !rhst(j,k,i)=scal(j,k,i)  !save T in rhst too

             enddo
!note: we need two T -> one goes to physical, the other one DOES NOT
           enddo
         enddo
!        d(psi)/dy comes from last substep in ten23, but if we are
!        in the very first substep we dont have it yet :(
!        Prepare d(psi)/dy 
!debug flag:  do it always
!         if (istep.eq.1.and.rkstep.eq.1) then
!             do i=pb,pe
!                call deryr2(nywk(0,0,i),ten23(0,0,i),mz)
!                chwk= d(psi)/dy= d(nywk)/dy -- F-F
!             enddo
!         endif 
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
!   ten12: RHS(Omega_y)(n-1)
!   ten13: RHS(phi)(n-1)
!   ten23: d(psi)/dy
!   rhst : T
!   drho : T
!*******************************************************************!

!-----------------(3) Obtain rho*u ----------------------------!
!
!         do i=pb,pe
!             do k=0,mz1
!                do j=0,my-1
!                   !rhou=mx+ikx*psi
!!for debugging eliminating psi contributtion
!                   !mxwk(j,k,i)=mxwk(j,k,i)+xalp(i-1)*nywk(j,k,i)    
!                   mxwk(j,k,i)=mxwk(j,k,i)    
!                   !rhow=mz+ikz*psi
!                   !mzwk(j,k,i)=mzwk(j,k,i)+xbet(k)*nywk(j,k,i)    
!                   mzwk(j,k,i)=mzwk(j,k,i)    
!                   !rhov=my+d(psi)/dy
!                   !mywk(j,k,i)=mywk(j,k,i)+ten23(j,k,i)    
!                   mywk(j,k,i)=mywk(j,k,i)    
!                enddo
!             enddo
!         enddo
!!---------- Adding modes 00 to rhou and rhow------------------------!
!---------- Only master computes 00  modes
        if (myid.eq.0) then
           do j=0,my-1  
              mxwk(j,0,1) = u00wk(j)
              mzwk(j,0,1) = w00wk(j)
           enddo
        endif
 

!debugging
        if (myid.eq.0) then
          if (istati.eq.1.and.rkstep.eq.1) then
          kkita=0d0
            do j=0,my-1           
               kkita=kkita+(hy(j+1)*(0.25-(u00wk(j)/
     .              (u00wk(my-1)-u00wk(0)))**2))

            enddo
            write(*,*) "dm using u00wk-->dm=",kkita 
!           write(*,*) "u00wk",(u00wk(j),j=0,my1)
           endif
         endif
!debug flag
!        if (myid.eq.0) then
!           write(*,*) "modes 00:"
!           write(*,*) "istep",istep,"rkstep",rkstep
!           write(*,*) "rhou(2)",mxwk(1,0,1),"rhou(my)",mxwk(my-1,0,1)
!           write(*,*) "rhow(2)",mzwk(1,0,1),"rhow(my)",mzwk(my-1,0,1)
!           write(*,*) "rhov(2)",mywk(1,0,1),"rhov(my)",mywk(my-1,0,1)
!           write(*,*) "dpsi(1)",ten23(0,0,1),"dpsi(my)",ten23(my-1,0,1)
!        endif

!       !Now build in u00/w00 the variable u00**/w00**
         if (myid.eq.0) then
           do j=0,my1
              u00wk(j)=u00wk(j)+dtxi*rf0u(j)
              w00wk(j)=w00wk(j)+dtxi*rf0w(j)
           enddo
        endif
!-------------------------------------------------------------------! 



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
!not computing temperature stuff
!         do i=pb,pe
!            call deryr2(rhst(0,0,i),ten13(0,0,i),mz)
!         enddo
         !calculate laplacian of T and save it in "nywk"
!         do i=pb,pe
!            call laplacian(rhst(0,0,i),nywk(0,0,i),mz,i)
!         enddo 


!write file debugging
!         if ((rkstep.eq.1).and.(istep.eq.1).and.(myid.eq.0)) then
!           fname=filstt(1:index(filstt,' ')-1)//'.kk'
!           open(69,file=fname,status='unknown',form='unformatted')
!           write(69)     (y(j), j=0,my1),
!     .                   (hy(j),j=1,my),                
!     .                (real(rhst(0,k,1)),k=0,mz1),
!     .                (real(nywk(0,k,1)),k=0,mz1)
!            write(*,*) "rhst(0,k,1)@rkstep=1,istep=1"
!            write(*,*) (real(rhst(0,k,1)),k=0,mz1)
!            write(*,*) "nywk(0,k,1)"
!            write(*,*) (real(nywk(0,k,1)),k=0,mz1)
!           close(69)
!         endif


!flag debug ACTHUNG
!      if (myid.eq.0) then 
!         write(*,*) "istep",istep,"rkstep",rkstep
!         write(*,*) "lap(T)00",(nywk(j,0,pb),j=256,258)
!      endif


!---------!Calculate dT/dz,dT/dx------------------------------------!
!Eliminating T related stuff
!         do i=pb,pe
!            do k=0,mz1
!               do j=0,my-1
!                  !dT/dz=i kz T
!                  ten23(j,k,i)=xbet(k)  *rhst(j,k,i)
!                  !dT/dx=i kx T
!                  ten12(j,k,i)=xalp(i-1)*rhst(j,k,i)
!               enddo
!            enddo
!        enddo

        !debugging flag
!        if (myid.eq.0) then
!           !write(*,*) "v(-10:10)", (mywk(j,0,pb),j=0,my1)
!           write(*,*) "rw(-10)", mzwk(0,0,pb),"rw(10)=",mzwk(my-1,0,pb)
!        endif


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

!------------aaf April 1st 2014, checked up to this point............!

!        !debugging flag
!        if (myid.eq.0) then
!           write(*,*) "u(-10)", mxwk(0,0,pb),"u(10)=",mxwk(my-1,0,pb)
!           write(*,*) "w(-10)", mzwk(0,0,pb),"w(10)=",mzwk(my-1,0,pb)
!           write(*,*) "ruu(-10)", nxwk(0,0,pb),"ruu(10)=",
!     .                            nxwk(my-1,0,pb)
!        endif

        if (rkstep.eq.3.and.ihist.eq.1.and.myid.eq.0) then
             H00u=0.
             H00w=0.
             do j=1,my-2
                H00u = H00u + rf0u(j)*(y(j+1)-y(j-1))
                H00w = H00w + rf0w(j)*(y(j+1)-y(j-1))
             enddo
             H00u = 0.25*H00u
             H00w = 0.25*H00w
         endif

!Computes Dt*gamma, Dt*xi before calculating this substep CFL(calculated
!within fouphysfou subroutine
         dtgamma = Deltat*gama(rkstep)
         dtxi    = Deltat*xi(rkstep)
         dtbeta  = Deltat*beta(rkstep)


!-------------(6a) Evolving Temperature ----------------------------!       
!Debug FLAG not evolvin the temperature
!ACHTUNG!
!          do i=pb,pe
!              call rkstepexp(scal(0,0,i),rhst(0,0,i), mz,i,dtgamma)
!          enddo

!!aaf debug flag
!        if (myid.eq.0) then
!           write(*,*) "check T(n+1):",(real(scal(j,0,pb)),j=256,258)
!        endif

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
            do j=0,my-1
               rf0u(j)= real(nxwk(j,0,1))
               rf0w(j)= real(nzwk(j,0,1))
            enddo
         endif  
!-------------------------------------------------------------------!




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
! Calculate all tensor Sij
          !S12, S13 and S23
          do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  !S12
                  ten12(j,k,i)=(ten12(j,k,i)+
     .                         xalp(i-1)*mywk(j,k,i))/2
                  !S13
                  ten13(j,k,i)=(xbet(k)*mxwk(j,k,i)+
     .                          xalp(i-1)*mzwk(j,k,i))/2
                  !S23
                  ten23(j,k,i)=(ten23(j,k,i)+
     .                         xbet(k)*mywk(j,k,i))/2
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
!   ten12: S12
!   ten13: S13
!   ten23: S23
!   drho : T                                                        
!   rhst : RHS(T) -->RIP buffer, we keep it until the end of substep
!*******************************************************************!
!WARNING:
! we need to enter tauij with T(n+1), but we need a buffer of size
! mgalz, and T(n+1) is in scal(size mz-1),so that before entering
! tauij we will save scal-->chwk2
! This way we dont need to do the PHYS-FOU of T(n+1)
!eliminating T stuff
!       do i=pb,pe
!           do k=0,mz1
!              do j=0,my-1
!                  chwk2(j,k,i)=scal(j,k,i)!now nxwk has T(n+1)
!              enddo
!           enddo
!       enddo
!WARNING
!------------
!chwk2=T(n+1)

!ACHTUNG!!
!For incompressible flow we don't need to do tauij
!
!!------------(7) Go to Phys and compute TAUij and Drho (and back)----!
!        call tauij(mxwk,mywk,mzwk,ten12,ten13,ten23,
!     .                  drho,chwk2,chwk,myid,rkstep)
!! outputs: 
!! mxwk = tau11 
!! mywk = tau22
!! mzwk = tau33
!! ten12= tau12
!! ten13= tau13
!! ten23= tau23
!! drho = drho
!!---------------------------------------------------------------------!
!For debugging and checking we avoid tauij for incompressible code
!instead:
       do i =pb,pe
          do k=0,mz1
             do j=0,my1
                 mxwk(j,k,i)=2*mxwk(j,k,i)
                 mywk(j,k,i)=2*mywk(j,k,i)
                 mzwk(j,k,i)=2*mzwk(j,k,i)
                 ten12(j,k,i)=2*ten12(j,k,i)
                 ten13(j,k,i)=2*ten13(j,k,i)
                 ten23(j,k,i)=2*ten23(j,k,i)
              enddo
           enddo
        enddo
!------------------------------------------------------------



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
            do j=0,my-1
                rf0u(j)=rf0u(j)+ire*real(mxwk(j,0,1)) !RHS(rhou00)
                rf0w(j)=rf0w(j)+ire*real(mzwk(j,0,1)) !RHS(rhow00)
            enddo
         endif  
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
!

!------------(8b) Obtain psi(n+1)-----------------------------------!
!    Calculate lap(psi) needed to obtain lap(psi(n+1))
!          do i=pb,pe
!             call laplacian(psi(0,0,i),ten12(0,0,i),mz,i)
!             !ten12=lap(psi)
!          enddo 
!flag debug ACTHUNG
!      if (myid.eq.0) then 
!         write(*,*) "drho",(drho(j,0,pb),j=256,258)
!        ! write(*,*) "lap(psi)",(ten12(j,0,pb),j=256,258)
!      endif
!    Obtain lap(psi(n+1))
!eliminating the updating of psi--Incompressible PSI=0 all day long
!          do i=pb,pe
!             do k=0,mz1
!                do j=0,my-1
!                   psi(j,k,i)=-drho(j,k,i)/dtbeta-
!     .             alpha(rkstep)/beta(rkstep)*ten12(j,k,i)
!                   !lap(psi(n+1))
!                enddo
!             enddo
!          enddo
!Calculate integral of(lap(psi(0,0)))
!         if (myid.eq.0) then !only master do this
!            bcpsi00=0d0 !reset counter
!            do j=1,my
!               bcpsi00=bcpsi00+psi(j,0,pb)*hy(j)
!            enddo
!!lappsi00 keeps the int(lap(PSI00))
!!This is the BC for lapPSI solver in mode00 at the TOP
!            bcpsi00=0.5*bcpsi00
!!debug flag
!            write(*,*) "bcpsi00",bcpsi00,"rkstep",rkstep

!        endif
!  Solve laplacian to find psi(n+1)
!         do i=pb,pe
!            do k=0,mz1
!               k1 = icx(k)
!               rK = bet2(k1)+alp2(i-1)
!!              call Lapvdv(psi(0,k,i),psi(0,k,i),ten23(0,k,i),rK)
!              call Lapvdvhom2(psi(0,k,i),psi(0,k,i),ten23(0,k,i),rK)
!              !Lapvdvhom2 solves 00 mode
!            enddo
!         enddo
!Now mod00 of lap(psi) has to be solve differently
!         if (myid.eq.0) then 
!debug:changed in order to find only the particular solution
!            call Lappsi00(psi(0,0,pb),psi(0,0,pb),ten23(0,0,pb),
!     .                    bcpsi00)
!!debug flag: ACHTUNG
!            write(*,*) "psi00(i+1)",(psi(j,0,pb),j=256,258)
!         endif

! Now:
! ten23 = d(psi(n+1)/dy) --> will be recycled next substep
! psi   = psi(n+1)

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
          do i=pb,pe
              call rkstepexp(vor(0,0,i),ten12(0,0,i), mz,i,dtgamma)
              call rkstepexp(phi(0,0,i),ten13(0,0,i), mz,i,dtgamma)
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
!
       

!rkstep=3
          if (rkstep.eq.3) then
            if (mod(istep,ncfl).eq.0) then
              if (myid.eq.0) then 
                write(*,*) 'he entrado con ', istep
              !nothing to do here now. 
              endif        
            endif  
         endif

! ----------------------------------------------------------------------!
!             this is the block for the 00 mode
!                   Only master does it
! ----------------------------------------------------------------------!
!        call mpi_barrier(MPI_COMM_WORLD,ierr) 
        if (myid.eq.0) then
! evolve modes 00 por rhou, rhow
           do j=0,my1
             u00wk(j)=u00wk(j)+dtgamma*rf0u(j)
             w00wk(j)=w00wk(j)+dtgamma*rf0w(j)
           enddo
!Disabling rk00exp--> make it simpler.
!           call rk00exp(u00wk,rf0u,dtgamma)
!           call rk00exp(w00wk,rf0w,dtgamma)
!debug flag 
!        write(*,*) "u00",(u00(j),j=0,my-1)
!compute mass
           massu=0.
           massw=0.
           do j=0,my-1
              !massu = massu+u00(j)*trp(j)
              !massw = massw+w00(j)*trp(j)
              massu = massu+u00wk(j)*hy(j+1)
              massw = massw+w00wk(j)*hy(j+1)
           enddo

!write file debugging
!         if ((rkstep.eq.1).and.(istep.eq.1)) then
!           fname=filstt(1:index(filstt,' ')-1)//'.kk'
!           open(69,file=fname,status='unknown',form='unformatted')
!           rewind(69)
!           write(69)     (y(j), j=0,my1),
!     .                   (hy(j),j=1,my),                
!     .                (u00wk(j),j=0,my1),
!     .                (w00wk(j),j=0,my1)
!           close(69)
!         endif

        endif   ! 00 modes

 10   continue  !!! End of R-K sub-step


! ----------------------------------------------------------------------!
!   Wx0 ans Wz0 are the box averaged vorticities at the wall            !
! ----------------------------------------------------------------------!

c        !     write history record     !

         if(myid.eq.0) then
!calculate eps

            if (ihist.eq.1) then
            reynota=0.
            !calculating int(ep)
            do j=1,my
               reynota=reynota+ep(j)*hy(j)
            enddo
             
               if (istep.eq.1) then
                   ctm2=0.
                   ttm =0.
               endif
324            format(17(d22.14))
               write(29,324) time,reynota,dm,ener,
!               write(29,324) time,-1.*Wz0,reynota,dm,ener,
!               write(29,324) time,-1.*Wz0,WzL,sqrt(reynota),ener,
     .                       massu,massw,MPI_WTIME()+iter_time,
     .                       commtimer-ctm2
               call flush(29)
325            format(i5,12(d14.6))
               write(*,325)istep,time,reynota,dm,Deltat,
!               write(*,325)istep,time,-1.*Wz0,reynota,dm,Deltat,
!               write(*,325)istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,
     .                     massu,massw,H00u,H00w,MPI_WTIME()+iter_time,
     .                    transtimer-ttm,commtimer-ctm2



            endif
            ctm2 = commtimer
            ttm  = transtimer 
            
         endif

c          	! time:
         time=time+Deltat
!debug flag
        if (myid.eq.0) then
            write(*,*) "time", time
        endif


         if(istati.eq.1) istati = 0
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
         write(*,*) 'Restarting from the last file...'
         
!initnan subroutine OUT OF ORDER          
!         call initnan(vor,phi,u00,w00,rf0u,rf0w,ch00wk,w00wk,nxwk,nzwk,
!     .                nywk,spwk,mxwk,mywk,mzwk,chwk,ten12,ten13,
!     .                ten23,scal,sp,myid)
      endif
      

      endsubroutine

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
!       ten12   =      dT/dx          rhouv                        !
!       ten13   =      dT/dy          rhouw                        !
!       ten23   =      dT/dz          rhovw                        !
!       rhst    =        T            RHS(T)                       !
!                                                                  !
!******************************************************************!
      subroutine foup(up1,up2,up3,ten11,ten22,ten33,
     .           rhst,ten12,ten13,ten23,chwk,myid,rkstep)

      implicit none

      include "mpif.h"
      include "ctes3D"

c  -------------------------   commons ------------------------- 

      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 

      integer nacumsp,jsptot
      common /spectra/   nacumsp,jsptot(2*nspec+1)
      save   /spectra/
      
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      integer nimag,nstep,nhist,ihist,icfl,ncfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
      save   /timacc/

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/

      real*8 commtimer,transtimer,totaltimer
      common /timers/ commtimer, transtimer, totaltimer
      save   /timers/
      
      integer nanerror 
      common /cnan/ nanerror
      save   /cnan/  

c     ---------------------- Variables ----------------------
      
      complex*8 up1   (0:my-1,0:mgalz-1,pb:pe),
     .          up2   (0:my-1,0:mgalz-1,pb:pe),
     .          up3   (0:my-1,0:mgalz-1,pb:pe),
     .          ten11 (0:my-1,0:mgalz-1,pb:pe),
     .          ten22 (0:my-1,0:mgalz-1,pb:pe),
     .          ten33 (0:my-1,0:mgalz-1,pb:pe),
     .          ten12 (0:my-1,0:mgalz-1,pb:pe),
     .          ten13 (0:my-1,0:mgalz-1,pb:pe),
     .          ten23 (0:my-1,0:mgalz-1,pb:pe),
     .          rhst  (0:my-1,0:mgalz-1,pb:pe)
               
      

      real*4 chwk(*),uner(9),aa
!      real*4 sp(0:nz1,1:2*nspec+1,8,pb:pe)
      
      integer myid,rkstep,i,j,k,kk,ierr
      
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     
      
c     ---------------------- Programa ----------------------      

c     Before Fpp, computes statistics
!ACHTUNG: this is BULLSHIT NOW, WE DONT HAVE U,V,W ANYMORE!!we have
!rhou,rhov,rhow
!       
      if (istati.eq.1.and.rkstep.eq.1) then
         nacum   = nacum   +1
         nacumsp = nacumsp +1
         do kk=1,9
            uner(kk) = 0.
         enddo
         do i=pb,pe  
            call cstatis(up1(0,0,i),up2(0,0,i),up3(0,0,i),
     .                  ten11(0,0,i),ten22(0,0,i),ten33(0,0,i),
     .                   chwk,uner,i)
!     .                   chwk,sp(0,1,1,i),uner,i)
         enddo 
         if (ihist.eq.1) then
             call MPI_ALLREDUCE(uner,ener,9,MPI_REAL,MPI_SUM,
     .                          MPI_COMM_CALC,ierr)
         
           do i=1,9
               ener(i)=sqrt(abs(ener(i)))
            enddo

!!aaf------- momentum thickness calculation---------!
!! dm will keep the momentum thickness
            dm=0d0
            do j=1,my
               dm=dm+(hy(j)*(0.25-(um(j)/(um(my)-um(1)))**2))
            enddo
         endif 
      endif    
!

c ----- substract umax/2 to u00 to increase dt !!!!
!      if (myid.eq.0) then
!         do j=0,my-1
!           up1(j,0,1) = up1(j,0,1) - a0
!         enddo
!      endif
!  local transpose yz to zy, and adding zeros between positive
!  and negative wavenumbers.
      call localyz2zy(up1,up1,chwk) !rhou
      call localyz2zy(up2,up2,chwk) !rhov
      call localyz2zy(up3,up3,chwk) !rhow
!      !TEMPERATURE
!      !call localyz2zy(ten22,ten22,chwk) !Lap(T) 
!      !call localyz2zy(ten12,ten12,chwk) !dT/dx         
!      !call localyz2zy(ten13,ten13,chwk) !dT/dy  
!      !call localyz2zy(ten23,ten23,chwk) !dT/dz        
!      !call localyz2zy(rhst ,rhst ,chwk) !T      
       
! inverse transform (fou --> fis @z),does the compacting aswell
      do i = pb,pe
         call fourz(up1(0,0,i),1)    !rhou
         call fourz(up2(0,0,i),1)    !rhov
         call fourz(up3(0,0,i),1)    !rhow
!         !TEMPERATURE
!         !call fourz(ten22(0,0,i),1)  !lap(T)
!         !call fourz(ten12(0,0,i),1)  ! dT/dx
!         !call fourz(ten13(0,0,i),1)  ! dT/dy
!         !call fourz(ten23(0,0,i),1)  ! dT/dz
!         !call fourz(rhst (0,0,i),1)  ! T
      enddo

      nanerror = 0

!  change plane to line(1:mx,lb:le)
!  and distribute package of lines to processors
      call chpl2ln(up1,up1,chwk,myid)
      call chpl2ln(up2,up2,chwk,myid)
      call chpl2ln(up3,up3,chwk,myid)
!temperature
!      !call chpl2ln(ten22,ten22,chwk,myid)
!      !call chpl2ln(ten12,ten12,chwk,myid)
!      !call chpl2ln(ten13,ten13,chwk,myid)
!      !call chpl2ln(ten23,ten23,chwk,myid)
!      !call chpl2ln(rhst ,rhst ,chwk,myid)

!==========================================================================!
! before this call upi,ten22,ten12,ten13, ten23 and rhst
! are still lines in fourier space x (phys z).
      call phys1(up1,up2,up3,ten11,ten22,ten33,
     .     ten12,ten13,ten23,rhst,myid,rkstep)
!==========================================================================!

c  ---------- back to the yz planes an fft 

!      !call chln2pl(rhst,rhst,chwk,myid)   !RHS(T)
      call chln2pl(up1,up1,chwk,myid)     !u
      call chln2pl(up2,up2,chwk,myid)     !v
      call chln2pl(up3,up3,chwk,myid)     !w
      call chln2pl(ten11,ten11,chwk,myid) !rhouu
      call chln2pl(ten22,ten22,chwk,myid) !rhovv  
      call chln2pl(ten33,ten33,chwk,myid) !rhoww
      call chln2pl(ten12,ten12,chwk,myid) !rhouv  
      call chln2pl(ten13,ten13,chwk,myid) !rhouw
      call chln2pl(ten23,ten23,chwk,myid) !rhovw

!  convert to fourier Z now.
      do i=pb,pe
!         !call fourz(rhst(0,0,i),-1)   !RHS(T)
         call fourz(up1(0,0,i),-1)   ! u
         call fourz(up2(0,0,i),-1)   ! v
         call fourz(up3(0,0,i),-1)   ! w
         call fourz(ten11(0,0,i),-1)  !rhouu
         call fourz(ten22(0,0,i),-1)  !rhovv
         call fourz(ten33(0,0,i),-1)  !rhoww
         call fourz(ten12(0,0,i),-1)  !rhouv
         call fourz(ten13(0,0,i),-1)  !rhouw
         call fourz(ten23(0,0,i),-1)  !rhovw
      enddo
!transpose back
!      !call localzy2yz(rhst,rhst,chwk)
      call localzy2yz(up1,up1,chwk)
      call localzy2yz(up2,up2,chwk)
      call localzy2yz(up3,up3,chwk) 
      call localzy2yz(ten11,ten11,chwk) 
      call localzy2yz(ten22,ten22,chwk) 
      call localzy2yz(ten33,ten33,chwk) 
      call localzy2yz(ten12,ten12,chwk)   
      call localzy2yz(ten13,ten13,chwk)   
      call localzy2yz(ten23,ten23,chwk)   
           
      endsubroutine


!*********************************************************************!
!                                                                     !
!         transform to physical and do operations                     !
!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!                                                                     !
!                                                                     !
!*********************************************************************!
      subroutine phys1(up1,up2,up3,ten11,ten22,ten33,
     .           ten12,ten13,ten23,rhst,myid,rkstep)

      implicit none

      include "mpif.h"
      include "ctes3D"

c --------------------------- Commons --------------------------


      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 


      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save   /tem/
      
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/

      integer nimag,nstep,nhist,ihist,icfl,ncfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
      save   /timacc/

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/



      real*4 up1wk,up2wk,up3wk,ten11wk,
     .       ten22wk, ten33wk,ten12wk,ten13wk,ten23wk,
     .       rhstwk
      real*8 up1wk8,up2wk8,up3wk8,
     .       ten11wk8,ten22wk8,ten33wk8,ten12wk8,ten13wk8,
     .       ten23wk8,rhstwk8,tmp1wk8,tmp2wk8,tmp3wk8
      common /wkhvect/ up1wk (0:mgalx+1),up2wk(0:mgalx+1),
     .             up3wk(0:mgalx+1),ten11wk(0:mgalx+1),
     .             ten22wk(0:mgalx+1),ten33wk(0:mgalx+1),
     .             ten12wk(0:mgalx+1),ten13wk(0:mgalx+1),
     .             ten23wk(0:mgalx+1),rhstwk(0:mgalx+1),
     .             up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .             up3wk8(0:mgalx+1),
     .             ten11wk8(0:mgalx+1),ten22wk8(0:mgalx+1),
     .             ten33wk8(0:mgalx+1),
     .             ten12wk8(0:mgalx+1),ten13wk8(0:mgalx+1),
     .             ten23wk8(0:mgalx+1),rhstwk8(0:mgalx+1),
     .             tmp1wk8(0:mgalx+1),tmp2wk8(0:mgalx+1),
     .             tmp3wk8(0:mgalx+1)


      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
      real*8  cflx,cfly,cflz,hxalp,hzbet,hyy,cfl0,reigmx1,aa
      real*4  up1(0:mx-1,lb:le),up2(0:mx-1,lb:le), 
     .        up3(0:mx-1,lb:le),ten11(0:mx-1,lb:le),
     .        ten22(0:mx-1,lb:le),ten33(0:mx-1,lb:le),
     .        ten12(0:mx-1,lb:le),ten13(0:mx-1,lb:le),
     .        ten23(0:mx-1,lb:le),
     .        rhst (0:mx-1,lb:le)

      real*8 duma,dumb,dumc,duma1,dumb1,dumc1
      real*8 pi
      !aaf adding Peclet
      real*8 peclet,sigma

     

c     ------------------------- Indices -----------------------

      integer vec(2,3),vecrc(2,3),vecy(2)         ! i,j
      integer resulm(2,3)

!      ----------------------- Programa ------------------------      
!initialization
      !aaf Peclet
      peclet=re     !for debuggin consider peclet=reynodls
      sigma=0.7D0
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
      !max of 1/Dy
      hyy=hy(1)
      do j=2,my
         hyy=min(hyy,hy(j))
      enddo
      hyy=1./hyy
! computing viscous Dt:
! computing viscous Dt inverse:
! Dt_visc=1/(2 nu)*Dy_min^2
      hyy=2/Re*(max(hyy,max(hxalp,hzbet)))**2
!      hyy=2/Re*(1./min(minval(hy),1./hxalp,1./hzbet))**2
!      hyy=max(hyy,max(hxalp,hzbet))
!      hyy=2./Re*hyy**2 !viscous time inverse

c   -------------------------------------------------

c        Move everything to pPP, line by line 

      duma = 0d0
      dumb = 0d0
      dumc = 0d0
      do 10 j = lb,le
       
c     copy lines...
!         if (myid.eq.0) then
!            duma = duma-MPI_WTIME()
!         endif
         do i=0,mx-1
            up1wk(i) = up1(i,j) !rhou
         enddo
         do i=0,mx-1           
            up2wk(i) = up2(i,j) !rhov
         enddo
         do i=0,mx-1
            up3wk(i) = up3(i,j) !rhow
         enddo
!!         !temperature
!         do i=0,mx-1
!            ten22wk(i) = ten22(i,j) !lap(T)
!         enddo
!         do i=0,mx-1
!            ten12wk(i) = ten12(i,j) !dT/dx
!         enddo
!         do i=0,mx-1
!            ten13wk(i) = ten13(i,j) !dT/dy
!         enddo
!         do i=0,mx-1
!            ten23wk(i) = ten23(i,j) !dT/dz
!         enddo
!         do i=0,mx-1
!            rhstwk(i) = rhst(i,j)   !T
!         enddo
!        !registering the mpi_wtime...
         if (myid.eq.0) then
            duma = duma+MPI_WTIME()
         endif

! convert lines in fourier to lines in physical for x
! upiwk8 means real 8 format
         call fourx(up1wk,up1wk8,1)    !rhou
         call fourx(up2wk,up2wk8,1)    !rhov
         call fourx(up3wk,up3wk8,1)    !rhow
!         !temperature
!         call fourx(ten22wk,ten22wk8,1) ! lap(T)
!         call fourx(rhstwk ,rhstwk8 ,1) ! T
!         call fourx(ten12wk,ten12wk8,1) ! dT/dx
!         call fourx(ten13wk,ten13wk8,1) ! dT/dy
!         call fourx(ten23wk,ten23wk8,1) ! dT/dz
! ----------------------------------------------------------------------!
!                                                                       !
!         computing in physical domain ( line by line)                  !
!                                                                       !
!          up1wk = H3 = u.omega_2 - v.omega_1 (F-F-P)                   !
! ----------------------------------------------------------------------!
!         !Calculate u,v,w
!         do i=0,mgalx-1
!            !tmp1wk8(i)= up1wk8(i)*rhstwk8(i) !u for wk 
!!debug rho=1
!            tmp1wk8(i)= up1wk8(i) !u for wk
!            up1wk(i)  = up1wk8(i)           !u for output
!
!            !tmp2wk8(i)= up2wk8(i)*rhstwk8(i) !v
!            tmp2wk8(i)= up2wk8(i) !v
!            up2wk(i)  = up2wk8(i)           !v for output
!
!            !tmp3wk8(i)= up3wk8(i)*rhstwk8(i) !w for wk
!            tmp3wk8(i)= up3wk8(i) !w for wk
!            up3wk(i)  = up3wk8(i)           !w for output
!         enddo

!         !Calculate rhs(T)
!         do i=0,mgalx-1
!            rhstwk(i)=-tmp1wk8(i)*ten12wk8(i)-tmp2wk8(i)*ten13wk8(i)-
!     .      tmp3wk8(i)*ten23wk8(i)+ !-u_i*dT/dxi
!     .      sigma/peclet*(rhstwk8(i)**sigma)*(ten12wk8(i)**2+
!     .      ten13wk8(i)**2+ten23wk8(i)**2) + !(1/Pe)s*T^s*(dT/dxi dT/dxi)
!     .      (1.0D0/peclet)*rhstwk8(i)**(sigma+1)*ten22wk8(i)
!         enddo         

!         !Calculate rhouu(ten11),rhovv(ten22),rhoww(ten33)
!         !rhou*u_i*u_i=rhou_i*u_i
         do i=0,mgalx-1
            ten11wk(i)=up1wk8(i)*up1wk8(i) !rhouu ready for output
            ten22wk(i)=up2wk8(i)*up2wk8(i) !rhovv ready for output
            ten33wk(i)=up3wk8(i)*up3wk8(i) !rhoww ready for output
         enddo
         
!         !Calculate the cross products rhou_i*u_j
         do i=0,mgalx-1
            ten12wk(i)=up1wk8(i)*up2wk8(i) !rhouv ready for output
            ten13wk(i)=up1wk8(i)*up3wk8(i) !rhouw ready for output
            ten23wk(i)=up3wk8(i)*up2wk8(i) !rhovw ready for ouput
         enddo
!ACHTUNG
!Cheking---write uu vv ww and check if its right
!
! now lines of mgalx, each mgalz one height of y.
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
!debuggin 
!            if (jj.eq.my) then
!            ! write(*,*) "u phys","myid=",myid
!            ! write(*,*) (up1wk8(i),i=0,mgalx-1)
!             write(*,*) "u*u","myid=",myid
!             write(*,*) (up1wk8(i)*up1wk8(i),i=0,1)
!             write(*,*) "uu","myid=",myid
!             write(*,*) (ten11wk8(i),i=0,1)
!             write(*,*) "v phys","myid=",myid
!             write(*,*) (up2wk8(i),i=0,1)
!             write(*,*) "w phys","myid=",myid
!             write(*,*) (up3wk8(i),i=0,1)
!            endif
            do i=0,mgalx-1
               aa = up2wk8(i)
               uuv(jj) = uuv(jj) +aa*up1wk8(i)**2
               wwv(jj) = wwv(jj) +aa*up3wk8(i)**2
               vvv(jj) = vvv(jj) +aa**3
            enddo
         endif

c            !      estimates maximum time step     !
! still within loop for each line
! now tmpiwk8 keeps the velocities
         if (rkstep.eq.1.and.icfl.eq.1) then
            cfl0=0.
            do i = 0,mgalx-1
               cflx = max(cflx,abs(up1wk8(i)))
               cfl0 = max(cfl0,abs(up2wk8(i)))
               cflz = max(cflz,abs(up3wk8(i)))
!cflx and cflz have same hx,hz for all domain
!but cfly needs to account for hy. 
            enddo
            cfly=max(cfly,cfl0/hy(jj))
         endif

c--------------------- back to F-P-P  

!         call fourx(rhstwk ,up1wk8,-1)  ! RHS(T)
         call fourx(up2wk  ,up1wk8,-1)  ! v
         call fourx(up3wk  ,up1wk8,-1)  ! w
         call fourx(ten11wk,up1wk8,-1)  ! rhouu
         call fourx(ten22wk,up1wk8,-1)  ! rhovv
         call fourx(ten33wk,up1wk8,-1)  ! rhoww
         call fourx(ten12wk,up1wk8,-1)  ! rhouv
         call fourx(ten13wk,up1wk8,-1)  ! rhouw
         call fourx(ten23wk,up1wk8,-1)  ! rhovw
         call fourx(up1wk  ,up1wk8,-1)  ! u
c --------    back to lines. We throw away the high modes (Dealiasing)

         if (myid.eq.0) then
            dumc = dumc-MPI_WTIME()
         endif
!         do i=0,mx-1            
!            rhst(i,j) = rhstwk(i) !RHS(T)
!         enddo
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




10    continue

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
            write(*,*) 'Dt_visc = ',CFL/hyy
         endif
!compute Deltat and dtr
         Deltat=CFL/reigmx1
         dtr=Re/Deltat


         if (reigmx1.lt.1e-1) then
            write(*,*) 'UGGG',myid,ihist,lb,le,pb,pe
         endif

      endif

      endsubroutine
      
      
      
!c----------------------------------------------------------------------
!c    Advances the 00 modes (linear term semi-implicit, order 3)
!c----------------------------------------------------------------------
!
!      subroutine rk00(f,rhs,nl,rkn1,dalre,dtgamma,dalbe,
!     .                dtri,dtxi,rkstep,flag)
!      implicit none
!      include "ctes3D"
!
!c     ----------------------- IN & OUT -----------------
!
!      real*8  f(my),rhs(my),nl(my)
!      real*8  rkn1,dalbe,dtri,dtxi,dalre,dtgamma
!      integer rkstep,flag
!c     ----------------------- work --------------------
!
!      integer i,j
!      real(8),dimension(my)   ::  fwkh1,fwkh2,dfwk,dfwkh1,dfwkh2
!      real(8)                 ::  coef(2),AA(2,2)
!
!c ---------------------------  commons -----------------
!      real*8 prem3,dt21,dt22
!      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
!      save   /d2y/
!
!      real*8 prem1,dt12
!      common /cfdiff/ prem1(7,my),dt12(7,my)
!      save   /cfdiff/
!
!      real*8  Re,alp,bet,a0
!      real*8  y,hy,fmap,trp,mss
!      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
!      save   /fis/
!
!      real*8 wk1,fwk,df
!      common /worksp/ wk1(5,my),fwk(my),df(my)
!      save   /worksp/
!!--------- ----------- BC functions----------------------------
!      real(8) :: Cvtop,Cdyvtop,BC_top
!      real(8) :: Cvbot,Cdyvbot,BC_bot
!
!      real*8 temp
!! ----------------------------------------------------------
!
!
!!u00 -> flag=1
!!w00 -> flag=2
!
!      if (flag.eq.1) then  !u00
!!BC constants
!!     Cvtop*u00(my)+Cdyvtop*u00(my)=BC_top
!!     Cvbot*u00(1)+Cdyvbot*u00(1)=BC_bot
!!     
!!     TOP
!!    -----------
!!Dirichlet
!      Cvtop=1d0
!      Cdyvtop=0d0
!      BC_top=1d0
!!Von Neuman
!!      Cvtop=0d0
!!      Cdyvtop=1d0
!!      BC_top=0d0
!
!!     BOT
!!     ---------
!!Dirichlet
!      Cvbot=1d0
!      Cdyvbot=0d0
!      BC_bot=-1d0
!!Von Neuman
!!      Cvbot=0d0
!!     Cdyvbot=1d0
!!      BC_bot=0d0
!
!      elseif (flag.eq.2) then   !w00
!!BC constants
!!     Cvtop*u00(my)+Vdyvtop*u00(my)=BC_top
!!     Cvbot*u00(1)+Vdyvbot*u00(1)=BC_bot
!!     
!!     TOP
!!    -----------
!        Cvtop=1d0
!        Cdyvtop=0d0
!        BC_top=0d0
!
!!     BOT
!!     ---------
!        Cvbot=1d0
!        Cdyvbot=0d0
!        BC_bot=0d0
!      endif 
!!-----------------------------------------------__!
!
!      do j=1,my      !--- to real*8, fitf
!         df(j)=f(j)
!      enddo
!!aaf  First RK-step we, dont know rhs yet.
!      if (rkstep.eq.1) then
!!         call deryyr(df,nl,rkn1,dalre,dtgamma)
!        fwk(1) =  dt22(3,1)*df(1) + dt22(4,1)*df(2) + 
!     $            dt22(5,1)*df(3)
!        fwk(2) =  dt22(2,2)*df(1) + dt22(3,2)*df(2) + 
!     $            dt22(4,2)*df(3) + dt22(5,2)*df(4)
!        do j=3,my-2
!           temp = 0d0
!           do i =1,5
!              temp = temp + dt22(i,j)*df(i+j-3)
!           enddo
!           fwk(j) = temp
!        enddo
!        fwk(my1)=dt22(1,my1)*df(my-3) + dt22(2,my1)*df(my-2)+ 
!     $            dt22(3,my1)*df(my-1) + dt22(4,my1)*df(my  )
!        fwk(my)=  dt22(1,my  )*df(my-2) + dt22(2,my  )*df(my-1)+ 
!     $            dt22(3,my  )*df(my)
!! solve second derivative of u(rkstep) -> fwk
!! fwk=d2y2(u00)
!        call banbks(prem3,my,fwk)
!
!! rhs for step 1
!        do j=1,my
!           df(j) = -rkn1*(df(j)+dalre*fwk(j)+dtgamma*nl(j))
!        enddo
!
!      else
!!aaf   This is done for all rksteps except ONE
!         do j=1,my
!            df(j) = -rkn1*(rhs(j)+dtgamma*nl(j))
!         enddo
!
!      endif
!
!!---------------- Calculation of new u(rkstep)------------------------!
!!
!      if (rkstep.ne.3) then
!         do j=1,my
!            rhs(j)=dtri*df(j)+dtxi*nl(j)
!         enddo
!      endif
!!Prepare wk1 matrix: wk1=dt22-rk dt21
!      wk1(1,1)  = 0d0
!      wk1(2,1)  = 0d0
!      wk1(3,1)  = 1d0
!      wk1(4,1)  = 0d0  
!      wk1(5,1)  = 0d0  
!      wk1(1,my) = 0d0
!      wk1(2,my) = 0d0
!      wk1(3,my) = 1d0
!      wk1(4,my) = 0d0
!      wk1(5,my) = 0d0
!
!      do j=2,my-1
!        do i=1,5
!          wk1(i,j)=dt22(i,j)-rkn1*dt21(i,j)
!        enddo
!      enddo
!
!!numerical recipes codebase--
!!given nxn band diagonal matrix wk1, construct LU decomposition
!!used in conjuction with banbks to solve band-diagonal sets of equations
!      call bandec(wk1,my)
!! ----------------------------------------------------------------------!
!      call calcrhslap(df,fwk)  !prepares RHS multiplying by dt21
!!solves band diagonal linear equations Ax=b --returning b (as fwk)
!!            A      b
!      call banbks(wk1,my,fwk)
!
!! Now fwk has u(rkstep) particular
!! Now solving the homogeneus solutions
!      fwkh1(1) = 1d0
!      fwkh2(1) = 0d0
!      do j=2,my-1
!         fwkh1(j)=0d0
!         fwkh2(j)=0d0
!      enddo  
!      fwkh1(my) = 0d0
!      fwkh2(my) = 1d0
!! solving for uh1, uh2
!      call banbks(wk1,my,fwkh1)
!      call banbks(wk1,my,fwkh2) 
!        
!!FLAG
!!      write(*,*) "fwk(1),fwk(2)",fwk(1),fwk(2)
!!.................................................
!!     Caculating the derivative:
!!     multiplying by dt12, preparing the derivative calculation
!      call calcrhsd1(fwk,dfwk)     
!      call calcrhsd1(fwkh1,dfwkh1) 
!      call calcrhsd1(fwkh2,dfwkh2) 
!
!!FLAG
!!      write(*,*) "dfwk(1),dfwk(2)",dfwk(1),dfwk(2)
!!.................................................
!
!!     derivative in y
!      call banbks7(prem1,my,dfwk  )!d(up)/dy
!      call banbks7(prem1,my,dfwkh1)!d(uh1)/dy
!      call banbks7(prem1,my,dfwkh2)!d(uh2)/dy
!
!
!!     Building the A matrix
!!     part multiplied by coef(1) on the left for BC
!!     at bottom
!      AA(1,1)=Cvbot+Cdyvbot*dfwkh1(1)*fmap(1)
!!     part multiplied by coef(2) on the left for BC
!!     at bottom
!      AA(1,2)=Cdyvbot*dfwkh2(1)*fmap(1)
!!     same for top BC
!      AA(2,1)=Cdyvtop*dfwkh1(my)*fmap(my)
!      AA(2,2)=Cvtop + Cdyvtop*dfwkh2(my)*fmap(my)
!!     RHS of linear sistem A · coef' = B
!      coef(1)=BC_bot-Cdyvbot*dfwk( 1)*fmap(1 )
!      coef(2)=BC_top-Cdyvtop*dfwk(my)*fmap(my)
!!FLAG
!!      write(*,*) "dfwk(1),dfwk(2)",dfwk(1),dfwk(2)
!!      write(*,*) "coef(1),coef(2)",coef(1),coef(2)
!!...............................................      
!!     solve system
!      call gaussj(AA,2,2,coef,2,1)
!!     sum all solutions for u(rkstep) fwk
!
!      do j=1,my
!         fwk(j)=fwk(j)+ coef(1)*fwkh1(j) + coef(2)*fwkh2(j)
!      enddo
!! new u(rkstep) has been calculated (fwk)
!! preparation of rhs for next step
!      if (rkstep.ne.3) then  
!        do j=1,my
!            rhs(j)=rhs(j)+dalbe*fwk(j)
!        enddo
!      else
!         do j=1,my
!            rhs(j)=fwk(j)
!         enddo
!      endif
!
!      do j=1,my
!         f(j)=fwk(j)
!      enddo
!
!
!      endsubroutine
!
!
c----------------------------------------------------------------------
c    Advances the 00 modes (linear term explicit, order 3)
c----------------------------------------------------------------------

      subroutine rk00exp(f,rhs,dtgamma)
      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real*8  f(my),rhs(my)
      real*8  dtgamma
c     ----------------------- work --------------------

      integer i,j

c ---------------------------  commons -----------------

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/


! ----------------------------------------------------------
          do j=1,my 
            f(j) = f(j)+dtgamma*rhs(j) 
            !new u00/w00:
          enddo
! not evolving the BC

      endsubroutine
!




!*********************************************************************!
!       subroutine   cstati.
!        compute the statistics in FPF
!
!*********************************************************************!

      subroutine cstatis(u1c,u2c,u3c,o1c,o2c,o3c,
     .                   chwkc,uner,spplane) 
!     .                   chwkc,spaux,uner,spplane) 
      implicit none

      include "mpif.h"
      include "ctes3D"

c    ------------------------ Commons ----------------------------

      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 


      integer nacumsp,jsptot
      common /spectra/   nacumsp,jsptot(2*nspec+1)
      save   /spectra/
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/
      

      integer nimag,nstep,nhist,ihist,icfl,ncfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
      save   /timacc/

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/
     
      real*4 uner(9)
      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/

c    ----------------------- Variables ---------------------------

      complex*8 u1c(my,0:mz1),u2c(my,0:mz1),u3c(my,0:mz1),
     .          o1c(my,0:mz1),o2c(my,0:mz1),o3c(my,0:mz1),
     .        chwkc(my,0:mz1)
!      real*4 spaux(0:nz1,2*nspec+1,8)  
      integer j,k,iy,kk,spplane
      real*4 hyy

      real*8 aa
      complex*16 cc

c     --------------------- Programa -----------------------------

c     -----  spectra



!      if (spplane.eq.1)  then          ! i = 0
!         do iy=1,2*nspec+1
!            do kk = 0,mz1
!               k = icx(kk)
!               spaux(k,iy,1) = spaux(k,iy,1)+u1c(jsptot(iy),kk)*
!     &                        conjg(u1c(jsptot(iy),kk))
!               spaux(k,iy,2) = spaux(k,iy,2)+u2c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk))
!               spaux(k,iy,3) = spaux(k,iy,3)+u3c(jsptot(iy),kk)*
!     &                        conjg(u3c(jsptot(iy),kk))
!               spaux(k,iy,4) = spaux(k,iy,4)+real(u1c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk)))
!               spaux(k,iy,5) = spaux(k,iy,5)+imag(u1c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk)))
!               spaux(k,iy,6) = spaux(k,iy,6)+o1c(jsptot(iy),kk)*
!     &                        conjg(o1c(jsptot(iy),kk))
!               spaux(k,iy,7) = spaux(k,iy,7)+o2c(jsptot(iy),kk)*
!     &                        conjg(o2c(jsptot(iy),kk))
!               spaux(k,iy,8) = spaux(k,iy,8)+o3c(jsptot(iy),kk)*
!     &                        conjg(o3c(jsptot(iy),kk))
!            enddo
!         enddo
!      else
!         do iy=1,2*nspec+1
!            do kk = 0,mz1
!               k = icx(kk)
!               spaux(k,iy,1) = spaux(k,iy,1)+2.*u1c(jsptot(iy),kk)*
!     &                        conjg(u1c(jsptot(iy),kk))
!               spaux(k,iy,2) = spaux(k,iy,2)+2.*u2c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk))
!               spaux(k,iy,3) = spaux(k,iy,3)+2.*u3c(jsptot(iy),kk)*
!     &                        conjg(u3c(jsptot(iy),kk))
!               spaux(k,iy,4) = spaux(k,iy,4)+2.*real(u1c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk)))
!               spaux(k,iy,5) = spaux(k,iy,5)+2.*imag(u1c(jsptot(iy),kk)*
!     &                        conjg(u2c(jsptot(iy),kk)))     
!               spaux(k,iy,6) = spaux(k,iy,6)+2.*o1c(jsptot(iy),kk)*
!     &                        conjg(o1c(jsptot(iy),kk))
!               spaux(k,iy,7) = spaux(k,iy,7)+2.*o2c(jsptot(iy),kk)*
!     &                        conjg(o2c(jsptot(iy),kk))
!               spaux(k,iy,8) = spaux(k,iy,8)+2.*o3c(jsptot(iy),kk)*
!     &                        conjg(o3c(jsptot(iy),kk))
!            enddo
!         enddo
!      endif
!             
c                    ----- intensities ------------
c                    ----- & dissipation ------------

      call deryr2(u2c(1,0),chwkc(1,0),mz)!,chwkc)!,chwkr)
      if (spplane.eq.1) then             ! i=1
         do k=0,mz1                     
            do j=1,my          
c              intensities ----------------
                
               ener(1) = u1c(j,k)*conjg(u1c(j,k))
               up(j) = up(j) + ener(1)
        
               ener(2) = u2c(j,k)*conjg(u2c(j,k))
               vp(j) = vp(j) + ener(2)
               
               ener(3) = u3c(j,k)*conjg(u3c(j,k))
               wp(j) = wp(j) + ener(3)

               ener(4) = real(u1c(j,k)*conjg(u2c(j,k)))
               uvr(j)= uvr(j) + ener(4)

               ener(5) = real(u1c(j,k)*conjg(u3c(j,k)))
               uwr(j)= uwr(j) + ener(5)

               ener(6) = real(u3c(j,k)*conjg(u2c(j,k)))
               vwr(j)= vwr(j) + ener(6)

               ener(7) = o1c(j,k)*conjg(o1c(j,k))
               w1p(j)= w1p(j) + ener(7)

               ener(8) = o2c(j,k)*conjg(o2c(j,k))
               w2p(j)= w2p(j) + ener(8)

               ener(9) = o3c(j,k)*conjg(o3c(j,k))
               w3p(j)= w3p(j) + ener(9)

c              dissipation  ----------------

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
            enddo
         enddo
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

c           dissipation  ----------------

               aa = ( alp2(spplane-1) + bet2(k) ) *
     &              ( u1c(j,k)*conjg(u1c(j,k)) +
     &                u2c(j,k)*conjg(u2c(j,k)) +
     &                u3c(j,k)*conjg(u3c(j,k)) ) +
     &              chwkc(j,k)*conjg(chwkc(j,k))

               cc = o1c(j,k) + xbet(k)*u2c(j,k)
               aa = aa + cc*conjg(cc)
               cc = o3c(j,k) - xalp(spplane-1)*u2c(j,k)
               aa = aa + cc*conjg(cc)

               ep(j) = ep(j) + 2.*aa

cc ------------ add this plane energy
  
               hyy = hy(j)
               do kk = 1,9
                  uner(kk) = uner(kk) + ener(kk)*hyy
               enddo

             enddo   ! End y loop
          enddo
      endif

c ------------ means
      if (spplane.eq.1) then 
         do j=1,my
            um(j) = um(j)+u1c(j,0)
            vm(j) = vm(j)+u2c(j,0)
            wm(j) = wm(j)+u3c(j,0)
            w1m(j)= w1m(j)+o1c(j,0)
            w2m(j)= w2m(j)+o2c(j,0)
            w3m(j)= w3m(j)+o3c(j,0)
         enddo
      endif

      if (spplane.eq.1) then
         Wx0a=Wx0a+o1c(1,0)
         Wz0a=Wz0a+o3c(1,0)
      endif

c------------- compute vorticity & Re stress  at walls

       if (spplane.eq.1) then
          Wx0 = o1c(1,0)
          Wz0 = o3c(1,0)
          WxL = o1c(my,0)
          WzL = o3c(my,0)
       endif

      endsubroutine


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
      implicit none
      include "ctes3D"
      integer m,i,j,k,xpl

      real*4 u(2,my,m),rhs(2,my,m),tmp
      real*8 dtgamma

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/


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
      implicit none
      include "ctes3D"
      integer m,i,j,k,xpl

      real*4 u(2,my,m),du(2,my,m),tmp
      real*8 wk1(my),wk2(my)
      real*8 coef

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/


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
!!      drho    =        T            Drho                         !
!!       tnext   =      T(n+1)          -                           !
!!******************************************************************!
!      subroutine tauij(ten11,ten22,ten33,
!     .           ten12,ten13,ten23,drho,tnext,chwk,myid,rkstep)
!
!      implicit none
!
!      include "mpif.h"
!      include "ctes3D"
!
!c  -------------------------   commons ------------------------- 
!
!      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
!      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
!     .               pbeg(0:numerop-1),pend(0:numerop-1),
!     .               plcalc(nplanes),plsav(nplanes),
!     .               pb,pe,lb,le,mmp,mml,procs
!      save /point/ 
!
!      
!      real*4 Deltat,CFL,time,dtr
!      common /tem/ Deltat,CFL,time,dtr
!      save /tem/
!
!      real*8  Re,alp,bet,a0
!      real*8  y,hy,fmap,trp,mss
!      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
!     $            trp(0:my-1),mss(0:my-1)
!      save   /fis/
!      
!      integer nimag,nstep,nhist,ihist,icfl,ncfl
!      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
!      save   /timacc/
!
!      integer iax,icx
!      real*4 alp2,bet2,ralp,rbet
!      complex*8    xalp, xbet
!      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
!     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
!      save  /wave/
!
!
!      real*8 commtimer,transtimer,totaltimer
!      common /timers/ commtimer, transtimer, totaltimer
!      save   /timers/
!      
!      integer nanerror 
!      common /cnan/ nanerror
!      save   /cnan/  
!
!c     ---------------------- Variables ----------------------
!      
!      complex*8 ten11 (0:my-1,0:mgalz-1,pb:pe),
!     .          ten22 (0:my-1,0:mgalz-1,pb:pe),
!     .          ten33 (0:my-1,0:mgalz-1,pb:pe),
!     .          ten12 (0:my-1,0:mgalz-1,pb:pe),
!     .          ten13 (0:my-1,0:mgalz-1,pb:pe),
!     .          ten23 (0:my-1,0:mgalz-1,pb:pe),
!     .          tnext (0:my-1,0:mgalz-1,pb:pe),
!     .          drho  (0:my-1,0:mgalz-1,pb:pe)
!               
!      
!!      real*8 rf0u(0:my-1),rf0w(0:my-1),dum
!
!      real*4 chwk(*),uner(9),aa
!!      real*4 sp(0:nz1,1:2*nspec+1,8,pb:pe)
!      
!      integer myid,rkstep,i,j,k,kk,ierr
!      
!      integer MPI_GROUP_WORLD
!      integer MPI_GROUP_CALC,MPI_COMM_CALC
!      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
!      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
!     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
!      save /MPI_GROUPS/
!     
!      
!!  local transpose yz to zy, and adding zeros between positive
!!  and negative wavenumbers.
!      call localyz2zy(ten11,ten11,chwk) !S11
!      call localyz2zy(ten22,ten22,chwk) !S22
!      call localyz2zy(ten33,ten33,chwk) !S33
!      call localyz2zy(ten12,ten12,chwk) !S12         
!      call localyz2zy(ten13,ten13,chwk) !S13  
!      call localyz2zy(ten23,ten23,chwk) !S23        
!!      call localyz2zy(drho ,drho ,chwk) !T      
!!      call localyz2zy(tnext,tnext,chwk) !T(n+1)      
!       
!! inverse transform (fou --> fis @z),does the compacting aswell
!      do i = pb,pe
!         call fourz(ten11(0,0,i),1)    
!         call fourz(ten22(0,0,i),1)    
!         call fourz(ten33(0,0,i),1)    
!         call fourz(ten12(0,0,i),1)  
!         call fourz(ten13(0,0,i),1)  
!         call fourz(ten23(0,0,i),1)  
!!         call fourz(drho (0,0,i),1) 
!!         call fourz(tnext(0,0,i),1) 
!      enddo
!
!      nanerror = 0
!
!!  change plane to line(1:mx,lb:le)
!!  and distribute package of lines to processors
!      call chpl2ln(ten11,ten11,chwk,myid)
!      call chpl2ln(ten22,ten22,chwk,myid)
!      call chpl2ln(ten33,ten33,chwk,myid)
!      call chpl2ln(ten12,ten12,chwk,myid)
!      call chpl2ln(ten13,ten13,chwk,myid)
!      call chpl2ln(ten23,ten23,chwk,myid)
!!      call chpl2ln(drho ,drho ,chwk,myid)
!!      call chpl2ln(tnext,tnext,chwk,myid)
!
!
!!==========================================================================!
!! before this call ten11,ten22,ten33,ten12,ten13, ten23,drho and tnext
!! are still lines in fourier space x (phys z).
!      call tensor(ten11,ten22,ten33,ten12,ten13,ten23,
!     .           drho,tnext,myid,rkstep)
!!==========================================================================!
!! they are still lines in fourier x after hvect
!
!c  ---------- back to the yz planes and fft 
!!      call chln2pl(drho,drho,chwk,myid)   !drho 
!      call chln2pl(ten11,ten11,chwk,myid) !tau11
!      call chln2pl(ten22,ten22,chwk,myid) !tau22
!      call chln2pl(ten33,ten33,chwk,myid) !tau33
!      call chln2pl(ten12,ten12,chwk,myid) !tau12
!      call chln2pl(ten13,ten13,chwk,myid) !tau13
!      call chln2pl(ten23,ten23,chwk,myid) !tau23
!
!!  convert to fourier Z now.
!      do i=pb,pe
!!         call fourz(drho (0,0,i),-1)  !drho
!         call fourz(ten11(0,0,i),-1)  !tau11
!         call fourz(ten22(0,0,i),-1)  !tau22
!         call fourz(ten33(0,0,i),-1)  !tau33
!         call fourz(ten12(0,0,i),-1)  !tau12
!         call fourz(ten13(0,0,i),-1)  !tau13
!         call fourz(ten23(0,0,i),-1)  !tau23
!         
!      enddo
!!transpose back
!!      call localzy2yz(drho , drho,chwk) 
!      call localzy2yz(ten11,ten11,chwk) 
!      call localzy2yz(ten22,ten22,chwk) 
!      call localzy2yz(ten33,ten33,chwk) 
!      call localzy2yz(ten12,ten12,chwk)   
!      call localzy2yz(ten13,ten13,chwk)   
!      call localzy2yz(ten23,ten23,chwk)   
!           
!
!      endsubroutine
!


!!*********************************************************************!
!!                                                                     !
!!         transform to physical and do operations                     !
!!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!!                                                                     !
!!                                                                     !
!!*********************************************************************!
!      subroutine tensor(ten11,ten22,ten33,ten12,ten13,ten23,
!     .                 drho,tnext,myid,rkstep)
!
!      implicit none
!
!      include "mpif.h"
!      include "ctes3D"
!
!c --------------------------- Commons --------------------------
!
!
!      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
!      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
!     .               pbeg(0:numerop-1),pend(0:numerop-1),
!     .               plcalc(nplanes),plsav(nplanes),
!     .               pb,pe,lb,le,mmp,mml,procs
!      save /point/ 
!
!
!      real*4 Deltat,CFL,time,dtr
!      common /tem/ Deltat,CFL,time,dtr
!      save   /tem/
!      
!      real*8  Re,alp,bet,a0
!      real*8  y,hy,fmap,trp,mss
!      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
!     $            trp(0:my-1),mss(0:my-1)
!      save   /fis/
!
!    
!      integer MPI_GROUP_WORLD
!      integer MPI_GROUP_CALC,MPI_COMM_CALC
!      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
!      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
!     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
!      save /MPI_GROUPS/
!     
!
!c     --------------------------- Variables -------------------      
!
!      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
!      real*8  cflx,cfly,cflz,hxalp,hzbet,hyy,cfl0,reigmx1,aa
!      real*4  tnext(0:mx-1,lb:le),
!     .        ten11(0:mx-1,lb:le),
!     .        ten22(0:mx-1,lb:le),ten33(0:mx-1,lb:le),
!     .        ten12(0:mx-1,lb:le),ten13(0:mx-1,lb:le),
!     .        ten23(0:mx-1,lb:le),
!     .        drho (0:mx-1,lb:le)
!      real*8  sigma,const1
!      
!
!      real*4 up1wk,up2wk,up3wk,ten11wk,
!     .       ten22wk, ten33wk,ten12wk,ten13wk,ten23wk,
!     .       rhstwk
!      real*8 up1wk8,up2wk8,up3wk8,
!     .       ten11wk8,ten22wk8,ten33wk8,ten12wk8,ten13wk8,
!     .       ten23wk8,rhstwk8,tmp1wk8,tmp2wk8,tmp3wk8
!      common /wkhvect/ up1wk (0:mgalx+1),up2wk(0:mgalx+1),
!     .             up3wk (0:mgalx+1),ten11wk(0:mgalx+1),
!     .             ten22wk(0:mgalx+1),ten33wk(0:mgalx+1),
!     .             ten12wk(0:mgalx+1),ten13wk(0:mgalx+1),
!     .             ten23wk(0:mgalx+1),
!     .             up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
!     .             up3wk8(0:mgalx+1),
!     .             ten11wk8(0:mgalx+1),ten22wk8(0:mgalx+1),
!     .             ten33wk8(0:mgalx+1),
!     .             ten12wk8(0:mgalx+1),ten13wk8(0:mgalx+1),
!     .             ten23wk8(0:mgalx+1),rhstwk8(0:mgalx+1),
!     .             tmp1wk8(0:mgalx+1),tmp2wk8(0:mgalx+1),
!     .             tmp3wk8(0:mgalx+1)
!
!
!
!
!
!
!!      ----------------------- Programa ------------------------      
!      sigma=0.7D0
!      const1=2D0/3D0
!
!c        Move everything to pPP, line by line 
!
!      do 10 j = lb,le
!
!!         do i=0,mx-1
!!            up1wk(i) = drho(i,j)!T(n)
!!         enddo
!!         do i=0,mx-1           
!!            up2wk(i) = tnext(i,j) !T(n+1)
!!         enddo
!         do i=0,mx-1
!            ten11wk(i) = ten11(i,j) !S11
!         enddo
!         do i=0,mx-1
!            ten22wk(i) = ten22(i,j) !S22
!         enddo
!         do i=0,mx-1
!            ten33wk(i) = ten33(i,j) !S33
!         enddo
!         do i=0,mx-1
!            ten12wk(i) = ten12(i,j) !S12
!         enddo
!         do i=0,mx-1
!            ten13wk(i) = ten13(i,j) !S13
!         enddo
!         do i=0,mx-1
!            ten23wk(i) = ten23(i,j) !S23
!         enddo
!
!! convert lines in fourier to lines in physical for x
!! upiwk8 means real 8 format
!!         call fourx(up1wk,up1wk8,1)    !T
!!         call fourx(up2wk,up2wk8,1)    !T(n+1)
!         call fourx(ten11wk,ten11wk8,1) ! S11
!         call fourx(ten22wk,ten22wk8,1) ! S22
!         call fourx(ten33wk,ten33wk8,1) ! S33
!         call fourx(ten12wk,ten12wk8,1) ! S12
!         call fourx(ten13wk,ten13wk8,1) ! S13
!         call fourx(ten23wk,ten23wk8,1) ! S23
!! ----------------------------------------------------------------------!
!!                                                                       !
!!         computing in physical domain ( line by line)                  !
!!                                                                       !
!! ----------------------------------------------------------------------!
!         !Calculate drho=1/T(n+1)-1/T
!         do i=0, mgalx-1
!!            up1wk(i)  = 1.0D0/up2wk8(i)-1.0D0/up1wk8(i)           !drho ready
!!            up2wk8(i) =const1*(ten11wk8(i)+ten22wk8(i)+ten33wk8(i)) !2/3Skk
!             up2wk8(i) =0d0 !incompressible !2/3Skk
!         enddo
!         !S11,S22,S33 like: tau11=mu*(2S11-2/3Skk)
!         !Sij like        : tauij=2*mu*Sij
!         do i=0,mgalx-1
!!            ten11wk(i)=up1wk8(i)**sigma*(2*ten11wk8(i)
!!     .                                   -up2wk8(i))
!!            ten22wk(i)=up1wk8(i)**sigma*(2*ten22wk8(i)
!!     .                                   -up2wk8(i))
!!            ten33wk(i)=up1wk8(i)**sigma*(2*ten33wk8(i)
!!     .                                   -up2wk8(i))
!!            ten12wk(i)=2*up1wk8(i)**sigma*ten12wk8(i)
!!            ten13wk(i)=2*up1wk8(i)**sigma*ten13wk8(i)
!!            ten23wk(i)=2*up1wk8(i)**sigma*ten23wk8(i)
!!
!            ten11wk(i)=2*ten11wk8(i)-up2wk8(i)
!            ten22wk(i)=2*ten22wk8(i)-up2wk8(i)
!            ten33wk(i)=2*ten33wk8(i)-up2wk8(i)
!            ten12wk(i)=2*ten12wk8(i)
!            ten13wk(i)=2*ten13wk8(i)
!            ten23wk(i)=2*ten23wk8(i)
! 
!
!         enddo
!
!!========================================================!
!!IMPORTANT POINT: all quantities are now PHYSICAL
!!========================================================!
!c--------------------- back to F-P-P  
!
!!         call fourx(up1wk  ,up1wk8,-1)  ! drho
!         call fourx(ten11wk,up1wk8,-1)  ! tau11
!         call fourx(ten22wk,up1wk8,-1)  ! tau22
!         call fourx(ten33wk,up1wk8,-1)  ! tau33
!         call fourx(ten12wk,up1wk8,-1)  ! tau12
!         call fourx(ten13wk,up1wk8,-1)  ! tau13
!         call fourx(ten23wk,up1wk8,-1)  ! tau23
!
!c --------    back to lines. We throw away the high modes (Dealiasing)
!
!!         do i=0,mx-1
!!            drho(i,j) = up1wk(i) !drho
!!         enddo
!         do i=0,mx-1            
!            ten11(i,j) = ten11wk(i)
!         enddo
!         do i=0,mx-1            
!            ten22(i,j) = ten22wk(i)
!         enddo
!         do i=0,mx-1            
!            ten33(i,j) = ten33wk(i)
!         enddo
!         do i=0,mx-1            
!            ten12(i,j) = ten12wk(i)
!         enddo
!         do i=0,mx-1            
!            ten13(i,j) = ten13wk(i)
!         enddo
!         do i=0,mx-1            
!            ten23(i,j) = ten23wk(i)
!         enddo
!
!10    continue
!
!
!      endsubroutine
      

