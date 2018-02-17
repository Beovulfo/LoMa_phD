!**********************************************************************!
!                                                                      !
!      Da un paso en el tiempo. Runge - Kutta  SMR                     !
!                                                                      !
!       Resuelve:    Gt  + Hg = 1/Re G"                                !
!                    V"t + Hv = 1/Re V""                               !
!                                                                      !
!   input:                                                             !
!     vor: vorticidad segun la direccion y (n)                         !
!     phi: laplaciana velocidad segun la direccion y  (n)              !
!     vorwk: copia de vor para el calc. de los term. no lineales       !
!     phiwk: copia de phi para el calc. de los term. no lineales       !
!     hg: Hg                                                           !
!     hv: Hv                                                           !
!     uwk: u                                                           !
!     dvordy: area de trabajo                                          !
!     chwk:   area de trabajo                                          !
!     u00,w00: modos 00                                                !
!     sp: Espectros                                                    !
!                                                                      !
!                                                                      !
!  output:                                                             !
!     vor: vorticidad segun la direccion y (n+1)                       !
!     phi: laplaciana velocidad segun la direccion y  (n+1)            !
!                                                                      !
!   updated jjs 07/01/01                                               !
!   in jik form jca                                                    !
!   to CFdiff by of                                                    !
!   RK SMR and plane-lines SHC                                         !
!                                                                      !
!   Warning : hvect subroutine in real*8                               !
!             00 modes in real*8                                       !
!                                                                      !
!    JJ for 00 mode;   may/15/05                                       !
!                                                                      !
!                                                                      !
!**********************************************************************!

      subroutine cross1(vor,
     .                  phi,
     .                  u00,w00,
     .                  rf0u,rf0w,u00wk,w00wk,
     .                  hv,
     .                  hg,
     .                  phiwk,
     .                  spwk,
     .                  vorwk,
     .                  dvordy,
     .                  uwk,
     .                  chwk,
     .                  sp,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

      !------------- Parameters for SMR R-K ------------------------!
      
      real*8     gama(3), alpha(4), beta(3), ibeta(3), xi(3)
      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
      parameter (alpha=(/ 29d0/96d0, -3d0/40d0, 1d0/6d0, 29d0/96d0/))
      parameter (beta =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /))
      parameter (ibeta=(/ 160d0/37d0, 24d0/5d0,   6d0     /))
      parameter (xi   =(/ -17d0/60d0, -5d0/12d0,  0d0     /))
   
 
!     ---------------------   commons ------------------------------- 

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .        ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,
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
      character*70 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                              filinp,filout,filstt
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

      complex*8 phi(0:my-1,0:mz1,pb:pe),vor(0:my-1,0:mz1,pb:pe)

      complex*8 hg    (0:my-1,0:mgalz-1,pb:pe),
     .          vorwk (0:my-1,0:mgalz-1,pb:pe),
     .          hv    (0:my-1,0:mgalz-1,pb:pe),
     .          phiwk (0:my-1,0:mgalz-1,pb:pe),
     .          uwk   (0:my-1,0:mgalz-1,pb:pe),
     .          dvordy(0:my-1,0:mgalz-1,pb:pe),
     .          chwk  (0:my-1,0:mgalz-1,pb:pe)

      real*8 u00(0:my-1),w00(0:my-1),rf0u(0:*),u00wk(0:*),
     .       rf0w(0:*),w00wk(0:*)

      real*4 sp  (0:nz1,1:2*nspec+1,8,pb:pe),
     .       spwk(0:nz1,1:  nspec+1,8,pb:pe)
     
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,kk
      
      real*8  rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire
      
      real*8  iter_time,write_time
      complex*8 temp1,temp2
      real*4  temporal(3),comun(3)
      
      real*4  reynota,H00u,H00w,dumu,dumw,cteu,ctew
      real*4  massu,massw,massu1,massu2,massw1,massw2,dum
       
      complex*8,allocatable :: tempa(:,:),tempb(:,:)
      
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
         tempa(0,1)=0d0
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
            
            call escru(vor,phi,u00,w00,sp,spwk,myid)

            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
            endif
         ENDIF

! ------------------- finished writing image --------------------------!

      if (myid.gt.numerop-1) goto 30    ! only for save procs

      if(irun.eq.0) then ! this is done only for the first step
         irun = 1

c          !   calcula la v a partir de phi !

         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               call Lapvdv(phi(0,k,i),hg(0,k,i),hv(0,k,i),rK)
            enddo
         enddo

        !     prepara phiwk,vorwk,u00wk,w00wk !

         if (myid.eq.0) then
            do j=0,my1
              u00wk(j)=u00(j)
              w00wk(j)=w00(j)
            enddo
            call deryr(u00wk,rf0u)
            call deryr(w00wk,rf0w)
         endif

         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  vorwk(j,k,i)=vor(j,k,i)
                  phiwk(j,k,i)=phi(j,k,i)
               enddo
            enddo
         enddo
         
      endif ! end special first step
      
               !  Runge-Kutta third order  !
                

      do 10 rkstep=1,3

c--------------  dvordy      : d (vorwk) / dy -- F-F
         do i=pb,pe
            call deryr2(vorwk(0,0,i),chwk(0,0,i),mz)
         enddo
         
c ---- Computation of velocity and vorticity 

! ----------------------------------------------------------------------!
!   at this point
!   vor:    r.h.s of vorticity
!   phi:    r.h.s of phi
!   vorwk:  vorticity
!   phiwk:  phi 
!   dvordy: d (vorwk) / dy 
!   hg:     v
!   hv:     dv/dy
!   uwk:    empty 
!   
! ----------------------------------------------------------------------!


c---------------------- computes  u,w, ome1, ome3

         do i=pb,pe
            do k=0,mz1
               temp1 = tempa(k,i)
               temp2 = tempb(k,i)
               do j=0,my-1
                  ! velocities

                  uwk(j,k,i)  = hv(j,k,i)*temp1-vorwk(j,k,i)*temp2

                  hv(j,k,i)   = hv(j,k,i)*temp2+vorwk(j,k,i)*temp1

                  ! vorticities
                  dvordy(j,k,i)= chwk(j,k,i)*temp2-phiwk(j,k,i)*temp1
     
                  phiwk(j,k,i) = chwk(j,k,i)*temp1+phiwk(j,k,i)*temp2

               enddo
            enddo
         enddo

! ----------------------------------------------------------------------!
!   at this point
!   vor:   r.h.s of vorticity
!   phi:   r.h.s of phi
!   vorwk: vor_y
!   phiwk: vor_x  
!   dvordy:vor_z  
!   uwk:   u   
!   hg:    v
!   hv:    w   
!
! ----------------------------------------------------------------------!

!---------- Only master computes 0  mode of vorticity

        if (myid.eq.0) then
           do j=0,my-1  
              dvordy(j,0,1) = -rf0u(j)
              phiwk (j,0,1) = rf0w (j)
              uwk   (j,0,1) = u00wk(j)
              hv    (j,0,1) = w00wk(j)
           enddo

        endif

        
        call hvhg(uwk,hg,hv,phiwk,vorwk,dvordy,
     .          rf0u,rf0w,chwk,sp,myid,rkstep)
     
     
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



! ----------------------------------------------------------------------!
!  at this point:  hv and hg are the outs
!  anything else has been destroyed 
! ----------------------------------------------------------------------!

c------- advances in time : nonlinear terms explicitly

c------- After computes CFL, we need the second derivative of 
c------- vor and phi

         rkn1    = Re*ibeta(rkstep)/deltat 
         dalre   = deltat*alpha(rkstep)*ire
         dtgamma = deltat*gama(rkstep)
         dalbe   = 1d0+alpha(rkstep+1)*ibeta(rkstep)
         dtri    = deltat*alpha(rkstep+1)*ire
         dtxi    = deltat*xi(rkstep)
       

         if (rkstep.eq.1.and.icfl.eq.1) then
            do i=pb,pe
              call deryyr2(vor(0,0,i),vorwk(0,0,i),hg(0,0,i),mz,i
     .                    ,rkn1,dalre,dtgamma)
              call deryyr2(phi(0,0,i),phiwk(0,0,i),hv(0,0,i),mz,i
     .                    ,rkn1,dalre,dtgamma)
            enddo
         else ! para todos los pasos con icfl 0 
            do i=pb,pe
               do k=0,mz1
                  do j=0,my-1
                     phiwk(j,k,i)=-rkn1*(phi(j,k,i)+dtgamma*hv(j,k,i))
                     vorwk(j,k,i)=-rkn1*(vor(j,k,i)+dtgamma*hg(j,k,i))
                  enddo
               enddo
            enddo
         endif
        
c          ! implicit time step for    !
c          ! viscous terms.            !


         if (rkstep.ne.3) then  
            do i=pb,pe
               rK = alp2(i-1)                                           
               call visc(vorwk(0,0,i),vor(0,0,i),                       
     .             phi(0,0,i),hg(0,0,i),hv(0,0,i),phiwk(0,0,i),         
     .             rK+rkn1,rK,dalbe,dtri,dtxi,1)                         
               do k=1,nz1                                               
                  k1 = icx(k)                                           
                  rK = bet2(k1)+alp2(i-1)                               
                  kk = mz-k                                         
                  call visc(vorwk(0,k,i),vor(0,k,i),                    
     .               phi(0,k,i),hg(0,k,i),hv(0,k,i),phiwk(0,k,i),       
     .               rK+rkn1,rK,dalbe,dtri,dtxi,1)                       
                  call visc(vorwk(0,kk,i),vor(0,kk,i),                   
     .               phi(0,kk,i),hg(0,kk,i),hv(0,kk,i),phiwk(0,kk,i),   
     .               rK+rkn1,rK,dalbe,dtri,dtxi,0)                       
               enddo                                                    
            enddo
C       icfl will change in the next step
         else       ! rkstep eq to 3
            if (mod(istep,ncfl).eq.0) then
               if (myid.eq.0) write(*,*)'he entrado con ',istep
               do i=pb,pe
                  rK = alp2(i-1)
                  call visc(vorwk(0,0,i),vor(0,0,i),
     .                 phi(0,0,i),hg(0,0,i),hv(0,0,i),phiwk(0,0,i),
     .                 rK+rkn1,rK,0d0,0d0,0d0,1)
                  do k=1,nz1
                     k1 = icx(k)
                     rK = bet2(k1)+alp2(i-1)
                     kk = mz-k
                     call visc(vorwk(0,k,i),vor(0,k,i),
     .                   phi(0,k,i),hg(0,k,i),hv(0,k,i),phiwk(0,k,i),
     .                   rK+rkn1,rK,0d0,0d0,0d0,1)
                     call visc(vorwk(0,kk,i),vor(0,kk,i),
     .                  phi(0,kk,i),hg(0,kk,i),hv(0,kk,i),phiwk(0,kk,i),
     .                  rK+rkn1,rK,0d0,0d0,0d0,0)
                  enddo
               enddo
            else   ! icfl =0 and rkstep equal to 3
               do i=pb,pe
                  rK = alp2(i-1)
                  call visc(vorwk(0,0,i),vor(0,0,i),
     .                phi(0,0,i),hg(0,0,i),hv(0,0,i),phiwk(0,0,i),
     .                rK+rkn1,rK,dalbe,dtri,0d0,1)
                  do k=1,nz1
                     k1 = icx(k)
                     rK = bet2(k1)+alp2(i-1)
                     kk = mz-k
                     call visc(vorwk(0,k,i),vor(0,k,i),
     .                  phi(0,k,i),hg(0,k,i),hv(0,k,i),phiwk(0,k,i),
     .                  rK+rkn1,rK,dalbe,dtri,0d0,1)
                     call visc(vorwk(0,kk,i),vor(0,kk,i),
     .                  phi(0,kk,i),hg(0,kk,i),hv(0,kk,i),phiwk(0,kk,i),
     .                  rK+rkn1,rK,dalbe,dtri,0d0,0)
                  enddo
               enddo
            endif
         endif



! ----------------------------------------------------------------------!
!
!             this is the block for the 00 mode
!                   Only master does it
!
! ----------------------------------------------------------------------!
        if (myid.eq.0) then
        
           call rk00(u00wk,u00,rf0u,rkn1,
     &           dalre,dtgamma,dalbe,dtri,dtxi,rkstep)
     
           call rk00(w00wk,w00,rf0w,rkn1,
     &           dalre,dtgamma,dalbe,dtri,dtxi,rkstep)
     
           call deryr(u00wk,rf0u)
           call deryr(w00wk,rf0w)
           
           massu=0.
           massw=0.
           do j=0,my-1
              massu = massu+u00wk(j)*trp(j)
              massw = massw+w00wk(j)*trp(j)
           enddo

         
        endif   ! 00 modes


 10   continue  !!! End of R-K sub-step


! ----------------------------------------------------------------------!
!   Wx0 ans Wz0 are the box averaged vorticities at the wall            !
! ----------------------------------------------------------------------!

c        !     write history record     !

         if(myid.eq.0) then

            reynota=0.5*re*re*(abs(wz0/re+uv0)+abs(wzl/re+uvL))

            if (ihist.eq.1) then
             
             
               if (istep.eq.1) then
                   ctm2=0.
                   ttm =0.
               endif
324            format(17(d22.14))
               write(29,324) time,-1.*Wz0,WzL,sqrt(reynota),ener,
     .                       massu,massw,MPI_WTIME()+iter_time,
     .                       commtimer-ctm2
               call flush(29)
325            format(i5,12(d14.6))
               write(*,325)istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,
     .                     massu,massw,H00u,H00w,MPI_WTIME()+iter_time,
     .                    transtimer-ttm,commtimer-ctm2

            endif
            ctm2 = commtimer
            ttm  = transtimer 
            
         endif

c          	! time:
         time=time+Deltat

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
      end if
      
      write(*,*) nanproblem

100   if (nanproblem.gt.0) then
         write(*,*) 'Corrupted Memory'
         write(*,*) 'Restarting from the last file...'
         
         
         call initnan(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,hv,hg,
     .                phiwk,spwk,vorwk,dvordy,uwk,chwk,sp,myid)
      endif
      

      endsubroutine

!*******************************************************************
!                                                                  !
!      computes the forcing terms (convective terms)               !
!                                                                  !
!    input:                                                        !
!                                                                  !
!      up1,up2,up3: Velocities                                     !
!      wp1,wp2,wp3: Vorticities                                    !
!      sp:          Array for spectra                              !
!                                                                  !
!   output:                                                        !
!     hv (up3): non linear term for the phi equation               !
!     hg (up2): non linear term for vorticity equation             !
!     rf0u: non linear term for the evolution of Kx=Kz=0  u        !
!     rf0w: non linear term for the evolution of Kx=Kz=0  w        !
!                                                                  !
!                                                                  !
!    updated jjs 22/12/00     (INCOMPLETE)                         !
!    single  jjs  4/01/01                                          !
!    low storage : 24/01/01 (incomplete)                           !
!    plnes-lines SHC: 29/09/05                                     !
!                                                                  !
!******************************************************************!
      subroutine hvhg(up1,up2,up3,wp1,wp2,wp3,
     .           rf0u,rf0w,chwk,sp,myid,rkstep)

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
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,
     .                  istati,ntimes,nacum,nstart
      save /statis/

      real*8 commtimer,transtimer,totaltimer
      common /timers/ commtimer, transtimer, totaltimer
      save   /timers/
      
      integer nanerror 
      common /cnan/ nanerror
      save   /cnan/  

c     ---------------------- Variables ----------------------
      
      complex*8 up1 (0:my-1,0:mgalz-1,pb:pe),
     .          up2 (0:my-1,0:mgalz-1,pb:pe),
     .          up3 (0:my-1,0:mgalz-1,pb:pe),
     .          wp1 (0:my-1,0:mgalz-1,pb:pe),
     .          wp2 (0:my-1,0:mgalz-1,pb:pe),
     .          wp3 (0:my-1,0:mgalz-1,pb:pe)
      complex*8,allocatable:: wk2(:,:,:)
      
      real*8 rf0u(0:my-1),rf0w(0:my-1),dum
      real*4 chwk(*),uner(9),aa
      real*4 sp(0:nz1,1:2*nspec+1,8,pb:pe)
      
      integer myid,rkstep,i,j,k,kk,ierr
      
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     
      
c     ---------------------- Programa ----------------------      
    
      allocate(wk2(0:my-1,0:mgalz-1,pb:pe))

c     Before Fpp, computes statistics
       
      if (istati.eq.1.and.rkstep.eq.1) then
         nacum   = nacum   +1
         nacumsp = nacumsp +1
         do kk=1,9
            uner(kk) = 0.
         enddo
         do i=pb,pe  
            call cstatis(up1(0,0,i),up2(0,0,i),up3(0,0,i),
     .                   wp1(0,0,i),wp2(0,0,i),wp3(0,0,i),
     .                   chwk,sp(0,1,1,i),uner,i)
         enddo 
         if (ihist.eq.1) then
             call MPI_ALLREDUCE(uner,ener,9,MPI_REAL,MPI_SUM,
     .                          MPI_COMM_CALC,ierr)
         
           do i=1,9
               ener(i)=sqrt(abs(ener(i)))
            enddo
         endif 
      endif    


c ----- substract umax/2 to u00 to increase dt !!!!
      if (myid.eq.0) then
         do j=0,my-1
           up1(j,0,1) = up1(j,0,1) - a0
         enddo
      endif
!      allocate(wk2(0:my-1,0:mgalz-1,pb:pe))

      call localyz2zy(up1,up1,chwk)
      call localyz2zy(up2,up2,chwk)
      call localyz2zy(up3,up3,chwk)
      call localyz2zy(wp1,wp1,chwk)      
      call localyz2zy(wp2,wp2,chwk)  
      call localyz2zy(wp3,wp3,chwk)         

      do i = pb,pe
         call fourz(up1(0,0,i),1)    !u
         call fourz(up2(0,0,i),1)    !v
         call fourz(up3(0,0,i),1)    !w
         call fourz(wp1(0,0,i),1)    ! omega_1
         call fourz(wp2(0,0,i),1)    ! omega_2
         call fourz(wp3(0,0,i),1)    ! omega_3
      enddo
      nanerror = 0
!      call check_nan(up1,my*mgalz*mmp*2,nanerror)
!      call check_nan(up2,my*mgalz*mmp*2,nanerror)
!      call check_nan(up3,my*mgalz*mmp*2,nanerror)
!      call check_nan(wp1,my*mgalz*mmp*2,nanerror)
!      call check_nan(wp2,my*mgalz*mmp*2,nanerror)
!      call check_nan(wp3,my*mgalz*mmp*2,nanerror)

      if (nanerror.eq.1) then
         write(*,*) 'NAN Found in ',myid,' before changes'
      endif
      call chpl2ln(up1,up1,chwk,myid)
      call chpl2ln(up2,up2,chwk,myid)
      call chpl2ln(up3,up3,chwk,myid)
      call chpl2ln(wp1,wp1,chwk,myid)
      call chpl2ln(wp2,wp2,chwk,myid)
      call chpl2ln(wp3,wp3,chwk,myid)

      deallocate(wk2)

c ------ hvect

!      if (myid.eq.0) then
!        dum  = -MPI_WTIME()
!      endif
!
      call hvect(up1,up2,up3,wp1,wp2,wp3,myid,rkstep)   

!      if (myid.eq.0) then
!        dum  = dum+MPI_WTIME()
!        write(*,*) ' Tiempo total hvect',dum
!      endif

c  ---------- back to the yz planes an fft 

      call chln2pl(up1,up1,chwk,myid)
      call chln2pl(up2,up2,chwk,myid)
      call chln2pl(up3,up3,chwk,myid)   

!      call check_nan(up1,my*mgalz*mmp*2,nanerror)
!      call check_nan(up2,my*mgalz*mmp*2,nanerror)
!      call check_nan(up3,my*mgalz*mmp*2,nanerror)


      if (nanerror.eq.1) then
         write(*,*) 'NAN Found in ',myid,' after changes'
      endif
c      call chln2plin(up2,chwk,myid)

      do i=pb,pe
         call fourz(up1(0,0,i),-1)    ! H3
         call fourz(up2(0,0,i),-1)    ! H1
         call fourz(up3(0,0,i),-1)    ! H2
      enddo

      call localzy2yz(up1,up1,chwk)
      call localzy2yz(up2,up2,chwk)
      call localzy2yz(up3,up3,chwk) 
           
      if (myid.eq.0) then
         do j=0,my-1
            rf0w(j) =  real(up1(j,0,pb))
            rf0u(j) =  real(up2(j,0,pb))
         enddo
      endif

           !   up3 = - dH3/dx + dH1/dz    !
           !   o1  = dH1/dx + dH3/dz      !

      do i=pb,pe
         do k=0,mz1
            do j=0,my-1
               wp1(j,k,i) =  xalp(i-1)*up2(j,k,i)+
     .                       xbet(k  )*up1(j,k,i)

               up2(j,k,i) = -xalp(i-1)*up1(j,k,i)+ ! hg  out
     .                       xbet(k  )*up2(j,k,i)  !
            enddo
         enddo
      enddo

      do i=pb,pe
         call deryr2(wp1(0,0,i),wp1(0,0,i),mz)
      enddo

           !   wp2 = d^2 H2/dx^2 + d^H2/dz^2      !

      do i=pb,pe
         do k=0,mz1
            aa = alp2(i-1)+bet2(k)
            do j=0,my-1
               up3(j,k,i) = - up3(j,k,i)*aa - wp1(j,k,i)   ! hv out
            enddo
         enddo
      enddo

      endsubroutine


!*********************************************************************!
!                                                                     !
!         computes the forcing terms (convective terms)               !
!         LINES--LINES--LINES--LINES--LINES--LINES--LINES             !
!                                                                     !
!                                                                     !
!    updated jjs 22/12/00     (INCOMPLETE)                            !
!    single  jjs  4/01/01                                             !
!    low storage : 24/01/01 (incomplete)                              ! 
!    plnes-lines SHC: 29/09/05                                        !
!*********************************************************************!
      subroutine hvect(up1,up2,up3,wp1,wp2,wp3,myid,rkstep)

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
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,
     .                  istati,ntimes,nacum,nstart
      save /statis/
  
      real*4 up1wk,up2wk,up3wk,wp1wk,wp2wk,wp3wk
      real*8 up1wk8,up2wk8,up3wk8,wp1wk8,wp2wk8,wp3wk8
      common /wkhvect/ up1wk (0:mgalx+1),up2wk(0:mgalx+1),
     .                 up3wk (0:mgalx+1),wp1wk(0:mgalx+1),
     .                 wp2wk (0:mgalx+1),wp3wk(0:mgalx+1),
     .                 up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .                 up3wk8(0:mgalx+1),wp1wk8(0:mgalx+1),
     .                 wp2wk8(0:mgalx+1),wp3wk8(0:mgalx+1)
     
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE)
      real*8  cflx,cfly,cflz,hxalp,hzalp,hyy,cfl0,reigmx1,aa
      real*4  up1(0:mx-1,lb:le),up2(0:mx-1,lb:le), 
     .        up3(0:mx-1,lb:le),wp1(0:mx-1,lb:le),
     .        wp2(0:mx-1,lb:le),wp3(0:mx-1,lb:le)

      real*8 duma,dumb,dumc,duma1,dumb1,dumc1

c     ------------------------- Indices -----------------------

      integer vec(2,3),vecrc(2,3),vecy(2)         ! i,j
      integer resulm(2,3)

!      ----------------------- Programa ------------------------      

      cflx = 0.
      cfly = 0.
      cflz = 0.
      cfl0 = 0.

      hxalp=alp*mx*0.5
      hzalp=bet*mz*0.5
      uv0=0.
      uvL=0.

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
            up1wk(i) = up1(i,j)
         enddo
         do i=0,mx-1           
            up2wk(i) = up2(i,j)
         enddo
         do i=0,mx-1
            up3wk(i) = up3(i,j)
         enddo
         do i=0,mx-1
            wp1wk(i) = wp1(i,j)
         enddo
         do i=0,mx-1
            wp2wk(i) = wp2(i,j)
         enddo
         do i=0,mx-1
            wp3wk(i) = wp3(i,j)
         enddo
        
         if (myid.eq.0) then
            duma = duma+MPI_WTIME()
         endif

         call fourx(up1wk,up1wk8,1)    !u
         call fourx(up2wk,up2wk8,1)    !v
         call fourx(up3wk,up3wk8,1)    !w
         call fourx(wp1wk,wp1wk8,1)    ! omega_1
         call fourx(wp2wk,wp2wk8,1)    ! omega_2
         call fourx(wp3wk,wp3wk8,1)    ! omega_3

c    ------------------------------------------------
c        compute statistics in PPP
c    ----------- triple products 

         jj = (j-1)/mgalz +1
         if (jj.lt.0) jj=1
         if (istati.eq.1 .and. rkstep.eq.1) then
            do i=0,mgalx-1
               aa = up2wk8(i)
               uuv(jj) = uuv(jj) +aa*up1wk8(i)**2
               wwv(jj) = wwv(jj) +aa*up3wk8(i)**2
               vvv(jj) = vvv(jj) +aa**3
            enddo

         endif

c            !      estimates maximum time step     !

         if (rkstep.eq.1.and.icfl.eq.1) then
            cfl0 = 0.
            do i = 0,mgalx-1
               cflx = max(cflx,abs(up1wk8(i)))
               cfl0 = max(cfl0,abs(up2wk8(i)))
               cflz = max(cflz,abs(up3wk8(i)))
            enddo

            cfly = max(cfly,cfl0/hy(jj))
         endif

! ----------------------------------------------------------------------!
!                                                                       !
!         computes u X omega  line by line                              !
!                                                                       !
!          up2wk = H1 = v.omega_3 - w.omega_2 (F-F-P)                   !
!          up3wk = H2 = w.omega_1 - u.omega_3 (F-F-P)                   !
!          up1wk = H3 = u.omega_2 - v.omega_1 (F-F-P)                   !
! ----------------------------------------------------------------------!

!         if (myid.eq.0) then
!            dumb = dumb-MPI_WTIME()
!         endif
         do i=0,mgalx-1
            up2wk(i) = up2wk8(i)*wp3wk8(i)
     .                -up3wk8(i)*wp2wk8(i)
            up3wk(i) = up3wk8(i)*wp1wk8(i)
     .                -up1wk8(i)*wp3wk8(i)
            up1wk(i) = up1wk8(i)*wp2wk8(i)
     .                -up2wk8(i)*wp1wk8(i)
         enddo         

!         if (myid.eq.0) then
!            dumb = dumb+MPI_WTIME()
!         endif
c--------------------- back to F-P-P  

         call fourx(up2wk,up1wk8,-1)  ! H1
         call fourx(up3wk,up1wk8,-1)  ! H2
         call fourx(up1wk,up1wk8,-1)  ! H3

c --------    back to lines. We throw away the high nodes (Dealiasing)
         if (myid.eq.0) then
            dumc = dumc-MPI_WTIME()
         endif
         do i=0,mx-1
            up1(i,j) = up1wk(i)
         enddo
         do i=0,mx-1
            up2(i,j) = up2wk(i)
         enddo
         do i=0,mx-1            
            up3(i,j) = up3wk(i)
         enddo
         if (myid.eq.0) then
            dumc = dumc+MPI_WTIME()
         endif

10    continue

c-------------------- computes Deltat

! --------------------------- WARNING --------------------------------!
!                                                                     !
!  Cambiamos los valores de cfly y cflz a la mitad.                   !
!  En el hre.dat, necesariamente cfl = 1                              !
!                                                                     !
!---------------------------------------------------------------------!
            
         cflz = cflz/2
         cfly = cfly/2 
         if (myid.eq.-2) then
            write(*,*) 'primer bucle',duma
            write(*,*) 'segun  bucle',dumb
            write(*,*) 'tercer bucle',dumc
         endif
      if (rkstep.eq.1.and.icfl.eq.1) then 
         cflx = cflx*hxalp
         cflz = cflz*hzalp
         cfl0 = max(cflx,max(cfly,cflz))
         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL8,MPI_MAX,
     .                     MPI_COMM_CALC,ierr)

         Deltat=CFL/reigmx1
         dtr=Re/Deltat

         if (reigmx1.lt.1e-1) then
            write(*,*) 'UGGG',myid,ihist,lb,le,pb,pe
         endif




!         if(myid.ne.0) then
!            call MPI_SEND(cflx,1,MPI_REAL8,0,
!     &                 myid,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(cfly,1,MPI_REAL8,0,
!     &                    myid,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(cflz,1,MPI_REAL8,0,
!     &                 myid,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(vec,6,MPI_INTEGER,0,
!     .                   myid,MPI_COMM_WORLD,ierr)
!         else
!            duma1 = cflx
!            dumb1 = cfly
!            dumc1 = cflz
!            resulm(1,1) = vec(1,1)
!            resulm(2,1) = vec(2,1)
!            resulm(1,2) = vec(1,2)
!            resulm(2,2) = vec(2,2)
!            resulm(1,3) = vec(1,3)
!            resulm(2,3) = vec(2,3)                   
!            do iproc=1,numerop-1
!               call MPI_RECV(duma,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call MPI_RECV(dumb,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call MPI_RECV(dumc,1,MPI_REAL8,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call MPI_RECV(vecrc,6,MPI_INTEGER,iproc,
!     &                      MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               
!               if (duma.gt.duma1) then
!                   duma1 = duma
!                   resulm(1,1) = vecrc(1,1)
!                   resulm(2,1) = vecrc(2,1)
!               endif
!               
!               if (dumb.gt.dumb1) then
!                   dumb1 = dumb
!                   resulm(1,2) = vecrc(1,2)
!                   resulm(2,2) = vecrc(2,2)
!               endif
!               
!               if (dumc.gt.dumc1) then
!                   dumc1 = dumc
!                   resulm(1,3) = vecrc(1,3)
!                   resulm(2,3) = vecrc(2,3)
!               endif
!            
!            enddo
!            write(*,*)'max de u',duma1,'li',resulm(1,1),'x',resulm(2,1)
!            write(*,*)'max de v',dumb1,'li',resulm(1,2),'x',resulm(2,2)
!            write(*,*)'max de w',dumc1,'li',resulm(1,3),'x',resulm(2,3)
!         endif

      endif

      endsubroutine
      
      
      
c----------------------------------------------------------------------
c    Advances the 00 modes (linear term semi-implicit, order 3)
c----------------------------------------------------------------------

      subroutine rk00(f,rhs,nl,rkn1,
     &           dalre,dtgamma,dalbe,dtri,dtxi,rkstep)
      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real*8  f(my),rhs(my),nl(my)
      real*8  rkn1,dalbe,dtri,dtxi,dalre,dtgamma
      integer rkstep
c     ----------------------- work --------------------

      integer i,j
      real*8  massin, mass1,massu

c ---------------------------  commons -----------------
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/

      real*8 wk1,fwk,df
      common /worksp/ wk1(5,my),fwk(my),df(my)
      save   /worksp/


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      This used to be deryyr (explicit part of the linear term)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      mass1=0.          !! mass due to nonlinear terms
      do j=1,my
         mass1 = mass1 + nl(j)*trp(j)
      enddo

      do j=1,my
         nl(j) = nl(j)-mass1
      enddo


      do j=1,my      !--- to real*8, fitf
         df(j)=f(j)
      enddo

      if (rkstep.eq.1) then
c   ------------------------  !! have to do the linear term explicitly

         fwk(1) =  dt22(3,1)*df(1) + dt22(4,1)*df(2) +
     .             dt22(5,1)*df(3)
         fwk(2) =  dt22(2,2)*df(1) + dt22(3,2)*df(2) +
     .             dt22(4,2)*df(3) + dt22(5,2)*df(4)

         do j=3,my-2
            fwk(j)=dt22(1,j)*df(j-2)
            do i=2,5
               fwk(j) = fwk(j) + dt22(i,j)*df(i+j-3)
            enddo
         enddo

         fwk(my1)= dt22(1,my1)*df(my-3) + dt22(2,my1)*df(my-2)+
     .             dt22(3,my1)*df(my-1) + dt22(4,my1)*df(my  )
         fwk(my) = dt22(1,my )*df(my-2) + dt22(2,my  )*df(my-1)+
     .             dt22(3,my )*df(my)

         call banbks(prem3,my,fwk)

! ---- add linear and nonlinear terms  ----

         mass1=0.          !! mass due to viscous terms
         do j=1,my
            mass1 = mass1 + df(j)*mss(j)
         enddo

         do j=1,my
            fwk(j)= fwk(j)-mass1
            df(j)  =-rkn1*(df(j)+dalre*fwk(j)+dtgamma*nl(j))
         enddo

      else
c ------------------------ !! the linear term is already in the rhs

         do j=1,my
            df(j) = -rkn1*(rhs(j)+dtgamma*nl(j))
         enddo

      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      from here down it used to be Lapv1 (implicit linear term)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (rkstep.ne.3) then      !!!!  explicit rhs
         do j=1,my
            rhs(j) = dtri*df(j)+dtxi*nl(j)
         enddo
      endif

c     calculo viscoso (matriz del poisson y LU)

      wk1(1,1)=0d0
      wk1(2,1)=0d0
      wk1(3,1)=1d0
      wk1(4,1)=0d0
      wk1(5,1)=0d0
      wk1(1,my)=0d0
      wk1(2,my)=0d0
      wk1(3,my)=1d0
      wk1(4,my)=0d0
      wk1(5,my)=0d0

      do j=2,my1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rkn1*dt21(i,j)
         enddo
      enddo

      call bandec(wk1,my)

c     resuelvo para el rhs de verdad ----------

      fwk(1) = 0d0
      fwk(2) = dt21(2,2)*df(1)+dt21(3,2)*df(2)+
     &         dt21(4,2)*df(3)+dt21(5,2)*df(4)
      do j=3,my-2
         fwk(j) = dt21(1,j)*df(j-2)+dt21(2,j)*df(j-1) +
     &            dt21(3,j)*df(j  )+dt21(4,j)*df(j+1) +
     &            dt21(5,j)*df(j+2)
      enddo
      fwk(my1) = dt21(1,my1)*df(my-3)+dt21(2,my1)*df(my-2)+
     .           dt21(3,my1)*df(my1) +dt21(4,my1)*df(my)

      fwk(my)  = 0d0

      call banbks(wk1,my,fwk)

      massin=0.
      massu =0.
      do j=1,my
         massin=massin+trp(j)*df(j)
         massu =massu +trp(j)*fwk(j)
      enddo
      massin=-massin/rkn1

c     resuelvo para el vector uniforme  --------------------

      df(1) = 0d0
      df(2) = dt21(2,2)+dt21(3,2)+ dt21(4,2)+dt21(5,2)
      do j=3,my-2
         df(j)    = dt21(1,j)+dt21(2,j)+dt21(3,j)+
     &              dt21(4,j)+dt21(5,j)
      enddo
      df(my1)= dt21(1,my1)+dt21(2,my1)+dt21(3,my1)+dt21(4,my1)
      df(my)  = 0d0

      call banbks(wk1,my,df)

c     -------------  adjust mass -----------
      mass1 =0.
      do j=1,my
         mass1=mass1+trp(j)*df(j)
      enddo

      mass1= (massin-massu)/mass1
      do j=1,my
         fwk(j) = fwk(j)+mass1*df(j)
      enddo

c     -----------  update rhs -----------
      if (rkstep.ne.3) then
         do j=1,my
            rhs(j) = rhs(j)+dalbe*fwk(j)
         enddo
      else
         do j=1,my
            rhs(j) = fwk(j)
         enddo
      endif

      do j=1,my
         f(j) = fwk(j)
      enddo


      end


!*********************************************************************!
!       subroutine   cstati.
!        compute the statistics in FPF
!
!*********************************************************************!

      subroutine cstatis(u1c,u2c,u3c,o1c,o2c,o3c,
     .                   chwkc,spaux,uner,spplane) 
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
     .       ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,
     .                  istati,ntimes,nacum,nstart
      save /statis/

c    ----------------------- Variables ---------------------------

      complex*8 u1c(my,0:mz1),u2c(my,0:mz1),u3c(my,0:mz1),
     .          o1c(my,0:mz1),o2c(my,0:mz1),o3c(my,0:mz1),
     .        chwkc(my,0:mz1)
      real*4 spaux(0:nz1,2*nspec+1,8)  
      integer j,k,iy,kk,spplane
      real*4 hyy

      real*8 aa
      complex*16 cc

c     --------------------- Programa -----------------------------

c     -----  spectra



      if (spplane.eq.1)  then          ! i = 0
         do iy=1,2*nspec+1
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
         do iy=1,2*nspec+1
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
