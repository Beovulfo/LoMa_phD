
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
     .                  scal1,
     .                  scal2,
     .                  scal3,
     .                  scal,
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

      real*8  scadiff
      common /fis2/ scadiff(1:nscalars)
      save   /fis2/
      
c ------------------- Variables ----------------------------------------

      complex*8 phi(0:my-1,0:mz1,pb:pe),vor(0:my-1,0:mz1,pb:pe),
     .         scal(0:my-1,0:mz1,pb:pe,nscalars)

      complex*8 hg    (0:my-1,0:mgalz-1,pb:pe),
     .          vorwk (0:my-1,0:mgalz-1,pb:pe),
     .          vorrhs(0:my-1,0:mgalz-1,pb:pe),
     .          hv    (0:my-1,0:mgalz-1,pb:pe),
     .          phiwk (0:my-1,0:mgalz-1,pb:pe),
     .          phirhs(0:my-1,0:mgalz-1,pb:pe),
     .          uwk   (0:my-1,0:mgalz-1,pb:pe),
     .          dvordy(0:my-1,0:mgalz-1,pb:pe),
     .          chwk  (0:my-1,0:mgalz-1,pb:pe),
     .          scal1 (0:my-1,0:mgalz-1,pb:pe,nscalars),
     .          scal2 (0:my-1,0:mgalz-1,pb:pe,nscalars),
     .          scal3 (0:my-1,0:mgalz-1,pb:pe,nscalars),
     .          scalrhs(0:my-1,0:mgalz-1,pb:pe,nscalars)

      real*8 u00(0:my-1),w00(0:my-1),rf0u(0:*),u00wk(0:*),
     .       rf0w(0:*),w00wk(0:*)

      real*4 sp  (0:nz1,1:2*nspec+1,8,pb:pe),
     .       spwk(0:nz1,1:  nspec+1,8,pb:pe)
     
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,kk,iscal
      
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
!set to zero the rhs terms
      do i = pb,pe
         do k = 0,mgalz-1
             do j = 0,my-1
                  vorrhs(j,k,i) =0.0
                  phirhs(j,k,i) =0.0 
              enddo
          enddo
       enddo

!set to zero the rhs terms for scal
      do iscal=1,nscalars
        do i = pb,pe
           do k = 0,mgalz-1
               do j = 0,my-1
                  scalrhs(j,k,i,iscal)=0.0
               enddo
            enddo
         enddo
       enddo


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
            
            call escru(vor,phi,u00,w00,sp,spwk,myid)
!write scalar field
!CAUTION escruscal needs to be reviewed for NSCALARS
            do iscal=1,nscalars
               call escruscal(scal(0,0,pb,iscal),iscal,myid)
            enddo
!aaf Here we need to call a escru for scalar field

            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
            endif
         ENDIF

! ------------------- finished writing image --------------------------!

      if (myid.gt.numerop-1) goto 30    ! only for save procs

      if(irun.eq.0) then ! this is done only for the first step
         irun = 1

!     prepare u00wk,w00wk,phiwk and vorwk !
         if (myid.eq.0) then
            do j=0,my1
              u00wk(j)=u00(j)
              w00wk(j)=w00(j)
            enddo

            call deryr(u00wk,rf0u)
            call deryr(w00wk,rf0w)
         endif


      endif ! end special first step
      
               !  Runge-Kutta third order  !
                
!----------- RUNGE-KUTTA SUBSTEPS START----------------------------------

      do 10 rkstep=1,3

c-------------save vor and phi in vorwk/phiwk to work---------------
         do i=pb,pe
            do k=0,mz1
               do j=0,my-1
                  vorwk(j,k,i)=vor(j,k,i)
                  phiwk(j,k,i)=phi(j,k,i)
               enddo
            enddo
         enddo

c------------save scalars in wk spaces
         do iscal=1,nscalars
           do i=pb,pe
              do k=0,mz1
                 do j=0,my-1
!save scalar field into scal1
!aaf care with size of scal1 vs size scal
                    scal1(j,k,i,iscal)=scal(j,k,i,iscal)
                 enddo
              enddo
            enddo
         enddo


c--------------  dvordy      : chwk= d (vorwk) / dy -- F-F
         do i=pb,pe
            call deryr2(vorwk(0,0,i),chwk(0,0,i),mz)
         enddo
c---------------------------------------------------------

c----calculate derivative of scal: scal2=d(scal1)/dy -- F-F
         do iscal=1,nscalars
           do i=pb,pe
              call deryr2(scal1(0,0,i,iscal),scal2(0,0,i,iscal),mz)
           enddo
         enddo

c---------Calculate d(scal)/dx,d(scal)/dz------------------------F-F
         do iscal=1,nscalars
           do i=pb,pe
              do k=0,mz1
                 do j=0,my-1
                    !dc/dz=i kz c
                    scal3(j,k,i,iscal)=xbet(k)*scal1(j,k,i,iscal)
                    !dc/dx=i kx c
                    scal1(j,k,i,iscal)=xalp(i-1)*scal1(j,k,i,iscal)
                 enddo
              enddo
          enddo
        enddo

c          !   calcula la v a partir de phi !
         do i=pb,pe
            do k=0,mz1
               k1 = icx(k)
               rK = bet2(k1)+alp2(i-1)
               call Lapvdv(phi(0,k,i),hg(0,k,i),hv(0,k,i),rK)
            enddo
         enddo


c ---- Computation of velocity and vorticity 

! ----------------------------------------------------------------------!
!   at this point
!   vor:    vorticity
!   phi:    phi
!   vorwk:  vorticity
!   phiwk:  phi 
!   chwk: d (vorwk) / dy 
!   hg:     v
!   hv:     dv/dy
!   uwk:    empty 
!aaf
!   scal: scalar field
!   scal1: d(scal)/dx
!   scal2: d(scal)/dy
!   scal3: d(scal)/dz
!   
! ----------------------------------------------------------------------!


c---------------------- computes  u,w, ome1, ome3
! Solving eq system from:
! (1) Continuity eq,
! (2)Vor_y definition ---> u, w
! (3) phi=(-rot(vor))_y = dvor_z/dx - dvor_x/dz 
! (4) Div(vor)=0; ---> vor_x, vor_z

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
!   vor:   vorticity
!   phi:   phi
!   vorwk: vor_y
!   phiwk: vor_x  
!   dvordy:vor_z  
!   uwk:   u   
!   hg:    v
!   hv:    w   
!aaf
!   scal: scalar field
!   scal1: d(scal)/dx
!   scal2: d(scal)/dy
!   scal3: d(scal)/dz
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
      
! hvhg computes non-linear terms, change fourier to physical and
! change again from physical to fourier (making calculations that
! need physical expression in between). 

        call hvhg(uwk,hg,hv,phiwk,vorwk,dvordy,
     .          rf0u,rf0w,chwk,scal1,scal2,scal3,
     .          sp,myid,rkstep)
! outputs: hg,hv,rf0u,rf0w,scal1
!   respecting the ordering of data buffer in hvhg call

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

c------- After computes CFL (computed in hvect)
         dtgamma = Deltat*gama(rkstep)
         dtxi    = Deltat*xi(rkstep)
       
!debug
!          if (myid.eq.0) then
!             do k=0,mz-1
!                do j=0,my-1
!                 scal2(j,k,pb)=scal(j,k,pb)
!                enddo
!            enddo
!          endif
!         write(*,*) 'myid',myid,'scal1(0,0,pb)',scal1(0,0,pb)
!explicit scheme
          do i=pb,pe
!compute new non linear term, update phi, saves the nonlinear term.
              call rkstepexp(vor(0,0,i),vorrhs(0,0,i),hg(0,0,i),mz,i,
     .             ire,dtgamma,dtxi)
                 
              call rkstepexp(phi(0,0,i),phirhs(0,0,i),hv(0,0,i),mz,i,
     .             ire, dtgamma,dtxi)
          enddo
!scalar rkstepexp
          do iscal=1,nscalars
             do i=pb,pe
               call rkstepexp(scal(0,0,i,iscal),scalrhs(0,0,i,iscal),
     .         scal1(0,0,i,iscal),mz,i,scadiff(iscal),dtgamma,dtxi)
             enddo
           enddo
!rkstep=3
          if (rkstep.eq.3) then
            if (mod(istep,ncfl).eq.0) then
              if (myid.eq.0) then 
                write(*,*) 'he entrado con ', istep
!       write(*,*) 'CHECK scal EVO',scal(0,0,pb)-scal2(0,0,pb)
!       write(*,*) 'CHECK scal2',scal2(0,0,pb)
!       write(*,*) 'CHECK scal',scal(0,0,pb)
              !nothing to do here now. 
              endif        
            endif  
         endif

! ----------------------------------------------------------------------!
!
!             this is the block for the 00 mode
!                   Only master does it
!
! ----------------------------------------------------------------------!
!        call mpi_barrier(MPI_COMM_WORLD,ierr) 
        if (myid.eq.0) then

           call rk00exp(u00wk,u00,rf0u,dtgamma,dtxi,rkstep,ire,1)
           call rk00exp(w00wk,w00,rf0w,dtgamma,dtxi,rkstep,ire,2)

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
!calculate eps

            if (ihist.eq.1) then
            reynota=0.
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
         
         
         call initnan(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,hv,hg,
     .                phiwk,spwk,vorwk,dvordy,uwk,chwk,scal1,scal2,
     .                scal3,scal,sp,myid)
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
!      scal1, scal2, scal3: scalar derivatives
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
!    scalar added aaf: 26/01/2014
!                                                                  !
!******************************************************************!
      subroutine hvhg(up1,up2,up3,wp1,wp2,wp3,
     .           rf0u,rf0w,chwk,scal1,scal2,scal3,
     .           sp,myid,rkstep)

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
      
      complex*8 up1 (0:my-1,0:mgalz-1,pb:pe),
     .          up2 (0:my-1,0:mgalz-1,pb:pe),
     .          up3 (0:my-1,0:mgalz-1,pb:pe),
     .          wp1 (0:my-1,0:mgalz-1,pb:pe),
     .          wp2 (0:my-1,0:mgalz-1,pb:pe),
     .          wp3 (0:my-1,0:mgalz-1,pb:pe),
     .         scal1(0:my-1,0:mgalz-1,pb:pe,nscalars),
     .         scal2(0:my-1,0:mgalz-1,pb:pe,nscalars),
     .         scal3(0:my-1,0:mgalz-1,pb:pe,nscalars)
               
      
      real*8 rf0u(0:my-1),rf0w(0:my-1),dum
      real*4 chwk(*),uner(9),aa
      real*4 sp(0:nz1,1:2*nspec+1,8,pb:pe)
      
      integer myid,rkstep,i,j,k,kk,ierr,iscal
      
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     
      
c     ---------------------- Programa ----------------------      

c     Before Fpp, computes statistics
       
      if (istati.eq.1.and.rkstep.eq.1) then
         nacum   = nacum   +1
         nacumsp = nacumsp +1
         do kk=1,9
            uner(kk) = 0.
         enddo
         do i=pb,pe  
            call cstatis(up1(0,0,i),up2(0,0,i),up3(0,0,i),
     .                  wp1(0,0,i),wp2(0,0,i),wp3(0,0,i),
     .                   chwk,sp(0,1,1,i),uner,i)
         enddo 
         if (ihist.eq.1) then
             call MPI_ALLREDUCE(uner,ener,9,MPI_REAL,MPI_SUM,
     .                          MPI_COMM_CALC,ierr)
         
           do i=1,9
               ener(i)=sqrt(abs(ener(i)))
            enddo

!aaf------- momentum thickness calculation---------!
! dm will keep the momentum thickness
            dm=0d0
            do j=1,my
               dm=dm+(hy(j)*(0.25-(um(j)/(um(my)-um(1)))**2))
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
!  local transpose yz to zy, and adding zeros between positive
!  and negative wavenumbers.
      call localyz2zy(up1,up1,chwk)
      call localyz2zy(up2,up2,chwk)
      call localyz2zy(up3,up3,chwk)
      call localyz2zy(wp1,wp1,chwk)      
      call localyz2zy(wp2,wp2,chwk)  
      call localyz2zy(wp3,wp3,chwk)         
       
! inverse transform (fou --> fis @z),does the compacting aswell
      do i = pb,pe
         call fourz(up1(0,0,i),1)    !u
         call fourz(up2(0,0,i),1)    !v
         call fourz(up3(0,0,i),1)    !w
         call fourz(wp1(0,0,i),1)    ! omega_1
         call fourz(wp2(0,0,i),1)    ! omega_2
         call fourz(wp3(0,0,i),1)    ! omega_3
      enddo

      nanerror = 0

!  change plane to line(1:mx,lb:le)
!  and distribute package of lines to processors
      call chpl2ln(up1,up1,chwk,myid)
      call chpl2ln(up2,up2,chwk,myid)
      call chpl2ln(up3,up3,chwk,myid)
      call chpl2ln(wp1,wp1,chwk,myid)
      call chpl2ln(wp2,wp2,chwk,myid)
      call chpl2ln(wp3,wp3,chwk,myid)

!Transformation for scalars
!
      do iscal=1,nscalars
       !change yz to zy and prepare for fourz calling
       call localyz2zy(scal1(0,0,pb,iscal),scal1(0,0,pb,iscal),chwk)  
       call localyz2zy(scal2(0,0,pb,iscal),scal2(0,0,pb,iscal),chwk)  
       call localyz2zy(scal3(0,0,pb,iscal),scal3(0,0,pb,iscal),chwk)  
       !for each plane now fourz
       do i=pb,pe
         call fourz(scal1(0,0,i,iscal),1)  ! dscal/dx
         call fourz(scal2(0,0,i,iscal),1)  ! dscal/dy
         call fourz(scal3(0,0,i,iscal),1)  ! dscal/dz
       enddo
       !change from planes to lines 
      call chpl2ln(scal1(0,0,pb,iscal),scal1(0,0,pb,iscal),chwk,myid)
      call chpl2ln(scal2(0,0,pb,iscal),scal2(0,0,pb,iscal),chwk,myid)
      call chpl2ln(scal3(0,0,pb,iscal),scal3(0,0,pb,iscal),chwk,myid)

      enddo


!
! before this call upi,wpi are still lines in fourier space x (phys z).
      call hvect(up1,up2,up3,wp1,wp2,wp3,scal1,scal2,scal3,myid,rkstep)
! up1,up2,up3 are H1,H2,H3
! they are still lines in fourier x after hvect
! scal1 is NLscal in fourier x

!      if (myid.eq.0) then
!        dum  = dum+MPI_WTIME()
!        write(*,*) ' Tiempo total hvect',dum
!      endif

c  ---------- back to the yz planes an fft 

      call chln2pl(up1,up1,chwk,myid)
      call chln2pl(up2,up2,chwk,myid)
      call chln2pl(up3,up3,chwk,myid)   


      if (nanerror.eq.1) then
         write(*,*) 'NAN Found in ',myid,' after changes'
      endif

!  convert to fourier Z now.
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

!Non-Linear Scalar term transformation
      do iscal=1,nscalars
      !change lines to planes
        call chln2pl(scal1(0,0,pb,iscal),scal1(0,0,pb,iscal),chwk,myid)
        !fourz inverse transformation
        do i=pb,pe
            call fourz(scal1(0,0,i,iscal),-1)
        enddo
        !change zy to yz
        call localzy2yz(scal1(0,0,pb,iscal),scal1(0,0,pb,iscal),chwk)   !NLscal out
        !NLscal is ready now to go out
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
      subroutine hvect(up1,up2,up3,wp1,wp2,wp3,scal1,scal2,scal3,
     .                 myid,rkstep)

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
  
      real*4 up1wk,up2wk,up3wk,wp1wk,wp2wk,wp3wk,
     .       scal1wk,scal2wk,scal3wk
      real*8 up1wk8,up2wk8,up3wk8,wp1wk8,wp2wk8,wp3wk8,
     .       scal1wk8,scal2wk8,scal3wk8
      common /wkhvect/ up1wk (0:mgalx+1),up2wk(0:mgalx+1),
     .                 up3wk (0:mgalx+1),wp1wk(0:mgalx+1),
     .                 wp2wk (0:mgalx+1),wp3wk(0:mgalx+1),
     .                 up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .                 up3wk8(0:mgalx+1),wp1wk8(0:mgalx+1),
     .                 wp2wk8(0:mgalx+1),wp3wk8(0:mgalx+1),
     .                 scal1wk(0:mgalx+1),scal2wk(0:mgalx+1),
     .                 scal3wk(0:mgalx+1),
     .                 scal1wk8(0:mgalx+1),scal2wk8(0:mgalx+1),
     .                 scal3wk8(0:mgalx+1)
    
     
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     

c     --------------------------- Variables -------------------      

      integer myid,rkstep,ierr,i,j,jj,iproc,istat(MPI_STATUS_SIZE),iscal
      real*8  cflx,cfly,cflz,hxalp,hzbet,hyy,cfl0,reigmx1,aa
      real*4  up1(0:mx-1,lb:le),up2(0:mx-1,lb:le), 
     .        up3(0:mx-1,lb:le),wp1(0:mx-1,lb:le),
     .        wp2(0:mx-1,lb:le),wp3(0:mx-1,lb:le),
     .        scal1(0:mx-1,lb:le,nscalars),
     .        scal2(0:mx-1,lb:le,nscalars),
     .        scal3(0:mx-1,lb:le,nscalars)

      real*8 duma,dumb,dumc,duma1,dumb1,dumc1
      real*8 pi
      
      real*8  up1fis8(0:mgalx+1,lb:le),
     .        up2fis8(0:mgalx+1,lb:le),
     .        up3fis8(0:mgalx+1,lb:le)
      

c     ------------------------- Indices -----------------------

      integer vec(2,3),vecrc(2,3),vecy(2)         ! i,j
      integer resulm(2,3)

!      ----------------------- Programa ------------------------      
!initialization
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
! computing viscous Dt inverse:
! Dt_visc=1/(2 nu)*Dy_min^2
      hyy=2/Re*(max(1/minval(hy),max(hxalp,hzbet)))**2
      
c   -------------------------------------------------

c        Move everything to pPP, line by line 

      duma = 0d0
      dumb = 0d0
      dumc = 0d0
      do 10 j = lb,le

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

! convert lines in fourier to lines in physical for x
! upiwk8 means real 8 format
         call fourx(up1wk,up1wk8,1)    !u
         call fourx(up2wk,up2wk8,1)    !v
         call fourx(up3wk,up3wk8,1)    !w
         call fourx(wp1wk,wp1wk8,1)    ! omega_1
         call fourx(wp2wk,wp2wk8,1)    ! omega_2
         call fourx(wp3wk,wp3wk8,1)    ! omega_3

!
! now lines of mgalx, each mgalz one height of y.
!========================================================!
!IMPORTANT POINT: all quantities are now PHYSICAL
!========================================================!
! now save the velocity P-P-P in buffer
         do i=0,mgalx-1
            up1fis8(i,j) = up1wk8(i)
         enddo
         do i=0,mgalx-1
            up2fis8(i,j) = up2wk8(i)
         enddo
         do i=0,mgalx-1            
            up3fis8(i,j) = up3wk8(i)
         enddo
!
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
            do i=0,mgalx-1
               aa = up2wk8(i)
               uuv(jj) = uuv(jj) +aa*up1wk8(i)**2
               wwv(jj) = wwv(jj) +aa*up3wk8(i)**2
               vvv(jj) = vvv(jj) +aa**3
            enddo

         endif

c            !      estimates maximum time step     !
! still within loop for each line
         if (rkstep.eq.1.and.icfl.eq.1) then
            cfl0 = 0.
            do i = 0,mgalx-1
               cflx = max(cflx,abs(up1wk8(i)))
               cfl0 = max(cfl0,abs(up2wk8(i)))
               cflz = max(cflz,abs(up3wk8(i)))
!cflx and cflz have same hx,hz for all domain
!but cfly needs to account for hy. 
            enddo
            cfly=max(cfly,cfl0/hy(jj))
         endif

! ----------------------------------------------------------------------!
!                                                                       !
!         computes u X omega  line by line                              !
!                                                                       !
!          up2wk = H1 = v.omega_3 - w.omega_2 (F-F-P)                   !
!          up3wk = H2 = w.omega_1 - u.omega_3 (F-F-P)                   !
!          up1wk = H3 = u.omega_2 - v.omega_1 (F-F-P)                   !
! ----------------------------------------------------------------------!

         do i=0,mgalx-1
            up2wk(i) = up2wk8(i)*wp3wk8(i)
     .                -up3wk8(i)*wp2wk8(i)
            up3wk(i) = up3wk8(i)*wp1wk8(i)
     .                -up1wk8(i)*wp3wk8(i)
            up1wk(i) = up1wk8(i)*wp2wk8(i)
     .                -up2wk8(i)*wp1wk8(i)
         enddo         

c--------------------- back to F-P-P  

         call fourx(up2wk,up1wk8,-1)  ! H1
         call fourx(up3wk,up1wk8,-1)  ! H2
         call fourx(up1wk,up1wk8,-1)  ! H3

c --------    back to lines. We throw away the high modes (Dealiasing)

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

!-----Cambiamos los valores de cfly y cflz a la mitad
!      cflz=cflz/2
!      cfly=cfly/2
      
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

      endif

!===================================================================!
!NOW lets do SCALARS transformations
!===================================================================!
      do iscal=1,nscalars
        do 20 j = lb,le
      !save line into wk space
           do i=0,mx-1
              scal1wk(i) = scal1(i,j,iscal)
           enddo
           do i=0,mx-1
              scal2wk(i) = scal2(i,j,iscal)
           enddo
           do i=0,mx-1
              scal3wk(i) = scal3(i,j,iscal)
           enddo
 
           !go to phys
           call fourx(scal1wk,scal1wk8,1)    ! d(scal)/dx
           call fourx(scal2wk,scal2wk8,1)    ! d(scal)/dy
           call fourx(scal3wk,scal3wk8,1)    ! d(scal)/dz
           
           !calculate non-linear terms using buffers of velocities
           do i=0,mgalx-1
              scal1wk(i)=-up1fis8(i,j)*scal1wk8(i)
     .                   -up2fis8(i,j)*scal2wk8(i)
     .                   -up3fis8(i,j)*scal3wk8(i)
           enddo

           !back to fourier in X
           call fourx(scal1wk,scal1wk8,-1)  ! NLscal
           !save back the results into main buffer scal1
           do i=0,mx-1
              scal1(i,j,iscal)=scal1wk(i)
           enddo

20      continue
      enddo !end loop for each scalar


      endsubroutine
      
      
      
c----------------------------------------------------------------------
c    Advances the 00 modes (linear term semi-implicit, order 3)
c----------------------------------------------------------------------

      subroutine rk00(f,rhs,nl,rkn1,dalre,dtgamma,dalbe,
     .                dtri,dtxi,rkstep,flag)
      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real*8  f(my),rhs(my),nl(my)
      real*8  rkn1,dalbe,dtri,dtxi,dalre,dtgamma
      integer rkstep,flag
c     ----------------------- work --------------------

      integer i,j
      real(8),dimension(my)   ::  fwkh1,fwkh2,dfwk,dfwkh1,dfwkh2
      real(8)                 ::  coef(2),AA(2,2)

c ---------------------------  commons -----------------
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/

      real*8 wk1,fwk,df
      common /worksp/ wk1(5,my),fwk(my),df(my)
      save   /worksp/
!--------- ----------- BC functions----------------------------
      real(8) :: Cvtop,Cdyvtop,BC_top
      real(8) :: Cvbot,Cdyvbot,BC_bot

      real*8 temp
! ----------------------------------------------------------


!u00 -> flag=1
!w00 -> flag=2

      if (flag.eq.1) then  !u00
!BC constants
!     Cvtop*u00(my)+Cdyvtop*u00(my)=BC_top
!     Cvbot*u00(1)+Cdyvbot*u00(1)=BC_bot
!     
!     TOP
!    -----------
!Dirichlet
      Cvtop=1d0
      Cdyvtop=0d0
      BC_top=1d0
!Von Neuman
!      Cvtop=0d0
!      Cdyvtop=1d0
!      BC_top=0d0

!     BOT
!     ---------
!Dirichlet
      Cvbot=1d0
      Cdyvbot=0d0
      BC_bot=-1d0
!Von Neuman
!      Cvbot=0d0
!     Cdyvbot=1d0
!      BC_bot=0d0

      elseif (flag.eq.2) then   !w00
!BC constants
!     Cvtop*u00(my)+Vdyvtop*u00(my)=BC_top
!     Cvbot*u00(1)+Vdyvbot*u00(1)=BC_bot
!     
!     TOP
!    -----------
        Cvtop=1d0
        Cdyvtop=0d0
        BC_top=0d0

!     BOT
!     ---------
        Cvbot=1d0
        Cdyvbot=0d0
        BC_bot=0d0
      endif 
!-----------------------------------------------__!

      do j=1,my      !--- to real*8, fitf
         df(j)=f(j)
      enddo
!aaf  First RK-step we, dont know rhs yet.
      if (rkstep.eq.1) then
!         call deryyr(df,nl,rkn1,dalre,dtgamma)
        fwk(1) =  dt22(3,1)*df(1) + dt22(4,1)*df(2) + 
     $            dt22(5,1)*df(3)
        fwk(2) =  dt22(2,2)*df(1) + dt22(3,2)*df(2) + 
     $            dt22(4,2)*df(3) + dt22(5,2)*df(4)
        do j=3,my-2
           temp = 0d0
           do i =1,5
              temp = temp + dt22(i,j)*df(i+j-3)
           enddo
           fwk(j) = temp
        enddo
        fwk(my1)=dt22(1,my1)*df(my-3) + dt22(2,my1)*df(my-2)+ 
     $            dt22(3,my1)*df(my-1) + dt22(4,my1)*df(my  )
        fwk(my)=  dt22(1,my  )*df(my-2) + dt22(2,my  )*df(my-1)+ 
     $            dt22(3,my  )*df(my)
! solve second derivative of u(rkstep) -> fwk
! fwk=d2y2(u00)
        call banbks(prem3,my,fwk)

! rhs for step 1
        do j=1,my
           df(j) = -rkn1*(df(j)+dalre*fwk(j)+dtgamma*nl(j))
        enddo

      else
!aaf   This is done for all rksteps except ONE
         do j=1,my
            df(j) = -rkn1*(rhs(j)+dtgamma*nl(j))
         enddo

      endif

!---------------- Calculation of new u(rkstep)------------------------!
!
      if (rkstep.ne.3) then
         do j=1,my
            rhs(j)=dtri*df(j)+dtxi*nl(j)
         enddo
      endif
!Prepare wk1 matrix: wk1=dt22-rk dt21
      wk1(1,1)  = 0d0
      wk1(2,1)  = 0d0
      wk1(3,1)  = 1d0
      wk1(4,1)  = 0d0  
      wk1(5,1)  = 0d0  
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0

      do j=2,my-1
        do i=1,5
          wk1(i,j)=dt22(i,j)-rkn1*dt21(i,j)
        enddo
      enddo

!numerical recipes codebase--
!given nxn band diagonal matrix wk1, construct LU decomposition
!used in conjuction with banbks to solve band-diagonal sets of equations
      call bandec(wk1,my)
! ----------------------------------------------------------------------!
      call calcrhslap(df,fwk)  !prepares RHS multiplying by dt21
!solves band diagonal linear equations Ax=b --returning b (as fwk)
!            A      b
      call banbks(wk1,my,fwk)

! Now fwk has u(rkstep) particular
! Now solving the homogeneus solutions
      fwkh1(1) = 1d0
      fwkh2(1) = 0d0
      do j=2,my-1
         fwkh1(j)=0d0
         fwkh2(j)=0d0
      enddo  
      fwkh1(my) = 0d0
      fwkh2(my) = 1d0
! solving for uh1, uh2
      call banbks(wk1,my,fwkh1)
      call banbks(wk1,my,fwkh2) 
        
!FLAG
!      write(*,*) "fwk(1),fwk(2)",fwk(1),fwk(2)
!.................................................
!     Caculating the derivative:
!     multiplying by dt12, preparing the derivative calculation
      call calcrhsd1(fwk,dfwk)     
      call calcrhsd1(fwkh1,dfwkh1) 
      call calcrhsd1(fwkh2,dfwkh2) 

!FLAG
!      write(*,*) "dfwk(1),dfwk(2)",dfwk(1),dfwk(2)
!.................................................

!     derivative in y
      call banbks7(prem1,my,dfwk  )!d(up)/dy
      call banbks7(prem1,my,dfwkh1)!d(uh1)/dy
      call banbks7(prem1,my,dfwkh2)!d(uh2)/dy


!     Building the A matrix
!     part multiplied by coef(1) on the left for BC
!     at bottom
      AA(1,1)=Cvbot+Cdyvbot*dfwkh1(1)*fmap(1)
!     part multiplied by coef(2) on the left for BC
!     at bottom
      AA(1,2)=Cdyvbot*dfwkh2(1)*fmap(1)
!     same for top BC
      AA(2,1)=Cdyvtop*dfwkh1(my)*fmap(my)
      AA(2,2)=Cvtop + Cdyvtop*dfwkh2(my)*fmap(my)
!     RHS of linear sistem A · coef' = B
      coef(1)=BC_bot-Cdyvbot*dfwk( 1)*fmap(1 )
      coef(2)=BC_top-Cdyvtop*dfwk(my)*fmap(my)
!FLAG
!      write(*,*) "dfwk(1),dfwk(2)",dfwk(1),dfwk(2)
!      write(*,*) "coef(1),coef(2)",coef(1),coef(2)
!...............................................      
!     solve system
      call gaussj(AA,2,2,coef,2,1)
!     sum all solutions for u(rkstep) fwk

      do j=1,my
         fwk(j)=fwk(j)+ coef(1)*fwkh1(j) + coef(2)*fwkh2(j)
      enddo
! new u(rkstep) has been calculated (fwk)
! preparation of rhs for next step
      if (rkstep.ne.3) then  
        do j=1,my
            rhs(j)=rhs(j)+dalbe*fwk(j)
        enddo
      else
         do j=1,my
            rhs(j)=fwk(j)
         enddo
      endif

      do j=1,my
         f(j)=fwk(j)
      enddo


      endsubroutine


c----------------------------------------------------------------------
c    Advances the 00 modes (linear term explicit, order 3)
c----------------------------------------------------------------------

      subroutine rk00exp(f,rhs,nl,dtgamma,dtxi,rkstep,ire,flag)
      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real*8  f(my),rhs(my),nl(my)
      real*8  dtxi,dtgamma,ire
      integer rkstep,flag
c     ----------------------- work --------------------

      integer i,j

c ---------------------------  commons -----------------
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/

      real*8 wk1,fwk,df
      common /worksp/ wk1(5,my),fwk(my),df(my)
      save   /worksp/
!--------- ----------- BC functions----------------------------

      real*8 temp
! ----------------------------------------------------------

!f(j) has the u00 in this substep
      do j=1,my      !--- to real*8, fitf
         df(j)=f(j)
      enddo
!         call deryyr(df,nl,rkn1,dalre,dtgamma)
        fwk(1) =  dt22(3,1)*df(1) + dt22(4,1)*df(2) + 
     $            dt22(5,1)*df(3)
        fwk(2) =  dt22(2,2)*df(1) + dt22(3,2)*df(2) + 
     $            dt22(4,2)*df(3) + dt22(5,2)*df(4)
        do j=3,my-2
           temp = 0d0
           do i =1,5
              temp = temp + dt22(i,j)*df(i+j-3)
           enddo
           fwk(j) = temp
        enddo
        fwk(my1)=dt22(1,my1)*df(my-3) + dt22(2,my1)*df(my-2)+ 
     $            dt22(3,my1)*df(my-1) + dt22(4,my1)*df(my  )
        fwk(my)=  dt22(1,my  )*df(my-2) + dt22(2,my  )*df(my-1)+ 
     $            dt22(3,my  )*df(my)
! solve second derivative of u(rkstep) -> fwk
! fwk=d2y2(u00)
        call banbks(prem3,my,fwk)
       
!       if (flag.eq.1) then
       !u00   
       !BC - Dirichlet
         ! f(1)=-1d0;
         ! f(my)=1d0
!     !  else
       !w00
       !BC - Dirichlet
!         ! f(1)=0d0;
!         ! f(my)=0d0;
!     !  endif

        do j=1,my 
           df(j) = nl(j)+ire*fwk(j) !local buffer
           !new u00/w00:
           f(j)  = f(j)+dtgamma*df(j)+dtxi*rhs(j)
           rhs(j)= df(j) !nl for the next substep
        enddo
       


      endsubroutine
!




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
