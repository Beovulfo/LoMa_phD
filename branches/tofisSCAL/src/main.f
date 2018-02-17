!***********************************************************************!
!                                                                       !
!     Compute and write fields in phys space                            !
!                                                                       !
!                                                                       !
!     VERSION 1 -> Re_tau 2000, CFD, pl-ln                              !
!                                                                       !
!     SHC     15/12/05                                                  !
!     OFA     22/07/13  Simplified for toni, only u,v,w, ...            !
!                                                                       !
!***********************************************************************!
 
      program main
      
      implicit none 
      include "ctes3D"
      
      ! -------------------- Commons ----------------------------------!

      integer iinp,iswap
      character*80 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(8)
      save /ficheros/
      
      real*4  Re,alp,bet,a0
      real*8  y,fmap
      common /fis/ Re,alp,bet,a0,y(my),fmap(my)
      save   /fis/

      integer nvecy,plx
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr,nvecy,plx
      save /tem/

      integer      icx
      real*4       alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     &              icx(0:mz1)
      save  /wave/

      ! ----------------------- constants -----------------------------!
      
      integer nbuffou,nbuffis,i,j,k,nlines,nbuffbig,jj,buc,fj
!      integer nbuffou,nbuffis,i,j,k,isn,nlines,nbuffbig,jj,buc,fj
      real*8  u00(my),w00(my),du00(my),dw00(my)
!aaf added fj

      integer time1,time2,time3,time4,time5,timepk1,timepk2,irate
      real*8  timetot,timel1,timel2a,timel2b,timel3,timel2

      integer tb1,tb2,tb3,tb4,tb5,ind,nfiles,ifiles
      real*4 ttot1,ttot2,ttot3,dum
      
      real*4, allocatable:: vor(:,:,:),phi(:,:,:),      
     .        u(:),v(:),w(:),du(:,:),dudy(:),dvdy(:),dwdy(:),dvordy(:),
     .        wk1(:),wk(:,:),pr(:),ps(:)
!     % wkb(:),press(:,:),presr(:,:),pr(:),ps(:),bcpress(:,:,:),stok (:),
!pr and ps added to buffer by aaf
     
      character*80 fnamesta,fnameima
!Loop for reading nfiles instead of only one!
      
         nfiles=10;
!        nfiles=1;
!*********       
      do ifiles=1,nfiles 
!All within the loop            
      ! ------------------------ Program ------------------------------!

c      /*   initializes everything    */  

!      isn = 19
      ind=1

c--------------- initializes commons and things 

      write(*,*) 'Initalizing ...'      
      
      call initcr(ifiles,nfiles)  
      nbuffou = 2*mz*my*plx
      
      allocate(phi(2*my,0:mz1,plx))
      allocate(vor(2*my,0:mz1,plx))
      allocate(v(nbuffou))
     
 
      !----------- begin operations in yz planes ---------------------!
      
      do i=1,mx1+1,plx  !,mx1+1,plx ! !
         write(*,*) 'Plane: ',i,'of ', mx1+1
        ! ---       reads ix=constant plane
         do jj=1,plx
            read(iinp)((phi(j,k,jj),j=1,2*my),k=0,mz1)
         enddo 
         
         call swap(phi  ,vor,vor,i,plx,1)  ! u is phi

      enddo

      close (iinp)

      write(*,*) '-----------------------------------------------'

      !  Lines - lines - lines 
      
      ! Deallocate memory
      
      deallocate(vor,phi,v) !,wkb)

      ! Planos xz ! ! Planos xz ! ! Planos xz !
      
      nbuffou = mx*mz*nvecy
      nbuffis = (mbx+2)*mbz
      
      ! -------------------- allocates buffers -------------------------!

      allocate(v (nbuffou),pr(nbuffou))
      allocate(phi (mx,mz,1),ps(nbuffou))
      allocate(wk1(nbuffis))
      allocate(wk(0:mx1,0:mz1))

      do j=0,mz1
         do i=0,mx1
            wk(i,j)=1.0/(alp2(i) + bet2(j))
         enddo
      enddo
      wk(0,0)=0.

      do j=1,my
         if (mod(j,10)==0) write(*,*) ' Computing plane', j,'of',my


         call getswapyz2xz(v     ,pr,pr,j,1)

c ---    reads a ix=constant plane

      fj=0
      if (ind.le.nspec) then
!         write(*,*)'chequeando',ind,j,jspecy(ind)
         if (jspecy(ind)==j) then
            fj=1
            ind=ind+1
         endif
      endif
      
!      write(*,*) 'entrando escrphys',j,u00(j)
      call escru1phys(v,wk1,wk1,j,fj)
!      write(*,*) 'saliendo de escrphys',j,u00(j)

      enddo
      
      deallocate(v,pr,ps,phi,wk1,wk)
!      deallocate(vor,v,dvdy,pr,ps,phi,dvordy,wk1,wk2,wk3)

!********
      enddo !loop for each file 

      end

!***********************************************************************!
!     Subroutine initcr                                                 !
!                                                                       !
!     Initializes commons                                               !
!     Adapted, SHC 20-12-05                                             !
!                                                                       !
!***********************************************************************!

      subroutine initcr(ifiles,nfiles)
      
      implicit none 
      include "ctes3D"

!     ------------------------- COmmons --------------------------------!

      real*4  Re,alp,bet,a0
      real*8  y,fmap
      common /fis/ Re,alp,bet,a0,y(my),fmap(my)
      save   /fis/

      integer nvecy,plx
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr,nvecy,plx
      save /tem/
      
      integer iinp,iswap
      character*80 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(8)
      save /ficheros/
     
!      real*8 d21,d22,dt21
!      common /secder/ d21(my,5),d22(my,5),dt21(5,my)
!      save   /secder/
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

!      real*8 prem3,dt22
!      common /d2y/ prem3(my,5),dt22(5,my)
!      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(my,7),dt12(7,my)
      save   /cfdiff/
      
      integer icx
      real*4 alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     &              icx(0:mz1)
      save  /wave/
      
      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,prr,prs,pso,
     &       ptt,Wx0a,Wz0a,uvr,uwr,vwr
      integer nacum
      common /stat1/ um(my), vm(my), wm(my), up(my), vp(my), wp(my),
     &       w1m(my),w2m(my),w3m(my),w1p(my),w2p(my),w3p(my),uvr(my),
     &  uwr(my),vwr(my),prr(my),prs(my),pso(my),ptt(my),Wx0a,Wz0a,nacum
      save /stat1/


      real*8 prod,disp,turb,visd,pstr,pdif
      common /tensors/ prod(my,6),disp(my,6),turb(my,6),
     &                 visd(my,6),pstr(my,6),pdif(my,6)
      save   /tensors/      

      real*8 pret,vi2d,turr
      common /nten/ pret(my,6),vi2d(my,6),turr(my,6)
      save   /nten/

      

!     ----------------- Workspaces -------------------------------------!

      real*8, allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      integer idat(20),mybu,mzbu,i,j,iproc,jj,istart,ifor,idif,k,
     &        mxe,mye,mze,ifiles,nfiles
      
!      real*4 dat(20),zero,Delt,u00(*),w00(*),wk2(2*my)
      real*4 dat(20),zero,Delt,wk2(2*my)
      character*80 text
      
!     -----------------PROGRAMA ----------------------------------------!
      
      zero = 0e0
      
!     ---------------- reads in data  ----------------------------------!
      
      !open hre.dat if first file
      if (ifiles.eq.1) then
         open(19,file='hre.dat',status='old')
      endif


      
64    read(19,'(a)') text
      if(text(1:2).eq.'CC') goto 64
      read(text,*) (idat(i),i=1,2)

      nvecy = idat(1)
      plx   = idat(2)
      

65    read(19,'(a)') text
      if(text(1:2).eq.'CC') goto 65
      read(text,'(a80)') filswap

166   read(19,'(a)') text
      if(text(1:2).eq.'CC') goto 166
      read(text,'(a80)') filinp

167   read(19,'(a)') text
      if(text(1:2).eq.'CC') goto 167
      read(text,'(a80)') filout

      iinp=1
      iswap=2
      if (ifiles.eq.nfiles) then
         close(19)
      endif
      
      extensiones(1) = 'scal'
                       
      ! Read head of file
      
      write(*,*) filinp,filout
      
      open(iinp,file=filinp,status='old',form='unformatted')
      write(*,*) 'archivo abierto' 

      read(iinp) time,Re,alp,bet,a0,mxe,mye,mze,
     &             (y(j),fmap(j),j=1,my),(wk2(j),j=1,2*my)
     
      if(mx .ne.mxe . OR .
     &   my .ne.mye . OR .
     &   mz .ne.mze       ) then
         write(*,*) 'parameters changed'
         write(*,*) 'mx=',mx,mxe
         write(*,*) 'my=',my,mye
         write(*,*) 'mz=',mz,mze
         write(*,*) 'MX/Y/Z not matching MXE/MYE/MZE,TRY AGAIN :-)'
         stop
      endif
     
     
      write(*,*) 'inicializamos' 

c    ---------  initializes fast fourier transforms ----

      call cfti(mgalz)
      call rfti(mgalx)

c    -----------  write header for output -------------

      write(*,'(a7,f8.2,a8,f6.3,a8,f6.3)') '  Re =',Re,'channel'
      write(*,'(a7,f8.3,a8,f6.3,a8,f6.3)')
     .              'alp =',alp,'  bet =',bet
      write(*,*)

      write(*,'(a8,i5,a8,i5,a8,i5)')
     .              'mgalx =',mgalx,'mgalz =',mgalz,'my =',my
      write(*,'(a8,i5,a8,i5,a8,i5)')
     .              'mx =',mx,'mz =',mz
      write(*,*)

      write(*,'(a10,i5,a8,i5)')
     .        'y planes =',nvecy,' xplanes',plx
      write(*,*)
      write(*,'(a,a)') 'reading from:  ',filinp
      write(*,*)
      write(*,'(a,a)') '  swap in :  ',filswap
      write(*,*)
      write(*,'(a,a)') '  out in :  ',filout
!AAF----------------------------------------------------

!      allocate(d11(my,7),d12(my,7))
      allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))

      call derivadas(d11,d12,d21,d22)          

      deallocate(d11,d12,d21,d22)
      write(*,*) 'Fin precalculo de matrices...'

c ------------  compute y coordinates, pointers and
c ------------  modes values for FOURIER      

      do k=0,nz1
         xbet(k) = cmplx(zero,bet*k)
         icx(k) = k
      enddo

      do k=nz1+1,mz1
         xbet(k) = cmplx(zero ,-bet*(mz1+1-k))
      enddo

      do k=1,nz1
         icx(mz-k) = k
      enddo
      
      do i=0,mx1
         xalp(i) = cmplx(zero ,alp*i)
      enddo

      do i=0,mx1
         alp2(i) = -xalp(i)**2
      enddo

      do j=0,mz1
         bet2(j) = -xbet(j)**2
      enddo

      do j=1,my
         um(j)  = 0d0
         vm(j)  = 0d0
         wm(j)  = 0d0
         up(j)  = 0d0
         vp(j)  = 0d0
         wp(j)  = 0d0
         w1m(j) = 0d0
         w2m(j) = 0d0
         w3m(j) = 0d0
         w1p(j) = 0d0
         w2p(j) = 0d0
         w3p(j) = 0d0
         prr(j) = 0d0
         prs(j) = 0d0
         pso(j) = 0d0
         ptt(j) = 0d0
         uvr(j) = 0d0
         uwr(j) = 0d0
         vwr(j) = 0d0
      enddo
      do k=1,6
         do j=1,my
            prod(j,k)= 0d0
            disp(j,k)= 0d0
            turb(j,k)= 0d0
            visd(j,k)= 0d0
            pstr(j,k)= 0d0
            pdif(j,k)= 0d0
            pret(j,k)= 0d0
            turr(j,k)= 0d0
         enddo
      enddo

      end
