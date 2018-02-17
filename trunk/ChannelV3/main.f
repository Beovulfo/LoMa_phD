!***********************************************************************!
!                                                                       !
!      Resuelve las ecuaciones de N-S incompresibles para una           !
!    canal tridimensional con una formulacion basada en                 !
!    la vorticidad y la laplaciana de la velocidad en la direccion      !
!    normal a la pared, Kim, Moin y Moser (1987).                       !
!                                                                       !
!       tocado j.jimenez (1/92)                                         !
!       hundido j.c.a. (1/01)                                           !
!       reflotado, version plano-recta shc (8/04)                       !
!                                                                       !
!                                                                       !
!       head =(time,Re,alp,bet,a0,mx,my,mz)                             !
!                                                                       !
!     REFERENCIAS                                                       !
!       - Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statis         !
!         tics in fully developped channel flow at low Reynolds         !
!         numbers, J. Fluid Mech., 177, 133-166                         !
!                                                                       !
!     ATENCION                                                          !
!          ctes3D modificado para incluir los nuevos grupos de MPI      !
!                                                                       !
!                                                                       !
!***********************************************************************!
      Program main
 
      implicit none 
 
      include "mpif.h"
      include "ctes3D"
      
!     ----------------------- Commons ----------------------------------
      
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
      integer nacumsp,jsptot
      common /spectra/   nacumsp,jsptot(2*nspec+1)
      save   /spectra/

      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
                 
!      ----------------------- Variables -------------------------------

      integer master,myid,iproc,numprocs,istat(MPI_STATUS_SIZE),ierr  
        
      integer nbuffsize,nwkasize
      integer ihv,ihg,iphiwk,ivorwk,irf0u,irf0w,iu00wk,iw00wk,idvordy,
     .        iuwk,ichwk,nstr,nspsize,iflag
     
      integer itags,newtag,imess,i,comproc(0:numerop-1),
     .                      savproc(numerop:numtot-1)      
      
      real*8  u00(my),w00(my),val 
      real*4, allocatable:: vor(:),phi(:),wk(:),sp(:)

!     ------------------------ Program ---------------------------------
      
                           ! MPI GROUPS !     

      call MPI_INIT(ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,numprocs,ierr)
      call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
      
      do i=0,numerop-1
         comproc(i) = i
      enddo
      do i=numerop,numtot-1
         savproc(i) = i
      enddo

      if (myid.lt.numerop) then 
          call MPI_Group_incl(MPI_GROUP_WORLD,numerop,comproc,
     .                        MPI_GROUP_CALC,ierr)
          call MPI_Group_excl(MPI_GROUP_WORLD,numerop,comproc,
     .                        MPI_GROUP_SAVE,ierr)
      else
          call MPI_Group_incl(MPI_GROUP_WORLD,numerosa,savproc,
     .                        MPI_GROUP_SAVE,ierr)
     
          call MPI_Group_excl(MPI_GROUP_WORLD,numerosa,savproc,
     .                        MPI_GROUP_CALC,ierr)   
      endif

      call MPI_COMM_CREATE (MPI_COMM_WORLD,
     .                       MPI_GROUP_CALC,MPI_COMM_CALC,ierr)
     
      call MPI_COMM_CREATE (MPI_COMM_WORLD,
     .                       MPI_GROUP_SAVE,MPI_COMM_SAVE,ierr)

      if (numprocs.ne.numtot) then
         write(*,*) 'wrong number of processors',numprocs
         write(*,*) 'compiled with', numtot
         stop
      endif

! ------------- be careful with the numbers of planes

      if (mod(nplanes,numerop).ne.0) then
         write(*,*) 'the number of planes',nplanes
         write(*,*) 'must be multiple of processors,',numerop
         stop
      endif
      
      write(*,*) 'MPI enviroment initialized ...', myid

!--------------- initializes commons and things 

      call initcr(myid)
      
!--------------  allocates spectra and buffers
    
      nspsize   = (2*nspec+1)*(nz1+1)*8*(pe-pb+1)
      nbuffsize = max((pe-pb+1)*2*my*mgalz,(le-lb+1)*mx)
      nwkasize  = 7*nbuffsize + 8*my !+ 9*nbuffsize
      allocate(sp(nspsize))
             
      if (myid.lt.numerop) then     ! Comp. Proc.
      
         allocate(vor(2*my*mz*(pe-pb+1)))
         allocate(phi(2*my*mz*(pe-pb+1)))
         allocate(wk(nwkasize))

      else                          ! Save proc
          nwkasize = 4*mmp*my*mz+1000
          allocate(wk(nwkasize))
          
      endif  
       
             ! zero everything, fitf !
     
      nacumsp = 0
      do i=1,nspsize
         sp(i)=0.
      enddo

      if (myid.lt.numerop) then
      
         do i=1,2*my*mz*mmp 
            vor(i) = 0. 
            phi(i) = 0. 
         enddo              

         do i=1,my          
            u00(i) = 0d0     
            w00(i) = 0d0  
         enddo              


          call getfil(vor,phi,u00,w00,wk,myid)

      endif
      

          ! start time advancement  !
 
      irf0u   = 1
      irf0w   = irf0u   + 2 * my  ! Real*8
      iu00wk  = irf0w   + 2 * my
      iw00wk  = iu00wk  + 2 * my
      ihv     = iw00wk  + 2 * my
      ihg     = ihv     + nbuffsize
      iphiwk  = ihg     + nbuffsize
      ivorwk  = iphiwk  + nbuffsize
      idvordy = ivorwk  + nbuffsize
      iuwk    = idvordy + nbuffsize
      ichwk   = iuwk    + nbuffsize
      
      
      if (myid.lt.numerop) then

      call cross1(vor,
     .            phi,
     .            u00,w00,
     .            wk(irf0u),wk(irf0w),wk(iu00wk),wk(iw00wk),
     .            wk(ihv),
     .            wk(ihg),
     .            wk(iphiwk),
     .            wk(idvordy),
     .            wk(ivorwk),
     .            wk(idvordy),
     .            wk(iuwk),
     .            wk(ichwk),
     .            sp,myid)
      else
         call cross1(wk(1),
     .               wk(nwkasize/2),
     .               u00,w00,
     .               wk(1),wk(1),wk(1),wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               sp,myid)
      endif      

                ! finalize procedure      !
                
                
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)     

      master=0
      if (myid.ne.0) then

         itags=myid
         call MPI_SEND(1.,1,MPI_REAL,master,
     &              itags,MPI_COMM_WORLD,ierr)

      else

         do iproc=1,numprocs-1

            call MPI_RECV(val,1,MPI_REAL,MPI_ANY_SOURCE,
     &              MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

            imess=istat(MPI_TAG)
            write(*,*) 'process',imess,'over'
            newtag=100
            call MPI_SEND(1.,1,MPI_REAL,imess,
     &                newtag,MPI_COMM_WORLD,ierr)
         enddo

      endif

      if (myid.ne.0) then

         call MPI_RECV(val,1,MPI_REAL,master,
     &              MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         if(istat(MPI_TAG).eq.100) goto 200

      else

         write(*,*) 'process',master,'over'

      endif

200   call MPI_FINALIZE(ierr)

      endprogram

!***********************************************************************!
!                   get data field from file                            !
!                                             jjs  27/1/01              !
!                                             shc  05/1/05              !
!***********************************************************************!
      subroutine getfil(vor,phi,u00,w00,wk1,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

!  --------------------- commons ---------------------------------------

      real*4 Deltat,CFL,time,dtr 
      common /tem/ Deltat,CFL,time,dtr 
      save   /tem/
      
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
      
      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt
      save   /ficheros/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/
      
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/


      
!   -------------------- variables -------------------------------------
     
      integer myid,master,istat(MPI_STATUS_SIZE),ierr       
      real*4 vor(2*my,mz,*),phi(2*my,mz,*),wk1(*)
      real*8 u00(*),w00(*)
      
      real*4  Ree,alpe,bete,a0e

      integer i,j,mxe,mye,mze,ntotr,iproc,mym,mmy2
      

!   --------------------- Program -------------------------------------

      master = 0

!      ---------------  read from file and distribute

      if (myid.eq.master) then     

         open(iinp,file=filinp,status='unknown',form='unformatted')
         
         read(iinp) time,Ree,alpe,bete,a0e,mxe,mye,mze,
     .               (y(j),fmap(j),j=1,my),(wk1(j),j=1,2*my)

       
         write(*,*) 
         write(*,*) 'reading input file ...' 
         write(*,*) 'time=',time,Ree,alpe,bete,mxe,mye,mze
         write(*,*)

         write(*,*) 'valores para vhange',my*mgalz*mmp


         ntotr=4*mye*mze
!     ------------- 00 modes ------------------

         mym = min(my,mye)
         do j=1,mym
            u00(j) = wk1(2*j-1)
            w00(j) = wk1(2*j)
!            write(*,*) u00(j), w00(j) 
         enddo
        
              
         

!     ------------ 00 modes only for master 

         do iproc=1,numerop-1
            call MPI_SEND(mxe,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mye,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mze,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
         enddo

c         ------------  other modes,  master node -----
         
         write(*,*) 'master reads its data',ntotr
              

         do i=1,(pe-pb+1)
            read(iinp) (wk1(j),j=1,ntotr)
            call assign(wk1,vor(1,1,i),phi(1,1,i),my,mz,mye,mze,myid)
         enddo


!         ------------  other modes,  distribute to slaves --
 
         do iproc=1,numerop-1
           mmy2 = min(pend(iproc),mxe/2)-pbeg(iproc)+1
           write (*,*) mmy2,pend(iproc),mxe/2
           if (mmy2.gt.0) then
              do i=1,mmy2
                 read(iinp) (wk1(j),j=1,ntotr)
                  call MPI_SEND(wk1,ntotr,MPI_REAL,
     &                          iproc,iproc,MPI_COMM_WORLD,ierr)
              enddo
           write(*,*) 'master reads proc no',iproc,' data and send them'
           endif
         enddo

         close(iinp)

      else          ! -------- this is done by all slaves

         call MPI_RECV(mxe,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mye,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mze,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

          ntotr=4*mye*mze
c         ------------  receive other modes ----------
       
         mmy2 = min(pe,mxe/2)-pb+1
         if (mmy2.gt.0) then
            do i=1,mmy2
               call MPI_RECV(wk1,ntotr,MPI_REAL,master,
     &                       MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
               call assign(wk1,vor(1,1,i),phi(1,1,i),my,mz,mye,mze,myid)
            enddo
            write(*,*) 'proc no.',myid,'receives data from master'
         endif
      endif
      
      

      endsubroutine
!***********************************************************************!
!
!                   Initializes everything
!                                           
!***********************************************************************!

      subroutine initcr(myid)

      implicit none 
      include "mpif.h"
      include "ctes3D"
      
!     ----------------------- Commons ----------------------------------

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

      integer nimag,nstep,nhist,ihist,icfl,ncfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl,ncfl
      save   /timacc/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt
      save /ficheros/

      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
      integer myslice,myplane,plbypl,req
      common /MPI_datatype/myslice(0:numerop-1),myplane(0:numerop-1),
     .                     plbypl (0:numerop-1),req(0:2*numerop-1)
      save /MPI_datatype/

      integer nacumsp,jsptot
      common/spectra/   nacumsp,jsptot(2*nspec+1)
      save  /spectra/

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/

      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save   /mesh/

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/
      
! ---------------------- common workspaces -----------------------------
      real*8 wk1,phipr,phipi,vpr,vpi,dvpr,dvpi,phi1,v1,dv1,dfr,dfi,
     .       dphi1,dphi2,fwk1,fwk2
      common /worksp/ wk1(5,my),fwk1(my),fwk2(my),dphi1(my),dphi2(my),
     .                phipr(my),phipi(my),vpr(my),vpi(my),
     .    dvpr(my),dvpi(my),phi1(my),v1(my),dv1(my),dfr(my),dfi(my)



      real*4 up1wk,up2wk,up3wk,wp1wk,wp2wk,wp3wk
      real*8 up1wk8,up2wk8,up3wk8,wp1wk8,wp2wk8,wp3wk8
      common /wkhvect/ up1wk(0:mgalx+1),up2wk(0:mgalx+1),
     .                 up3wk(0:mgalx+1),wp1wk(0:mgalx+1),
     .                 wp2wk(0:mgalx+1),wp3wk(0:mgalx+1),
     .                 up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .                 up3wk8(0:mgalx+1),wp1wk8(0:mgalx+1),
     .                 wp2wk8(0:mgalx+1),wp3wk8(0:mgalx+1)

! ---------------------- variables -------------------------------------

      integer istat(MPI_STATUS_SIZE),ierr,myid,idat(30)
      integer i,j,iproc,k,ifor,blocks,elem,stride,mxe,mye,mze

      real*8, allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      integer,allocatable::aux1(:),aux2(:)

      real*4 dat(20),pi,Delt,zero,Ree,alpe,bete,a0e

      character*80 text

c ---------------------- Programa --------------------------------------

      Pi=4.*atan(1.)
      zero = 0.0
c                ! reads in data !

      if(myid.eq.0) then
         open(19,file='hre.dat',status='old')

965      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 965
         read(text,*) (dat(j),j=1,4)

966      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 966
         read(text,*) (idat(j),j=5,8)

964      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 964
         read(text,*) idat(1),dat(7)

967      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 967
         read(text,*) dat(5),dat(6), idat(10)

968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) (idat(j), j=11,14)

969      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 969
         read(text,*) idat(2),idat(3),dat(8),dat(9),dat(10),dat(11)

65       read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 65
         read(text,'(a70)') filout

166      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 166
         read(text,'(a70)') filinp

167      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 167
         read(text,'(a70)') filstt

         close(19)
         do iproc=1,numtot-1
            call MPI_SEND(dat,11,MPI_REAL,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(idat,19,MPI_INTEGER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
         enddo
      else
         call MPI_RECV(dat,11,MPI_REAL,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(idat,19,MPI_INTEGER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Re   =  dat(1)
      alp  = dat(2)
      bet  = dat(3)
      a0   = dat(4) 
      Deltat = dat(5)           ! Este parametro ya no sirve      
      cfl  = dat(6)
      gamma= dat(7)

      imesh=idat(1)
      nstep=idat(5)
      nimag=idat(6)
      nhist=idat(7)
      ncfl = idat(10)   
      ifor=idat(11)      ! Este parametro ya no sirve
      id22=idat(12)
      nstart=idat(13)
      ntimes=idat(14)
      
      iinp=2
      ispf=15
      istati=0
      isn=3
      iout=1
 
c ---------------  compute y coordinates, pointers and
c ---------------  modes values for FOURIER

      if (myid.eq.0) then
         open(iinp,file=filinp,status='unknown',form='unformatted')
         read(iinp) time,Ree,alpe,bete,a0e,mxe,mye,mze,
     .               (y(j),fmap(j),j=1,my)
         rewind(iinp)
         close(iinp)
         do iproc=1,numtot-1
            call MPI_SEND(y,my,MPI_REAL8,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(fmap,my,MPI_REAL8,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
         enddo
      else
         call MPI_RECV(y,my,MPI_REAL8,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(fmap,my,MPI_REAL8,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
      
      endif
        

!      call malla8(my,y,fmap)

      if (myid.lt.numerop) then
      
         call pointers_calc(pbeg,pend,lbeg,lend,procs,myid)
         pb=pbeg(myid)
         pe=pend(myid)
         lb=lbeg(myid)
         le=lend(myid)
         mmp = pe-pb+1
         mml = le-lb+1
      else
         call pointers_save(pbeg,pend)   
         pb=pbeg(myid-numerop)
         pe=pend(myid-numerop)         
         mmp = pe-pb+1
         lb = 0
         le = 3 
      endif
      
      if (myid.eq.numerop) then

         call MPI_SEND(pbeg,numerop,MPI_INTEGER,0,
     &                    0,MPI_COMM_WORLD,ierr)
         call MPI_SEND(pend,numerop,MPI_INTEGER,0,
     &                    0,MPI_COMM_WORLD,ierr)
                
      elseif (myid.eq.0) then
         allocate(aux1(0:numerop-1),aux2(0:numerop-1))
         
         call MPI_RECV(aux1,numerop,MPI_INTEGER,numerop,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     
         call MPI_RECV(aux2,numerop,MPI_INTEGER,numerop,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr) 
      endif
      
      
      if (myid.eq.0) then
         j=0
         do i=1,nplanes
            if (i.lt.pend(j)+1) then
                plcalc(i) = j
            else
                j=j+1
                plcalc(i) = j
            endif
         enddo
         j=0
         do i=1,nplanes
            if (i.lt.aux2(j)+1) then
                plsav(i) = j+numerop
            else
                j=j+1
                plsav(i) = j+numerop
            endif
         enddo
      endif
      
      call MPI_BCAST(filout,70,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plcalc,nplanes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plsav ,nplanes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      if (myid.gt.numerop) return
          
      

c    ------------  initializes fast fourier transforms and CFDiff ----
      call cfti(mgalz)
c      call initfft()
      call rfti(mgalx)

      write (*,*) 'precalculo matrices ...',myid,y(2),y(my-1)
      
      allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))
      
      call derivadas(d11,d12,d21,d22)      
            
      deallocate(d11,d12,d21,d22)
      
      write (*,*) 'FIN precalculo matrices ...',myid
      
!    --------------  coefficients for trapz averaging ------

      do j=2,my-1
         trp(j)= 0.25d0*(y(j+1)-y(j-1))
         hy(j) = (y(j+1)-y(j-1))/2.5d0/2d0
      enddo
      trp(1) = 0.25d0*(y(2) -y(1))
      trp(my)= 0.25d0*(y(my)-y(my-1))
      hy(1)  = (y(2)-y(1))/2.5d0
      hy(my) = (y(my)-y(my-1))/2.5d0

      call premass(trp,mss)                

c   --------------  prepare spectra -----

      do j=1,nspec
         jsptot(j)           =    jspecy(j)
         jsptot(2*nspec+2-j) = my-jspecy(j)+1
      enddo
      jsptot(nspec+1) = (my+1)/2

c ------------------ Re/dt --------------------------------
      dtr = Re/Delt

c ------------------ MPI Datatypes ------------------------

      blocks = le-lb+1
      elem   = 2*(pe-pb+1)
      stride = mx
      do iproc=0,numerop-1

         if (iproc.ne.myid) then

             call MPI_TYPE_VECTOR(blocks,elem,
     .                            stride,MPI_REAL,myslice(iproc),
     .                            ierr)

             call MPI_TYPE_COMMIT(myslice(iproc),ierr)

         endif
      enddo
      
c ---------------------

      blocks = (pe-pb+1)
      elem   = 2*(le-lb+1)
      stride = 2*my*mgalz

      do iproc=0,numerop-1

         if (iproc.ne.myid) then
             elem = 2* (lend(iproc)-lbeg(iproc)+1)
             call MPI_TYPE_VECTOR(blocks,elem,
     .                            stride,MPI_REAL,myplane(iproc),
     .                            ierr)

             call MPI_TYPE_COMMIT(myplane(iproc),ierr)

         endif
      enddo

c ----------------------

c ------ slide plane by plane.

      blocks = (le-lb+1)
      elem   = 2
      stride = mx

      do iproc=0,numerop-1

         if (iproc.ne.myid) then
             call MPI_TYPE_VECTOR(blocks,elem,
     .                            stride,MPI_REAL,plbypl(iproc),
     .                            ierr)

             call MPI_TYPE_COMMIT(plbypl(iproc),ierr)

         endif
      enddo


c --------------  initialize stats -------------

      do j=1,my
           um(j)  = 0.
           vm(j)  = 0.
           wm(j)  = 0.
           up(j)  = 0.
           vp(j)  = 0.
           wp(j)  = 0.
           uvr(j) = 0.
           uwr(j) = 0.
           vwr(j) = 0.
           w1m(j) = 0.
           w2m(j) = 0.
           w3m(j) = 0.
           w1p(j) = 0.
           w2p(j) = 0.
           w3p(j) = 0.
           ep(j)  = 0.
           uuv(j) = 0.
           wwv(j) = 0.
           vvv(j) = 0.
        enddo


c ------------------- Computes wavenumbers ---------------


      do k=0,nz1
         xbet(k) = cmplx(zero,bet*k)
         rbet(K) = bet*k
         icx(k) = k
      enddo

      do k=nz1+1,mz1
         xbet(k) = cmplx(zero ,-bet*(mz1+1-k))
         rbet(k) = -bet*(mz1+1-k)
      enddo

      do k=1,nz1
         icx(mz-k) = k
      enddo

      do i=0,mx1
         iax(2*i+1) = i
         iax(2*i+2) = i
      enddo

      do i=0,mx1
         xalp(i) = cmplx(zero ,alp*i)
         ralp(i) = alp*i
      enddo

      do i=0,mx1
         alp2(i) = -xalp(i)**2
      enddo

      do j=0,mz1
         bet2(j) = -xbet(j)**2
      enddo        

c --------------  write header for output -------------


      if(myid.eq.0) then
         write(*,'(a7,f8.2,a8,f6.3,a8,f6.3)')
     .                    '  Re =',Re,'channel'
         write(*,'(a7,f8.3,a8,f6.3,a8,f6.3)')
     .                    'alp =',alp,'  bet =',bet
         write(*,*)

         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mgalx =',mgalx,'mgalz =',mgalz,'my =',my
         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mx =',mx,'mz =',mz
         write(*,*)

         write(*,'(a10,i6,a9,i6,a9,i5)')
     .     'nstep =',nstep,'nimag =',nimag,'nhist =',nhist
         write(*,'(a10,e10.4,a8,f5.2)')
     .     'Delt =',Delt,'  CFL =',CFL
         write(*,*)
         write(*,'(a,a)')
     .     'reading from:  ',filinp
         write(*,*)
         write(*,'(a,a)')
     .     '  write in :  ',filout
      endif

      endsubroutine
       
!***********************************************************************!
!
!     Routine assign.
!     Distribute data from (mye,mze)-mesh to (my,mz)-mesh 
!
!          single    jjs 4/01/01, 
!          rewritten jjs 28/01/01
!          adapted   shc 5/09/04
!***********************************************************************!
      subroutine assign(work,vor,phi,my,mz,mye,mze)

      implicit none

      integer   mym,mzm,klen,kini1,kini2,i,k,k1,k2,my,mz,mye,mze
      real*4    vor(2*my,mz),phi(2*my,mz)
      real*4    work(2,2*mye,mze)
      
      
      mym = min(my,mye)
      mzm = min(mz,mze)

      klen = (mzm+1)/2
      kini1=mze - klen + 1
      kini2=mz - klen + 1
      
      do k=1,klen
         do i=1,2*my
            vor(i,k)=work(1,i,k)
            phi(i,k)=work(2,i,k)
         enddo
      enddo

      do k=1,klen-1
         k1 = k + kini1
         k2 = k + kini2
         do i=1,2*my
            vor(i,k2)=work(1,i,k1)
            phi(i,k2)=work(2,i,k1)
         enddo
      enddo

      endsubroutine


!***********************************************************************!
!
!     Routine initNaN
!     Restart the run in case of NAN
!
!          adapted   shc 30/05/05
!***********************************************************************!


      subroutine initnan(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,hv,hg,
     .                  phiwk,spwk,vorwk,dvordy,uwk,chwk,sp,myid)
      include "ctes3D"

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
     
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
            
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
     
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,nanerror,kk
      
!
      write(*,*) ' entre aqui'
      
!     Reinitialize everything

      
      do i=pb,pe
         do k=0,mgalz-1
            do j=0,my-1
               hg    (j,k,i) = 0.0
               vorwk (j,k,i) = 0.0
               hv    (j,k,i) = 0.0
               phiwk (j,k,i) = 0.0
               uwk   (j,k,i) = 0.0
               dvordy(j,k,i) = 0.0
               chwk  (j,k,i) = 0.0
            enddo
         enddo
      enddo
      do i=pb,pe
         do k=0,mz-1
            do j=0,my-1 
               phi(j,k,i) = 0.0
               vor(j,k,i) = 0.0
            enddo
         enddo
      enddo
      do j=0,my-1 
         u00     (j) = 0d0
         w00     (j) = 0d0
         rf0u    (j) = 0d0
         rf0w    (j) = 0d0
         u00wk   (j) = 0d0
         rf0w    (j) = 0d0
      enddo      
      
      do i=pb,pe
         do kk=1,7
            do j=1,2*nspec+1
               do k=0,nz-1
                  sp  (k,j,kk,i) = 0.0
               enddo
            enddo
         enddo
      enddo
      
!     Statisitics

      do j=1,my
           um(j)  = 0.
           vm(j)  = 0.
           wm(j)  = 0.
           up(j)  = 0.
           vp(j)  = 0.
           wp(j)  = 0.
           uvr(j) = 0.
           uwr(j) = 0.
           vwr(j) = 0.
           w1m(j) = 0.
           w2m(j) = 0.
           w3m(j) = 0.
           w1p(j) = 0.
           w2p(j) = 0.
           w3p(j) = 0.
           ep(j)  = 0.
           uuv(j) = 0.
           wwv(j) = 0.
           vvv(j) = 0.
      enddo
      
!     Call assign

      call getfil(vor,phi,u00,w00,chwk,myid)
      
      
      
!     call Cross Again

      if (myid.lt.numerop) then

      call cross1(vor,
     .            phi,
     .            u00,w00,
     .            rf0u,rf0w,u00wk,w00wk,
     .            hv,
     .            hg,
     .            phiwk,
     .            spwk,
     .            vorwk,
     .            dvordy,
     .            uwk,
     .            chwk,
     .            sp,myid)

      else
         call cross1(vor,
     .               phi,
     .               u00,w00,
     .               vor,vor,vor,vor,
     .               vor,
     .               vor,
     .               vor,
     .               vor,
     .               vor,
     .               vor,
     .               vor,
     .               vor,
     .               sp,myid)
      endif
      
      endsubroutine
      
      subroutine check_nan(u,n,flag)
      implicit none
      
      real*4 u(n)
      integer j,i,n,flag
      
      do j=1,n
         if (u(j).ne.u(j)) flag = 1
      enddo
      end 


