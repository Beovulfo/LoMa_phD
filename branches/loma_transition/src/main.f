***********************************************************************!
!                      L o M A    D N S                                 !
!-----------------------------------------------------------------------!
!  Solves the incompressible N-S equations for a 3D box domain with     !
!  periodic BC @ "x,z" and general BC's for "y".                        !
!  Spectral method for x,z Compact Finite Differences for y.            ! 
!  Formulation based on y-Vorticity and Laplacian of v, Kim,Moin and    !
!  Moser (1987). Based on LISO code. EXPLICIT SOLVER                    !
!                                                                       ! 
!       LISO:                                                           !
!       --------------------------------------------                    !
!       tocado j.jimenez (1/92)                                         !
!       hundido j.c.a. (1/01)                                           !
!       reflotado, version plano-recta shc (8/04)                       !
!                                                                       !
!       LOMA:                                                           !
!       ----------------------------------------                        !
!       changed to work with Robin BC ofa & aaf (06/2013)               !
!                                                                       ! 
!              head =(time,Re,alp,bet,a0,mx,my,mz)                      !
!                                                                       !
!     REFERENCES                                                        !
!       - Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statis         !
!         tics in fully developped channel flow at low Reynolds         !
!         numbers, J. Fluid Mech., 177, 133-166                         !
!                                                                       !
!       - Spalart,R., Moser,R. and Rogers,M.(1990) Spectral Methods     !
!         for the Navier-Stokes Equations with One Infinite and         !
!         two Periodic Directions,Journal of Comp. Physics 96, 297-324  !
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
                 
!      ----------------------- Variables -------------------------------

      integer master,myid,iproc,numprocs,istat(MPI_STATUS_SIZE),ierr  
        
      integer nbuffsize,nwkasize
      integer inxwk,inywk,inzwk,
     .        imxwk,imywk,imzwk,
     .        iten12,iten13,iten23,
     .        irf0u,irf0w,iu00wk,iw00wk,
     .        irf0v,iv00wk,
     .        ichwk,ichwk2,irhst,idrho,nstr,nspsize,iflag
     
      integer itags,newtag,imess,i,comproc(0:numerop-1),
     .                      savproc(numerop:numtot-1)      
      
      real*8  u00(my),w00(my)
      real*4  val 
      real*8  v00(my) 
!      !real*4, allocatable:: vor(:),phi(:),wk(:),sp(:),
      real*4, allocatable:: vor(:),phi(:),wk(:),
     .                      psi(:),scal(:)

!     ------------------------ Program ---------------------------------
      
                           ! MPI GROUPS !     

      call MPI_INIT(ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,numprocs,ierr)
      call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
!     distribute groups comp/save      
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
      nbuffsize = max((pe-pb+1)*2*my*mgalz,(le-lb+1)*mx)
      nwkasize  = 7*nbuffsize + 12*my + 6*nbuffsize
!added v00 and v00wk,rf0v
!      nwkasize  = 7*nbuffsize + 8*my + 6*nbuffsize
!      allocate(sp(nspsize))
             
      if (myid.lt.numerop) then     ! Comp. Proc.
      
         allocate(vor(2*my*mz*(pe-pb+1)))
         allocate(phi(2*my*mz*(pe-pb+1)))
         allocate(psi(2*my*mz*(pe-pb+1)))
         allocate(scal(2*my*mz*(pe-pb+1)))
         allocate(wk(nwkasize))

      else                          ! Save proc
          nwkasize = 4*mmp*my*mz + 1000
          !nwkasize = 6*mmp*my*mz+1000
          allocate(wk(nwkasize))
          
      endif  
       
             ! make zero everything  
     

      if (myid.lt.numerop) then
      
         do i=1,2*my*mz*mmp 
            vor(i) = 0. 
            phi(i) = 0. 
            psi(i) = 0. 
            scal(i) = 0. 
         enddo              

         do i=1,my          
            u00(i) = 0d0     
            w00(i) = 0d0  
            v00(i) = 0d0  
         enddo              

! load data from file filinp (stated on hre.dat)
!adding scalar with  getscal
          call getfil(vor,phi,u00,w00,wk,myid)
!v00 not included yet,initially can be zero.
          call getscal(scal,wk,1,myid)
      endif
      

          ! start time advancement  !
 
      irf0u   = 1
      irf0w   = irf0u   + 2 * my  ! Real*8
      irf0v   = irf0w   + 2 * my  ! Real*8
      iu00wk  = irf0v   + 2 * my
      iw00wk  = iu00wk  + 2 * my  
      iv00wk  = iw00wk  + 2 * my  
      imxwk   = iv00wk  + 2 * my
      imywk   = imxwk   + nbuffsize
      imzwk   = imywk   + nbuffsize
      inxwk   = imzwk   + nbuffsize
      inywk   = inxwk   + nbuffsize
      inzwk   = inywk   + nbuffsize
      irhst   = inzwk   + nbuffsize
      idrho   = irhst   + nbuffsize
      iten12  = idrho   + nbuffsize
      iten13  = iten12  + nbuffsize
      iten23  = iten13  + nbuffsize
      ichwk2  = iten23  + nbuffsize
      ichwk   = ichwk2  + nbuffsize
 


      
      
      if (myid.lt.numerop) then

      call cross1(vor,
     .            phi,
     .            psi,
     .            scal,
     .            u00,w00,v00,
     .            wk(irf0u),wk(irf0w),wk(irf0v),
     .            wk(iu00wk),wk(iw00wk),wk(iv00wk),
     .            wk(imxwk),
     .            wk(imywk),
     .            wk(imzwk),
     .            wk(inxwk),
     .            wk(inywk),
     .            wk(inzwk),
     .            wk(iten12),
     .            wk(iten13),
     .            wk(iten23),
     .            wk(irhst),
     .            wk(idrho),
     .            wk(ichwk2),
     .            wk(ichwk),
     .            myid)
!     .            sp,myid)
      else
         call cross1(wk(1),
     .               wk(nwkasize/2),
     .               wk(1),
     .               wk(1),
     .               u00,w00,v00,
     .               wk(1),wk(1),wk(1),
     .               wk(1),wk(1),wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
!     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               wk(1),
     .               myid)
!     .               sp,myid)
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

!-------------------------------------------------!
!deallocate
      if (myid.lt.numerop) then     ! Comp. Proc.
         deallocate(vor)
         deallocate(phi)
         deallocate(psi)
         deallocate(scal)
         deallocate(wk)

      else                          ! Save proc
         deallocate(wk)
          
      endif  
!........................................! 

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
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save   /ficheros/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),
     .             fmap(my),trp(my),mss(my)
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
      
      real*8 fmape(my) 
      real*4  Ree,alpe,bete,a0e

      integer i,j,mxe,mye,mze,ntotr,iproc,mym,mmy2
      

!   --------------------- Program -------------------------------------

      master = 0

!      ---------------  read from file and distribute

      if (myid.eq.master) then     

         open(iinp,file=filinp,status='unknown',form='unformatted')
         
         read(iinp) time,Ree,alpe,bete,a0e,mxe,mye,mze,
     .               (y(j),fmape(j),j=1,my),(wk1(j),j=1,2*my)

         write(*,*) 
         write(*,*) 'reading input file ...' 
         write(*,*) 'time=',time,'Ree=',Ree,'alpe=',alpe,
     .   'bete=',bete,'mxe=',mxe,'mye=',mye,'mze =', mze
         write(*,*)

         write(*,*) 'valores para vhange',my*mgalz*mmp


         ntotr=4*mye*mze
!     ------------- 00 modes ------------------

         mym = min(my,mye)
         do j=1,mym
            u00(j) = wk1(2*j-1)
            w00(j) = wk1(2*j)
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
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save   /ficheros/


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


      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12,endder
      common /cfdiff/ prem1(7,my),dt12(7,my),endder(my)
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

!   !--------------------- variables -------------------------------------
      integer ndat, nidat,iscal
      parameter(ndat=5,nidat=7)

      integer istat(MPI_STATUS_SIZE),ierr,myid,idat(nidat)
      integer i,j,iproc,k,ifor,blocks,elem,stride,mxe,mye,mze

      real*8, allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      integer,allocatable::aux1(:),aux2(:)

      real*4 dat(ndat),pi,Delt,zero,Ree,alpe,bete,a0e

      character*80 text

c ---------------------- Programa --------------------------------------

      Pi=4.*atan(1.)
      zero = 0.0
c                ! reads input data from hre.dat !

      if(myid.eq.0) then
         open(19,file='hre.dat',status='old')

965      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 965
         read(text,*) (dat(j),j=1,4)
!         dat=[Re,alp,bet,-u]

966      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 966
         read(text,*) (idat(j),j=1,3)
!        idat(5:7)=[nstep,nimag,nhist]


967      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 967
         read(text,*) dat(5), idat(4)
!        dat(6)=CFL, idat(10)=ncfl

968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) (idat(j), j=5,7)
!        idat(12)=first output file
!        idat(13)=stats(0=no stats)
!        idat(14)=#steps between stats 

65       read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 65
         read(text,'(a70)') filout

166      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 166
         read(text,'(a70)') filinp

167      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 167
         read(text,'(a70)') filstt

168      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 168
         read(text,'(a70)') filscainp



         close(19)
         do iproc=1,numtot-1
            call MPI_SEND(dat,ndat,MPI_REAL,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(idat,nidat,MPI_INTEGER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
         enddo
      else
         call MPI_RECV(dat,ndat,MPI_REAL,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(idat,nidat,MPI_INTEGER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Re   =  dat(1)
      alp  = dat(2)
      bet  = dat(3)
      a0   = dat(4) 
      Deltat = 1e-6     ! this data is not used 
      cfl  = dat(5)
      nstep=idat(1)
      nimag=idat(2)
      nhist=idat(3)
      ncfl = idat(4)   
      id22=idat(5)
      
      iinp=2
      ispf=15 !aaf ??
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
!  save in plcalc(i) the processor "j" that is
!  taking care of the plane "i".
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
      call cfti(mgalz) !complex Fourier Transform init
c      call initfft()
      call rfti(mgalx) !real Fourier Transform init

      write (*,*) 'precalculo matrices ...',myid
      
      allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))
      
      call derivadas(d11,d12,d21,d22)      
            
      deallocate(d11,d12,d21,d22)
      
      write (*,*) 'FIN precalculo matrices ...',myid
      
!    --------------  coefficients for trapz averaging ------

      do j=2,my-1
         trp(j)= 0.25d0*(y(j+1)-y(j-1))
         hy(j) = (y(j+1)-y(j-1))/2d0
      enddo
      trp(1) = 0.25d0*(y(2) -y(1))
      trp(my)= 0.25d0*(y(my)-y(my-1))
      hy(1)  = (y(2)-y(1))/2d0
      hy(my) = (y(my)-y(my-1))/2d0

! mss not used anymore
      call premass(trp,mss)!prep matrix for 2nd derivative mean          

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
     .     'Delt =',Deltat,'  CFL =',CFL
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

      
      subroutine check_nan(u,n,flag)
      implicit none
!      
      real*4 u(n)
      integer j,i,n,flag
!      
      do j=1,n
         if (u(j).ne.u(j)) flag = 1
      enddo
      end 
!
!


!***********************************************************************!
!                   get scalar field from file                          !
!                   aaf 2014/01                                         !
!***********************************************************************!
      subroutine getscal(scal,wk1,iscal,myid)

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
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
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
      real*4 scal(2*my,mz,*),wk1(*)
      
      real*8 fmape(my),ydumb(my)
      real*4 timedumb
      real*4 Ree,alpe,bete,a0e

      integer i,j,mxe,mye,mze,ntotr,iproc,mym,mmy2,iscal
      character*3 ext1
      character*84 fnamesca 

!   --------------------- Program -------------------------------------

      master = 0

!      ---------------  read from file and distribute

      if (myid.eq.master) then     
         !write iscal (001,002...) into ext1 variable
         write(ext1,'(i3.3)') iscal
         !create scalar input file name:ex:finpsca001.dat
         fnamesca=filscainp(1:index(filscainp,' ')-1)//ext1//'.dat'

         open(iinp,file=fnamesca,status='unknown',form='unformatted')
!      for scalars the data should match the one given by main files 
         read(iinp) timedumb,Ree,alpe,bete,a0e,mxe,mye,mze,
     .               (ydumb(j),fmape(j),j=1,my),(wk1(j),j=1,2*my)
         
         ntotr=2*mye*mze

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
            call assignscal(wk1,scal(1,1,i),mye,mze)
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
           write(*,*) 'master reads proc n',iproc,' scalar sending it'
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


          ntotr=2*mye*mze
c         ------------  receive other modes ----------
       
         mmy2 = min(pe,mxe/2)-pb+1
         if (mmy2.gt.0) then
            do i=1,mmy2
               call MPI_RECV(wk1,ntotr,MPI_REAL,master,
     &                       MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
               call assignscal(wk1,scal(1,1,i),mye,mze)
            enddo
            write(*,*) 'proc no.',myid,'receives scalar from master'
         endif
      endif
      
      

      endsubroutine
!

!***********************************************************************!
!
!     Routine assignscal                                                ! 
!     save data from wk1(1:ntotr) to matrix form                        !
!***********************************************************************!
      subroutine assignscal(work,scal,my,mz)

      implicit none

      integer   j,k,my,mz
      real*4    scal(2*my,mz)
      real*4    work(2*my,mz)
      
      
      do k=1,mz
         do j=1,2*my
            scal(j,k)=work(j,k)
         enddo
      enddo


      endsubroutine



