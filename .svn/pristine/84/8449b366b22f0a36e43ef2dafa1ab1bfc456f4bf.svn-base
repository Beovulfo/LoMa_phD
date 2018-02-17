***********************************************************************!
!                      L o M A    D N S                                 !
!-----------------------------------------------------------------------!
!  Solves the Low Mach N-S equations for a 3D box domain with           !
!  periodic BC @ "x,z" and general BC's for "y".                        !
!  Spectral method for x,z Compact Finite Differences for y.            ! 
!  Formulation based on y-Vorticity and Laplacian of v, Kim,Moin and    !
!  Moser (1987). Based on LISO code.
!                                                                       !
!  Solves Turbulent Mixing layer with variable density, using
!  explicit time-step scheme.                                           !
!                                                                       ! 
!       LISO:                                                           !
!       --------------------------------------------                    !
!       tocado j.jimenez (1/92)                                         !
!       hundido j.c.a. (1/01)                                           !
!       reflotado, version plano-recta shc (8/04)                       !
!                                                                       !
!       LOMA:                                                           !
!       ----------------------------------------                        !
!       a. almagro[aaf] (2012-)
!       Explicit                                                        !
!       changed to work with Robin BC ofa & aaf (06/2013)               !
!       eliminated commons by modules - aaf(06/2014)                    !
!
!!       LOMAHZ:                                                           !
!       ----------------------------------------                        !
!       a. almagro[aaf] (2015-)
!       Reactive ML
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
      use spectra
      use wave
      use point
      use fis
      use MPI_GROUPS
      use combustion
      use prec
 
      implicit none 
 
      include 'mpif.h'
      include 'ctes3D'
      
!      ----------------------- Variables -------------------------------

      integer master,myid,iproc,numprocs,istat(MPI_STATUS_SIZE),ierr  
        
      integer nbuffsize,nwkasize,nbuffsizephys
      integer inxwk,inywk,inzwk,
     .        imxwk,imywk,imzwk,
     .        iten12,iten13,iten23,
     .        irf0u,irf0w,iu00wk,iw00wk,
     .        irf0v,iv00wk,
     .        ichwk,iTphys,irhst,idrho,nstr,nspsize,iflag,
     .        irhsz,ilapz,imfz,idhache,idzeta
      integer ivor,iphi,ipsi,iscal
     
      integer itags,newtag,imess,i,comproc(0:numerop-1),
     .                      savproc(numerop:numtot-1)      
      
      real(realprec) :: u00(my),w00(my)
      real(4) :: val 
      real(realprec) :: val8 
      real(realprec) ::  v00(my) 
      real(realprec), allocatable:: vor(:),phi(:),wk(:),
     .                      psi(:),scal(:), mfz(:)
      real(4), allocatable :: sp(:)

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
      nspsize   = (nspec+1)*(nz1+1)*12*(pe-pb+1)
      nbuffsize =max((pe-pb+1)*2*my*mgalz,(le-lb+1)*mx)
       !complex double prec
      !nbuffsize = max((pe-pb+1)*2*my*mgalz,(le-lb+1)*mx)
      !nbuffsizephys will keep things in physical space and Real(8) (x2)
      nbuffsizephys =(le-lb+1)*mgalx
      nwkasize  = 14*nbuffsize + 6*my + 3*nbuffsizephys
      !nwkasize  = 7*nbuffsize + 12*my + 8*nbuffsize
      !nwkasize  = 7*nbuffsize + 12*my + 6*nbuffsize
!added v00 and v00wk,rf0v
!      nwkasize  = 7*nbuffsize + 8*my + 6*nbuffsize
      allocate(sp(nspsize))
             
      if (myid.lt.numerop) then     ! Comp. Proc.
      
         allocate(vor(2*my*mz*(pe-pb+1)))
         allocate(phi(2*my*mz*(pe-pb+1)))
         allocate(psi(2*my*mz*(pe-pb+1)))
         allocate(scal(2*my*mz*(pe-pb+1)))
         allocate(mfz(2*my*mz*(pe-pb+1)))
         allocate(wk(nwkasize))

      else                          ! Save proc
          nwkasize = 5*2*mmp*my*mz+1000
          !nwkasize = 8*mmp*my*mz 
          allocate(wk(nwkasize))
          
      endif  
      !index for SAVE PROCS
       ivor = 1
       iphi = ivor + 2*mmp*my*mz
       ipsi = iphi + 2*mmp*my*mz
       iscal = ipsi + 2*mmp*my*mz
       imfz  = iscal + 2*mmp*my*mz

       
             ! make zero everything  
     
      nacumsp = 0
      do i=1,nspsize
         sp(i)=0.
      enddo


      if (myid.lt.numerop) then
      
         do i=1,2*my*mz*mmp 
            vor(i) = 0. 
            phi(i) = 0. 
            psi(i) = 0. 
            scal(i) = 0. 
            mfz(i) = 0. 
         enddo              
         do i=1,my          
            u00(i) = 0d0     
            w00(i) = 0d0  
            v00(i) = 0d0  
         enddo              
         !clean wk
         do i = 1,nwkasize
              wk(i) = 0.0
         enddo

          call getfil(vor,phi,psi,scal,mfz,u00,v00,w00,wk,myid)!pending
          !call getfil(vor,phi,psi,scal,u00,v00,w00,wk,myid)
      endif
      
          ! start time advancement  !
      irf0u   = 1
      irf0w   = irf0u   +  my  ! Real*8
      irf0v   = irf0w   +  my  ! Real*8
      iu00wk  = irf0v   +  my
      iw00wk  = iu00wk  +  my  
      iv00wk  = iw00wk  +  my  
      imxwk   = iv00wk  +  my
      imywk   = imxwk   + nbuffsize
      imzwk   = imywk   + nbuffsize
      inxwk   = imzwk   + nbuffsize
      inywk   = inxwk   + nbuffsize
      inzwk   = inywk   + nbuffsize
      irhst   = inzwk   + nbuffsize
      irhsz   = irhst   + nbuffsize
      ilapz   = irhsz   + nbuffsize
      idrho   = ilapz   + nbuffsize
      iten12  = idrho   + nbuffsize
      iten13  = iten12  + nbuffsize
      iten23  = iten13  + nbuffsize
      iTphys   = iten23  + nbuffsize
      idhache = iTphys  + nbuffsizephys
      idzeta  = idhache + nbuffsizephys
      ichwk   = idzeta  + nbuffsizephys
      
      
      if (myid.lt.numerop) then

      call cross1(vor,
     .            phi,
     .            psi,
     .            scal,
     .            mfz,
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
     .            wk(irhsz),
     .            wk(ilapz),
     .            wk(idrho),
     .            wk(iTphys),
     .            wk(idhache),
     .            wk(idzeta),
     .            wk(ichwk), !same position chwk and spwk
     .            wk(ichwk),
     .            sp,myid)
      else
         call cross1(wk(ivor),
     .               wk(iphi),
     .               wk(ipsi),
     .               wk(iscal),
     .               wk(imfz),
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

!!-------------------------------------------------!
!deallocate
      if (myid.lt.numerop) then     ! Comp. Proc.
         deallocate(vor)
         deallocate(phi)
         deallocate(psi)
         deallocate(scal)
         deallocate(mfz)
         deallocate(wk)
      else                          ! Save proc
         deallocate(wk)
      endif  
!!........................................! 
200   call MPI_FINALIZE(ierr)


      endprogram

!***********************************************************************!
!                   get data field from file                            !
!                                             jjs  27/1/01              !
!                                             shc  05/1/05              !
!***********************************************************************!
      subroutine getfil(vor,phi,psi,scal,mfz,u00,v00,w00,wk1,myid)
      use tem
      use point
      use ficheros
      use fis
      use MPI_GROUPS
      use combustion
      use prec
   

      implicit none

      include "mpif.h"
      include "ctes3D"

!   -------------------- variables -------------------------------------
     
      integer myid,master,istat(MPI_STATUS_SIZE),ierr       
      real(realprec) ::vor(2*my,mz,*),phi(2*my,mz,*),
     .          psi(2*my,mz,*),scal(2*my,mz,*),mfz(2*my,mz,*)
      real(4) :: wk1(*)
      real(realprec):: u00(*),v00(*),w00(*)
      !real(4) :: wk00(3*my)
      
      real(realprec) :: fmape(my) 
      real(4) ::  Ree,alpe,bete,a0e

      integer i,j,mxe,mye,mze,ntotr,iproc,mym,mmy2
      

!   --------------------- Program -------------------------------------

      master = 0

!      ---------------  read from file and distribute

      if (myid.eq.master) then     

         open(iinp,file=filinp,status='unknown',form='unformatted')
         
         read(iinp) time,Ree,alpe,bete,a0e,mxe,mye,mze,
     .               (y(j),fmape(j),j=1,my),(wk1(j),j=1,3*my)

         write(*,*) 
         write(*,*) 'reading input file ...' 
         write(*,*) 'time=',time,'Ree=',Ree,'alpe=',alpe,
     .   'bete=',bete,'mxe=',mxe,'mye=',mye,'mze =', mze
     .   ,"my=",my
         write(*,*)

         write(*,*) 'valores para vhange',my*mgalz*mmp


         !ntotr=4*mye*mze
         ntotr=10*mye*mze !vor phi psi scal mfz x2
!     ------------- 00 modes ------------------

         mym = min(my,mye)
         do j=1,mym
            u00(j) = wk1(3*j-2)
            v00(j) = wk1(3*j-1)
            w00(j) = wk1(3*j)
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
              
         !write(*,*) "pb=",pb,"pe=",pe
         do i=1,(pe-pb+1)
            read(iinp) (wk1(j),j=1,ntotr)
            !write(*,*) "iplane=",i
            call assign(wk1,vor(1,1,i),phi(1,1,i),psi(1,1,i),
     .           scal(1,1,i),mfz(1,1,i),my,mz,mye,mze,myid,i)
         enddo
         write(*,*) 'data read by MASTER'


!         ------------  other modes,  distribute to slaves --
 
         do iproc=1,numerop-1
           mmy2 = min(pend(iproc),mxe/2)-pbeg(iproc)+1
           write (*,*) mmy2,pend(iproc),mxe/2
           if (mmy2.gt.0) then
              do i=1,mmy2
                 read(iinp) (wk1(j),j=1,ntotr)
                 call MPI_SEND(wk1,ntotr,MPI_REAL,
     &                       iproc,iproc,MPI_COMM_WORLD,ierr)
              enddo
           write(*,*) 'master reads proc no',iproc,' data and send them'
           endif
           write(*,*) 'mmy2=',mmy2
         enddo

         close(iinp)

      else          ! -------- this is done by all slaves

         call MPI_RECV(mxe,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mye,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mze,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

          ntotr=10*mye*mze
c         ------------  receive other modes ----------
       
         mmy2 = min(pe,mxe/2)-pb+1
         if (mmy2.gt.0) then
            do i=1,mmy2
               call MPI_RECV(wk1,ntotr,MPI_REAL,master,
     &                       MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
               call assign(wk1,vor(1,1,i),phi(1,1,i),psi(1,1,i),
     .                      scal(1,1,i),mfz(1,1,i),my,mz,mye,mze,myid,i)
            enddo
            write(*,*) 'proc no.',myid,'receives data from master'
         endif
      endif
      
      

      end 
!***********************************************************************!
!
!                   Initializes everything
!                                           
!***********************************************************************!
      subroutine initcr(myid)
      use timacc
      use fis
      use tem
      use ficheros
      use point
      use MPI_datatype
      use matrices
      use diag
      use wave
      use wkhvect
      use statis
      use spectra
      use combustion

      implicit none 
      include 'mpif.h'
      include 'ctes3D'
      
      
!   !--------------------- variables -------------------------------------
      integer ndat, nidat
      parameter(ndat=12,nidat=7)

      integer istat(MPI_STATUS_SIZE),ierr,myid,idat(nidat)
      integer i,j,iproc,k,ifor,blocks,elem,stride,mxe,mye,mze

      real(8),allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      integer,allocatable::aux1(:),aux2(:)

      real(4) dat(ndat),Delt,zero,Ree,alpe,bete,a0e
      real(8) :: Zf0
      real(8) :: pi

      character*100 text

c ---------------------- Programa --------------------------------------

      pi=4.*atan(1.)
      zero = 0.0
c                ! reads input data from hre.dat !

      if(myid.eq.0) then
         open(19,file='hre.dat',status='old')

964      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 964
         read(text,*) (dat(j),j=1,4)
!         dat=[Re,alp,bet,-u]
965      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 965
         read(text,*) (dat(j),j=6,12)


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
      Deltat = 1E-5 
      cfl  = dat(5)
      sigma = dat(6)
      Pr  = dat(7)
      peclet = Pr*Re
      !combustion part
      S   = dat(8)
      Lf  = dat(9)
      gam = dat(10)
      betha=dat(11)
      T0f=dat(12)
      
      !additional calculations
      Zs = 1.0/(1.0+S)
      Zsbar = 1/(S/Lf+1)
      iLm = (S/Lf+1)/(S+1)
      Hf = T0f-1-gam*(S+1)/S
      !Diffusitivy time-step factor:
      ! the maximum diffusivity factor of the 3 equations implied
      ! Momentum:nu/nu_a=>T^(sigma+1)/Re
      ! H equation:nu/(nu_a*Pr) => T^(sigma+1)/Pe
      ! Z equation: nu/(nu_a*Pr)*(S/Lf+1)/(S+1) => T^(sigma+1)/Pe*iLm
      if (myid.eq.0) then
        diffcoef=0d0;
        diffcoef=max(max(1d0/Re,1d0/peclet),iLm/peclet)
        write(*,*) "ire=,",1/Re,"ipeclet=",1/peclet, "iLm=",iLm
        write(*,*) 'Diffusivity coefficient for time-step = ',diffcoef
      !!else
      !  !write(*,*) 'Diffusivity coefficient for time-step = ',diffcoef
      endif
      
      !Calculate integration constant for Tsmooth
      !Zf0 = 1d0 !value of Z on fuel stream
      !cteintT = gam*(S+1d0)/S-gam/((1d0-Zs)*Zs)*smootht(Zf0)
      
      nstep=idat(1)
      nimag=idat(2)
      nhist=idat(3)
      ncfl = idat(4)   
      id22=idat(5)
      ista=idat(5) !use initially same number for stats and imag
      nstart = idat(6)
      ntimes = idat(7)
      
      iinp=2
      ispf=15 !aaf ??
      istati=0 !stats accumulation FLAG
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
            call MPI_SEND(diffcoef,1,MPI_REAL8,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)

         enddo
         !print
         write(*,*) "first part of image read by master"
      else
         call MPI_RECV(y,my,MPI_REAL8,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(fmap,my,MPI_REAL8,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(diffcoef,1,MPI_REAL8,0,
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
         !debug
         write(*,*) "pointers_calc myid,pb,pe,mmp,mml:",
     .               myid, pb,pe,mmp,mml
      else
         call pointers_save(pbeg,pend)   
         pb=pbeg(myid-numerop)
         pe=pend(myid-numerop)         
         mmp = pe-pb+1
         !debug
         !write(*,*) "pointers_save info---"
         !write(*,*) myid, pb,pe,mmp
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
      
      call MPI_BCAST(filout,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
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

c ---------preparing jsptot for spectra comput....
      do j=1,nspec
         jsptot(j)           =  jspecy(j)
         !jsptot(2*nspec+2-j) = my-jspecy(j)+1
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


c ------------------- Computes wavenumbers ---------------


      do k=0,nz1
         xbet(k) = cmplx(zero,bet*k)
         !rbet(K) = bet*k
         icx(k) = k
      enddo

      do k=nz1+1,mz1
         xbet(k) = cmplx(zero ,-bet*(mz1+1-k))
         !rbet(k) = -bet*(mz1+1-k)
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
         !ralp(i) = alp*i
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
     .                    '  Re =',Re,'momentum'
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
         write(*,'(a10,e12.4,a8,f5.2)')
     .     'Delt =',Deltat,'  CFL =',CFL
         write(*,*)
         write(*,'(a,a)')
     .     'reading from:  ',filinp
         write(*,*)
         write(*,'(a,a)')
     .     '  write in :  ',filout
      endif

      end
       
!***********************************************************************!
!
!     Routine assign.
!     Distribute data from (mye,mze)-mesh to (my,mz)-mesh 
!
!          single    jjs 4/01/01, 
!          rewritten jjs 28/01/01
!          adapted   shc 5/09/04
!***********************************************************************!
      subroutine assign(work,vor,phi,psi,scal,mfz,
     .                  my,mz,mye,mze,myid,planex)
      use prec

      implicit none

      integer   mym,mzm,klen,kini1,kini2,j,k,k1,k2,my,mz,mye,mze
      real(realprec)    vor(2*my,mz),phi(2*my,mz)
      real(realprec)    psi(2*my,mz),scal(2*my,mz),mfz(2*my,mz)
      real(4)    work(5,2*mye,mze)
      real(realprec)    zero
      integer planex
      integer   myid
      integer flag
      
      zero = 0.0
      flag = 0
      
      mym = min(my,mye)
      mzm = min(mz,mze)

      klen = (mzm+1)/2
      kini1=mze - klen + 1
      kini2=mz - klen + 1
      !write(*,*) "assign for plane=",planex,"(",myid,")"
      
      do k=1,klen
         do j=1,2*my
            vor(j,k)=work(1,j,k)
            phi(j,k)=work(2,j,k)
            psi(j,k)=work(3,j,k)
            scal(j,k)=work(4,j,k)
            mfz(j,k)=work(5,j,k)
         enddo
      enddo

      do k=1,klen-1
         k1 = k + kini1
         k2 = k + kini2
         do j=1,2*my
            vor(j,k2) = work(1,j,k1)
            phi(j,k2) = work(2,j,k1)
            psi(j,k2) = work(3,j,k1)
            scal(j,k2)= work(4,j,k1)
            mfz(j,k2) = work(5,j,k1)
         enddo
      enddo

!      !SET IMAGINARY PART ON MODES (0,kz) TO ZERO
!      if (myid.eq.0) then
!         if (planex.eq.1) then
!            write(*,*) "Setting imaginary parts of modes (0 , kz)
!     .                   to zero."
!            do k=1,mz
!               do j=2,2*my,2 !Imaginary parts are even
!                       vor(j,k)  = 0.0
!                       phi(j,k)  = 0.0
!                       psi(j,k)  = 0.0
!                       scal(j,k) = 0.0
!               enddo
!            enddo
!         endif
!      endif



      end

!----------------------------------------------------------------------! 
      subroutine check_nan(u,n,flag)
      use prec
      implicit none
!      
      real(realprec) u(n)
      integer j,i,n,flag
!      
      do j=1,n
         if (u(j).ne.u(j)) flag = 1
      enddo
      end 


!!***********************************************************************!
!!                   get scalar field from file                          !
!!                   aaf 2014/01                                         !
!!***********************************************************************!
!      subroutine getscal(scal,wk1,iscal,myid)
!      use tem
!      use point
!      use ficheros
!      use fis
!      !use MPI_GROUPS
!
!      implicit none
!
!      include "mpif.h"
!      include "ctes3D"
!
!     
!!   -------------------- variables -------------------------------------
!     
!      integer myid,master,istat(MPI_STATUS_SIZE),ierr       
!      real(4) ::  scal(2*my,mz,*),wk1(*)
!      
!      real(8) :: fmape(my),ydumb(my)
!      real(4) :: timedumb
!      real(4) ::  Ree,alpe,bete,a0e
!
!      integer i,j,mxe,mye,mze,ntotr,iproc,mym,mmy2,iscal
!      character(3) ext1
!      character(100) fnamesca 
!
!!   --------------------- Program -------------------------------------
!
!      master = 0
!
!!      ---------------  read from file and distribute
!
!      if (myid.eq.master) then     
!         !write iscal (001,002...) into ext1 variable
!         write(ext1,'(i3.3)') iscal
!         !create scalar input file name:ex:finpsca001.dat
!         fnamesca=filscainp(1:index(filscainp,' ')-1)//ext1//'.dat'
!
!         open(iinp,file=fnamesca,status='unknown',form='unformatted')
!!      for scalars the data should match the one given by main files 
!         read(iinp) timedumb,Ree,alpe,bete,a0e,mxe,mye,mze,
!     .               (ydumb(j),fmape(j),j=1,my),(wk1(j),j=1,2*my)
!         
!         ntotr=2*mye*mze
!
!!     ------------ 00 modes only for master 
!
!         do iproc=1,numerop-1
!            call MPI_SEND(mxe,1,MPI_INTEGER,iproc,
!     &                iproc,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(mye,1,MPI_INTEGER,iproc,
!     &                iproc,MPI_COMM_WORLD,ierr)
!            call MPI_SEND(mze,1,MPI_INTEGER,iproc,
!     &                iproc,MPI_COMM_WORLD,ierr)
!         enddo
!
!
!c         ------------  other modes,  master node -----
!         
!         write(*,*) 'master reads its data',ntotr
!
!         do i=1,(pe-pb+1)
!            read(iinp) (wk1(j),j=1,ntotr)
!            call assignscal(wk1,scal(1,1,i),mye,mze)
!         enddo
!
!
!!         ------------  other modes,  distribute to slaves --
! 
!         do iproc=1,numerop-1
!           mmy2 = min(pend(iproc),mxe/2)-pbeg(iproc)+1
!           write (*,*) mmy2,pend(iproc),mxe/2
!           if (mmy2.gt.0) then
!              do i=1,mmy2
!                 read(iinp) (wk1(j),j=1,ntotr)
!                  call MPI_SEND(wk1,ntotr,MPI_REAL,
!     &                          iproc,iproc,MPI_COMM_WORLD,ierr)
!              enddo
!           write(*,*) 'master reads proc n',iproc,' scalar sending it'
!           endif
!         enddo
!
!         close(iinp)
!
!      else          ! -------- this is done by all slaves
!
!         call MPI_RECV(mxe,1,MPI_INTEGER,master,
!     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!         call MPI_RECV(mye,1,MPI_INTEGER,master,
!     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!         call MPI_RECV(mze,1,MPI_INTEGER,master,
!     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!
!
!          ntotr=2*mye*mze
!c         ------------  receive other modes ----------
!       
!         mmy2 = min(pe,mxe/2)-pb+1
!         if (mmy2.gt.0) then
!            do i=1,mmy2
!               call MPI_RECV(wk1,ntotr,MPI_REAL,master,
!     &                       MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
!               call assignscal(wk1,scal(1,1,i),mye,mze)
!            enddo
!            write(*,*) 'proc no.',myid,'receives scalar from master'
!         endif
!      endif
!      
!      
!
!      end
!!
!
!!***********************************************************************!
!!
!!     Routine assignscal                                                ! 
!!     save data from wk1(1:ntotr) to matrix form                        !
!!***********************************************************************!
!      subroutine assignscal(work,scal,my,mz)
!
!      implicit none
!
!      integer   j,k,my,mz
!      real*4    scal(2*my,mz)
!      real*4    work(2*my,mz)
!      
!      
!      do k=1,mz
!         do j=1,2*my
!            scal(j,k)=work(j,k)
!         enddo
!      enddo
!
!
!      end
!
!
