***********************************************************************!
!                      L o M A    D N S                                 !
!-----------------------------------------------------------------------!
! CREATES INITIAL CONDITIONS FOR LOMACTE
!***********************************************************************!
      Program main
      use inparams
      use wave
      use fis
      use tem
      use ficheros
!      use ifport 
 
      implicit none 
 
      include 'ctes3D'
      
!----------------------- Variables -------------------------------
      integer myid,i,j,k,jj
      real*4  kh,amp,phs1,phs2,ky
      real*8  nh,rK,k0
      real*4  urms,vrms,wrms
      complex*8  box1, box2
      integer k1 
      real*4  ywin(my)
      complex*8,allocatable :: 
     .                       ufou(:,:,:),vfou(:,:,:),wfou(:,:,:),
     .                       vor(:,:,:),phi(:,:,:), psi(:,:,:),
     .                       scal(:,:,:),rhofou(:,:,:),hzero(:,:,:),
     .                       aux1(:,:,:),aux2(:,:,:),
     .                       auxu(:,:,:),auxv(:,:,:),auxw(:,:,:)
      real*8 euu,evv,eww
      real*8,dimension(my) :: rhou00,rhov00,rhow00,rho00,T00
      integer :: iseed
      integer, dimension(:), allocatable :: aseed
      integer,dimension(8) :: dtseed


      !! ----- Set up random_number_number seed portably -----
      CALL RANDOM_SEED(size=iseed)
      ALLOCATE(aseed(1:iseed))
      CALL RANDOM_SEED(get=aseed)
      CALL DATE_AND_TIME(values=dtseed)
      aseed(iseed)=dtseed(8) 
      aseed(1)=dtseed(8)*dtseed(7)*dtseed(6)
      CALL RANDOM_SEED(put=aseed)
      DEALLOCATE(aseed)
      !! ----- Done setting up random_number_number seed -----



!--------------- initializes derivadas,fouriers,modes and read mesh
      myid = 0 
      call initcr(myid)

      write(*,*) "nhmin = ",nhmin, "swin = ", swin
!     Allocate !
!     For FOU - PHYS
      allocate(ufou(my,mgalz,nplanes))
      allocate(vfou(my,mgalz,nplanes))
      allocate(wfou(my,mgalz,nplanes))
      allocate(rhofou(my,mgalz,nplanes))
      allocate(aux1(my,mgalz,nplanes))
      allocate(aux2(my,mgalz,nplanes))
      allocate(auxu(my,mgalz,nplanes))
      allocate(auxv(my,mgalz,nplanes))
      allocate(auxw(my,mgalz,nplanes))
      !Fourier only --> writing
      allocate(vor(my,mz,nplanes))
      allocate(phi(my,mz,nplanes))
      allocate(psi(my,mz,nplanes))
      allocate(hzero(my,mz,nplanes))
      allocate(scal(my,mz,nplanes))

!     Initialization
      !set all velocities to zero
      do i = 1, nplanes
         do k = 1, mgalz
            do j = 1,my
                ufou(j,k,i) = 0.0
                vfou(j,k,i) = 0.0
                wfou(j,k,i) = 0.0
                rhofou(j,k,i) = 0.0
            enddo
         enddo
      enddo
      do i = 1, nplanes
         do k = 1, mz
            do j = 1,my
                vor(j,k,i) = 0.0
                phi(j,k,i) = 0.0
                psi(j,k,i) = 0.0
                hzero(j,k,i) = 0.0
                scal(j,k,i) = 0.0
            enddo
         enddo
      enddo



! Creating y window FILTER (gaussian type)
      !swin = 4.0*dm0 !filter width on y
      do j =1,my
          ywin(j) = exp(-0.5*(y(j)/swin)**2)
      enddo


!----------------------------------------------!
! Input control parameters for spectrum input  !              
!  on hre.dat                                  !
!----------------------------------------------!

! initialization of energy boxes
      euu = 0d0
      evv = 0d0
      eww = 0d0

      k0 = nhmax*alp
      write(*,*) "k0 = ", k0

!filling fourier modes of U
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5
             nh = (icx(k-1)**2 + i**2)**0.5 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    !ky   = 2*kh/nymax*jj
                    !ky   = jj/(2.0*pi)
                    call random_number(amp)
                    amp = 2.0*(0.5 - amp)
                    call random_number(phs1)
                    call random_number(phs2)
                    !y loop
                    do j=1,my
                        ufou(j,k,i) = ufou(j,k,i)+
     .                    amp*exp(2.0*pi*cmplx(0.0,1.0)*phs1)*
     .                    ywin(j)*
!     .                   sin(ky*y(j)+2.0*pi*phs2)*ywin(j)*
     .                   (kh/k0)**2*exp(-(kh/k0)**2)
                    enddo
                enddo
!----_VIP
                !Dont let imaginary part on modes 0,kz
                if (i.eq.1) then
                   do j=1,my
                        ufou(j,k,i) = real(ufou(j,k,i)) 
                   enddo
                endif

             !Compute integral of u^2
                !EXCLUDING MODE 00...
             endif
         enddo
      enddo


!filling fourier modes of V
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    !ky   = jj/(2.0*pi)
                    !ky   = kh/khref*jj
                    call random_number(amp)
                    amp = 2.0*(0.5 - amp)
                    call random_number(phs1)
                    call random_number(phs2)
                    !y loop
                    do j=1,my
                        vfou(j,k,i) = vfou(j,k,i)+
     .                    amp*exp(2.0d0*pi*cmplx(0.0,1.)*phs1)*
     .                    ywin(j)*
!     .                   sin(ky*y(j)+2.0*pi*phs2)*ywin(j)*
     .                   (kh/k0)**2*exp(-(kh/k0)**2)
                    enddo
                enddo
!----_VIP
                !Dont let imaginary part on modes 0,kz
                if (i.eq.1) then
                   do j=1,my
                        vfou(j,k,i) = real(vfou(j,k,i)) 
                   enddo
                endif

             !Compute integral of v^2
             endif
         enddo
      enddo

!filling fourier modes of W
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    !ky   = jj/(2.0*pi)
                    call random_number(amp)
                    amp = 2.0*(0.5 - amp)
                    call random_number(phs1)
                    call random_number(phs2)
                    !y loop
                    do j=1,my
                        wfou(j,k,i) = wfou(j,k,i)+
     .                    amp*exp(2.0d0*pi*cmplx(0.0,1.0)*phs1)*
     .                    ywin(j)*
!     .                   sin(ky*y(j)+2.0*pi*phs2)*ywin(j)*
     .                   (kh/k0)**2*exp(-(kh/k0)**2)
                    enddo
                enddo
!----_VIP
                !Dont let imaginary part on modes 0,kz
                if (i.eq.1) then
                   do j=1,my
                        wfou(j,k,i) = real(wfou(j,k,i)) 
                   enddo
                endif

             endif
         enddo
      enddo

!Eliminateing 2D modes
      do i = 1,nplanes
         do j=1,my
            ufou(j,1,i) = 0.0
            vfou(j,1,i) = 0.0
            wfou(j,1,i) = 0.0
         enddo
      enddo

      do k = 1,mz
         do j=1,my
            ufou(j,k,1) = 0.0
            vfou(j,k,1) = 0.0
            wfou(j,k,1) = 0.0
         enddo
      enddo


      !URMS at mid
      j=my/2+1
      write(*,*) "plane at mid = ",j 
      urms = 0.0
      vrms = 0.0
      wrms = 0.0
      do i=1,nplanes
         do k = 1,mz
           if (i.eq.1) then
               if (k.ne.1) then
                   urms = urms + abs(ufou(j,k,i)*conjg(ufou(j,k,i)))
               endif
           else
               urms = urms + 2.0*abs(ufou(j,k,i)*conjg(ufou(j,k,i)))
           endif
         enddo
      enddo

      do i=1,nplanes
         do k = 1,mz
           if (i.eq.1) then
               if (k.ne.1) then
                   vrms = vrms + abs(vfou(j,k,i)*conjg(vfou(j,k,i)))
               endif
           else
               vrms = vrms + 2.0*abs(vfou(j,k,i)*conjg(vfou(j,k,i)))
           endif
         enddo
      enddo

      do i=1,nplanes
         do k = 1,mz
           if (i.eq.1) then
               if (k.ne.1) then
                   wrms = wrms + abs(wfou(j,k,i)*conjg(wfou(j,k,i)))
               endif
           else
               wrms = wrms + 2.0*abs(wfou(j,k,i)*conjg(wfou(j,k,i)))
           endif
         enddo
      enddo

      urms = urms**0.5
      vrms = vrms**0.5
      wrms = wrms**0.5
 
 
      
      write(*,*) "URMS input before normalizing was: ", urms
      write(*,*) "VRMS input before normalizing was: ", vrms
      write(*,*) "WRMS input before normalizing was: ", wrms
      
      !normalizing u with energy Euu
      do i = 1,nplanes
         do k = 1,mz
            do j=1,my
                 ufou(j,k,i) = ufou(j,k,i)/urms * euutarget
                 vfou(j,k,i) = vfou(j,k,i)/vrms * evvtarget
                 wfou(j,k,i) = wfou(j,k,i)/wrms * ewwtarget
            enddo
         enddo
      enddo


      write(*,*) "URMS target: ", euutarget
      write(*,*) "VRMS target: ", evvtarget
      write(*,*) "WRMS target: ", ewwtarget


      ! Adding 00 modes
      !replicate same shape function as Temperature when s=2.0
      
      do j=1,my
          rho00(j) = (1d0 + (2.0-1d0)/(2.0+1d0)*tanh(-y(j)/(2d0*dm0)))
          T00(j)   = 1d0/rho00(j) 
      enddo
      !s can be used to multiply scalar max
      do j=1,my
          scal(j,1,1)=s*(T00(j)-T00(1))/(T00(my)-T00(1))
         !scal from 0 to 1 with tanh
          !scal(j,1,1)=(0.5d0 + 0.5*tanh(-y(j)/(2d0*dm0))) 
          !without perturbation on rho/T we only need to specify 00 mode
          ufou(j,1,1) = u0*tanh(-y(j)/(2d0*dm0))
          vfou(j,1,1) = 0.0
          wfou(j,1,1) = 0.0
      enddo

       !show scalar shape
       write(*,*) "Scal(min) = ",real(scal(1,1,1))
       write(*,*) "Scal(max) = ",real(scal(my,1,1))

!Save velocities for avoiding fourier back
      do i = 1, nplanes
         do k = 1, mz
            do j = 1,my
                auxu(j,k,i) = ufou(j,k,i)
                auxv(j,k,i) = vfou(j,k,i)
                auxw(j,k,i) = wfou(j,k,i)
            enddo
         enddo
      enddo

!================================================================!
!going to phys 

      write(*,*) "starting local tranposing...."
      call localyz2zy(ufou,ufou,aux1)
      call localyz2zy(vfou,vfou,aux1)
      call localyz2zy(wfou,wfou,aux1)
      
      write(*,*) "Transposed yz2zy variables."

      do i = 1, nplanes
         call fourz(ufou(1,1,i),1) !U
         call fourz(vfou(1,1,i),1) !V
         call fourz(wfou(1,1,i),1) !W
      enddo

      write(*,*) "Fourz transformation made."

      call chpl2ln(ufou,ufou,aux1,myid)
      call chpl2ln(vfou,vfou,aux1,myid)
      call chpl2ln(wfou,wfou,aux1,myid)

      write(*,*) "Changed pl2ln."


!==============================================!
      call physscal(ufou,vfou,wfou)
!==============================================!
!Save velocities for avoiding fourier back
      do i = 1, nplanes
         do k = 1, mz
            do j = 1,my
                ufou(j,k,i)= auxu(j,k,i)
                vfou(j,k,i)= auxv(j,k,i) 
                wfou(j,k,i)= auxw(j,k,i)
            enddo
         enddo
      enddo




!================================================================!

!---------------------------------------------------------------!
! We have now:
!       - ufou = RHOU
!       - vfou = RHOV
!       - wfou = RHOW
!       - rhofou = T (fou-phys-fou)
!       - scal = T
!---------------------------------------------------------------!
! Checking modes 00
      open(15,file = 'fields00.txt',status='unknown')
324   format(4(d22.14))
      write(15,324) (y(j),real(ufou(j,1,1)),real(vfou(j,1,1)),
     .              real(scal(j,1,1)),j=1,my)
      close(15)

      open(11,file = 'ywindow.txt',status='unknown')
326   format(2(d22.14))
      write(11,326) (y(j),ywin(j),j=1,my)
      close(11)


! VOR y = curl (rhoui) y
      do i = 1, nplanes
         do k = 1,mz
            do j = 1,my
               vor(j,k,i) = -xalp(i-1)*wfou(j,k,i)+xbet(k-1)
     .                        *ufou(j,k,i)
            enddo
         enddo
      enddo

! PHI = lap(my)
! my (aux1) = rhov - d(psi)/dy
      do i =1, nplanes
         do k = 1, mz
            do j = 1,my
               aux1(j,k,i) = vfou(j,k,i) 
            enddo
         enddo
      enddo
! also needed second derivative of my
      do i = 1, nplanes
          call deryyr2(aux1(1,1,i),aux2(1,1,i),mz)
      enddo

!now we can calculate laplacian of my 
      do i=1,nplanes
         do k=1,mz
            k1 = icx(k-1)
            rK = bet2(k1)+alp2(i-1)
            do j = 1, my
               phi(j,k,i) = -rK*aux1(j,k,i)+aux2(j,k,i)
            enddo
         enddo
      enddo

!=========PREPARING WRITING=================================!
!
! Save mode 00
      do j = 1, my
           rhou00(j) = real(ufou(j,1,1))
           !rhov00(j) = real(vfou(j,1,1))
           rhov00(j) = 0d0
           !rhow00(j) = real(wfou(j,1,1))
           rhow00(j) = 0d0
      enddo

      do i = 1, nplanes
         do k = 1, mz
            do j = 1,my
                psi(j,k,i) = 0.0
            enddo
         enddo
      enddo

! Now write the data
      call escru(vor,phi,psi,scal,rhou00,rhow00,rhov00,myid)
      call escru_reactml(vor,phi,psi,hzero,scal,
     .           rhou00,rhow00,rhov00,myid)

! INPUT FILE WRITTEN!
      write(*,*) "Input file written at: ", filout

! Checking modes 55
      open(17,file = 'fields55vorphi.txt',status='unknown')
      write(17,324) (y(j),real(vor(j,6,6)),real(phi(j,6,6)),
     .           real(psi(j,6,6)),real(scal(j,6,6)),j=1,my)
      close(17)



      deallocate(ufou,vfou,wfou,rhofou,vor,phi,psi,scal) 
      deallocate(aux1,aux2) 
      

      endprogram

     

!***********************************************************************!
!
!                   Initializes everything
!                                           
!***********************************************************************!
      subroutine initcr(myid)
      use inparams
      use timacc
      use fis
      use tem
      use ficheros
      use matrices
      use wave
      use wkhvect
      use diag
      use statis
      use point

      implicit none 
      include 'ctes3D'
      
      
!   !--------------------- variables -------------------------------------
      integer ndat0,ndat, nidat
      parameter(ndat0=3,ndat=6,nidat=7)

      integer myid,idat(nidat)
      integer i,j,k,mxe,mye,mze

      real(8),allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      integer,allocatable::aux1(:),aux2(:)

      real(4) dat(ndat),Delt,zero,Ree,alpe,bete,a0e
      real(8) dat0(ndat0)

      character*100 text

c ---------------------- Programa --------------------------------------

      zero = 0.0
c                ! reads input data from hre.dat !

         open(19,file='hre.dat',status='old')

965      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 965
         read(text,*) (dat(j),j=1,4)
!         dat=[Re,alp,bet,-u]

966      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 966
         read(text,*) dm0, u0, s


967      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 967
         read(text,*) euutarget, evvtarget, ewwtarget


968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) nymax, nhmin,nhmax, swin
!        idat(5)= nymax (number of y modes)
!        idat(6)= nhmin (minimum horizontal wave number)
!        idat(7)= nhmax (maximum horizontal wave number)

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


      Re   =  dat(1)
      alp  = dat(2)
      bet  = dat(3)
      a0   = dat(4) 
      Deltat = 1e-6     ! this data is not used 
      
      iinp=2
      ispf=15 !aaf ??
      istati=0 !stats accumulation FLAG
      isn=3
      iout=1
 
!        READ MESH

         write(*,*) "Reading mesh..."
         open(iinp,file=filinp,status='old')
         read(iinp,*) (y(j),j=1,my)
         rewind(iinp)
         close(iinp)
!need pointers for fourier

         call pointers_calc(pbeg,pend,lbeg,lend,procs,myid)
         pb=pbeg(myid)
         pe=pend(myid)
         lb=lbeg(myid)
         le=lend(myid)
         mmp = pe-pb+1
         mml = le-lb+1
         write(*,*) "pb = ", pb, " pe = ", pe, "lb = ", lb, "le = ", le
       


          
c    ------------  initializes fast fourier transforms and CFDiff ----
      call cfti(mgalz) !complex Fourier Transform init
c      call initfft()
      call rfti(mgalx) !real Fourier Transform init

      write (*,*) 'precalculo matrices ...'
      
      allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))
      
      call derivadas(d11,d12,d21,d22)      
            
      deallocate(d11,d12,d21,d22)
      
      write (*,*) 'FIN precalculo matrices ...'
      
!    --------------  coefficients for trapz averaging ------

      do j=2,my-1
         trp(j)= 0.25d0*(y(j+1)-y(j-1))
         hy(j) = (y(j+1)-y(j-1))/2d0
      enddo
      trp(1) = 0.25d0*(y(2) -y(1))
      trp(my)= 0.25d0*(y(my)-y(my-1))
      hy(1)  = (y(2)-y(1))/2d0
      hy(my) = (y(my)-y(my-1))/2d0


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


         write(*,'(a7,f8.2,a8,f6.3,a8,f6.3)')
     .                    '  Re = ',Re,'Momentum'
         write(*,'(a7,f8.3,a8,f6.3,a8,f6.3)')
     .                    'alp = ',alp,'  bet = ',bet
         write(*,*)

         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mgalx = ',mgalx,' mgalz = ',mgalz,' my = ',my
         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mx = ' ,mx,' mz = ',mz
         write(*,*)
         write(*,'(a8,f8.3,a8,f8.3,a8,f8.3)')
     .           'dm0 = ', dm0, ' umax = ', u0, ' s = ', s
         write(*,*)
         write(*,'(a,a)')
     .     'reading from:  ',filinp
         write(*,*)
         write(*,'(a,a)')
     .     '  write in :  ',filout

      end
       


! ===============================================
! Last step for physical change
!  PHYS - GOEST TO PHYS and computes RHOUi
! ===============================================
      subroutine physscal(u,v,w)
      use point
      use fis
      use wkhvect

      implicit none
      include "ctes3D"

!   -----------------Variables--------------------

      integer i,j, jj
      real(4) :: u(0:mx-1,lb:le),v(0:mx-1,lb:le),
     .           w(0:mx-1,lb:le) 
      real(4) maxup,maxvp,maxwp
     
      !check perturbations induced
      maxup = 0.0
      maxvp = 0.0
      maxwp = 0.0

      do j = lb, le
! copy lines
         do i = 0,mx-1
            up1wk(i) = u(i,j)
         enddo

         do i = 0,mx-1
            up2wk(i) = v(i,j)
         enddo

         do i = 0,mx-1
            up3wk(i) = w(i,j)
         enddo
         
         call fourx(up1wk,up1wk8,1) 
         call fourx(up2wk,up2wk8,1) 
         call fourx(up3wk,up3wk8,1) 

!computing in physcal domain.....
         do i = 0, mgalx-1
              maxup = max(maxup,abs(up1wk8(i)))
              maxvp = max(maxvp,abs(up2wk8(i)))
              maxwp = max(maxwp,abs(up3wk8(i)))
         enddo

      enddo
      
      !print max velocitites
      write(*,*) "MAXIMUM UP = ",maxup
      write(*,*) "MAXIMUM VP = ",maxvp
      write(*,*) "MAXIMUM WP = ",maxwp

      end subroutine


