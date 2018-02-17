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
      use combustion

      implicit none 
      include 'ctes3D'
!----------------------- Variables -------------------------------
      integer :: iseed
      integer, dimension(:), allocatable :: aseed
      integer,dimension(8) :: dtseed
      integer myid,i,j,k,jj
      real*4  kh,kt,amp,phs1,phs2,ky
      real*8  nh,rK,k0
      real*4  urms,vrms,wrms
      complex*8  box1, box2
      integer k1 
      real*8  ywin(my)
      !real*4  ywin(my)
      complex*8,allocatable :: 
     .                       ufou(:,:,:),vfou(:,:,:),wfou(:,:,:),
     .                       vor(:,:,:),phi(:,:,:), psi(:,:,:),
     .                       scal(:,:,:),mfz(:,:,:),rhofou(:,:,:),
     .                       aux1(:,:,:),aux2(:,:,:),
     .                       auxu(:,:,:),auxv(:,:,:),auxw(:,:,:)
      real*8 euu,evv,eww
      real*8,dimension(my) :: rhou00,rhov00,rhow00,rho00,T00,Z00,YO,H00
      real*8 :: yf,CYO
      !! ----- Set up random_number seed portably -----
      CALL RANDOM_SEED(size=iseed)
      ALLOCATE(aseed(1:iseed))
      CALL RANDOM_SEED(get=aseed)
      CALL DATE_AND_TIME(values=dtseed)
      aseed(iseed)=dtseed(8) 
      aseed(1)=dtseed(8)*dtseed(7)*dtseed(6)
      CALL RANDOM_SEED(put=aseed)
      DEALLOCATE(aseed)
      !! ----- Done setting up random_number seed -----


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
      allocate(scal(my,mz,nplanes))
      allocate(mfz(my,mz,nplanes))

       
!     Initialization
      Zs=0.2
      Lf=1.0d0
      gam=3.73
      betha=20
      S = 1d0/Zs-1d0
      Hf = -gam*(S+1d0)/S
      yf = 1.4 !depends on Zs
      cteintT = gam*(S+1d0)/S-gam/((1d0-Zs)*Zs)*smooth(1d0)
      write(*,*) "smooth1=",smooth(1d0)

      do j=1,my
        T00(j) = 1d0+gam*exp(-(y(j)-yf)**2*0.25)
        !write(*,*) "T00(j)=",T00(j)
      enddo
      !cteint for YO
      CYO = 1d0-(-1d0/Zs)*smooth(0d0)
      write(*,*) "CYO=,",CYO
      write(*,*) "cteintT=,",cteintT
      !calculate Z00
      do j=1,my
         Z00(j) = 0.5d0*(1d0+tanh(-y(j)/2d0))
      enddo 
      !calculate YO
      do j=1,my
         YO(j) = CYO+(-1d0/Zs)*smooth(Z00(j))
      enddo
      !calculate H
      do j=1,my
         H00(j) = T00(j)-1d0+(gam*(S+1d0)/S*(YO(j)-1d0))
         !write(*,*) H00(j)
      enddo
      
      
       

      
      !set all velocities to zero
      do i = 1, nplanes
         do k = 1, mgalz
            do j = 1,my
                ufou(j,k,i) = 0.0
                vfou(j,k,i) = 0.0
                wfou(j,k,i) = 0.0
                rhofou(j,k,i) = 0.0
                aux1(j,k,i) = 0.0
                aux2(j,k,i) = 0.0
            enddo
         enddo
      enddo

      do i = 1, nplanes
         do k = 1, mz
            do j = 1,my
                vor(j,k,i) = 0.0
                phi(j,k,i) = 0.0
                psi(j,k,i) = 0.0
                scal(j,k,i) = 0.0
                mfz(j,k,i) = 0.0
            enddo
         enddo
      enddo

! Creating y window FILTER (gaussian type)
      !swin = 4.0*dm0 !filter width on y
      do j =1,my
          ywin(j) = exp(-0.5d0*(y(j)/swin)**2)
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

!filling fourier modes of U'
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    ky   = 2*kh/nymax*jj
                    kt =(kh**2+ky**2)**0.5d0
                    !ky   = jj/(2.0*pi)
                    call random_number(amp)
                    amp = 2.0d0*(0.5d0 - amp)
                    call random_number(phs1)
                    !call random_number(phs2)
                    !y loop
                    do j=1,my
                        ufou(j,k,i) = ufou(j,k,i)+
     .                  amp*exp(cmplx(0.0,1.0)*(ky*y(j)+2.0*pi*phs1))*
     .                  ywin(j)*(kt/k0)**2*exp(-(kt/k0)**2)
!     .                  ywin(j)*(kh/k0)**2*exp(-(kh/k0)**2)
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


!filling fourier modes of V'
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    ky   = 2*kh/nymax*jj
                    kt =(kh**2+ky**2)**0.5d0
                    !ky   = jj/(2.0*pi)
                    call random_number(amp)
                    amp = 2.0d0*(0.5d0 - amp)
                    call random_number(phs1)
                    !call random_number(phs2)
                    !y loop
                    do j=1,my
                        vfou(j,k,i) = vfou(j,k,i)+
     .                  amp*exp(cmplx(0.0,1.0)*(ky*y(j)+2.0*pi*phs1))*
     .                  ywin(j)*(kt/k0)**2*exp(-(kt/k0)**2)
!     .                  ywin(j)*(kh/k0)**2*exp(-(kh/k0)**2)
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

!filling fourier modes of W'
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
             if (nh.ge.nhmin) then
                do jj = 1, nymax
                    ky   = 2*kh/nymax*jj
                    kt =(kh**2+ky**2)**0.5d0
                    !ky   = jj/(2.0*pi)
                    call random_number(amp)
                    amp = 2.0d0*(0.5d0 - amp)
                    call random_number(phs1)
                    !call random_number(phs2)
                    !y loop
                    do j=1,my
                        wfou(j,k,i) = wfou(j,k,i)+
     .                  amp*exp(cmplx(0.0,1.0)*(ky*y(j)+2.0*pi*phs1))*
     .                  ywin(j)*(kt/k0)**2*exp(-(kt/k0)**2)
!     .                  ywin(j)*(kh/k0)**2*exp(-(kh/k0)**2)
                    enddo
                enddo
!
                !Dont let imaginary part on modes 0,kz
                if (i.eq.1) then
                   do j=1,my
                        wfou(j,k,i) = real(wfou(j,k,i)) 
                   enddo
                endif

             endif
         enddo
      enddo



!-------------MAKE SOLENOIDAL: curl (vec(u))
!      auxu = dw/dy-dv/dz
!      auxv = du/dz-dw/dx
!      auxw = dv/dx-du/dy
!.................................
      do i = 1,nplanes
         call deryr2(wfou(1,1,i),aux1(1,1,i),mz)
      enddo 
      do i = 1,nplanes
         call deryr2(ufou(1,1,i),aux2(1,1,i),mz)
      enddo 
      do i= 1,nplanes
         do k=1,mz
            do j=1,my
               auxu(j,k,i)=aux1(j,k,i)-xbet(k-1)*vfou(j,k,i)          
               auxv(j,k,i)=xbet(k-1)*ufou(j,k,i)-xalp(i-1)*wfou(j,k,i)  
               auxw(j,k,i)=xalp(i-1)*vfou(j,k,i)-aux2(j,k,i)  
            enddo
         enddo
      enddo
!
!-------Ensuring that each mode is properly scaled for
!       spectrum shape
!
      do i=1,nplanes
         do k = 1,mz
             kh = (alp2(i-1)+bet2(k-1))**0.5d0
             nh = (icx(k-1)**2 + i**2)**0.5d0 
             !Rectangle Filter for too big and too small wavenumbers
              do j=1,my
                        ufou(j,k,i) = auxu(j,k,i)
!     .                   (kh/k0)**2*exp(-(kh/k0)**2)
!     .                    ywin(j)*

                        vfou(j,k,i) = auxv(j,k,i)
!     .                   (kh/k0)**2*exp(-(kh/k0)**2)

                        wfou(j,k,i) = auxw(j,k,i)
!     .                   (kh/k0)**2*exp(-(kh/k0)**2)
               enddo
          enddo
      enddo
!!!Eliminateing 2D modes
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



!------------------

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
               vrms = vrms + 2.0d0*abs(vfou(j,k,i)*conjg(vfou(j,k,i)))
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
               wrms = wrms + 2.0d0*abs(wfou(j,k,i)*conjg(wfou(j,k,i)))
           endif
         enddo
      enddo

      urms = urms**0.5d0
      vrms = vrms**0.5d0
      wrms = wrms**0.5d0
 
 
      
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
      do j=1,my
          rhofou(j,1,1) = rho00(j)
          !without perturbation on rho/T we only need to specify 00 mode
          ufou(j,1,1) = u0*tanh(-y(j)/(2d0*dm0)) !u00
          vfou(j,1,1) = 0.0
          wfou(j,1,1) = 0.0
      enddo

       !show Density ratio
       !write(*,*) "Density ratio s12 = ",real(rho00(1)/rho00(my))
       !write(*,*) "Density ratio s21 = ",real(rho00(my)/rho00(1))

!================================================================!

!---------------------------------------------------------------!
! We have now:
!       - ufou = RHOU
!       - vfou = RHOV
!       - wfou = RHOW
!       - rhofou = T (fou-phys-fou)
!       - scal = T
!---------------------------------------------------------------!
      open(11,file = 'ywindow.txt',status='unknown')
326   format(2(d22.14))
      write(11,326) (y(j),ywin(j),j=1,my)
      close(11)

!1! Checking modes 55
!1      open(13,file = 'fields55.txt',status='unknown')
!1328   format(5(d22.14))
!1      write(13,324) (y(j),real(ufou(j,6,6)),real(vfou(j,6,6)),
!1     .           real(wfou(j,6,6)),real(rho00(j)),j=1,my)
!1      close(13)

      !Computing VOR-PHI-PSI 
      ! Save drhov/dy on aux1
      do i = 1,nplanes
         call deryr2(vfou(1,1,i),aux1(1,1,i),mz)
      enddo 
     
      box1=0.0
      !Save div(rhoui) on aux2
      do i = 1,nplanes
         do k = 1,mz
            do j = 1,my
               aux2(j,k,i) =xalp(i-1)*ufou(j,k,i)+
     .          xbet(k-1)*wfou(j,k,i)+aux1(j,k,i)
               !check sum
               box1 = box1 + aux2(j,k,i)
            enddo
          enddo
       enddo
      
       !Check
       write(*,*) "Checking divergence matrix:"
       write(*,*) "sum(div(rhou_i)) = ", box1





!     Solve laplacian to find PSI lap(psi) = div (rho u) 
      do i=1,nplanes
         do k=1,mz
            k1 = icx(k-1)
            rK = bet2(k1)+alp2(i-1)
            !call Lapvdv(aux2(1,k,i),psi(1,k,i),aux1(1,k,i),rK)
            call Lapvdvhom(aux2(1,k,i),psi(1,k,i),aux1(1,k,i),rK)
         enddo
      enddo

      !Check dPSI/dy
      box2 = 0.0
      do i = 1,nplanes
         do k = 1,mz
            do j = 1,my
               box2 = box2 + psi(j,k,i)
            enddo
          enddo
       enddo

       !Check
       write(*,*) "Checking psi matrix:"
       write(*,*) "sum(psi(ij)) = ", box2



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
               aux1(j,k,i) = vfou(j,k,i) - aux1(j,k,i)
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

!If some perturbations are introduced we need to pass rhofou
! containing T to scal ready to write
! SCAL
      do i = 1, nplanes
         do k = 1, mz
            do j = 1, my
                scal(j,k,i) = 0.0 !put all zeros
                mfz(j,k,i) = 0.0 !put all zeros
            enddo
         enddo
      enddo
        do j=1,my
            !H and Z
            scal(j,1,1)=H00(j)
            mfz(j,1,1)=Z00(j)
         enddo
      ! Save mode 00
         do j = 1, my
           rho00(j) = 1d0/Temper(H00(j),Z00(j))
           rhou00(j) =rho00(j)*u0*tanh(-y(j)/(2d0*dm0)) 
           !rhov00(j) = real(vfou(j,1,1))
           rhov00(j) = 0d0
           !rhow00(j) = real(wfou(j,1,1))
           rhow00(j) = 0d0
         enddo
! Checking modes 00
      open(15,file = 'fields00.txt',status='unknown')
324   format(5(d22.14))
      write(15,324) (y(j),rhou00(j),real(scal(j,1,1)),
     .           real(mfz(j,1,1)),  
     .           rho00(j),j=1,my)
      close(15)



!=========PREPARING WRITING=================================!
!

! Now write the data
      call escru(vor,phi,psi,scal,mfz,rhou00,rhow00,rhov00,myid)

! INPUT FILE WRITTEN!
      write(*,*) "Input file written at: ", filout

! Checking modes 55
      open(17,file = 'fields55vorphi.txt',status='unknown')
      write(17,324) (y(j),real(vor(j,6,6)),real(phi(j,6,6)),
     .           real(psi(j,6,6)),real(scal(j,6,6)),j=1,my)
      close(17)



      deallocate(ufou,vfou,wfou,vor,phi,psi,scal) 
      deallocate(aux1,aux2,auxu,auxv,auxw) 
      

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
         read(text,*) dm0, u0, sden

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
     .           'dm0 = ', dm0, ' umax = ', u0, ' s = ', sden
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
      subroutine phys(u,v,w,rho)
      use point
      use fis
      use wkhvect

      implicit none
      include "ctes3D"

!   -----------------Variables--------------------

      integer i,j, jj
      real(4) :: u(0:mx-1,lb:le),v(0:mx-1,lb:le),
     .           w(0:mx-1,lb:le),rho(0:mx-1,lb:le) 
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

         do i = 0,mx-1
            rhstwk(i) = rho(i,j)
         enddo
         
         call fourx(up1wk,up1wk8,1) 
         call fourx(up2wk,up2wk8,1) 
         call fourx(up3wk,up3wk8,1) 
         call fourx(rhstwk,rhstwk8,1) 

!computing in physcal domain.....
         do i = 0, mgalx-1
              maxup = max(maxup,abs(up1wk8(i)))
              maxvp = max(maxvp,abs(up2wk8(i)))
              maxwp = max(maxwp,abs(up3wk8(i)))

              up1wk(i) = up1wk8(i)*rhstwk8(i)
              up2wk(i) = up2wk8(i)*rhstwk8(i)
              up3wk(i) = up3wk8(i)*rhstwk8(i)
              !rhstwk(i) = 1./rhstwk8(i) !get T instead of rho
         enddo

         call fourx(up1wk,up1wk8,-1) 
         call fourx(up2wk,up1wk8,-1) 
         call fourx(up3wk,up1wk8,-1) 
         !call fourx(rhstwk,up1wk8,-1) 

        ! do i=0,mx-1
        !    rho(i,j) = rhstwk(i)
        ! enddo

         do i=0,mx-1
            u(i,j) = up1wk(i)
         enddo

         do i=0,mx-1
            v(i,j) = up2wk(i)
         enddo

         do i=0,mx-1
            w(i,j) = up3wk(i)
         enddo

      enddo
      
      !print max velocitites
      write(*,*) "MAXIMUM UP = ",maxup
      write(*,*) "MAXIMUM VP = ",maxvp
      write(*,*) "MAXIMUM WP = ",maxwp

      end subroutine


