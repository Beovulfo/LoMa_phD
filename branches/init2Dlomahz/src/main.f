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
      use base

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
      complex(8)::  box1, box2
      integer k1 
      real*8  ywin(my)
      real(8) :: upmax,scaling
      !real*4  ywin(my)
      complex*8,allocatable :: 
     .                       ufou(:,:,:),vfou(:,:,:),wfou(:,:,:),
     .                       vor(:,:,:),phi(:,:,:), psi(:,:,:),
     .                       scal(:,:,:),mfz(:,:,:),rhofou(:,:,:),
     .                       aux1(:,:,:),aux2(:,:,:),
     .                       hfou(:,:,:),zfou(:,:,:),
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

      write(*,*) "gamma = ",gam, "Zs = ", Zs
      write(*,*) " LF = ",Lf,    "betha = ", betha

!     Allocate !
!     For FOU - PHYS
      allocate(ufou(my,mgalz,nplanes))
      allocate(vfou(my,mgalz,nplanes))
      allocate(wfou(my,mgalz,nplanes))
      !allocate(rhofou(my,mgalz,nplanes))
      allocate(hfou(my,mgalz,nplanes))
      allocate(zfou(my,mgalz,nplanes))
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

       
      do j=1,my
         Z00(j) = Z0(j)
         H00(j) = H0(j)
         !T00(j) = 1/rho0(j)
         !Z00(j) = 0.5d0*(1d0+tanh(-y(j)/2d0))
      enddo 
      
      !set all velocities to zero
      do i = 1, nplanes
         do k = 1, mgalz
            do j = 1,my
                ufou(j,k,i) = 0.0
                vfou(j,k,i) = 0.0
                wfou(j,k,i) = 0.0
                !rhofou(j,k,i) = 0.0
                hfou(j,k,i) = 0.0
                zfou(j,k,i) = 0.0
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
      !do j =1,my
      !    ywin(j) = exp(-(y(j)/swin)**2)
      !enddo


!----------------------------------------------!
! Input control parameters for spectrum input  !              
!  on hre.dat                                  !
!----------------------------------------------!

! initialization of energy boxes
      euu = 0d0
      evv = 0d0
      eww = 0d0

      !DEBUG
      !k0 = nhmax*alp
      !write(*,*) "k0 = ", k0
      !only MODE 1,0 perturbation if pflag = 1 
      pflag = 1;
      if (pflag == 1) then
          !Calculate upp
          !u=i/a0*dv/dy
          call deryr(vpp,dvdyr)  
          call deryr(vppi,dvdyi)  
          do j=1,my
            upC(j) =-1.0/a0*dvdyi(j) + ju*1.0/a0*dvdyr(j)
            write(*,*) "j=",j,"upC(j)=",upC(j)
          enddo
          !Find the maximum of upC
          upmax=0.0
          do j=1,my
             upmax=max(upmax,abs(upC(j)))
          enddo
          scaling = 1.0/upmax
          write(*,*) "upmax before normalizing was:",upmax

          do j=1,my
             upC(j)=upC(j)*scaling
          enddo
          
          do j=1,my
            vpC(j) =scaling*(vpp(j)+ju*vppi(j))
          enddo
          upmax=0.0
          do j=1,my
             upmax=max(upmax,abs(vpC(j)))
            !write(*,*) "j=",j,"vpC(j)=",vpC(j)
          enddo
          write(*,*) "vpmax after normalizing was:",upmax
          write(*,*) "scaling:",scaling
          !Sanity check: continuity equation
          box1=0.0
          do j=1,my
            !box1=box1+real(upC(j))
            box1=box1+ju*a0*upC(j)+scaling*(dvdyr(j)+ju*dvdyi(j))
          enddo
          write(*,*) "Sum of Continuity eq.=",box1

          !Need derivative of base flow
          call deryr(Z0,dZ0dy)
          call deryr(H0,dH0dy)
          !compute zp and hp
          do j=1,my
            zpC(j) = -1.0/(ju*a0*(Um0(j)-ccomplex))*vpC(j)*dZ0dy(j)
            hpC(j) = -1.0/(ju*a0*(Um0(j)-ccomplex))*vpC(j)*dH0dy(j)
          enddo
          
     
        
         !Mode 1,0 is the subharmonic of max growth
         !Lx should be such that alp = alp_max/2.0
         !Adding perturbation of eigenfunction 
         do j=1,my
            ufou(j,1,2) = upC(j)*A1
!!     .                 A1*exp(ju*2.0*pi*phi1)

            vfou(j,1,2) = vpC(j)*A1
!!!     .                 A1*exp(cmplx(0.0,1.0)*2.0*pi*phi1)
            scal(j,1,2) = hpC(j)*A1
!!!     .                 A1*exp(cmplx(0.0,1.0)*2.0*pi*phi1)
            mfz(j,1,2) = zpC(j)*A1
!!     .                 A1*exp(cmplx(0.0,1.0)*2.0*pi*phi1)
         enddo
         !Mode 2,0 is the max growth alpha
         do j=1,my
            ufou(j,1,3) = A2*upC(j)
            vfou(j,1,3) = A2*vpC(j)
            scal(j,1,3) = A2*hpC(j)
             mfz(j,1,3) = A2*zpC(j)
         enddo
      endif
!----_VIP

      ! Adding 00 modes
      do j=1,my
          scal(j,1,1) = H0(j)
          mfz(j,1,1) = Z0(j)
          ufou(j,1,1) = Um0(j)!u00
          !vfou(j,1,1) = 0.0 !The stability problem was
                            !solved for planar perturbations
          vfou(j,1,1) = V0(j)
          wfou(j,1,1) = 0.0
      enddo


      ! Save H and Z before going to phys
      do i = 1, nplanes
         do k = 1, mz
            do j = 1, my
                hfou(j,k,i) = scal(j,k,i)
                zfou(j,k,i) = mfz(j,k,i)
            enddo
         enddo
      enddo
!!================================================================!
!!going to phys 
!
      write(*,*) "starting local tranposing...."
      call localyz2zy(ufou,ufou,aux1)
      call localyz2zy(vfou,vfou,aux1)
      call localyz2zy(wfou,wfou,aux1)
      call localyz2zy(hfou,hfou,aux1)
      call localyz2zy(zfou,zfou,aux1)
!      
      write(*,*) "Transposed yz2zy variables."
!
      do i = 1, nplanes
         call fourz(ufou(1,1,i),1) !U
         call fourz(vfou(1,1,i),1) !V
         call fourz(wfou(1,1,i),1) !W
         call fourz(hfou(1,1,i),1) !rho
         call fourz(zfou(1,1,i),1) !rho
      enddo
!
      write(*,*) "Fourz transformation made."
!
      call chpl2ln(ufou,ufou,aux1,myid)
      call chpl2ln(vfou,vfou,aux1,myid)
      call chpl2ln(wfou,wfou,aux1,myid)
      call chpl2ln(hfou,hfou,aux1,myid)
      call chpl2ln(zfou,zfou,aux1,myid)
!
      write(*,*) "Changed pl2ln."
!
!
!!==============================================!
      call phys(ufou,vfou,wfou,hfou,zfou)
!!==============================================!
!
      write(*,*) "rhoui calculated..."
!
      call chln2pl(ufou,ufou,aux1,myid)
      call chln2pl(vfou,vfou,aux1,myid)
      call chln2pl(wfou,wfou,aux1,myid)
!
      write(*,*) "Changed ln2pl."
!
      do i=1,nplanes
         call fourz(ufou(1,1,i),-1) !U
         call fourz(vfou(1,1,i),-1) !V
         call fourz(wfou(1,1,i),-1) !W
      enddo 
!
      write(*,*) "Fourz inverse transformation done."
!
      call localzy2yz(ufou,ufou,aux1)
      call localzy2yz(vfou,vfou,aux1)
      call localzy2yz(wfou,wfou,aux1)
!
      write(*,*) "Local zy2yz transposing made."
 



!================================================================!

!---------------------------------------------------------------!
! We have now:
!       - ufou = RHOU
!       - vfou = RHOV
!       - wfou = RHOW
!---------------------------------------------------------------!

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

      ! Save mode 00
      do j = 1, my
           rhou00(j) = ufou(j,1,1)
           rhov00(j) = vfou(j,1,1)
           rhow00(j) = wfou(j,1,1)
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
      use base
      use combustion

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
         read(text,*) (dat(j),j=1,6)
!         dat=[Re,alp,bet,alp0,c]
! c is the eigenvalue 
!966      read(19,'(a)') text
!         if(text(1:2).eq.'CC') goto 966
!         read(text,*) nymax, nhmin,nhmax, swin
!966      read(19,'(a)') text
!         if(text(1:2).eq.'CC') goto 966
!         read(text,*) dm0, u0, sden

967      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 967
!         read(text,*) euutarget, evvtarget, ewwtarget
         read(text,*) A1, phi1, A2

968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) gam, S, Lf, betha
!        idat(5)= nymax (number of y modes)
!        idat(6)= nhmin (minimum horizontal wave number)
!        idat(7)= nhmax (maximum horizontal wave number)

65       read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 65
         read(text,'(a70)') filout

166      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 166
         read(text,'(a70)') filinp

         close(19)


      Re   =  dat(1)
      alp  = dat(2)
      bet  = dat(3)
      a0   = dat(4) 
      creal    = dat(5) 
      cimag    = dat(6) 
      ccomplex = creal+ju*cimag
      Deltat = 1e-6     ! this data is not used 
!     Initialization for combustion
      Zs = 1/(1+S)
      S = 1d0/Zs-1d0
      Hf = -gam*(S+1d0)/S
      
      iinp=2
      ista = 4
      ispf=15 !aaf ??
      istati=0 !stats accumulation FLAG
      isn=3
      iout=1
 
!        READ MESH

         write(*,*) "Reading mesh and mean profiles.."
         open(iinp,file=filinp,status='old')
         read(iinp,*) (y(j),Um0(j),V0(j),H0(j),Z0(j),rho0(j),
     .       vpp(j),vppi(j),j=1,my)
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
     .     'reading from:  ',filinp, 'and', filstt
         write(*,*)
         write(*,'(a,a)')
     .     '  write in :  ',filout

      end
       


! ===============================================
! Last step for physical change
!  PHYS - GOEST TO PHYS and computes RHOUi
! ===============================================
      subroutine phys(u,v,w,hff,zf)
      use point
      use fis
      use wkhvect
      use combustion

      implicit none
      include "ctes3D"

!   -----------------Variables--------------------

      integer i,j, jj
      real(4) :: u(0:mx-1,lb:le),v(0:mx-1,lb:le),
     .           w(0:mx-1,lb:le),hff(0:mx-1,lb:le),
     .           zf(0:mx-1,lb:le) 
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
            rhstwk(i) = hff(i,j)
         enddo

         do i = 0,mx-1
            ten11wk(i) = zf(i,j)
         enddo
         
         call fourx(up1wk,up1wk8,1) 
         call fourx(up2wk,up2wk8,1) 
         call fourx(up3wk,up3wk8,1) 
         call fourx(rhstwk,rhstwk8,1) 
         call fourx(ten11wk,ten11wk8,1) 

!computing in physcal domain.....
         do i = 0, mgalx-1
              maxup = max(maxup,abs(up1wk8(i)))
              maxvp = max(maxvp,abs(up2wk8(i)))
              maxwp = max(maxwp,abs(up3wk8(i)))
              !Calculate density from h z
              rhstwk8(i) =1d0/Temper(rhstwk8(i),ten11wk8(i))
              up1wk(i) = up1wk8(i)*rhstwk8(i)
              up2wk(i) = up2wk8(i)*rhstwk8(i)
              up3wk(i) = up3wk8(i)*rhstwk8(i)
         enddo

         call fourx(up1wk,up1wk8,-1) 
         call fourx(up2wk,up1wk8,-1) 
         call fourx(up3wk,up1wk8,-1) 

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


