!***********************************************************************!
!                                                                       !
!               writes a intermediate solution                          !
!    single jjs  4/01/01                                                !
!    plane  lines shc 01/03/05                                          !
!    Compressible Flow AAF 10/07/14                                     !
!                                                                       !
!                                                                       !
!***********************************************************************!

      subroutine escrphys(u,v,w,dudy,dwdy,scal,mfz,dscal,dmfz,
     %                    wk,wk1,wk1r,wk2,wk2r,wk3,wk3r,wk4,wk4r,
     %                    wk5,wk5r,wk6,wk6r,wk7,wk7r,wkT,wkTr,
     %                    ufou,ufour,vfou,vfour,wfou,wfour,
     %                    u00,v00,w00,du00,dw00,j,fj)
      ! Esta subrutina escribe en disco un plano yz, un yx y los
      ! marcados en ctes3d para j. 
      use spectra
      use wave
      use point
      use fis
      use combustion
      use ficheros
      use tem


      implicit none
      include 'ctes3D'
      
!     --------------------------------------------------------
      integer,intent(in)::j,fj
      integer plyx,plyz !planes for writing
            
      complex*8,dimension(0:mx1,0:mz1)::u,v,w,dudy,dwdy,scal,wk1,
     .                    ufou,vfou,wfou,mfz,dmfz,
     .                    dscal,wk2,wk3,wk4,wk5,wk6,wk7,wkT
      real*8    u00(my),v00(my),w00(my),du00(my),dw00(my)
      real*4 wk1r(mgalx+2 ,mgalz),wk2r(mgalx+2,mgalz),
     .       wk3r(mgalx+2 ,mgalz),
     .       wk4r(mgalx+2, mgalz),
     .       wk5r(mgalx+2 ,mgalz),
     .       wk6r(mgalx+2 ,mgalz),
     .       wk7r(mgalx+2 ,mgalz),
     .       wkTr(mgalx+2 ,mgalz),
     .       ufour(mgalx+2,mgalz),
     .       vfour(mgalx+2,mgalz),
     .       wfour(mgalx+2,mgalz),
     .       wk(0:mx1,0:mz1)
      integer i,k,kk,nacum,jj
      complex*8 xa,xb
      character*85 fname
      character*5  extp(27)
      character*5  exten
      ! ----------------------------------------------------

      !Compute jspecy
       do jj=1,nspec
          jspecy(jj) = ispec1 + jj -1
       enddo

      kk=0; nacum=1
      plyx=48  !plyx RP
      plyz=1  !plyz we anto to write
      extp=(/ '.Zfyz','.Zfyx','.Zfxz','.Hfyz','.Hfyx','.Hfxz',
     &       '.upyz','.upyx','.upxz','.vpyz','.vpyx','.vpxz',
     &       '.wpyz','.wpyx','.wpxz','.Tfyz','.Tfyx','.Tfxz',
     &       '.o1yz','.o1yx','.o1xz','.o2yz','.o2yx','.o2xz',
     &       '.o3yz','.o3yx','.o3xz'/)

      ! Velocidades       

      if(j.eq.1) then
         do k=1,27
            exten = extp(k)
            fname= filout(1:index(filout,' ')-1)//extp(k)
            if (exten(4:5).eq.'yz') then
                  nacum = plyz
                  !write(*,*) "plyz = ",nacum
            elseif (exten(4:5).eq.'yx') then
                  nacum = plyx
                  !write(*,*) "plyx = ",nacum
            endif
            open (1000+k,file=fname,status='unknown',form='unformatted')
            write(1000+k) time,Re,alp,bet,mgalx,my,mgalz,nspec,nacum 
           write(1000+k)(jspecy(jj),jj=1,nspec),(y(jj),fmap(jj),jj=1,my)
         enddo
      endif
!      write(*,*)'despues',j
 
      ! u

      do k=0,mz1
         do i=0,mx1
            !Save rho u into wk1
            wk1(i,k) = u(i,k) 
            !save scal in wk2
            wkT(i,k) = scal(i,k)
            wk7(i,k) = mfz(i,k)
         enddo
      enddo

      wk1r(1,1) = u00(j) !rhou00
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1) !rhou
      call fourxz(wkT,wkT,1,1) !scal (H)
      call fourxz(wk7,wk7,1,1) !scal (Z)

      !writing "Z"
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                i=1,mgalx),k=1,mgalz)
 
      !writing "H"
      kk=kk+1
      write(1000+kk) (wkTr(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wkTr(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wkTr(i,k),
     .                i=1,mgalx),k=1,mgalz)
 

      !Calculate T
      do k=1,mgalz
         do i=1,mgalx
             wkTr(i,k) = Temper(wkTr(i,k),wk7r(i,k))
         enddo
      enddo

      !Save u 
      do k=1,mgalz
         do i=1,mgalx
            wk1r(i,k) = wk1r(i,k)*wkTr(i,k) !now wk1r keeps U phys
            ufour(i,k) = wk1r(i,k) !Copy used to go back to FOU
         enddo
      enddo

      !writing "Z"
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                i=1,mgalx),k=1,mgalz)
! ... change it to phys & multiply by T in order take rho out
      !writing "U"
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),
     .                i=1,mgalx),k=1,mgalz)

      call fourxz(ufou,ufour,-1,1) !phys--> fou 
      !U is U FOU
      !wk1r is U PHYS

      ! v

      do k=0,mz1
         do i=0,mx1
            wk2(i,k) = v(i,k) !rhov already calculated
         enddo
      enddo

      wk2r(1,1) = v00(j) !rhov00
      wk2r(2,1) = 0d0
      call fourxz(wk2,wk2,1,1) !rhov

      !Save v
      do k=1,mgalz
         do i=1,mgalx
            wk2r(i,k) = wk2r(i,k)*wkTr(i,k) 
            vfour(i,k)= wk2r(i,k) 
         enddo
      enddo
      !Now write V
      kk=kk+1
      write(1000+kk) (wk2r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk2r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk2r(i,k),i=1,mgalx),
     .                                               k=1,mgalz)
      !Save v on FOURIER 
      call fourxz(vfou,vfour,-1,1) !phys--> fou

      ! w

      do k=0,mz1
         do i=0,mx1
            wk3(i,k) = w(i,k) !rhow
         enddo
      enddo

      wk3r(1,1) = w00(j)
      wk3r(2,1) = 0d0
      call fourxz(wk3,wk3,1,1)

      do k=1,mgalz
         do i=1,mgalx
            wk3r(i,k) = wk3r(i,k)*wkTr(i,k) !w 
            wfour(i,k) = wk3r(i,k) !w 
         enddo
      enddo
      !Write W
      kk=kk+1
      write(1000+kk) (wk3r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk3r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk3r(i,k),i=1,mgalx),
     .                                               k=1,mgalz)
      !Save w 
      call fourxz(wfou,wfour,-1,1) !phys--> fou

      ! T

      kk=kk+1
      write(1000+kk) (wkTr(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wkTr(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wkTr(i,k),i=1,mgalx),k=1,mgalz)      


!=======================================================!
! VORTICITY
!=======================================================!
!VARIABLES :
!----------
!wk1r,wk2r,wk3r    => PHYSICAL VELOCITIES (u,v,w)
!u,v,w             => FOURIER VELOCITIES (U,V,W)
!wkTr              => PHYSICAL TEMPERATURE (T)
!du/dy,dw/dy,dscal/dy => FOURIER DERIVATIVES
!=======================================================!

      ! VOR X

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            !Save in buffer all quantities we need to Transform to FIS
            wk4(i,k) = dscal(i,k)!dT/dy
            wk5(i,k) = dwdy(i,k) !drhow/dy
            wk6(i,k) = xb*vfou(i,k) !ikz*v 
         enddo
      enddo

      call fourxz(wk4,wk4,1,1) !dT/dy

      wk5r(1,1) = dw00(j) !drhow/dy(0,0)
      wk5r(2,1) = 0d0
      call fourxz(wk5,wk5,1,1) !drhow/dy

      call fourxz(wk6,wk6,1,1) !ikz*v
      do k=1,mgalz
         do i=1,mgalx
           !VOR X = T*(d(rhow)/dy-w*drho/dy)-dv/dz
           !VOR X = T*(d(rhow)/dy+w*(1/(T^2))*dT/dy-dv/dz
           wk7r(i,k) = wkTr(i,k)* wk5r(i,k)+wk3r(i,k)*wk4r(i,k)/
     .                 wkTr(i,k)  - wk6r(i,k) 
         enddo
      enddo

      !writing "vor X"
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)

      ! VOR Y

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            !VORY = du/dz-dw/dx
            wk7(i,k) = xb*ufou(i,k) - xa*wfou(i,k)   
         enddo
      enddo

      call fourxz(wk7,wk7,1,1) !Convert to PHYS 

      !writing "vor Y"
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)


      ! VOR Z

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            !Save in buffer all quantities we need to Transform to FIS
            !wk4(i,k) = dscal(i,k)!dT/dy
            wk5(i,k) = dudy(i,k) !drhou/dy
            wk6(i,k) = xa*vfou(i,k) !ikx*v 
! calculate it on FOU and...
         enddo
      enddo

      wk5r(1,1) = du00(j) !drhou/dy(0,0)
      wk5r(2,1) = 0d0
      call fourxz(wk5,wk5,1,1) !drhou/dy
      call fourxz(wk6,wk6,1,1) !ikx*v
      do k=1,mgalz
         do i=1,mgalx
           !VOR Z =- T*d(rhou)/dy-u*1/T*dT/dy+dv/dx
           wk7r(i,k) = -wkTr(i,k)* wk5r(i,k)-wk1r(i,k)*wk4r(i,k)/
     .                 wkTr(i,k)**2  + wk6r(i,k) 
         enddo
      enddo

      !writing "vor Z"
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)




      end  

