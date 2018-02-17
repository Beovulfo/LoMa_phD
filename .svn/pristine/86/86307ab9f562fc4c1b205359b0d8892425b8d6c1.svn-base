!***********************************************************************!
!                                                                       !
!               writes a intermediate solution                          !
!     AAF - PQR May 2017                                                                       !
!                                                                       !
!***********************************************************************!

      subroutine escrphys(u,v,w,dudy,dwdy,dvdy,mfz,dscal,scal,dmfz,
     %                    wk,wk1,wk1r,wk2,wk2r,wk3,wk3r,wk4,wk4r,
     %                    wk5,wk5r,wk6,wk6r,wk7,wk7r,wkT,wkTr,
     %                    wk8,wk8r,wk9,wk9r,
     %                    ufou,ufour,vfou,vfour,wfou,wfour,
     %                    u00,v00,w00,du00,dv00,dw00,j,fj)
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
     .                    ufou,vfou,wfou,mfz,dmfz,dvdy,
     .                   dscal,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wkT
      real*8    u00(my),v00(my),w00(my),du00(my),dv00(my),dw00(my)
      real*4 wk1r(mgalx+2 ,mgalz),wk2r(mgalx+2,mgalz),
     .       wk3r(mgalx+2 ,mgalz),
     .       wk4r(mgalx+2, mgalz),
     .       wk5r(mgalx+2 ,mgalz),
     .       wk6r(mgalx+2 ,mgalz),
     .       wk7r(mgalx+2 ,mgalz),
     .       wk8r(mgalx+2 ,mgalz),
     .       wk9r(mgalx+2 ,mgalz),
     .       wkTr(mgalx+2 ,mgalz),
     .       ufour(mgalx+2,mgalz),
     .       vfour(mgalx+2,mgalz),
     .       wfour(mgalx+2,mgalz),
     .       wk(0:mx1,0:mz1)
      integer i,k,kk,nacum,jj
      real*4 duma,dumb
      complex*8 xa,xb
      character*85 fname
      character*5  extp(18)
      character*5  exten
      ! ----------------------------------------------------

      !Compute jspecy
       do jj=1,nspec
          jspecy(jj) = ispec1 + jj -1
       enddo

      kk=0; nacum=1
      plyx=48  !plyx RP
      plyz=1  !plyz we anto to write

      extp=(/'.ppyz','.ppyx','.ppxz','.qqyz','.qqyx','.qqxz',
     &       '.rryz','.rryx','.rrxz','.qsyz','.qsyx','.qsxz',
     &       '.rsyz','.rsyx','.rsxz','.qwyz','.qwyx','.qwxz'/)


      if(j.eq.1) then
         do k=1,18
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
      write(*,*)'plane',j
 
      ! MAKE ALL TRANSFORMATIONS NEEDED TO HAVE JACOBIAN MATRIX
      ! COMPONENTS
      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            !Save  rhou rhov rhow
            wk1(i,k) = u(i,k)!ru 
            wk2(i,k) = v(i,k)!rv
            wk3(i,k) = w(i,k)!rw
            !Save drhou_i/dy
            wk4(i,k) = dudy(i,k)!drhoudy 
            wk5(i,k) = dvdy(i,k)!rhovdy
            wk6(i,k) = dwdy(i,k)!rdhowdy
            !Save H and Z
            wkT(i,k) = scal(i,k)!H
            wk7(i,k) = mfz(i,k)!Z
            !wk7(i,k) = dscal(i,k)!dscal
            !and their derivatives
            wk8(i,k) = dscal(i,k)!dH/dy
            wk9(i,k) = dmfz(i,k)!dZ/dy
         enddo
      enddo

      wk1r(1,1) = u00(j) !rhou00
      wk1r(2,1) = 0d0
      wk2r(1,1) = v00(j) !rhov00
      wk2r(2,1) = 0d0
      wk3r(1,1) = w00(j) !rhow00
      wk3r(2,1) = 0d0
      wk4r(1,1) = du00(j) !drhou00
      wk4r(2,1) = 0d0
      wk5r(1,1) = dv00(j) !drhov00
      wk5r(2,1) = 0d0
      wk6r(1,1) = dw00(j) !drhow00
      wk6r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1) !rhou
      call fourxz(wk2,wk2,1,1) !rhov
      call fourxz(wk3,wk3,1,1) !rhow
      call fourxz(wk4,wk4,1,1) !drhoudy
      call fourxz(wk5,wk5,1,1) !drhovdy
      call fourxz(wk6,wk6,1,1) !drhowdy
      call fourxz(wkT,wkT,1,1) !H
      call fourxz(wk7,wk7,1,1) !Z
      call fourxz(wk8,wk8,1,1) !dH/dy
      call fourxz(wk9,wk9,1,1) !dZ/dy
      !Calculate T
      do k=1,mgalz
         do i=1,mgalx
           duma=wkTr(i,k)!H
           dumb=wk7r(i,k)!Z 
           wkTr(i,k) = Temper(duma,dumb) !Temp
           !dTdy = dT/dH dH/dy + dZb/dy*dT/dZb
           wk7r(i,k) = wk8r(i,k) + wk9r(i,k)*dTdZ(dumb)!dT/dy
         enddo
      enddo

            !Calculate derivatives of y with chain rule
      do k=1,mgalz
         do i=1,mgalx
            wk4r(i,k) = wk4r(i,k)*wkTr(i,k)!Tdrhoudy
     .       + wk1r(i,k)*wk7r(i,k)!rhou*dTdy
            wk5r(i,k) = wk5r(i,k)*wkTr(i,k)!Tdrhovdy
     .       + wk2r(i,k)*wk7r(i,k)!rhov*dTdy
            wk6r(i,k) = wk6r(i,k)*wkTr(i,k)!Tdrhowdy
     .       + wk3r(i,k)*wk7r(i,k)!rhow*dTdy
            ufour(i,k) = wk1r(i,k)*wkTr(i,k)!u
            vfour(i,k) = wk2r(i,k)*wkTr(i,k)!v
            wfour(i,k) = wk3r(i,k)*wkTr(i,k)!w
         enddo
      enddo
      !NOW PHYS BUFFERS:
      !---------------------------
      ! ufour <- u
      ! vfour <- v
      ! wfour <- w
      ! wk4r  <- du/dy
      ! wk5r  <- dv/dy
      ! wk6r  <- dw/dy 
      !----------------------------------
      !With these buffers we can create any invariant
      !AVaILABLE BUFFERS:
      !wk1r,wk2r,wk3r,scal,wkTr,wk7
! ... change it to phys & multiply by T in order take rho out
      call fourxz(ufou,ufour,-1,1) !phys--> fou 
      call fourxz(vfou,vfour,-1,1) !phys--> fou 
      call fourxz(wfou,wfour,-1,1) !phys--> fou 
      !U is U FOU
      !wk1r is U PHYS
      !NOW
      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            !Compute dudx,dudz,dvdx...
            wk1(i,k) = xb*ufou(i,k)!dudz 
            wk2(i,k) = xb*vfou(i,k)!dvdz
            wk3(i,k) = xb*wfou(i,k)!dwdz
            ufou(i,k)= xa*ufou(i,k)!dudx
            vfou(i,k)= xa*vfou(i,k)!dvdx
            wfou(i,k)= xa*wfou(i,k)!dwdx
         enddo
      enddo
      !NOW GO TO PHYS AGAIN
      call fourxz(wk1,wk1,1,1) !dudz
      call fourxz(wk2,wk2,1,1) !dvdz
      call fourxz(wk3,wk3,1,1) !dwdz
      call fourxz(ufou,ufou,1,1) !dudx
      call fourxz(vfou,vfou,1,1) !dvdx
      call fourxz(wfou,wfou,1,1) !dwdx


      !writing "pp"
      do k=1,mgalz
         do i=1,mgalx
             wk7r(i,k) = -ufour(i,k)-wk5r(i,k)-wk3r(i,k)
         enddo
      enddo
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                i=1,mgalx),k=1,mgalz)


      !computing "qq"
      do k=1,mgalz
         do i=1,mgalx
            wk7r(i,k) = -0.5*(ufour(i,k)**2+wk5r(i,k)**2+
     .                        wk3r(i,k)**2) -
!    du/dy dv/dx + du/dz dw/dx + dv/dz dw/dy
     .      (wk4r(i,k)*vfour(i,k)+wk1r(i,k)*wfour(i,k)+
     .      wk2r(i,k)*wk6r(i,k))  
         enddo
      enddo
      !Now write qq
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),i=1,mgalx),
     .                                               k=1,mgalz)

      !RR
      do k=1,mgalz
         do i=1,mgalx
            wk7r(i,k) = -1/3*(ufour(i,k)**3+wk5r(i,k)**3+
     .                        wk3r(i,k)**3) -
!    du/dy dv/dx du/dx + du/dy dv/dy dv/dx + du/dy dv/dz dw/dx
!    du/dz dw/dx du/dx + dv/dx du/dz dw/dy + du/dz dw/dz dw/dx
!    dv/dy dv/dz dw/dy + dv/dz dw/dz dw/dy 
! --------------------------------------------------+
!    du/dx (du/dy dv/dx + du/dz dw/dx) 
     .    (ufour(i,k)*(wk4r(i,k)*vfour(i,k)+wk1r(i,k)*wfour(i,k))+ 
!    dv/dy (du/dy dv/dx + dv/dz dw/dy) 
     .     wk5r(i,k) *(wk4r(i,k)*vfour(i,k)+wk2r(i,k)*wk6r(i,k))+ 
!    dw/dz (du/dz dw/dx + dv/dz dw/dy) 
     .     wk3r(i,k) *(wk1r(i,k)*wfour(i,k)+wk2r(i,k)*wk6r(i,k))+ 
!    du/dy dv/dz dw/dx + dv/dx du/dz dw/dy) 
     . wk4r(i,k)*wk2r(i,k)*wfour(i,k)+vfour(i,k)*wk1r(i,k)*wk6r(i,k))
         enddo
      enddo
 
      !Write RR
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),i=1,mgalx),
     .                                               k=1,mgalz)

      ! QS
      !Computing QS
      !QS = -1/2(Sii^2)-(S12^2+S13^2+S23^2)
      do k=1,mgalz
         do i=1,mgalx
            wk7r(i,k) = -0.5*(ufour(i,k)**2+wk5r(i,k)**2+
     .                        wk3r(i,k)**2) -
!    S12 = 1/2*(du/dy + dv/dx)
     .    ((1/2*(wk4r(i,k)+vfour(i,k)))**2  + 
!    S23 = 1/2*(dv/dz + dw/dy)
     .     (1/2*(wk2r(i,k)+ wk6r(i,k)))**2  + 
!    S13 = 1/2*(du/dz + dw/dx)
     .     (1/2*(wk1r(i,k)+wfour(i,k)))**2) 
         enddo
      enddo
 
      !Writing QS
      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)

      ! RS
      ! Same as R but with Sij instead of Aij, and is symmetric
      do k=1,mgalz
         do i=1,mgalx
            wk7r(i,k) = -1/3*(ufour(i,k)**3+wk5r(i,k)**3+
     .                        wk3r(i,k)**3) -
! --------------------------------------------------+
!     !S11 (S12**2 + S13**2) 
     . (ufour(i,k)*((1/2*(wk4r(i,k)+vfour(i,k)))**2+
     .              (1/2*(wk1r(i,k)+wfour(i,k)))**2)+ 
!     !S22 (S12**2  + S23**2) 
     .  wk5r(i,k)*((1/2*(wk4r(i,k)+vfour(i,k)))**2+
     .             (1/2*(wk2r(i,k)+ wk6r(i,k)))**2)+ 
!     !S33 (S13**2  + S23**2) 
     .  wk3r(i,k)*((1/2*(wk1r(i,k)+wfour(i,k)))**2+
     .             (1/2*(wk2r(i,k)+ wk6r(i,k)))**2)+ 
!     !2*(S12 S23 S13) 
     .     (wk4r(i,k)+vfour(i,k))*1/2*(wk2r(i,k)+wk6r(i,k))*
     . 1/2*(wk1r(i,k)+wfour(i,k)))
         enddo
      enddo
 


      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)


      ! QW
      !QW = 1/4*Omega_i*Omega_i
      !computing "qq"
      do k=1,mgalz
         do i=1,mgalx
            wk7r(i,k) = 0.25d0*(
             !Vor_x = dw/dy - dv/dz
     .      (wk6r(i,k)-wk2r(i,k))**2+
             !Vor_y = -dw/dx + du/dz
     .      (-wfour(i,k)+wk1r(i,k))**2+
             !Vor_z = dv/dx - du/dy
     .      (vfour(i,k)-wk4r(i,k))**2)
         enddo
      enddo


      kk=kk+1
      write(1000+kk) (wk7r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk7r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk7r(i,k),
     .                               i=1,mgalx),k=1,mgalz)



      end  

