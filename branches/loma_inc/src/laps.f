!=======================================================================
!                                                                      !
!                 	Viscous step solvers Package                   !
!                                                                      !
!======================================================================!

!   General Boundary conditions: 
!   k^2*C_vvv(i)*v(0)   + k*C_dyv(i)*dvdy(0) + C_phi(i)*phi(0) = BC(i), i=1,2
!   k^2*C_vvv(i)*v(1)   + k*C_dyv(i)*dvdy(1) + C_phi(i)*phi(1) = BC(i), i=3,4
!   k  *C_ome(1)*ome(0) +   C_dyo(1)*domedy(0)                 = BC(5)
!   k  *C_ome(2)*ome(1) +   C_dyo(2)*domedy(1)                 = BC(6)
!------------------------------------------------------------------------

!     Additional subroutines (before where functions)

!------------------------------------------------------------------------
      subroutine calcrhslap(f,rhslap)
!     Calculates rhs of laplacian
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhslap
!     ----------------------- Commons  ---------------------------------

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/
      
      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/ 
            
!--------------------------------------------------------------------------! 

      rhslap(1) = 0d0
!      rhslap(1) = dt21(3,1)*f(1)+dt21(4,1)*f(2)+dt21(5,1)*f(3)
      rhslap(2) = dt21(2,2)*f(1)+dt21(3,2)*f(2)+dt21(4,2)*f(3)+
     .            dt21(5,2)*f(4) 

      do j=3,my-2
         rhslap(j) = dt21(1,j)*f(j-2)    
         do i=2,5
            rhslap(j) = rhslap(j)+dt21(i,j)*f(j-3+i)
         enddo
      enddo

      rhslap(my1)=dt21(1,my1)*f(my-3)+dt21(2,my1)*f(my-2)+
     $             dt21(3,my1)*f(my1)+dt21(4,my1)*f(my)
      rhslap(my)  = 0d0
!      rhslap(my)  = dt21(1,my)*f(my-2)+dt21(2,my)*f(my-1)+
!     .              dt21(3,my)*f(my)
      end subroutine calcrhslap


!------------------------------------------------------------------------------

      subroutine calcrhsd1(f,rhsd1)
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhsd1
!     ----------------------- Commons  ---------------------------------

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/
      
      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/ 
!      ----------------------------------------------------------------
        

      rhsd1(1   ) = dt12(4,1   )*f(1   ) + dt12(5,1   )*f(2) 
     $            + dt12(6,1   )*f(3   ) + dt12(7,1   )*f(4)
                     
      rhsd1(2   ) = dt12(3,2   )*f(1   ) + dt12(4,2   )*f(2) 
     $            + dt12(5,2   )*f(3   ) + dt12(6,2   )*f(4) 
     $            + dt12(7,2   )*f(5   )

      rhsd1(3   ) = dt12(2,3   )*f(1   ) + dt12(3,3   )*f(2)  
     $            + dt12(4,3   )*f(3   ) + dt12(5,3   )*f(4) 
     $            + dt12(6,3   )*f(5   ) + dt12(7,3   )*f(6)

      do j=4,my-3
         rhsd1(j)=dt12(1,j)*f(j-3)
         do i=2,7
            rhsd1(j) = rhsd1(j) + dt12(i,j)*f(i+j-4)
         enddo
      enddo
    
      rhsd1(my-2) = dt12(1,my-2)*f(my-5) + dt12(2,my-2)*f(my-4)  
     $             + dt12(3,my-2)*f(my-3) + dt12(4,my-2)*f(my-2)
     $             + dt12(5,my-2)*f(my-1) + dt12(6,my-2)*f(my)    
  
      rhsd1(my-1) = dt12(1,my-1)*f(my-4) + dt12(2,my-1)*f(my-3) 
     $             + dt12(3,my-1)*f(my-2) + dt12(4,my-1)*f(my-1) 
     $             + dt12(5,my-1)*f(my)
          
      rhsd1(my)   = dt12(1,my  )*f(my-3) + dt12(2,my  )*f(my-2) 
     $             + dt12(3,my  )*f(my-1) + dt12(4,my  )*f(my)

      end subroutine calcrhsd1

!***********************************************************************!
!                                                                       !
!        LAPVDV                                                         !
!                                                                       !
!         solve the problem     v'' - rK v = phi                        !
!         form a matrix   ->         {v}*[A] = {f}                      ! 
!         General boundary conditions are allowed:                      !
!    C_vvv(i)*v(0)+C_dyv(i)*dvdy(0)=BC(i), i=1    !
!    C_vvv(i)*v(1)+C_dyv(i)*dvdy(1)=BC(i), i=3    !
!                                                                       !
!                                                                       !
!        Using Compact Finite Difference                                !
!  input:                                                               !
!     phi: vector phi.                                                  !
!       n: number of terms                                              !
!      rK: wavenumber, k=sqrt(rk)
! output:                                                               !
!       v: solution in the plane fourier x fourier                      !
!    dvdy: dvndy in the same plane                                      !
!                                                                       !
!***********************************************************************!
      subroutine Lapvdv(phi,v,dvdy,rK)

      implicit none
      include "ctes3D"

      integer i,j
      real(4),dimension(2,my) :: phi, v, dvdy
      real(8)                 :: rK
      real(4)                 :: zero

      real(8)                 :: coef(4),AA(2,2)

!--------- ----------- BC functions----------------------------
      real(8) :: Cvtop,Cdyvtop,BC_top
      real(8) :: Cvbot,Cdyvbot,BC_bot

! common functions with VISC
      real(8) ::  C_vvv_top1,C_dyv_top1,BC_phi_top1
      real(8) ::  C_vvv_bot1,C_dyv_bot1,BC_phi_bot1

! ----------------------- workspace ---------------------------

      real(8),dimension(my)   ::  vpr, vpi, v1, v2
      real(8),dimension(my)   :: dvpr,dvpi,dv1,dv2

!------------------ commons ------------------------------------
 
      
      real*8 wk1,fwk1,fwk2,fwk,dphi1,dphi2
      common /worksp/wk1(5,my),fwk1(my),fwk2(my),fwk(my),dphi1(my),
     .               dphi2(my)
      save /worksp/


      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/
      
! ----------------------- Program ---------------------------
!     TOP
!     -----------
!      Cvtop=1d0   !couette
!      Cdyvtop=0d0
!      BC_top=0d0
!
! Same BC as VISC 1
      Cvtop=C_vvv_top1(rK) !ml
      Cdyvtop=C_dyv_top1(rK)
      BC_top=BC_phi_top1(rK)

!     BOT
!     ---------
!      Cvbot=1d0 !couette
!      Cdyvbot=0d0
!      BC_bot=0d0

      Cvbot=C_vvv_top1(rK) 
      Cdyvbot=C_dyv_bot1(rK)
      BC_bot=BC_phi_top1(rK)


      zero = 0d0

      if(rK.eq.0) then
        do j=1,my
          v   (1,j) =  zero
          v   (2,j) =  zero
          dvdy(1,j) =  zero
          dvdy(2,j) =  zero
        enddo
       return
      endif

! calculate v (resolve poisson):

      do j=2,my-1
       do i=1,5
          wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
       enddo
      enddo

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

      call bandec(wk1,my)

! resolve v (real coefficients):

      do j=1,my
       v1(j)=phi(1,j)
       v2(j)=phi(2,j)
      enddo

      call calcrhslap(v1,vpr)
      call calcrhslap(v2,vpi) 

      v1(1) = 1d0
      v2(1) = 0d0
      do j=2,my-1
         v1(j)=zero
         v2(j)=zero
      enddo  
      v1(my) = 0d0
      v2(my) = 1d0

      call banbks(wk1,my,vpr)     
      call banbks(wk1,my,vpi)
      call banbks(wk1,my,v1)     
      call banbks(wk1,my,v2)
! Now we have:
! vpr = Real part of vp (particular) (after solving laplacian of phi)
! vpi = Imaginary part of vp (after solving laplacian of phi)
! v1  = vh1 (homogeneus 1)
! v2  = vh2 (homogeneus 2)

      call calcrhsd1(vpr,dvpr)
      call calcrhsd1(vpi, dvpi)
      call calcrhsd1(v1 , dv1 )
      call calcrhsd1(v2 , dv2 )
!preparing for obtaining derivative y of vpr,vpi,v1,v2.   

      call banbks7(prem1,my,dvpr)
      call banbks7(prem1,my,dvpi)
      call banbks7(prem1,my,dv1)
      call banbks7(prem1,my,dv2)

!-----------------------------------------------------------------------
!  Boundary conditions: Only BOT1 and TOP1 
      AA(1,1) = Cvbot + Cdyvbot*dv1(1)*fmap(1)
      AA(1,2) = Cdyvbot*dv2(1)*fmap(1)
      AA(2,1) = Cdyvtop*dv1(my)*fmap(my)
      AA(2,2) = Cvtop + Cdyvtop*dv2(my)*fmap(my)
! REAL:
      coef(1) = BC_bot-Cdyvbot*dvpr( 1)*fmap( 1)
      coef(2) = BC_top-Cdyvtop*dvpr(my)*fmap(my)
! IMAG:
      coef(3) = BC_bot-Cdyvbot*dvpi( 1)*fmap( 1)
      coef(4) = BC_top-Cdyvtop*dvpi(my)*fmap(my)


      call gaussj(AA,2,2,coef,2,2)

      do j=1,my
         v(1,j) = vpr(j) + coef(1)*v1(j) + coef(2)*v2(j)
         v(2,j) = vpi(j) + coef(3)*v1(j) + coef(4)*v2(j)
         dvdy(1,j) =fmap(j)*(dvpr(j) + coef(1)*dv1(j) + coef(2)*dv2(j))
         dvdy(2,j) =fmap(j)*(dvpi(j) + coef(3)*dv1(j) + coef(4)*dv2(j))
      enddo 

      end subroutine Lapvdv

!***********************************************************************!
!                                                                       !
!        LAPPSI00 -Solve d2(psi)/dy^2= phi                              !
!                                                                       !
!         solve the problem     v'' = phi                        !
!         form a matrix   ->       {v}*[A] = {f}                      ! 
!         General boundary conditions are allowed:                      !
!    C_vvv(i)*v(0)+C_dyv(i)*dvdy(0)=BCdyv00 BOT   !
!    C_vvv(i)*v(1)+C_dyv(i)*dvdy(1)=-BCdyv00 TOP  !
!                                                                       !
!        Using Compact Finite Difference                                !
!  input:                                                               !
!     phi: vector phi.                                                  !
!       n: number of terms                                              !
! output:                                                               !
!       v: solution in the plane fourier x fourier                      !
!    dvdy: dvndy in the same plane                                      !
!                                                                       !
!***********************************************************************!
      subroutine Lappsi00(phi,v,dvdy,bcdyv00)

      implicit none
      include "ctes3D"

      integer i,j
      real(4),dimension(2,my) :: phi, v, dvdy
      real(8)                 :: rK
      real(8)                 :: zero
      real(4)                 :: bcdyv00
      real(8)                 :: coef(2),AA(2,2)

! ----------------------- workspace ---------------------------

      real(8),dimension(my)   ::  vpr, vpi, v1, v2
      real(8),dimension(my)   :: dvpr,dvpi,dv1,dv2

!------------------ commons ------------------------------------
 
      
      real*8 wk1,fwk1,fwk2,fwk,dphi1,dphi2
      common /worksp/wk1(5,my),fwk1(my),fwk2(my),fwk(my),dphi1(my),
     .               dphi2(my)
      save /worksp/


      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/
      

! ----------------------- Program ---------------------------
! Remember how CFD scheme works:
! [dt21]*v''=[dt22]f  (1)
! So initial eq:
!v''-rk v = p multiplied byu [dt21]:
![dt21]v''-rk[dt21]v=[dt21]f
! now using (1):
!([dt22]-rk[dt21])v=[dt21]f
!Mode 00 rk=0
      zero = 0d0
      rk=zero
!     TOP
!     -----------
!      Cvtop= zero
!      Cdyvtop= 1d0
!      BC_top= bcdyv00
!     BOT
!     ---------
!      Cvbot=zero 
!      Cdyvbot=1d0
!      BC_bot=-bcdyv00


! calculate v (resolve poisson):
      do j=2,my-1
       do i=1,5
          !wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
          wk1(i,j)=dt22(i,j)
       enddo
      enddo

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

      call bandec(wk1,my)

! resolve v (real coefficients):

      do j=1,my
       v1(j)=phi(1,j)
      enddo

      call calcrhslap(v1,vpr)

      v1(1) = 1d0
      v2(1) = 0d0
      do j=2,my-1
         v1(j)=zero
         v2(j)=zero
      enddo  
      v1(my) = 0d0
      v2(my) = 1d0

      call banbks(wk1,my,vpr)     
      call banbks(wk1,my,v1)     
      call banbks(wk1,my,v2)
!Theoretical solution of hom. solutions
       !v1=0.5*(1+1/y(1)*y);
       !v2=0.5*(1+1/y(my)*y);
! Now we have:
! vpr = Real part of vp (particular) (after solving laplacian of phi)
! vpi = Imaginary part of vp (after solving laplacian of phi)
! v1  = vh1 (homogeneus 1)
! v2  = vh2 (homogeneus 2)


      call calcrhsd1(vpr,dvpr)
      call calcrhsd1(v1 , dv1 )
      call calcrhsd1(v2 , dv2 )
!preparing for obtaining derivative y of vpr,vpi,v1,v2.   

      call banbks7(prem1,my,dvpr)
      call banbks7(prem1,my,dv1)
      call banbks7(prem1,my,dv2)

!-----------------------------------------------------------------------
!trying to put only the particular solution

!  Boundary conditions + Set value at top = zero
!(1) f'=BCbot
!(2) f'=BCtop
!second equation PSI@TOP=0
!mode00 is only  Real
      !Bot
!      AA(1,1) = dv1(1)*fmap(1)
!      AA(1,2) = dv2(1)*fmap(1)
!      !Top
!      AA(2,1) = v1(my)
!      AA(2,2) = v2(my)1

!      coef(1) = -dvpr(1)*fmap(1)-bcdyv00 !bottom negative
!      coef(2) = -vpr(my) + 0d0

!      write(*,*) "bcdyv00",bcdyv00
!debug 
!      write(*,*) "AA",AA(1,1),AA(1,2),AA(2,1),AA(2,2)
!      write(*,*) "coefs",coef(1),coef(2)
!Solve 
!      call gaussj(AA,2,2,coef,1,2)

!      write(*,*) "coefs",coef(1),coef(2)

!Imaginary part Solution
      do j=1,my
!         v(1,j) = vpr(j) + coef(1)*v1(j) + coef(2)*v2(j)
!         dvdy(1,j) =fmap(j)*(dvpr(j) + coef(1)*dv1(j) + 
!     .              coef(2)*dv2(j))
         v(1,j)= vpr(j) !only particular solution
         dvdy(1,j)=fmap(j)*dvpr(j) !only particular solution
         v(2,j) = zero
         dvdy(2,j) =zero
      enddo 
!      !debug
!      write(*,*) "dvdy(1,1)(BOT)",dvdy(1,1)

      end subroutine Lappsi00











!***********************************************************************!
!                                                                       !
!         resuelve el problema  v'' - rK v = phi                        !
!         forma matricial ->         {v}*[A] = {f}                      ! 
!    Atchung!!! Sienmpre suponemos c.c. homogeneas                      !
!                                                                       !
!        DIFERENCIAS FINITAS COMPACTAS                                  !
!  input:                                                               !
!     phi: vector phi.                                                  !
!       n: numero de terminos                                           !
!      rK: constante independiente.                                     !
! output:                                                               !
!       v: solucion en el plano fourier x fourier                       !
!    dvdy: dvndy en el mismo plano anterior                             !
!                                                                       !
!***********************************************************************!
      subroutine Lapvdvhom(phi,v,dvdy,rK)

      implicit none
      include "ctes3D"
            
      integer i,j
      real*4 phi(2,my),v(2,my),dvdy(2,my)
      real*8 rK
      real*4 zero
     
c ------------------ commons ------------------------------------
 
      
      real*8 wk1,fwk1,fwk2,fwk,dphi1,dphi2
      common /worksp/wk1(5,my),fwk1(my),fwk2(my),fwk(my),dphi1(my),
     .               dphi2(my)
      save /worksp/


      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/
      
c ----------------------- Programa ---------------------------
      
      zero = 0.
      
      if(rK.eq.0) then
         do j=1,my
            v   (1,j) =  zero
            v   (2,j) =  zero
            dvdy(1,j) =  zero
            dvdy(2,j) =  zero
         enddo
         return
      endif
      
c      calculo v (resuelvo el poisson):

      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
         enddo
      enddo

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

      call bandec(wk1,my)
      
c     resuelvo  v (coeficientes reales):

      do j=1,my
         dphi1(j)=phi(1,j)
         dphi2(j)=phi(2,j)
      enddo

      fwk1(1) = zero
      fwk2(1) = zero
      
      fwk1(2) = dt21(2,2)*dphi1(1)+dt21(3,2)*dphi1(2)+
     .          dt21(4,2)*dphi1(3)+dt21(5,2)*dphi1(4)
      fwk2(2) = dt21(2,2)*dphi2(1)+dt21(3,2)*dphi2(2)+
     .          dt21(4,2)*dphi2(3)+dt21(5,2)*dphi2(4)    
      
      do j=3,my-2
         fwk1(j)=dt21(1,j)*dphi1(j-2)
         fwk2(j)=dt21(1,j)*dphi2(j-2)    
         do i=2,5
            fwk1(j)=fwk1(j)+dt21(i,j)*dphi1(j-3+i)
            fwk2(j)=fwk2(j)+dt21(i,j)*dphi2(j-3+i)       
         enddo
      enddo

      fwk1(my-1)=dt21(1,my-1)*dphi1(my-3)+dt21(2,my-1)*dphi1(my-2)+
     .           dt21(3,my-1)*dphi1(my-1)+dt21(4,my-1)*dphi1(my)
      fwk2(my-1)=dt21(1,my-1)*dphi2(my-3)+dt21(2,my-1)*dphi2(my-2)+
     .           dt21(3,my-1)*dphi2(my-1)+dt21(4,my-1)*dphi2(my)  
        
      fwk1(my) = zero
      fwk2(my) = zero      

      call banbks(wk1,my,fwk1)     
      call banbks(wk1,my,fwk2)
      
c     guardo la v:

      do j=1,my
         v(1,j) = fwk1(j)
         v(2,j) = fwk2(j)
      enddo


      fwk(1) = dt12(4,1)*fwk1(1)+dt12(5,1)*fwk1(2) + 
     .         dt12(6,1)*fwk1(3)+dt12(7,1)*fwk1(4)
      fwk(2) = dt12(3,2)*fwk1(1)+dt12(4,2)*fwk1(2) +
     .         dt12(5,2)*fwk1(3)+dt12(6,2)*fwk1(4) +
     .         dt12(7,2)*fwk1(5)
      fwk(3) = dt12(2,3)*fwk1(1)+dt12(3,3)*fwk1(2) +
     .         dt12(4,3)*fwk1(3)+dt12(5,3)*fwk1(4) +
     .         dt12(6,3)*fwk1(5)+dt12(7,3)*fwk1(6)

      do j=4,my-3
         fwk(j) = dt12(1,j)*fwk1(j-3)
         do i=2,7
            fwk(j) = fwk(j) + dt12(i,j)*fwk1(i+j-4)
         enddo
      enddo
      fwk(my-2) = dt12(1,my-2)*fwk1(my-5)+dt12(2,my-2)*fwk1(my-4)+
     .            dt12(3,my-2)*fwk1(my-3)+dt12(4,my-2)*fwk1(my-2)+
     .            dt12(5,my-2)*fwk1(my-1)+dt12(6,my-2)*fwk1(my)

      fwk(my-1) = dt12(1,my-1)*fwk1(my-4)+dt12(2,my-1)*fwk1(my-3)+
     .            dt12(3,my-1)*fwk1(my-2)+dt12(4,my-1)*fwk1(my-1)+
     .            dt12(5,my-1)*fwk1(my)

      fwk(my)   = dt12(1,my  )*fwk1(my-3)+dt12(2,my  )*fwk1(my-2)+
     .            dt12(3,my  )*fwk1(my-1)+dt12(4,my  )*fwk1(my) 
  
      call banbks7(prem1,my,fwk)

      do j=1,my
         dvdy(1,j)=fwk(j)*fmap(j)
      enddo

c--------------------------------------------------------------------------
c   Imaginario:

      fwk(1)= dt12(4,1)*fwk2(1)+dt12(5,1)*fwk2(2) +
     .        dt12(6,1)*fwk2(3)+dt12(7,1)*fwk2(4)
      fwk(2)= dt12(3,2)*fwk2(1)+dt12(4,2)*fwk2(2) +
     .        dt12(5,2)*fwk2(3)+dt12(6,2)*fwk2(4) +
     .        dt12(7,2)*fwk2(5)
      fwk(3)= dt12(2,3)*fwk2(1)+dt12(3,3)*fwk2(2) +
     .        dt12(4,3)*fwk2(3)+dt12(5,3)*fwk2(4) +
     .        dt12(6,3)*fwk2(5)+dt12(7,3)*fwk2(6)

      do j=4,my-3
         fwk(j)=dt12(1,j)*fwk2(j-3)
         do i=2,7
            fwk(j)=fwk(j) + dt12(i,j)*fwk2(i+j-4)
         enddo
      enddo
      fwk(my-2)=dt12(1,my-2)*fwk2(my-5)+dt12(2,my-2)*fwk2(my-4)+
     .          dt12(3,my-2)*fwk2(my-3)+dt12(4,my-2)*fwk2(my-2)+
     .          dt12(5,my-2)*fwk2(my-1)+dt12(6,my-2)*fwk2(my)

      fwk(my-1)=dt12(1,my-1)*fwk2(my-4)+dt12(2,my-1)*fwk2(my-3)+
     .          dt12(3,my-1)*fwk2(my-2)+dt12(4,my-1)*fwk2(my-1)+
     .          dt12(5,my-1)*fwk2(my)
      fwk(my)=  dt12(1,my  )*fwk2(my-3)+dt12(2,my  )*fwk2(my-2)+
     .          dt12(3,my  )*fwk2(my-1)+dt12(4,my  )*fwk2(my)

      call banbks7(prem1,my,fwk)

      do j=1,my
         dvdy(2,j)=fwk(j)*fmap(j)
      enddo

      end
      
!***********************************************************************!
!                                                                       !
!         resuelve el problema  v'' - rK v = phi                        !
!         forma matricial ->         {v}*[A] = {f}                      ! 
!    Atchung!!! Sienmpre suponemos c.c. homogeneas                      !
!                                                                       !
!        DIFERENCIAS FINITAS COMPACTAS                                  !
!  input:                                                               !
!     phi: vector phi.                                                  !
!       n: numero de terminos                                           !
!      rK: constante independiente.                                     !
! output:                                                               !
!       v: solucion en el plano fourier x fourier                       !
!    dvdy: dvndy en el mismo plano anterior                             !
!                                                                       !
!***********************************************************************!
      subroutine Lapvdvhom2(phi,v,dvdy,rK)

      implicit none
      include "ctes3D"
            
      integer i,j
      real*4 phi(2,my),v(2,my),dvdy(2,my)
      real*8 rK
      real*4 zero
     
c ------------------ commons ------------------------------------
 
      
      real*8 wk1,fwk1,fwk2,fwk,dphi1,dphi2
      common /worksp/wk1(5,my),fwk1(my),fwk2(my),fwk(my),dphi1(my),
     .               dphi2(my)
      save /worksp/


      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/
      
c ----------------------- Programa ---------------------------
      
      zero = 0.
!Same as Lapvdhom but solving 00 mode     
!      if(rK.eq.0) then
!         do j=1,my
!            v   (1,j) =  zero
!            v   (2,j) =  zero
!            dvdy(1,j) =  zero
!            dvdy(2,j) =  zero
!         enddo
!         return
!      endif
      
c      calculo v (resuelvo el poisson):

      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
         enddo
      enddo

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

      call bandec(wk1,my)
      
c     resuelvo  v (coeficientes reales):

      do j=1,my
         dphi1(j)=phi(1,j)
         dphi2(j)=phi(2,j)
      enddo

      fwk1(1) = zero
      fwk2(1) = zero
      
      fwk1(2) = dt21(2,2)*dphi1(1)+dt21(3,2)*dphi1(2)+
     .          dt21(4,2)*dphi1(3)+dt21(5,2)*dphi1(4)
      fwk2(2) = dt21(2,2)*dphi2(1)+dt21(3,2)*dphi2(2)+
     .          dt21(4,2)*dphi2(3)+dt21(5,2)*dphi2(4)    
      
      do j=3,my-2
         fwk1(j)=dt21(1,j)*dphi1(j-2)
         fwk2(j)=dt21(1,j)*dphi2(j-2)    
         do i=2,5
            fwk1(j)=fwk1(j)+dt21(i,j)*dphi1(j-3+i)
            fwk2(j)=fwk2(j)+dt21(i,j)*dphi2(j-3+i)       
         enddo
      enddo

      fwk1(my-1)=dt21(1,my-1)*dphi1(my-3)+dt21(2,my-1)*dphi1(my-2)+
     .           dt21(3,my-1)*dphi1(my-1)+dt21(4,my-1)*dphi1(my)
      fwk2(my-1)=dt21(1,my-1)*dphi2(my-3)+dt21(2,my-1)*dphi2(my-2)+
     .           dt21(3,my-1)*dphi2(my-1)+dt21(4,my-1)*dphi2(my)  
        
      fwk1(my) = zero
      fwk2(my) = zero      

      call banbks(wk1,my,fwk1)     
      call banbks(wk1,my,fwk2)
      
c     guardo la v:

      do j=1,my
         v(1,j) = fwk1(j)
         v(2,j) = fwk2(j)
      enddo


      fwk(1) = dt12(4,1)*fwk1(1)+dt12(5,1)*fwk1(2) + 
     .         dt12(6,1)*fwk1(3)+dt12(7,1)*fwk1(4)
      fwk(2) = dt12(3,2)*fwk1(1)+dt12(4,2)*fwk1(2) +
     .         dt12(5,2)*fwk1(3)+dt12(6,2)*fwk1(4) +
     .         dt12(7,2)*fwk1(5)
      fwk(3) = dt12(2,3)*fwk1(1)+dt12(3,3)*fwk1(2) +
     .         dt12(4,3)*fwk1(3)+dt12(5,3)*fwk1(4) +
     .         dt12(6,3)*fwk1(5)+dt12(7,3)*fwk1(6)

      do j=4,my-3
         fwk(j) = dt12(1,j)*fwk1(j-3)
         do i=2,7
            fwk(j) = fwk(j) + dt12(i,j)*fwk1(i+j-4)
         enddo
      enddo
      fwk(my-2) = dt12(1,my-2)*fwk1(my-5)+dt12(2,my-2)*fwk1(my-4)+
     .            dt12(3,my-2)*fwk1(my-3)+dt12(4,my-2)*fwk1(my-2)+
     .            dt12(5,my-2)*fwk1(my-1)+dt12(6,my-2)*fwk1(my)

      fwk(my-1) = dt12(1,my-1)*fwk1(my-4)+dt12(2,my-1)*fwk1(my-3)+
     .            dt12(3,my-1)*fwk1(my-2)+dt12(4,my-1)*fwk1(my-1)+
     .            dt12(5,my-1)*fwk1(my)

      fwk(my)   = dt12(1,my  )*fwk1(my-3)+dt12(2,my  )*fwk1(my-2)+
     .            dt12(3,my  )*fwk1(my-1)+dt12(4,my  )*fwk1(my) 
  
      call banbks7(prem1,my,fwk)

      do j=1,my
         dvdy(1,j)=fwk(j)*fmap(j)
      enddo

c--------------------------------------------------------------------------
c   Imaginario:

      fwk(1)= dt12(4,1)*fwk2(1)+dt12(5,1)*fwk2(2) +
     .        dt12(6,1)*fwk2(3)+dt12(7,1)*fwk2(4)
      fwk(2)= dt12(3,2)*fwk2(1)+dt12(4,2)*fwk2(2) +
     .        dt12(5,2)*fwk2(3)+dt12(6,2)*fwk2(4) +
     .        dt12(7,2)*fwk2(5)
      fwk(3)= dt12(2,3)*fwk2(1)+dt12(3,3)*fwk2(2) +
     .        dt12(4,3)*fwk2(3)+dt12(5,3)*fwk2(4) +
     .        dt12(6,3)*fwk2(5)+dt12(7,3)*fwk2(6)

      do j=4,my-3
         fwk(j)=dt12(1,j)*fwk2(j-3)
         do i=2,7
            fwk(j)=fwk(j) + dt12(i,j)*fwk2(i+j-4)
         enddo
      enddo
      fwk(my-2)=dt12(1,my-2)*fwk2(my-5)+dt12(2,my-2)*fwk2(my-4)+
     .          dt12(3,my-2)*fwk2(my-3)+dt12(4,my-2)*fwk2(my-2)+
     .          dt12(5,my-2)*fwk2(my-1)+dt12(6,my-2)*fwk2(my)

      fwk(my-1)=dt12(1,my-1)*fwk2(my-4)+dt12(2,my-1)*fwk2(my-3)+
     .          dt12(3,my-1)*fwk2(my-2)+dt12(4,my-1)*fwk2(my-1)+
     .          dt12(5,my-1)*fwk2(my)
      fwk(my)=  dt12(1,my  )*fwk2(my-3)+dt12(2,my  )*fwk2(my-2)+
     .          dt12(3,my  )*fwk2(my-1)+dt12(4,my  )*fwk2(my)

      call banbks7(prem1,my,fwk)

      do j=1,my
         dvdy(2,j)=fwk(j)*fmap(j)
      enddo

      end
      








 
    
!***********************************************************************!
!           LAPV1                                                       !
!    THIS SUBROUTINE IS NOT USED                                        !
!           solves   u'' - rK u = f                                     !
!                                                                       !
!                    con u(-1) = u(+1)= 0                               !
!aaf   added homogeneus solutions for general BC                        !
!aaf   similar to Lapvdv subroutine                                     !
!                                                                       !
!      Using Compact Finite Difference                                  !
!  input:                                                               !
!            f : forzante                                               !
!            rK: constant                                               !
!    dalbe,dtri: constant                                               !
!           ind: Logical FLAG                                           !
!                                                                       !
! output:                                                               !
!       f: Solution of system                                           !
!                                                                       !
!       if dalbe=dtri=0.0                                               !
!       u: solution of the system                                       !
!       if dalbe ~=0                                                    !
!       u: solution + some things prepared for SMR                      !
!                                                                       !
!***********************************************************************!

      subroutine Lapv1(f,u,rK,dalbe,dtri,ind)
 
      implicit none
      include "ctes3D" 
!    ----------------------- IN & OUT ----------------- 
    
      real(8),dimension(my),intent(inout) :: f, u
      real(8),    intent(in)              :: rk,dalbe,dtri
      integer(4), intent(in)              :: ind

!     ----------------------- Commons -------------------
     
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .            trp(0:my-1),mss(0:my-1)
      save   /fis/
!aaf  if possible: trp, mss, etc... should be removed from
!     /fis/ common
 
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

!aaf      real*8 wk1,fwk,df
!      common /worksp/ wk1(5,my),fwk(my),df(my)
!      save  /worksp/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/

!--------- ----------- BC functions----------------------------
      real(8) :: C_phi_top1, C_vvv_top1,C_dyv_top1,BC_phi_top1
      real(8) :: C_phi_top2, C_vvv_top2,C_dyv_top2,BC_phi_top2
      real(8) :: C_phi_bot1, C_vvv_bot1,C_dyv_bot1,BC_phi_bot1
      real(8) :: C_phi_bot2, C_vvv_bot2,C_dyv_bot2,BC_phi_bot2
      real(8) :: C_ome_top,  C_dyo_top,BC_ome_top
      real(8) :: C_ome_bot,  C_dyo_bot,BC_ome_bot
!     ----------------------- work --------------------     

      integer(4)              ::  i,j
      real(8),dimension(my)   ::  fwk,df,fwkh1,fwkh2,dfwk,dfwkh1,dfwkh2
      real(8)                 ::  coef(2),AA(2,2)
      real(8)                 ::  wk1(5,my)


!     ----------------------- Program ------------------      
!     calculate v (set up poisson) 
!     Construct LHS
!     not dealing with boundary here
      do j=2,my-1
        do i=1,5
          wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
        enddo
      enddo

!     boundary here (wall)
      wk1(1,1)  = 0d0
      wk1(2,1)  = 0d0
      wk1(3,1)  = 1d0
      wk1(4,1)  = 0d0  !aaf
      wk1(5,1)  = 0d0  !aaf
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0

!numerical recipes codebase--
!given nxn band diagonal matrix wk1, construct LU decomposition
!used in conjuction with banbks to solve band-diagonal sets of equations
      call bandec(wk1,my)

!     resolve v (real coefficients):

      do j=1,my
         df(j)=f(j)
      enddo

! ----------------------------------------------------------------------!
      call calcrhslap(df,fwk)
!     construct RHS
!-----------------------------------------------------------------------!

!given arrays (from bandec) and rhs, 
!solves band diagonal linear equations Ax=b --returning b (as fwk)
!            A      b
      call banbks(wk1,my,fwk)

!    TO BE COMPLETED
!-----------------------------------------------------------------------!
!     solve homogeneous probl x2
!     Creating homogeneus vectors (BC1-fwkh1 and BC2-fwkh2)
      fwkh1(1) = 1d0
      fwkh2(1) = 0d0
      do j=2,my-1
         fwkh1(j)=0d0
         fwkh2(j)=0d0
      enddo  
      fwkh1(my) = 0d0
      fwkh2(my) = 1d0

!     Caculating the derivative:
!     preparing rhs 
      call calcrhsd1(fwk,dfwk)
      call calcrhsd1(fwkh1,dfwkh1)
      call calcrhsd1(fwkh2,dfwkh2)

!     derivative in y
      call banbks7(prem1,my,dfwk)
      call banbks7(prem1,my,dfwkh1)
      call banbks7(prem1,my,dfwkh2)

!     Building the A matrix
!     part multiplied by coef(1) on the left for BC
!     at bottom
      AA(1,1)=C_vvv_bot1(rk)+C_dyv_bot1(rk)*dfwkh1(1)*fmap(1)

!     part multiplied by coef(2) on the left for BC
!     at bottom
      AA(1,2)=C_dyv_bot1(rk)*dfwkh2(1)*fmap(1)

!     same for top BC
      AA(2,1)=C_dyv_top1(rk)*dfwkh1(my)*fmap(my)
      AA(2,2)=C_vvv_top1(rk) + C_dyv_top1(rk)*dfwkh2(my)*fmap(my)

!     RHS of linear sistem A · coef' = B
      coef(1)=BC_phi_bot1(rk)-C_dyv_bot1(rk)*dfwk( 1)*fmap(1 )
      coef(2)=BC_phi_top1(rk)-C_dyv_top1(rk)*dfwk(my)*fmap(my)

!     solve system
      call gaussj(AA,2,2,coef,2,1) !check with oscar
    
!     sum all solutions for fwk
      do j=1,my
         fwk(j)=fwk(j)+ coef(1)*fwkh1(j) + coef(2)*fwkh2(j)
      enddo

!aaf should be check
!     here, fwk = u00, df = rhs for u00
!special step
!----------------------------------------------------------------------!

      if (ind.eq.1) then
        do j=1,my
           u(j) = dalbe*fwk(j)+dtri*df(j)
        enddo
      else
        do j=1,my
           u(j) = fwk(j)
        enddo
      endif

      do j=1,my
         f(j) = fwk(j)
      enddo

      endsubroutine

!########################################################################!
!#               Subroutine Visc                                         !
!#                                                                       ! 
!#          Standalone routine--not in wavespace module                  !
!#        currently necessary for inhomogeneous operations               !
!#                                                                       !
!#      Solves the problems                                              !
!#                                                                       !
!#      The omega equation                                               !  
!#                           vor'' - rK1 vor = g                         !
!#      AND                                                              !
!#                                                                       !
!#      The phi/v equation                                               ! 
!#                           phi'' - rK1 phi = f                         !
!#                           v''   - rK2 v   = phi                       !
!#                                                                       !
!#      with GENERAL boundary condtions:                                 ! 
!#  C_phi(i)*phi(0)+k^2*C_vvv(1)*v(0)+k*C_dyv(i)*dvdy(0) = BC(i), i=1,2  !
!#  C_phi(i)*phi(1)+k^2*C_vvv(1)*v(1)+k*C_dyv(i)*dvdy(1) = BC(i), i=3,4  !
!#                 k*C_ome(1)*ome(0) + C_dyo(1)*domedy(0)= BC(5)         !
!#                 k*C_ome(2)*ome(1) + C_dyo(2)*domedy(1)= BC(6)         !
!#                                                                       !
!#      Using:  Compact Finite Difference                                !
!#                                                                       !
!#  input:                                                               !
!#           vorwk: forcing for omega (omega*)                           !
!#           phiwk: forcing for phi   (phi*)                             !
!#           v    : NL terms for omega                                   !
!#           dvdy : NL terms for vor                                     !
!#            rK2 : wavenumber, k=sqrt(rk1)                              !
!#            rK1 : constant (wavenumber + Re/(dt beta))                 !
!#      dalbe,dtri: constants                                            !
!#                                                                       !
!# output:                                                               !
!#           vorwk: vor                                                  !
!#           phiwk: phi                                                  !
!#           v    : v                                                    !
!#           dvdy : dvdy                                                 !
!#           vor  : vor + dt*alpha/Re*nabla vor + chi*dt*NL (dable~=0)   !
!#                  vor                                     (dable==0)   !
!#           phi  : phi + dt*alpha/Re*nabla phi + chi*dt*NL (dable~=0)   !
!#                  phi                                     (dable==0)   !
!#                                                                       !
!#                                                                       !
!########################################################################!
      subroutine visc(vorwk, vor, phi, v, dvdy, phiwk, rk1, rk2, dalbe,
     $                dtri, dtxi,flag)
    
      implicit none 
      include "ctes3D"
!aaf variables set to real 4
!                           i/o     o     o   i/o  i/o   i/o
      real(4),dimension(2,my) :: vorwk,  vor,  phi,   v, dvdy,phiwk
      real(8)                 :: rk1,rk2,dalbe,dtri,dtxi !coefficients

         
!      !These blocks are used as workspace / are private to routine
      real(8)                 :: zero
      real(8),dimension(5,my) :: wk1
!      real(8),dimension(5,my) :: wk1,wk2
      real(8),dimension(my)   :: phipr,phipi,phi1,phi2,
     .                           vpr,vpi,v1,v2,v3,v4
      real(8),dimension(my)   :: dvpr,dvpi,dv1,dv2,dv3,dv4
      integer(4)              :: n,i,j,flag      !indicies
      real(8)                 :: coef(8), AA(4,4)


!     ----------------------- Commons  ---------------------------------

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/
      
      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/ 
       
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .            trp(0:my-1),mss(0:my-1)
      save   /fis/  

      
!--------- ----------- BC functions----------------------------
      real(8) ::  C_phi_top1, C_vvv_top1,C_dyv_top1,BC_phi_top1
      real(8) ::  C_phi_top2, C_vvv_top2,C_dyv_top2,BC_phi_top2
      real(8) ::  C_phi_bot1, C_vvv_bot1,C_dyv_bot1,BC_phi_bot1
      real(8) ::  C_phi_bot2, C_vvv_bot2,C_dyv_bot2,BC_phi_bot2
      real(8) ::  C_ome_top,  C_dyo_top,BC_ome_top
      real(8) ::  C_ome_bot,  C_dyo_bot,BC_ome_bot


!     ! DO WE NEED ALL THESE COMMONS?? aalmagro     
!      real*8 wk1,phipr,phipi,vpr,vpi,dvpr,dvpi,phi1,v1,dv1,dfr,dfi,
!     .       dphi1,dphi2,fwk1,fwk2
!      common /worksp/ wk1(5,my),fwk1(my),fwk2(my),dphi1(my),dphi2(my),
!     .                phipr(my),phipi(my),vpr(my),vpi(my),
!     .    dvpr(my),dvpi(my),phi1(my),v1(my),dv1(my),dfr(my),dfi(my)
            
!     ----------------------- Program ----------------------------------
! build wk1 operator
      zero = 0d0

      wk1(1,1)  = 0d0
      wk1(2,1)  = 0d0
      wk1(3,1)  = 1d0
      wk1(4,1)  = 0d0
      wk1(5,1)  = 0d0
      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk1*dt21(i,j)
         enddo
      enddo
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0

      call bandec(wk1,my)


!     First, solve omega (only needs wk1 with rk1, homogeneous solutions are used 
!                         later for phi-v bi-laplacian).       

      do j=1,my  ! copy forcing (vor*) into vor
         vor(1,j) = vorwk(1,j)
         vor(2,j) = vorwk(2,j) 
      enddo 
      
      do j=1,my
         phi1(j)=vorwk(1,j)
         phi2(j)=vorwk(2,j)
      enddo

      call calcrhslap(phi1,phipr)
      call calcrhslap(phi2,phipi)

      phi1(1 ) = 1d0
      phi2(1 ) = 0d0
      do j=2,my-1
         phi1(j) = zero
         phi2(j) = zero
      enddo
      phi1(my) = 0d0 
      phi2(my) = 1d0 
      
      call banbks(wk1,my,phipr)  ! omega particular solution  real 
      call banbks(wk1,my,phipi)  ! omega particular solution  imag 
      call banbks(wk1,my,phi1)   ! omega homogeneous solution 1
      call banbks(wk1,my,phi2)   ! omega homogeneous solution 2

!     Derivative of rhs
     
      call calcrhsd1(phipr,dvpr)
      call calcrhsd1(phipi,dvpi)
      call calcrhsd1(phi1,dv1)
      call calcrhsd1(phi2,dv2)

      call banbks7(prem1,my,dvpr)! d/dy omega particular solution  re
      call banbks7(prem1,my,dvpi)! d/dy omega particular solution  im
      call banbks7(prem1,my,dv1) ! d/dy omega homogeneous solution 1
      call banbks7(prem1,my,dv2) ! d/dy omega homogeneous solution 2

! boundary conditions for omega: 
! k*C_ome(1)*ome(0) + C_dyo(1)*dome(0) = BC(5)
! k*C_ome(2)*ome(1) + C_dyo(2)*dome(1) = BC(6)
      AA(1,1) = C_ome_bot(rk2) + C_dyo_bot(rk2)*dv1(1)*fmap(1)
      AA(1,2) = C_dyo_bot(rk2)*dv2(1)*fmap(1)
      AA(2,1) = C_dyo_top(rk2)*dv1(my)*fmap(my)
      AA(2,2) = C_ome_top(rk2) + C_dyo_top(rk2)*dv2(my)*fmap(my)
!  REAL:
      coef(1) = BC_ome_bot(rk2)-C_dyo_bot(rk2)*dvpr( 1)*fmap( 1)
      coef(2) = BC_ome_top(rk2)-C_dyo_top(rk2)*dvpr(my)*fmap(my)
!  IMAG:
      coef(5) = BC_ome_bot(rk2)-C_dyo_bot(rk2)*dvpi( 1)*fmap( 1)
      coef(6) = BC_ome_top(rk2)-C_dyo_top(rk2)*dvpi(my)*fmap(my)
 
      call gaussj(AA,2,4,coef,2,2)

! keep omega
      do j=1,my
        vorwk(1,j) = phipr(j) + coef(1)*phi1(j) + coef(2)*phi2(j)
        vorwk(2,j) = phipi(j) + coef(5)*phi1(j) + coef(6)*phi2(j)
      enddo
!     ! return solution omega + some stuff
      if (dalbe.ne.0d0) then
! vor is the rhs(istep), that is equal to LHS of equation
! therefore it gives the rhs for the next istep
        do j=1,my
           vor(1,j) = dalbe*vorwk(1,j)+dtri*vor(1,j)+dtxi*v(1,j)
           vor(2,j) = dalbe*vorwk(2,j)+dtri*vor(2,j)+dtxi*v(2,j)
        enddo
      else   ! just return omega
        do j=1,my
          vor(1,j) = vorwk(1,j)
          vor(2,j) = vorwk(2,j)
        enddo
      endif
            
!     Second, solve phi & v. Homogeneous solutions have already been computed (phi1,phi2). 
!
!     PHI: 
      do j=1,my  ! copy forcing (phi*) into phi
         phi(1,j) = phiwk(1,j)
         phi(2,j) = phiwk(2,j) 
      enddo 

!aaf  Here we are changing from real 4 to real8
      do j=1,my
         vpr(j)=phiwk(1,j)
         vpi(j)=phiwk(2,j)
      enddo

      call calcrhslap(vpr,phipr)
      call calcrhslap(vpi,phipi)

      call banbks(wk1,my,phipr)      
      call banbks(wk1,my,phipi)
!
!     V: 

      wk1(1,1)  = 0d0
      wk1(2,1)  = 0d0
      wk1(3,1)  = 1d0
      wk1(4,1)  = 0d0
      wk1(5,1)  = 0d0
      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk2*dt21(i,j)
         enddo
      enddo
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0
    
      call bandec(wk1,my)

      call calcrhslap(phipr,vpr)
      call calcrhslap(phipi,vpi)
      call calcrhslap(phi1, v1 )
      call calcrhslap(phi2, v2 )
!FLAG
!      open(unit=356,file='test_vpr.dat',form='formatted')
!      do j= 1,my
!       write(356,*) vpr(j),vpi(j)
!      enddo
      

      v3(1) = 1d0
      v4(1) = 0d0 
      do j=2,my-1
         v3(j) = 0d0 
         v4(j) = 0d0 
      enddo
      v3(my) = 0d0 
      v4(my) = 1d0 


      call banbks(wk1,my,vpr) 
      call banbks(wk1,my,vpi) 
      call banbks(wk1,my,v1) 
      call banbks(wk1,my,v2) 
      call banbks(wk1,my,v3) 
      call banbks(wk1,my,v4) 
!FLAG
!      do j= 1,my
!       write(356,*) vpr(j),vpi(j)
!      enddo
!      close(356)
      
!
!     dVdy: 
      call calcrhsd1(vpr,dvpr)      
      call calcrhsd1(vpi,dvpi)      
      call calcrhsd1(v1 ,dv1 )      
      call calcrhsd1(v2 ,dv2 )      
      call calcrhsd1(v3 ,dv3 )      
      call calcrhsd1(v4 ,dv4 )      

      call banbks7(prem1,my,dvpr)
      call banbks7(prem1,my,dvpi)
      call banbks7(prem1,my,dv1)
      call banbks7(prem1,my,dv2)
      call banbks7(prem1,my,dv3)
      call banbks7(prem1,my,dv4)

! boundary conditions for phi/v: 
! C_phi(1)*phi(0) + k^2*C_vvv(1)*v(0) + k*C_dyv(1)*dv(0) = BC(1) -> AA(:,1)
! C_phi(2)*phi(0) + k^2*C_vvv(2)*v(0) + k*C_dyv(2)*dv(0) = BC(2) -> AA(:,2)
! C_phi(3)*phi(1) + k^2*C_vvv(3)*v(1) + k*C_dyv(3)*dv(1) = BC(3) -> AA(:,3)
! C_phi(4)*phi(1) + k^2*C_vvv(4)*v(1) + k*C_dyv(4)*dv(1) = BC(4) -> AA(:,4)
!aaf term for vh1/vh2 is missing?
      AA(1,1) = C_phi_bot1(rk2) +    C_dyv_bot1(rk2)*dv1(1)*fmap(1)
      AA(1,2) =                      C_dyv_bot1(rk2)*dv2(1)*fmap(1)
      AA(1,3) = C_vvv_bot1(rk2) +    C_dyv_bot1(rk2)*dv3(1)*fmap(1)
      AA(1,4) =                      C_dyv_bot1(rk2)*dv4(1)*fmap(1)

      AA(2,1) = C_phi_bot2(rk2) +    C_dyv_bot2(rk2)*dv1(1)*fmap(1)    
      AA(2,2) =                      C_dyv_bot2(rk2)*dv2(1)*fmap(1)
      AA(2,3) = C_vvv_bot2(rk2) +    C_dyv_bot2(rk2)*dv3(1)*fmap(1)
      AA(2,4) =                      C_dyv_bot2(rk2)*dv4(1)*fmap(1)

      AA(3,1) =                      C_dyv_top1(rk2)*dv1(my)*fmap(my)
      AA(3,2) = C_phi_top1(rk2) +    C_dyv_top1(rk2)*dv2(my)*fmap(my)
      AA(3,3) =                      C_dyv_top1(rk2)*dv3(my)*fmap(my)
      AA(3,4) = C_vvv_top1(rk2) +    C_dyv_top1(rk2)*dv4(my)*fmap(my)

      AA(4,1) =                      C_dyv_top2(rk2)*dv1(my)*fmap(my)
      AA(4,2) = C_phi_top2(rk2) +    C_dyv_top2(rk2)*dv2(my)*fmap(my)
      AA(4,3) =                      C_dyv_top2(rk2)*dv3(my)*fmap(my)
      AA(4,4) = C_vvv_top2(rk2) +    C_dyv_top2(rk2)*dv4(my)*fmap(my)

!  REAL:
      coef(1) = BC_phi_bot1(rk2) - C_dyv_bot1(rk2)*dvpr( 1)*fmap( 1)
      coef(2) = BC_phi_bot2(rk2) - C_dyv_bot2(rk2)*dvpr( 1)*fmap( 1)
      coef(3) = BC_phi_top1(rk2) - C_dyv_top1(rk2)*dvpr(my)*fmap(my)
      coef(4) = BC_phi_top2(rk2) - C_dyv_top2(rk2)*dvpr(my)*fmap(my)
!  IMAG:
      coef(5) = BC_phi_bot1(rk2) - C_dyv_bot1(rk2)*dvpi( 1)*fmap( 1)
      coef(6) = BC_phi_bot2(rk2) - C_dyv_bot2(rk2)*dvpi( 1)*fmap( 1)
      coef(7) = BC_phi_top1(rk2) - C_dyv_top1(rk2)*dvpi(my)*fmap(my)
      coef(8) = BC_phi_top2(rk2) - C_dyv_top2(rk2)*dvpi(my)*fmap(my)

      call gaussj(AA,4,4,coef,2,2)
!     ! keep phi
      do j=1,my
         phiwk(1,j) = phipr(j) + coef(1)*phi1(j) + coef(2)*phi2(j)
         phiwk(2,j) = phipi(j) + coef(5)*phi1(j) + coef(6)*phi2(j)
      enddo 
! return solution phi + some stuff
      if (dalbe.ne.0d0) then
         do j=1,my
           phi(1,j) = dalbe*phiwk(1,j)+dtri*phi(1,j)+dtxi*dvdy(1,j)
           phi(2,j) = dalbe*phiwk(2,j)+dtri*phi(2,j)+dtxi*dvdy(2,j)
         enddo
      else  
!just return phi
        do j=1,my
           phi(1,j) = phiwk(1,j)
           phi(2,j) = phiwk(2,j)
        enddo
      endif
    
!special step -- cfl changing next step
      if (rk2.ne.0d0) then
         do j=1,my 
            dvdy(1,j)=(dvpr(j)+coef(1)*dv1(j)+coef(2)*dv2(j)+
     .                 coef(3)*dv3(j)+coef(4)*dv4(j))*fmap(j)
            dvdy(2,j)=(dvpi(j)+coef(5)*dv1(j)+coef(6)*dv2(j)+
     .                 coef(7)*dv3(j)+coef(8)*dv4(j))*fmap(j)
               v(1,j)=  vpr(j)+coef(1)* v1(j)+coef(2)* v2(j)+
     .                 coef(3)* v3(j)+coef(4)* v4(j)
               v(2,j)=  vpi(j)+coef(5)* v1(j)+coef(6)* v2(j)+
     .                 coef(7)* v3(j)+coef(8)* v4(j)
         enddo 
      else
         do j=1,my
            dvdy(1,j) = zero
            dvdy(2,j) = zero
               v(1,j) = zero  
               v(2,j) = zero          
         enddo
      endif
    
      end subroutine visc



!-----------------------------------------------------------------------------!
! FUNCTIONS FOR BOUNDARY CONDITIONS
!
!     Declaring constansts for general BC's - ROBIN                                  
!     k^2*C_vvv(i)*v(0)   + k*C_dyv(i)*dvdy(0) + C_phi(i)*phi(0) = BC(i), i=1,2      
!     k^2*C_vvv(i)*v(1)   + k*C_dyv(i)*dvdy(1) + C_phi(i)*phi(1) = BC(i), i=3,4      
!     k  *C_ome(1)*ome(0) +   C_dyo(1)*domedy(0)                 = BC(5)             
!     k  *C_ome(2)*ome(1) +   C_dyo(2)*domedy(1)                 = BC(6)             

!   NOTE   : If only one C_phi is not zero, that eq should be i=2,4
!            (for compatibility with Lapvdv)
!   NOTE(2): The wavenumber k enters in the BC to ensure dimensional
!            correctness between the different terms 
!   NOTE(3): C_phi,C_vvv,C_dyv,C_ome,C_dyo and BC are real, so no coupling 
!            between real and imaginary parts is allowed
                                                                       
!      real(8) C_phi(4),C_vvv(4),C_dyv(4),C_ome(2),C_dyo(2),BC(6)                    
!      common /ctesbc/ C_phi,C_vvv,C_dyv,C_ome,C_dyo,BC                              
!    !     phi/v:   bot    top                                                       
!      data C_phi/   1, 0,   1, 0/
!      data C_vvv /  0, 1,   0, 1/
!      data C_dyv /  0, 1,   0, 1/
!    !     ome                  bot top
!      data C_ome /                1, 1/
!      data C_dyo /                0, 0/
!!    !     rhs
!      data BC    /  0, 0,  0, 0, 0,  0/
!      save   /ctesbc/
! --------------------------------------------------------------------------------
!
! PHI/V TOP EQ 1 - BC SHARED BY VISC and LAPVDV
!     k^2*C_vvv(i)*v(0)   + k*C_dyv(i)*dvdy(0) + C_phi(i)*phi(0) = BC(i), i=1,2      
      function C_phi_top1(rk2)
      real(8) :: rk2
      real(8) :: C_phi_top1
      C_phi_top1=0d0 !mixing layer
!      C_phi_top1=0d0 !Channel
      end function

      function C_vvv_top1(rk2)
      real(8) :: rk2
      real(8) :: C_vvv_top1
      C_vvv_top1=sqrt(rk2) !mixing layer
!      C_vvv_top1=1d0
      end function

      function C_dyv_top1(rk2)
      real(8) :: rk2
      real(8) :: C_dyv_top1
!      C_dyv_top1=0d0 !m.l.test
      C_dyv_top1=1d0 !m.l.
!      C_dyv_top1=0d0  !channel
      end function

      function BC_phi_top1(rk2)
      real(8) :: rk2
      real(8) ::  BC_phi_top1
      BC_phi_top1=0d0 !m.l.
      end function
! PHI TOP EQ 2
!     k^2*C_vvv(i)*v(0)   + k*C_dyv(i)*dvdy(0) + C_phi(i)*phi(0) = BC(i), i=1,2      
      function C_phi_top2(rk2)
      real(8) :: rk2
      real(8) :: C_phi_top2
      C_phi_top2=1d0 !m.l.
      end function

      function C_vvv_top2(rk2)
      real(8) :: rk2
      real(8) :: C_vvv_top2
      C_vvv_top2=0d0 !m.l.
!      C_vvv_top2=0d0  !channel
      end function

      function C_dyv_top2(rk2)
      real(8) :: rk2
      real(8) :: C_dyv_top2
      C_dyv_top2=0d0 !m-l
!      C_dyv_top2=1  !channel
      end function

      function BC_phi_top2(rk2)
      real(8) :: rk2
      real(8) ::  BC_phi_top2
      BC_phi_top2=0d0
      end function
!
!------------------------------------------------------------------------------

! PHI/V BOT EQ 1 --BC shared by VISC and LAPVDV
!    C_vvv(i)*v(1)   + C_dyv(i)*dvdy(1) + C_phi(i)*phi(1) = BC(i), i=3,4      
  
      function C_phi_bot1(rk2)
      real(8) :: rk2
      real(8) :: C_phi_bot1
      C_phi_bot1=0d0 !m-l
      end function

      function C_vvv_bot1(rk2)
      real(8) :: rk2
      real(8) :: C_vvv_bot1
      C_vvv_bot1=sqrt(rk2) !m-l
!      C_vvv_bot1=1d0
      end function

      function C_dyv_bot1(rk2)
      real(8) :: rk2
      real(8) :: C_dyv_bot1
!      C_dyv_bot1=0d0 !m-l test
      C_dyv_bot1=-1d0 !m-l
!     If y<0 is different to y>0
      end function

      function BC_phi_bot1(rk2)
      real(8) :: rk2
      real(8) ::  BC_phi_bot1
      BC_phi_bot1=0d0 !m-l
      end function
!

! PHI BOT EQ 2
!     k^2*C_vvv(i)*v(1)   + k*C_dyv(i)*dvdy(1) + C_phi(i)*phi(1) = BC(i), i=3,4      
      function C_phi_bot2(rk2)
      real(8) :: rk2
      real(8) :: C_phi_bot2
      C_phi_bot2=1d0 !m-l
!       C_phi_bot2=0d0  !channel
      end function

      function C_vvv_bot2(rk2)
      real(8) :: rk2
      real(8) :: C_vvv_bot2
      C_vvv_bot2=0d0
      end function

      function C_dyv_bot2(rk2)
      real(8) :: rk2
      real(8) :: C_dyv_bot2
      C_dyv_bot2=0d0 !m-l
!      C_dyv_bot2=1d0
      end function

      function BC_phi_bot2(rk2)
      real(8) :: rk2
      real(8) ::  BC_phi_bot2
      BC_phi_bot2=0d0
      end function
!


!
! OMEGA TOP FUNCTION
!     k  *C_ome(1)*ome(0) +   C_dyo(1)*domedy(0)                 = BC(5)             
      function C_ome_top(rk2)
      real(8) :: rk2
      real(8) :: C_ome_top
      C_ome_top=1d0
      end function 

      function C_dyo_top(rk2)
      real(8) :: rk2
      real(8) :: C_dyo_top
      C_dyo_top=0d0
      end function 

      function BC_ome_top(rk2)
      real(8) :: rk2
      real(8) :: BC_ome_top
      BC_ome_top=0d0
      end function 
! OMEGA BOT FUNCTION      
!     k  *C_ome(2)*ome(1) +   C_dyo(2)*domedy(1)                 = BC(6)             
      function C_ome_bot(rk2)
      real(8) :: rk2
      real(8) :: C_ome_bot
      C_ome_bot=1d0
      end function 

      function C_dyo_bot(rk2)
      real(8) :: rk2
      real(8) :: C_dyo_bot
      C_dyo_bot=0d0
      end function 

      function BC_ome_bot(rk2)
      real(8) :: rk2
      real(8) :: BC_ome_bot
      BC_ome_bot=0d0
      end function 
       
 
