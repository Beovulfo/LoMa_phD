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
      use matrices
      
!     Calculates rhs of laplacian
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhslap
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

!------------------------------------------------------------------------
      subroutine calcrhsfree(f,rhslap)
      use matrices
!     Calculates dt21*rhs without 0 @ extremes
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhslap
            
!--------------------------------------------------------------------------! 

      rhslap(1) = dt21(3,1)*f(1)+dt21(4,1)*f(2)+dt21(5,1)*f(3)
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
      rhslap(my)  = dt21(1,my)*f(my-2)+dt21(2,my)*f(my-1)+
     .              dt21(3,my)*f(my)

      end subroutine calcrhsfree




!------------------------------------------------------------------------------

      subroutine calcrhsd1(f,rhsd1)
      use matrices
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhsd1
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


!LAPVDV classic
!***********************************************************************!
!                                                                       !
!        LAPVDV                                                         !
!                                                                       !
!         solve the problem     v'' - rK v = phi                        !
!         form a matrix   ->         {v}*[A] = {f}                      ! 
!         General boundary conditions are allowed:                      !
!    C_phi(i)*phi(0)+k^2*C_vvv(i)*v(0)+k*C_dyv(i)*dvdy(0)=BC(i), i=1    !
!    C_phi(i)*phi(1)+k^2*C_vvv(i)*v(1)+k*C_dyv(i)*dvdy(1)=BC(i), i=3    !
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
      use worksp
      use fis
      use matrices

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


      zero = 0.

      if(rK.eq.0d0) then
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
      use fis
      use worksp
      use matrices

      implicit none
      include "ctes3D"

c ---------------Variables------------------------------------------- 
      integer i,j
      real(4),dimension(2,my) :: phi,v,dvdy
      real(8) :: rK
      real(4) :: zero

c ----------------------- Programa ---------------------------
      
      zero = 0.
      
      if(rK.eq.0d0) then
         do j=1,my
            v(1,j) =  zero
            v(2,j) =  zero
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
! aaf 2014-05
! Integrates (once) du in order to give u 
! all in real 8
!***********************************************************************!
      subroutine inty8(du,u)
      use worksp
      use fis
      use matrices


      implicit none
      include "ctes3D"

      integer i,j
      real(8),dimension(my) :: u, du
      real(8)                 :: zero, suma


! ----------------------- workspace ---------------------------

      real(8),dimension(my)   ::  f
      real(8),dimension(7,my)   :: wk2

!-----------------Program------------------------------------!

      zero = 0d0
! calculate v (solve d2vdy2=phi):
      !do j=1,my
      do j=1,my
       do i=1,7
          wk2(i,j)=dt12(i,j)
       enddo
      enddo
      !BC @ j=1 --> f=0.0
!      wk2(1,1) = zero
!      wk2(2,1) = zero
!      wk2(3,1) = zero
!      wk2(4,1) = 1.0d0
!      wk2(5,1) = zero
!      wk2(6,1) = zero
!      wk2(7,1) = zero
!     
      !create bandwitch matrix 
      call bandec7(wk2,my)

      !save the rhs into f
      do j=1,my
       f(j)=du(j)/fmap(j)
      enddo

      !prepare rhs multiplying by dt11
      vpr(1   ) =   dt11(4,1   )*f(1   ) + dt11(5,1   )*f(2) 
     $            + dt11(6,1   )*f(3   ) + dt11(7,1   )*f(4)
      !setting value of the integral
      !setting value of the integral
      !vpr(1   ) =   0d0 
!                     
      vpr(2   ) =   dt11(3,2   )*f(1   ) + dt11(4,2   )*f(2) 
     $            + dt11(5,2   )*f(3   ) + dt11(6,2   )*f(4) 
     $            + dt11(7,2   )*f(5   )

      vpr(3   ) =   dt11(2,3   )*f(1   ) + dt11(3,3   )*f(2)  
     $            + dt11(4,3   )*f(3   ) + dt11(5,3   )*f(4) 
     $            + dt11(6,3   )*f(5   ) + dt11(7,3   )*f(6)

      do j=4,my-3
         vpr(j)=dt11(1,j)*f(j-3)
         do i=2,7
            vpr(j) = vpr(j) + dt11(i,j)*f(i+j-4)
         enddo
      enddo
    
      vpr(my-2) =    dt11(1,my-2)*f(my-5) + dt11(2,my-2)*f(my-4)  
     $             + dt11(3,my-2)*f(my-3) + dt11(4,my-2)*f(my-2)
     $             + dt11(5,my-2)*f(my-1) + dt11(6,my-2)*f(my)    
  
      vpr(my-1) =    dt11(1,my-1)*f(my-4) + dt11(2,my-1)*f(my-3) 
     $             + dt11(3,my-1)*f(my-2) + dt11(4,my-1)*f(my-1) 
     $             + dt11(5,my-1)*f(my)
          
      vpr(my)   =    dt11(1,my  )*f(my-3) + dt11(2,my  )*f(my-2) 
     $             + dt11(3,my  )*f(my-1) + dt11(4,my  )*f(my)

!Inverting dt12 and premultiplying by dt11*f gives integration
      !call banbks7(wk2,my,vpr)     
      call banbks7(wk2,my,vpr)     
!Now we put zero on y=-Ly and make rest half of rhov00 @ Ly
       do j=1,my
         u(j)  = vpr(j) - vpr(1)
         !u(j)  = vpr(j) - (vpr(1)+vpr(my))/2.0d0
         !u(j)  = vpr(j) - vpr(1) - (vpr(my)-vpr(1))/2.0d0
      enddo 
!This way we make -DM/2 @ y=-Ly, and DM/2 @ y = Ly
      

      end subroutine inty8
!---------------------------------------------------------------! 


!***********************************************************************!
      subroutine solvepsi(phi,v,dvdy,rK)
      use worksp
      use fis
      use matrices

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

      Cvbot=C_vvv_bot1(rK) 
      Cdyvbot=C_dyv_bot1(rK)
      BC_bot=BC_phi_bot1(rK)


      zero = 0.

      if(rK.eq.0d0) then
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

      end subroutine solvepsi
 



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
       
 


