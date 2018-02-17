!=======================================================================!
!                                                                       !
!                 	Viscous step solvers Package                    !
!                           O.F. 2003                                   !
!                           S.H. may 2005                               !
!                                                                       !
!=======================================================================!

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
      subroutine Lapvdv(phi,v,dvdy,rK)

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
!      Solves the problems  v'' - rK1 v = phi                           !
!                with b.c.  v(-1)=0 v(+1)=0                             !
!      AND                                                              ! 
!                           phi'' - rK1 phi = f                         !
!                           v''   - rK2 v   = phi                       !
!                                                                       !
!                    bc:  v(-1)=v(1)=dvdy(1)=dvdy(-1)= 0                !
!                                                                       !
!                                                                       !
!      DIFERENCIAS FINITAS COMPACTAS                                    !
!                                                                       !
!  input:                                                               !
!           phi: field phi.                                             !
!            rK: constant                                               !
!    dalbe,dtri: constant                                               !
!                                                                       !
! output:                                                               !
!       phi: Solution of system                                         !
!                                                                       !
!       if dalbe=dtri=0.0                                               !
!       v: solution of the system                                       !
!       if dalbe ~=0                                                    !
!       v: solution + some things prepared for SMR                      !
!                                                                       !
!                                                                       !
!                                                                       !
!***********************************************************************!
      subroutine visc(vor,vvor,phi,v,dvdy,f,rk1,rk2,dalbe,dtri,dtxi)

      implicit none
      include "ctes3D"
      
c     ----------------------- IN & OUT ---------------------------------
      ! LAP 
      
      real*4  vor(2,my),vvor(2,my)
      
      ! BILAP
      
      real*4  phi(2,my),v(2,my),dvdy(2,my),f(2,my)        
      real*8  rk1,rk2,dalbe,dtri,dtxi

c     ----------------------- Commons  ---------------------------------

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/
      
      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save   /cfdiff/ 
       
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/  
           
      real*8 wk1,phipr,phipi,vpr,vpi,dvpr,dvpi,phi1,v1,dv1,dfr,dfi,
     .       dphi1,dphi2,fwk1,fwk2
      common /worksp/ wk1(5,my),fwk1(my),fwk2(my),dphi1(my),dphi2(my),
     .                phipr(my),phipi(my),vpr(my),vpi(my),
     .    dvpr(my),dvpi(my),phi1(my),v1(my),dv1(my),dfr(my),dfi(my)
            
      
c     ----------------------- work    ----------------------------------

      integer n,i,j
      real*8 zero,det,Ar,Ai,Br,Bi

c     ----------------------- Program ----------------------------------

      zero = 0d0

      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk1*dt21(i,j)
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
      
!     Laplaciano       

      call bandec(wk1,my)
      
      do j=1,my
         dphi1(j)=vor(1,j)
         dphi2(j)=vor(2,j)
      enddo

      fwk1(1) = zero
      fwk2(1) = zero
      fwk1(2) = dt21(2,2)*dphi1(1)+dt21(3,2)*dphi1(2)+
     .          dt21(4,2)*dphi1(3)+dt21(5,2)*dphi1(4)
      fwk2(2) = dt21(2,2)*dphi2(1)+dt21(3,2)*dphi2(2)+
     .          dt21(4,2)*dphi2(3)+dt21(5,2)*dphi2(4)   
      do j=3,my-2
         fwk1(j) = dt21(1,j)*dphi1(j-2)    
         fwk2(j) = dt21(1,j)*dphi2(j-2)  
         do i=2,5
            fwk1(j) = fwk1(j)+dt21(i,j)*dphi1(j-3+i)
            fwk2(j) = fwk2(j)+dt21(i,j)*dphi2(j-3+i)    
         enddo
      enddo
      fwk1(my-1)=dt21(1,my-1)*dphi1(my-3)+dt21(2,my-1)*dphi1(my-2)+
     .           dt21(3,my-1)*dphi1(my-1)+dt21(4,my-1)*dphi1(my  )
      fwk2(my-1)=dt21(1,my-1)*dphi2(my-3)+dt21(2,my-1)*dphi2(my-2)+
     .           dt21(3,my-1)*dphi2(my-1)+dt21(4,my-1)*dphi2(my  )
          
      fwk1(my)  = zero
      fwk2(my)  = zero
      
      call banbks(wk1,my,fwk1) 
      call banbks(wk1,my,fwk2)
c		 	guardo la v:
      do j=1,my
         vor(1,j) = fwk1(j)
         vor(2,j) = fwk2(j)
      enddo

      if (dalbe.ne.0d0) then
         do j=1,my
            vvor(1,j) = dalbe*fwk1(j)+dtri*dphi1(j)+dtxi*v(1,j)
            vvor(2,j) = dalbe*fwk2(j)+dtri*dphi2(j)+dtxi*v(2,j)
         enddo
      else
         do j=1,my
            vvor(1,j) = fwk1(j)
            vvor(2,j) = fwk2(j)
         enddo
      endif
      
      
!     BILAP                   !      
      

c      phip'' - rk1*phip = f     phip(1) = phip(-1) = 0 
           
      do j=1,my
         dfr(j)=f(1,j)
         dfi(j)=f(2,j)
      enddo

      phipr(1) = 0d0
      phipi(1) = 0d0
      
      phipr(2) = dt21(2,2)*dfr(1)+dt21(3,2)*dfr(2)+
     .           dt21(4,2)*dfr(3)+dt21(5,2)*dfr(4)
      phipi(2) = dt21(2,2)*dfi(1)+dt21(3,2)*dfi(2)+
     .           dt21(4,2)*dfi(3)+dt21(5,2)*dfi(4)  
     
      do j=3,my-2
         phipr(j)=dt21(1,j)*dfr(j-2)
         phipi(j)=dt21(1,j)*dfi(j-2)    
         do i=2,5
            phipr(j)=phipr(j)+dt21(i,j)*dfr(j-3+i)
            phipi(j)=phipi(j)+dt21(i,j)*dfi(j-3+i)       
         enddo
      enddo
      phipr(my-1) = dt21(1,my-1)*dfr(my-3)+dt21(2,my-1)*dfr(my-2)+
     .              dt21(3,my-1)*dfr(my-1)+dt21(4,my-1)*dfr(my)
      phipi(my-1) = dt21(1,my-1)*dfi(my-3)+dt21(2,my-1)*dfi(my-2)+
     .              dt21(3,my-1)*dfi(my-1)+dt21(4,my-1)*dfi(my)      

      phipr(my)=0d0
      phipi(my)=0d0 
      
      call banbks(wk1,my,phipr)      
      call banbks(wk1,my,phipi)
      
c     phi1'' - rk1*phi1 = 0     phi1(1) =0  phi1(-1)=1 

      phi1(1)=1d0
      do j=2,my
         phi1(j)=0d0
      enddo

      call banbks(wk1,my,phi1)

c     AHORA RESUELVO LAS VELOCIDADES
c     de nuevo la matriz wk1 es la misma para las 3 velocidades:

      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk2*dt21(i,j)
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

c      condiciones de contorno:

      vpr(1) = 0d0
      vpi(1) = 0d0
      v1(1)  = 0d0
      
      vpr(2) = dt21(2,2)*phipr(1)+dt21(3,2)*phipr(2)+
     .         dt21(4,2)*phipr(3)+dt21(5,2)*phipr(4)
      vpi(2) = dt21(2,2)*phipi(1)+dt21(3,2)*phipi(2)+
     .         dt21(4,2)*phipi(3)+dt21(5,2)*phipi(4)  
      v1(2)  = dt21(2,2)*phi1 (1)+dt21(3,2)*phi1 (2)+
     .         dt21(4,2)*phi1 (3)+dt21(5,2)*phi1 (4)      
     
      do j=3,my-2
         vpr(j) = dt21(1,j)*phipr(j-2)
         vpi(j) = dt21(1,j)*phipi(j-2)
         v1 (j) = dt21(1,j)*phi1 (j-2)
         do i=2,5
            vpr(j) = vpr(j)+dt21(i,j)*phipr(j-3+i)
      	    vpi(j) = vpi(j)+dt21(i,j)*phipi(j-3+i)
            v1 (j) = v1(j) +dt21(i,j)*phi1 (j-3+i)
         enddo
      enddo
      vpr(my-1)= dt21(1,my-1)*phipr(my-3)+dt21(2,my-1)*phipr(my-2)+
     .           dt21(3,my-1)*phipr(my-1)+dt21(4,my-1)*phipr(my)
      vpi(my-1)= dt21(1,my-1)*phipi(my-3)+dt21(2,my-1)*phipi(my-2)+
     .           dt21(3,my-1)*phipi(my-1)+dt21(4,my-1)*phipi(my)
      v1(my-1) = dt21(1,my-1)*phi1 (my-3)+dt21(2,my-1)*phi1 (my-2)+
     .           dt21(3,my-1)*phi1 (my-1)+dt21(4,my-1)*phi1 (my)
     
      vpr(my) = 0d0      
      vpi(my) = 0d0
      v1(my)  = 0d0 
      
      call banbks(wk1,my,vpr)
      call banbks(wk1,my,vpi)
      call banbks(wk1,my,v1)
     
c     REAL

      dvpr(1)= dt12(4,1)*vpr(1)+dt12(5,1)*vpr(2)+
     .         dt12(6,1)*vpr(3)+dt12(7,1)*vpr(4)
      dvpi(1)= dt12(4,1)*vpi(1)+dt12(5,1)*vpi(2)+
     .         dt12(6,1)*vpi(3)+dt12(7,1)*vpi(4)
      dv1(1)=  dt12(4,1)* v1(1)+dt12(5,1)* v1(2)+
     .         dt12(6,1)* v1(3)+dt12(7,1)* v1(4)
                   
      dvpr(2)= dt12(3,2)*vpr(1)+dt12(4,2)*vpr(2)+
     .         dt12(5,2)*vpr(3)+dt12(6,2)*vpr(4)+
     .         dt12(7,2)*vpr(5)
      dvpi(2)= dt12(3,2)*vpi(1)+dt12(4,2)*vpi(2)+
     .         dt12(5,2)*vpi(3)+dt12(6,2)*vpi(4)+
     .         dt12(7,2)*vpi(5) 
      dv1(2)=  dt12(3,2)*v1 (1)+dt12(4,2)*v1 (2)+
     .         dt12(5,2)*v1 (3)+dt12(6,2)*v1 (4)+
     .         dt12(7,2)*v1 (5)        
         
      dvpr(3)= dt12(2,3)*vpr(1)+dt12(3,3)*vpr(2)+
     .         dt12(4,3)*vpr(3)+dt12(5,3)*vpr(4)+
     .         dt12(6,3)*vpr(5)+dt12(7,3)*vpr(6)     
      dvpi(3)= dt12(2,3)*vpi(1)+dt12(3,3)*vpi(2)+
     .         dt12(4,3)*vpi(3)+dt12(5,3)*vpi(4)+
     .         dt12(6,3)*vpi(5)+dt12(7,3)*vpi(6)
      dv1(3) = dt12(2,3)*v1 (1)+dt12(3,3)*v1 (2)+
     .         dt12(4,3)*v1 (3)+dt12(5,3)*v1 (4)+
     .         dt12(6,3)*v1 (5)+dt12(7,3)*v1 (6)        

      do j=4,my-3
         dvpr(j)=dt12(1,j)*vpr(j-3)
         dvpi(j)=dt12(1,j)*vpi(j-3)
         dv1(j) =dt12(1,j) *v1(j-3)        
         do i=2,7
            dvpr(j) = dvpr(j) + dt12(i,j)*vpr(i+j-4)
            dvpi(j) = dvpi(j) + dt12(i,j)*vpi(i+j-4)
            dv1(j)  = dv1 (j) + dt12(i,j)*v1 (i+j-4)
         enddo
      enddo
      
      dvpr(my-2) = dt12(1,my-2)*vpr(my-5)+dt12(2,my-2)*vpr(my-4)+
     .             dt12(3,my-2)*vpr(my-3)+dt12(4,my-2)*vpr(my-2)+
     .             dt12(5,my-2)*vpr(my-1)+dt12(6,my-2)*vpr(my)    
      dvpi(my-2) = dt12(1,my-2)*vpi(my-5)+dt12(2,my-2)*vpi(my-4)+
     .             dt12(3,my-2)*vpi(my-3)+dt12(4,my-2)*vpi(my-2)+
     .             dt12(5,my-2)*vpi(my-1)+dt12(6,my-2)*vpi(my)
      dv1(my-2)  = dt12(1,my-2)*v1 (my-5)+dt12(2,my-2)*v1 (my-4)+
     .             dt12(3,my-2)*v1 (my-3)+dt12(4,my-2)*v1 (my-2)+
     .             dt12(5,my-2)*v1 (my-1)+dt12(6,my-2)*v1 (my)
     
      dvpr(my-1) = dt12(1,my-1)*vpr(my-4)+dt12(2,my-1)*vpr(my-3)+
     .             dt12(3,my-1)*vpr(my-2)+dt12(4,my-1)*vpr(my-1)+
     .             dt12(5,my-1)*vpr(my)
      dvpi(my-1) = dt12(1,my-1)*vpi(my-4)+dt12(2,my-1)*vpi(my-3)+
     .             dt12(3,my-1)*vpi(my-2)+dt12(4,my-1)*vpi(my-1)+
     .             dt12(5,my-1)*vpi(my) 
      dv1(my-1)  = dt12(1,my-1)*v1 (my-4)+dt12(2,my-1)*v1 (my-3)+
     .             dt12(3,my-1)*v1 (my-2)+dt12(4,my-1)*v1 (my-1)+
     .             dt12(5,my-1)*v1 (my)  
            
      dvpr(my)   = dt12(1,my  )*vpr(my-3)+dt12(2,my  )*vpr(my-2)+
     .             dt12(3,my  )*vpr(my-1)+dt12(4,my  )*vpr(my)
      dvpi(my)   = dt12(1,my  )*vpi(my-3)+dt12(2,my  )*vpi(my-2)+
     .             dt12(3,my  )*vpi(my-1)+dt12(4,my  )*vpi(my)
      dv1(my)    = dt12(1,my  )*v1 (my-3)+dt12(2,my  )*v1 (my-2)+
     .             dt12(3,my  )*v1 (my-1)+dt12(4,my  )*v1 (my)        


      call banbks7(prem1,my,dvpr)
      call banbks7(prem1,my,dvpi)
      call banbks7(prem1,my,dv1)



c     CALCULO LOS COEFICIENTES A y B IMPONIENDO dvdy(1)=bctdv, dvdy(-1)=bcbdv

      det =  -dv1(1)*dv1 (1) + dv1(my )*dv1 (my)

      Ar  = -(dv1(my)*dvpr(my) - dv1(1 )*dvpr(1))/det
      Br  = (-dv1(1 )*dvpr(my) + dv1(my)*dvpr(1))/det
      Ai  = -(dv1(my)*dvpi(my) - dv1(1 )*dvpi(1))/det
      Bi  = (-dv1(1 )*dvpi(my) + dv1(my)*dvpi(1))/det
      
      do j=1,my
         phipr(j) = phipr(j)+Ar*phi1(j)+ Br*phi1(my+1-j)    
         phipi(j) = phipi(j)+Ai*phi1(j)+ Bi*phi1(my+1-j)    
         f(1,j)   = phipr(j)                                
         f(2,j)   = phipi(j)                                
      enddo

      if (dalbe.ne.0d0) then
         do j=1,my
            phi(1,j) = dalbe*phipr(j)+dtri*dfr(j)+dtxi*dvdy(1,j)
            phi(2,j) = dalbe*phipi(j)+dtri*dfi(j)+dtxi*dvdy(2,j)
         enddo
      else
         do j=1,my
            phi(1,j) = f(1,j)
            phi(2,j) = f(2,j)
         enddo
      endif
      
      if (rk2.ne.0d0) then
         do j=1,my
            v(1,j) =    vpr(j)  +Ar*v1(j)  + Br*v1 (my+1-j)        
            v(2,j) =    vpi(j)  +Ai*v1(j)  + Bi*v1 (my+1-j)         
            dvdy(1,j) = (dvpr(j)+Ar*dv1(j) - Br*dv1(my+1-j))*fmap(j)
            dvdy(2,j) = (dvpi(j)+Ai*dv1(j) - Bi*dv1(my+1-j))*fmap(j)
         enddo
      else
         do j=1,my
            dvdy(1,j) = zero
            dvdy(2,j) = zero
            v(1,j)    = zero  
            v(2,j)    = zero          
         enddo
      endif

      end


!***********************************************************************!
!                                                                       !
!           solves   u'' - rK u = f                                     !
!                                                                       !
!                    con u(-1) = u(+1)= 0                               !
!                                                                       !
!      DIFERENCIAS FINITAS COMPACTAS                                    !
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

c     ----------------------- IN & OUT ----------------- 
    
      real*4  f(my),u(my)
      real*8  rk,dalbe,dtri
      integer ind
      
c     ----------------------- Commons -------------------
      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/


      real*8 wk1,fwk,df
      common /worksp/ wk1(5,my),fwk(my),df(my)
      save  /worksp/
c     ----------------------- work --------------------     
            
      integer i,j

c     ----------------------- Program ------------------      


c     calculo v (resuelvo el poisson):

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

c      resuelvo  v (coeficientes reales):

      do j=1,my
         df(j)=f(j)
      enddo

      fwk(1) = 0d0
      fwk(2) = dt21(2,2)*df(1)+dt21(3,2)*df(2)+
     .         dt21(4,2)*df(3)+dt21(5,2)*df(4)
      do j=3,my-2
         fwk(j)    = dt21(1,j)*df(j-2)
         do i=2,5
            fwk(j) = fwk(j)+dt21(i,j)*df(j-3+i)
         enddo
      enddo
      fwk(my-1)= dt21(1,my-1)*df(my-3)+dt21(2,my-1)*df(my-2)+
     .           dt21(3,my-1)*df(my-1)+dt21(4,my-1)*df(my)
     
      fwk(my)  = 0d0

      call banbks(wk1,my,fwk)
      

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
          
      end
