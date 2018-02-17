c!******************************************************************!
c!                                                                  !
c!           resuelve el problema  v'' - rK v = phi                 !
c!         forma matricial ->         {v}*[A] = {f}                 ! 
c!    Atchung!!! Sienmpre suponemos c.c. homogeneas                 !
c!                                                                  !
c!        DIFERENCIAS FINITAS COMPACTAS                             !
c!  input:                                                          !
c!     phi: vector phi.                                             !
c!       n: numero de terminos                                      !
c!      rK: constante independiente.                                !
c!      wk: area de trabajo de 9*n minimo (en real*8 !!!)           ! 
c! output:                                                          !
c!       v: solucion en el plano fourier x fourier                  !
c!    dvdy: dvndy en el mismo plano anterior                        !
c!                                                                  !
c!******************************************************************!
      subroutine Lapvdv(phi,v,dvdy,rK)

      implicit none
      include "ctes3D"
      
      integer i,j
      real*4 phi(2,my),v(2,my),dvdy(2,my)
      real*8 rK
      real*4 zero
     
c ------------------ commons ------------------------------------
 
      
      real*8 wk1,fwk1,fwk2,fwk,dphi1,dphi2
      common /worksp/wk1(my,5),fwk1(my),fwk2(my),fwk(my),dphi1(my),
     .               dphi2(my)
      save /worksp/


      real*8  Re,alp,bet,a0
      real*8  y,fmap
      common /fis/ Re,alp,bet,a0,y(my),fmap(my)
      save   /fis/
      
      real*8 d21,d22,dt21
      common/secder/ d21(my,5),d22(my,5),dt21(5,my)
      save  /secder/

      real*8 prem1,dt12
      common /cfdiff/ prem1(my,7),dt12(7,my)
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

      do i=1,5
         do j=2,my-1
            wk1(j,i)=d22(j,i)-rk*d21(j,i)
         enddo
      enddo

      wk1(1,1)=0d0
      wk1(1,2)=0d0
      wk1(1,3)=1d0
      wk1(1,4)=0d0
      wk1(1,5)=0d0
      wk1(my,1)=0d0
      wk1(my,2)=0d0
      wk1(my,3)=1d0
      wk1(my,4)=0d0
      wk1(my,5)=0d0

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


c!******************************************************************!
c!                                                                  !
c!           resuelve el problema  v'' - rK v = phi                 !
c!                                                                  !
c!                                 con v(-1)=0 v(+1)=0              !
c!	     DIFERENCIAS FINITAS COMPACTAS	                    !
c!  input:                                                          !
c!     phi: vector phi.                                             !
c!       n: numero de terminos 				            !
c!      rK: constante independiente.                                !
c!      wk: area de trabajo de 9*n minimo  (en real*8 !!!)	    ! 
c! output:                                                          !
c!       v: solucion en el plano fourier x fourier 	            !
c!                                                                  !
c!******************************************************************!
      subroutine Lapv(phi,v,rk,dalbe,dtri)

      implicit none
      include "ctes3D"
      
c     ----------------------- IN & OUT -----------------

      real*4  phi(2,my),v(2,my)
      real*8  rk,dalbe,dtri

c     ----------------------- Commons  -----------------

      real*8  d21,d22,dt21
      common  /secder/ d21(my,5),d22(my,5),dt21(5,my)
      save    /secder/
  
      real*8 wk1,fwk1,fwk2,dphi1,dphi2
      common /worksp/ wk1(my,5),fwk1(my),fwk2(my),dphi1(my),dphi2(my)
      save /worksp/
c     ----------------------- work    ------------------

      integer n,i,j
      real*8  zero

c     ----------------------- Program ------------------

      zero = 0d0
      do i=1,5
         do j=2,my-1
            wk1(j,i)=d22(j,i)-rk*d21(j,i)
         enddo
      enddo


      wk1(1,1)=0d0
      wk1(1,2)=0d0
      wk1(1,3)=1d0
      wk1(1,4)=0d0
      wk1(1,5)=0d0
      wk1(my,1)=0d0
      wk1(my,2)=0d0
      wk1(my,3)=1d0
      wk1(my,4)=0d0
      wk1(my,5)=0d0

      call bandec(wk1,my)
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
         phi(1,j) = fwk1(j)
         phi(2,j) = fwk2(j)
      enddo

      if (dalbe.ne.0d0) then
         do j=1,my
            v(1,j) = dalbe*fwk1(j)+dtri*dphi1(j)
            v(2,j) = dalbe*fwk2(j)+dtri*dphi2(j)
         enddo
      else
         do j=1,my
            v(1,j) = fwk1(j)
            v(2,j) = fwk2(j)
         enddo
      endif
      
      end


c!******************************************************************!
c!                                                                  !
c!           resuelve el problema                                   !
c!                                                                  !
c!           phi'' - rK1 phi = f                                    !
c!           v''   - rK2 v   = phi                                  !
c!                                                                  !
c!	 cc:  v(-1)=v(1)=dvdy(1)=dvdy(-1)=0                          !
c!                                                                  !
c!                                                                  !
c!                                                                  !
c!	     DIFERENCIAS FINITAS COMPACTAS			     !
c!      output: phi,v,dvdy vectores de 2*my elementos (re,im)       !
c!	 	 phi,v,dvdy puedem sobreescribir el forzante f       !
c!                                                                  !
c!      workspace:                                                  !
c!              vec1,vec2,vec3    vector de 2*my \                  !
c!              vec4,vec5,vec6    vector de my   |> real*4          !
c!              vec7,vec8,vec9    vector de my   /                  !
c!		 work  vector de 9*my (en real*8)		     !
c!                                                                  !
c!                                                                  !
c!******************************************************************!
      subroutine visc(var,bcb,bct,rk)
      
      implicit none
      include "ctes3D"

c     ----------------------- IN & OUT -----------------

      real*4    var(2,my)
      complex*8 bcb,bct
      real*8    rk

c     ----------------------- Commons  -------------------

      real*8  d21,d22,dt21
      common  /secder/ d21(my,5),d22(my,5),dt21(5,my)
      save    /secder/

      real*8 prem1,dt12
      common /cfdiff/ prem1(my,7),dt12(7,my)
      save   /cfdiff/

      real*4  Re,alp,bet,a0
      real*8  y,fmap
      common /fis/ Re,alp,bet,a0,y(my),fmap(my)
      save   /fis/

      real*8 wk1,vpr,vpi,dvpr,dvpi,v1,dv1,dfr,dfi
      common /worksp/ wk1(my,5),vpr(my),vpi(my),
     .    dvpr(my),dvpi(my),v1(my),dv1(my),dfr(my),dfi(my)

c     ----------------------- work    --------------------

      integer i,j,jj
      real*8 zero,det,Ar,Ai,Br,Bi,rbct,rbcb,ibct,ibcb

c     ----------------------- Program ------------------

     
      zero = 0d0
      rbcb = real(bcb)   ! 0d0 !
      rbct = real(bct)   ! 0d0 !
      ibcb = aimag(bcb)  ! 0d0 !
      ibct = aimag(bct)  ! 0d0 !

c     RESUELVO LA BILAPLACIANA EN 2 ETAPAS:
c     PRIMERO RESUELVO LAS phi DE LOS TRES PROBLEMAS (particular+2homog)
c     preparo la matriz wk1, que es la misma para todos los poisson phi:

      do i=1,5
         do j=2,my-1
            wk1(j,i)=d22(j,i)-rk*d21(j,i)
         enddo
      enddo
      wk1(1,1)=0d0
      wk1(1,2)=0d0
      wk1(1,3)=1d0
      wk1(1,4)=0d0
      wk1(1,5)=0d0
      wk1(my,1)=0d0
      wk1(my,2)=0d0
      wk1(my,3)=1d0
      wk1(my,4)=0d0
      wk1(my,5)=0d0
      call bandec(wk1,my)

c      phip'' - rk1*phip = f     phip(1) = phip(-1) = 0 
           
      do j=1,my
         dfr(j)=var(1,j)
         dfi(j)=var(2,j)
      enddo

      vpr(1) = 0d0
      vpi(1) = 0d0
      
      vpr(2) = dt21(2,2)*dfr(1)+dt21(3,2)*dfr(2)+
     .         dt21(4,2)*dfr(3)+dt21(5,2)*dfr(4)
      vpi(2) = dt21(2,2)*dfi(1)+dt21(3,2)*dfi(2)+
     .         dt21(4,2)*dfi(3)+dt21(5,2)*dfi(4)  
     
      do j=3,my-2
         vpr(j)=dt21(1,j)*dfr(j-2)
         vpi(j)=dt21(1,j)*dfi(j-2)    
         do i=2,5
            vpr(j)=vpr(j)+dt21(i,j)*dfr(j-3+i)
            vpi(j)=vpi(j)+dt21(i,j)*dfi(j-3+i)       
         enddo
      enddo
      vpr(my-1) = dt21(1,my-1)*dfr(my-3)+dt21(2,my-1)*dfr(my-2)+
     .            dt21(3,my-1)*dfr(my-1)+dt21(4,my-1)*dfr(my)
      vpi(my-1) = dt21(1,my-1)*dfi(my-3)+dt21(2,my-1)*dfi(my-2)+
     .            dt21(3,my-1)*dfi(my-1)+dt21(4,my-1)*dfi(my)      

      vpr(my)=0d0
      vpi(my)=0d0 
      
      call banbks(wk1,my,vpr)      
      call banbks(wk1,my,vpi)
      
c     phi1'' - rk1*phi1 = 0     phi1(1) =0  phi1(-1)=1 

      v1(1)=1d0
      do j=2,my
         v1(j)=0d0
      enddo

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


      det =  1d0/(-dv1(1)*dv1 (1) + dv1(my )*dv1 (my))

      Ar =(-dv1(my)*(dvpr(my)-rbct/fmap(my))
     $    + dv1(1) *(dvpr(1) -rbcb/fmap(1)))*det
      Br =(-dv1(1) *(dvpr(my)-rbct/fmap(my))
     $    + dv1(my)*(dvpr(1) -rbcb/fmap(1)))*det

      Ai =(-dv1(my)*(dvpi(my)-ibct/fmap(my))
     $    + dv1(1) *(dvpi(1) -ibcb/fmap(1)))*det
      Bi =(-dv1(1) *(dvpi(my)-ibct/fmap(my))
     $    + dv1(my)*(dvpi(1) -ibcb/fmap(1)))*det
     

      jj=my
      
      do j=1,my
         var(1,j)  =    vpr(j) + Ar*v1(j) + Br*v1(jj)
         var(2,j)  =    vpi(j) + Ai*v1(j) + Bi*v1(jj)
         jj=jj-1
      enddo

      end
