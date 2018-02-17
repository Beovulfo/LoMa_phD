!-----------------------------------------------------------------------!
!	PAQUETE DE SUBRUTINAS PARA DIFERENCIAS FINITAS COMPACTAS	            !
!	O.F, nov 2002				                                             !
!		               		                                             !
!	updated dic 2002  		                                             !
!	touched ene 2003  ---------> this subroutines don't pass worksp      !
!   d11,d12,d21,d22 & solver real*8  		                              !
!-----------------------------------------------------------------------!


!-----------------------------------------------------------------------
! pred1m60:                     
!                   
! Genera los coeficientes alpha,beta,delta y
!      a0,a1,a2 y a3 para el esquema de  
! diferencias finitas compactas con malla uniforme:
!                      
! delta*(du(j+2)+du(j-2))beta*(du(j+2)+du(j-2))+
!      alpha*(du(j+1)+du(j-1))+du(j) =
!  = a1*(u(j+1)-u(j-1)) + a2(u(j+2)-u(j-2)) + a3(u(j+2)-u(j-3))
!          
!    En la frontera del dominio:
!    alpha*du(2) + du(1) =
!     = a1*u(1) + a2*u(2) + a3*u(3) + a4*u(4) + a5*u(5)
!                      
!    alpha*(du(3)+du(1)) + du(2) = 
!     = a1*(u(3)-u(1)) 
!                
!-----------------------------------------------------------------------
      subroutine pred1m6o(my,den,ben,aln,one,alp,bep,dep,
     ,                       a3n,a2n,a1n, a0,a1p,a2p,a3p)
      implicit none

      integer my
      real*8 alpha,beta,delta,aaa,bbb,ccc
      real*8 den(my),ben(my),aln(my),one(my),alp(my),bep(my),dep(my)
      real*8 a3n(my),a2n(my),a1n(my), a0(my),a1p(my),a2p(my),a3p(my)

      real*8 h,pi
      integer j

      pi = 4d0*atan(1d0)

      do j=1,my
         one(j)=1d0
      enddo

      h=2d0/(my-1)


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!               centro del dominio
!               IMPONGO LOS COEFICIENTES A PARTIR DE CALCULOS EXTERNOS


      aaa =  0.56926818603810d0
      bbb =  0.31997939759919d0
      ccc =  0.03751544855390d0
      alpha = -3d0/4d0  +(33d0/20d0  )*ccc
     .                  +(239d0/180d0)*aaa + (82d0/45d0 )*bbb
      beta  =  3d0/10d0 +(24d0/25d0  )*ccc
     .                  -(88d0/225d0 )*aaa + (34d0/225d0)*bbb
      delta = -1d0/2d1  +(39d0/100d0 )*ccc
     .                  +(19d0/300d0 )*aaa + (2d0/75d0  )*bbb

      do j=4,my-3
         den(j)=delta
         ben(j)=beta
         aln(j)=alpha
         alp(j)=alpha
         bep(j)=beta
         dep(j)=delta
         a3n(j)=-ccc/h
         a2n(j)=-bbb/h
         a1n(j)=-aaa/h
         a0(j)=0d0
         a1p(j)= aaa/h
         a2p(j)= bbb/h
         a3p(j)= ccc/h
      enddo


c*************************************************************************
c*************************************************************************
c       fronteras del dominio
c       ESTANDAR:SOLO IMPONGO ECUACIONES DE ORDEN
c       USANDO TODA LA INFORMACION DISPONIBLE:

      j=1
      den(j) = 0d0              !   0                 !a30n
      ben(j) = 0d0              !   0                 !a20n
      aln(j) = 0d0              !   0                 !a10n
      alp(j) = 9d0              !   0                 !a10p
      bep(j) = alp(j)           !   0                 !a20p
      dep(j) = 1.d0             !   0                 !a30p

      a3n(j) = 0d0              !   0                 !b30n
      a2n(j) = 0d0              !   0                 !b20n
      a1n(j) = 0d0              !   0                 !b10n
      a0 (j) = -11d0/3d0/h      !   -11d0/6d0/h       !b00
      a1p(j) = -9d0/h           !   18d0/6d0/h        !b10p
      a2p(j) = -a1p(j)          !   -9d0/6d0/h        !b20p
      a3p(j) = -a0(j)           !   2d0/6d0/h         !b30p

      j=2
      den(j) = 0d0              ! 0
      ben(j) = 0d0              ! 0
      aln(j) = 1d0/16d0         !0.07968949d0           ! 0
      alp(j) = 9d0/4d0          !2.043726096d0            ! 0
      bep(j) = 1d0              ! 0
      dep(j) = aln(j)           ! 0

      a3n(j) = 0d0              ! 0
      a2n(j) = 0d0              ! 0
      a1n(j) = -25d0/96d0/h     !-0.3091292240d0/h          !-6d0/24d0/h
      a0 (j) = -5d0/3d0/h       !-1.483312091d0/h           !-20d0/24d0/h
      a1p(j) = 0d0              !  36d0/24d0/h
      a2p(j) = -a0(j)           ! -12d0/24d0/h
      a3p(j) = -a1n(j)          !   2d0/24d0/h


      j=3
      den(j) =  0d0                     !   0
      ben(j) =  2.538032824d-2            !   0
      aln(j) =  0.347040189d0             !   0
      alp(j) =  1d0                 !   0
      bep(j) =  aln(j)                  !   0
      dep(j) =  ben(j)                  !   0

      a3n(j) =  0d0                     !   0
      a2n(j) =  -9.841716629d-2/h        !    6d0/120d0/h
      a1n(j) =  -0.5830132425d0/h          !  -60d0/120d0/h
      a0 (j) =  -0.5037154692d0/h         !  -40d0/120d0/h
      a1p(j) =  -a0(j)               !  120d0/120d0/h
      a2p(j) =  -a1n(j)                   !  -30d0/120d0/h
      a3p(j) =  -a2n(j)                     !    4d0/120d0/h

* -----------------  the other end --------------

      j=my
      den(j) =   dep(1)
      ben(j) =   bep(1)
      aln(j) =   alp(1)
      alp(j) =   aln(1)
      bep(j) =   ben(1)
      dep(j) =   den(1)
      a3n(j) =   -a3p(1)
      a2n(j) =   -a2p(1)
      a1n(j) =   -a1p(1)
      a0 (j) =   -a0(1)
      a1p(j) =   -a1n(1)
      a2p(j) =   -a2n(1)
      a3p(j) =   -a3n(1)

      j=my-1
      den(j) =   dep(2)
      ben(j) =   bep(2)
      aln(j) =   alp(2)
      alp(j) =   aln(2)
      bep(j) =   ben(2)
      dep(j) =   den(2)
      a3n(j) =   -a3p(2)
      a2n(j) =   -a2p(2)
      a1n(j) =   -a1p(2)
      a0 (j) =   -a0(2)
      a1p(j) =   -a1n(2)
      a2p(j) =   -a2n(2)
      a3p(j) =   -a3n(2)

      j=my-2
      den(j) =    dep(3)
      ben(j) =    bep(3)
      aln(j) =    alp(3)
      alp(j) =    aln(3)
      bep(j) =    ben(3)
      dep(j) =    den(3)
      a3n(j) =    -a3p(3)
      a2n(j) =    -a2p(3)
      a1n(j) =    -a1p(3)
      a0 (j) =    -a0(3)
      a1p(j) =    -a1n(3)
      a2p(j) =    -a2n(3)
      a3p(j) =    -a3n(3)

      return
      end


c***************************************************************************
c*	prederiv2:							
c*								
c*	Genera los coeficientes alpha,beta,a0,a1 y a2 para el esquema de 
c*	diferencias finitas compactas:					
c*								
c*	beta*(d2u(j+2)+d2u(j-2))+alpha*(d2u(j+1)+d2u(j-1))+d2u(j) =		
c*	 = a0*u(j)+a1*(u(j+1)+u(j-1)) + a2(u(j+2)+u(j-2))	
c*		
c*						
c*				
c***************************************************************************
      subroutine prederiv2(my,ben,aln,aaa,alp,bep,a2n,a1n,a0,a1p,a2p,x)

      implicit none

      integer my,n
      real*8 x(my)
      real*8 ben(my),aln(my),aaa(my),alp(my),bep(my)
      real*8 a2n(my),a1n(my),a0(my),a1p(my),a2p(my)
      real*8 c(9,9),b(9)
   
      integer j

      do j=1,my
      aaa(j)=1d0
      enddo

c		Ahora calculo los coeficientes del centro del dominio
  
      do j=3,my-2
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=0d0
      c(1,4)=0d0
      c(1,5)=1d0
      c(1,6)=1d0
      c(1,7)=1d0
      c(1,8)=1d0
      c(1,9)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=0d0
      c(2,5)=(x(j-2)-x(j))
      c(2,6)=(x(j-1)-x(j))
      c(2,7)=0d0
      c(2,8)=(x(j+1)-x(j))
      c(2,9)=(x(j+2)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=2d0
      c(3,4)=2d0
      c(3,5)=-(x(j-2)-x(j))**2
      c(3,6)=-(x(j-1)-x(j))**2
      c(3,7)=0d0
      c(3,8)=-(x(j+1)-x(j))**2
      c(3,9)=-(x(j+2)-x(j))**2
      b(3)=-2d0

      do n=3,8
      c(n+1,1)=n*(n-1)*(x(j-2)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,3)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,4)=n*(n-1)*(x(j+2)-x(j))**(n-2)
      c(n+1,5)=-(x(j-2)-x(j))**n
      c(n+1,6)=-(x(j-1)-x(j))**n
      c(n+1,7)=0d0
      c(n+1,8)=-(x(j+1)-x(j))**n
      c(n+1,9)=-(x(j+2)-x(j))**n
      b(n+1)=0d0
      enddo
      
c		call LSLRG(9,c,9,b,1,sol)
      call gaussj(c,9,9,b,1,1)
      

      ben(j)=b(1)
      aln(j)=b(2)
      alp(j)=b(3)
      bep(j)=b(4)
      a2n(j)=b(5)
      a1n(j)=b(6)
      a0(j)=b(7)
      a1p(j)=b(8)
      a2p(j)=b(9)


      enddo


c			fronteras del dominio  

c*************************************************************************
c*************************************************************************
      j=1
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
c      c(1,5)=1d0
c      c(1,6)=1d0
      b(1)=0d0
      
      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j+1)-x(j))
      c(2,4)=(x(j+2)-x(j))
c      c(2,5)=(x(j+3)-x(j))
c      c(2,6)=(x(j+4)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j+1)-x(j))**2
      c(3,4)=-(x(j+2)-x(j))**2
c      c(3,5)=-(x(j+3)-x(j))**2
c      c(3,6)=-(x(j+4)-x(j))**2
      b(3)=-2d0

c      do n=3,5
      do n=3,3
      c(n+1,1)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,2)=0d0
      c(n+1,3)=-(x(j+1)-x(j))**n
      c(n+1,4)=-(x(j+2)-x(j))**n
c      c(n+1,5)=-(x(j+3)-x(j))**n
c      c(n+1,6)=-(x(j+4)-x(j))**n
      b(n+1)=0d0      
      enddo

c		call LSLRG(5,c,9,b,1,sol)
c      call gaussj(c,6,9,b,1,1)
      call gaussj(c,4,9,b,1,1)

      ben(j)=0d0
      aln(j)=0d0
      alp(j)=b(1)
      bep(j)=0d0
c      a2n(j)=b(6)
c      a1n(j)=b(5)
      a2n(j)=0d0
      a1n(j)=0d0
      a0(j)=b(2)
      a1p(j)=b(3)
      a2p(j)=b(4)
*************************************************************************
c************************************************************************
      j=2
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=0d0
      c(2,5)=(x(j+1)-x(j))
      b(2)=0d0
      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=0d0
      c(3,5)=-(x(j+1)-x(j))**2
      b(3)=-2d0

      do n=3,4
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,3)=-(x(j-1)-x(j))**n
      c(n+1,4)=0d0
      c(n+1,5)=-(x(j+1)-x(j))**n
      b(n+1)=0d0
      enddo

      call gaussj(c,5,9,b,1,1)
      
      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=b(2)
      bep(j)=0d0
      a2n(j)=0d0
      a1n(j)=b(3)
      a0(j)=b(4)
      a1p(j)=b(5)
      a2p(j)=0d0

c*************************************************************************
c*************************************************************************
      j=my-1
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=0d0
      c(2,5)=(x(j+1)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=0d0
      c(3,5)=-(x(j+1)-x(j))**2
      b(3)=-2d0

      do n=3,4
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,3)=-(x(j-1)-x(j))**n
      c(n+1,4)=0d0
      c(n+1,5)=-(x(j+1)-x(j))**n
      b(n+1)=0d0
      enddo
      
      call gaussj(c,5,9,b,1,1)
      
      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=b(2)
      bep(j)=0d0
      a2n(j)=0d0
      a1n(j)=b(3)
      a0(j)=b(4)
      a1p(j)=b(5)
      a2p(j)=0d0
      
c*************************************************************************
c*************************************************************************	
      j=my
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
c      c(1,5)=1d0
c      c(1,6)=1d0
      b(1)=0d0
      
      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=(x(j-2)-x(j))
c      c(2,5)=(x(j-3)-x(j))
c      c(2,6)=(x(j-4)-x(j))
      b(2)=0d0
      
      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=-(x(j-2)-x(j))**2
c      c(3,5)=-(x(j-3)-x(j))**2
c      c(3,6)=-(x(j-4)-x(j))**2
      b(3)=-2d0

c      do n=3,5
      do n=3,3
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=0d0
      c(n+1,3)=-(x(j-1)-x(j))**n
      c(n+1,4)=-(x(j-2)-x(j))**n
c      c(n+1,5)=-(x(j-3)-x(j))**n
c      c(n+1,6)=-(x(j-4)-x(j))**n
      b(n+1)=0d0
      enddo
      
c		call LSLRG(5,c,9,b,1,sol)
c      call gaussj(c,6,9,b,1,1)
      call gaussj(c,4,9,b,1,1)
      
      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=0d0
      bep(j)=0d0
      a2n(j)=b(4)
      a1n(j)=b(3)
      a0(j)=b(2)
c      a1p(j)=b(5)
c      a2p(j)=b(6)
      a1p(j)=0d0
      a2p(j)=0d0
      
      return
      end

c***************************************************************************
c*	deryr2 : Campos, todos
c*      deryyr2: Campos 2 derivada
c* 	deryr  : 1 dimensión, 1 derivada real*4
c*      deryyr : 2 diemsnsiones, 2 derivadas.
c*	
c*	Calculan la derivada primera y segunda
c*	u,du: vector dato y salida, (du puede sobrescribir u)
c*	m: numero de vectores
c*	n: tamanno de los vectores
c* 
c***************************************************************************
      subroutine deryr2(u,du,m)
      implicit none
      include "ctes3D"

      real*4 u(2,my,m),du(2,my,m)


      real*4  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .             trp(0:my-1),mss(0:my-1)
      save   /fis/
      
      real*8 prem1,dt12
      common /cfdiff/ prem1(my,7),dt12(7,my)
      save   /cfdiff/  
 
      real*8  wk1(my),wk2(my)
      integer m,i,j,k

c		!!!  calculo el termino independiente
      do k=1,m
      
         wk1(1)=dt12(4,1)*u(1,1,k)+dt12(5,1)*u(1,2,k)+
     .          dt12(6,1)*u(1,3,k)+dt12(7,1)*u(1,4,k)

         wk2(1)=dt12(4,1)*u(2,1,k)+dt12(5,1)*u(2,2,k)+
     .          dt12(6,1)*u(2,3,k)+dt12(7,1)*u(2,4,k)
     
         wk1(2)=dt12(3,2)*u(1,1,k)+dt12(4,2)*u(1,2,k)+
     .          dt12(5,2)*u(1,3,k)+dt12(6,2)*u(1,4,k)+
     .          dt12(7,2)*u(1,5,k)

         wk2(2)=dt12(3,2)*u(2,1,k)+dt12(4,2)*u(2,2,k)+
     .          dt12(5,2)*u(2,3,k)+dt12(6,2)*u(2,4,k)+
     .          dt12(7,2)*u(2,5,k)

         wk1(3)=dt12(2,3)*u(1,1,k)+dt12(3,3)*u(1,2,k)+
     .          dt12(4,3)*u(1,3,k)+dt12(5,3)*u(1,4,k)+
     .          dt12(6,3)*u(1,5,k)+dt12(7,3)*u(1,6,k)

         wk2(3)=dt12(2,3)*u(2,1,k)+dt12(3,3)*u(2,2,k)+
     .	       dt12(4,3)*u(2,3,k)+dt12(5,3)*u(2,4,k)+
     .          dt12(6,3)*u(2,5,k)+dt12(7,3)*u(2,6,k)

         do j=4,my-3
            wk1(j)=dt12(1,j)*u(1,j-3,k)
            wk2(j)=dt12(1,j)*u(2,j-3,k)
            do i=2,7
               wk1(j)=wk1(j) + dt12(i,j)*u(1,i+j-4,k)
               wk2(j)=wk2(j) + dt12(i,j)*u(2,i+j-4,k)
            enddo
         enddo

         wk1(my-2)=dt12(1,my-2)*u(1,my-5,k)+dt12(2,my-2)*u(1,my-4,k)+
     .             dt12(3,my-2)*u(1,my-3,k)+dt12(4,my-2)*u(1,my-2,k)+
     .             dt12(5,my-2)*u(1,my-1,k)+dt12(6,my-2)*u(1,my,k)

         wk2(my-2)=dt12(1,my-2)*u(2,my-5,k)+dt12(2,my-2)*u(2,my-4,k)+
     .             dt12(3,my-2)*u(2,my-3,k)+dt12(4,my-2)*u(2,my-2,k)+
     .             dt12(5,my-2)*u(2,my-1,k)+dt12(6,my-2)*u(2,my,k)

         wk1(my-1)=dt12(1,my-1)*u(1,my-4,k)+dt12(2,my-1)*u(1,my-3,k)+
     .             dt12(3,my-1)*u(1,my-2,k)+dt12(4,my-1)*u(1,my-1,k)+
     .             dt12(5,my-1)*u(1,my,k)

         wk2(my-1)=dt12(1,my-1)*u(2,my-4,k)+dt12(2,my-1)*u(2,my-3,k)+
     .             dt12(3,my-1)*u(2,my-2,k)+dt12(4,my-1)*u(2,my-1,k)+
     .             dt12(5,my-1)*u(2,my,k)

         wk1(my)  =dt12(1,my  )*u(1,my-3,k)+dt12(2,my )*u(1,my-2,k) +
     .             dt12(3,my  )*u(1,my-1,k)+dt12(4,my )*u(1,my,k)
 
         wk2(my)  =dt12(1,my  )*u(2,my-3,k)+dt12(2,my )*u(2,my-2,k) +
     .             dt12(3,my  )*u(2,my-1,k)+dt12(4,my )*u(2,my,k)


         call banbks7(prem1,my,wk1)

         do j=1,my
            du(1,j,k)=wk1(j)*fmap(j)
         enddo
c              backsubstitution, storage and MAPPING

         call banbks7(prem1,my,wk2)

         do j=1,my
            du(2,j,k)=wk2(j)*fmap(j)
         enddo

      enddo    

      endsubroutine

c------------------------------------------------------
c                 Funcion
c------------------------------------------------------


      subroutine deryr(u,du)

      implicit none
      include "ctes3D"
      integer i,j

      real*4 u(my),du(my)
      real*8 fwk(my)

      real*8 prem1,dt12
      common /cfdiff/ prem1(my,7),dt12(7,my)
      save   /cfdiff/  

      real*4  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     .             trp(0:my-1),mss(0:my-1)
      save   /fis/

c-------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 casi-pentadiagonal!!


      fwk(1)= 	dt12(4,1)*u(1)+dt12(5,1)*u(2) +
     .        	dt12(6,1)*u(3)+dt12(7,1)*u(4)
     
      fwk(2)= 	dt12(3,2)*u(1)+dt12(4,2)*u(2) +
     .        	dt12(5,2)*u(3)+dt12(6,2)*u(4) +
     .        	dt12(7,2)*u(5)
     
      fwk(3)= 	dt12(2,3)*u(1)+dt12(3,3)*u(2) +
     .        	dt12(4,3)*u(3)+dt12(5,3)*u(4) +
     .        	dt12(6,3)*u(5)+dt12(7,3)*u(6)

      do j=4,my-3
         fwk(j)=dt12(1,j)*u(j-3)
         do i=2,7
            fwk(j)=fwk(j) + dt12(i,j)*u(i+j-4)
         enddo
      enddo
      fwk(my-2)=dt12(1,my-2)*u(my-5)+dt12(2,my-2)*u(my-4)+
     .          dt12(3,my-2)*u(my-3)+dt12(4,my-2)*u(my-2)+
     .          dt12(5,my-2)*u(my-1)+dt12(6,my-2)*u(my)

      fwk(my-1)=dt12(1,my-1)*u(my-4)+dt12(2,my-1)*u(my-3)+
     .          dt12(3,my-1)*u(my-2)+dt12(4,my-1)*u(my-1)+
     .          dt12(5,my-1)*u(my)
      fwk(my)=  dt12(1,my  )*u(my-3)+dt12(2,my  )*u(my-2)+
     .          dt12(3,my  )*u(my-1)+dt12(4,my  )*u(my)


      call banbks7(prem1,my,fwk)

      do j=1,my
         du(j)=fwk(j)*fmap(j)
      enddo


      end 

c----------------------------------------------------------------------^M
c                 Funcion
c------------------------------------------------------------------------^M


      subroutine deryyr(u,ddu)
      implicit none
      include "ctes3D"

      integer i,j

      real*4 u(my),ddu(my),nl(my),dalre,dtgamma
      real*8 fwk(my),rkn1

      real*8 prem3,dt22
      common /d2y/ prem3(my,5),dt22(5,my)
      save   /d2y/

c------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 pentadiagonal!!!

      fwk(1) =  dt22(3,1)*u(1) + dt22(4,1)*u(2) + 
     .          dt22(5,1)*u(3)
      fwk(2) =  dt22(2,2)*u(1) + dt22(3,2)*u(2) +
     .          dt22(4,2)*u(3) + dt22(5,2)*u(4)
     
      do j=3,my-2
         fwk(j)=dt22(1,j)*u(j-2)
         do i=2,5
            fwk(j) = fwk(j) + dt22(i,j)*u(i+j-3)
         enddo
      enddo
      
      fwk(my-1)=dt22(1,my-1)*u(my-3) + dt22(2,my-1)*u(my-2)+
     .          dt22(3,my-1)*u(my-1) + dt22(4,my-1)*u(my  )  
      fwk(my)=  dt22(1,my  )*u(my-2) + dt22(2,my  )*u(my-1)+
     .          dt22(3,my  )*u(my)     

      call banbks(prem3,my,fwk)

      do j=1,my
         ddu(j) = fwk(j)
      enddo

      end

c----------------------------------------------------------------------
c                 Funcion deryyr2
c     Crea el termino derecho para SMR cuando va a cambiar el deltat
c----------------------------------------------------------------------

      subroutine deryyr2(u,du,m)

      implicit none
      include "ctes3D"
      integer m,i,j,k,xpl

      real*4 u(2,my,m),du(2,my,m)
      real*8 wk1(my),wk2(my),rkn1
      real*4 dalre,dtgamma,coef

      real*8 prem3,dt22
      common /d2y/ prem3(my,5),dt22(5,my)
      save   /d2y/

      integer icx
      real*4 alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              icx(0:mz1)
      save  /wave/

      do k=1,m 

         wk1(1) = dt22(3,1)*u(1,1,k)+dt22(4,1)*u(1,2,k)+
     .            dt22(5,1)*u(1,3,k)

         wk2(1) = dt22(3,1)*u(2,1,k)+dt22(4,1)*u(2,2,k)+
     .            dt22(5,1)*u(2,3,k)
         
         wk1(2) = dt22(2,2)*u(1,1,k)+dt22(3,2)*u(1,2,k)+
     .            dt22(4,2)*u(1,3,k)+dt22(5,2)*u(1,4,k)
 
         wk2(2) = dt22(2,2)*u(2,1,k)+dt22(3,2)*u(2,2,k)+
     .            dt22(4,2)*u(2,3,k)+dt22(5,2)*u(2,4,k)

         do j=3,my-2
            wk1(j)= dt22(1,j)*u(1,j-2,k)
            wk2(j)= dt22(1,j)*u(2,j-2,k)
            do i=2,5
               wk1(j) = wk1(j) + dt22(i,j)*u(1,i+j-3,k)
               wk2(j) = wk2(j) + dt22(i,j)*u(2,i+j-3,k)
            enddo
         enddo

         wk1(my-1)=dt22(1,my-1)*u(1,my-3,k)+dt22(2,my-1)*u(1,my-2,k)+
     .             dt22(3,my-1)*u(1,my-1,k)+dt22(4,my-1)*u(1,my  ,k)
         wk2(my-1)=dt22(1,my-1)*u(2,my-3,k)+dt22(2,my-1)*u(2,my-2,k)+
     .             dt22(3,my-1)*u(2,my-1,k)+dt22(4,my-1)*u(2,my  ,k)

         wk1(my)=  dt22(1,my  )*u(1,my-2,k)+dt22(2,my  )*u(1,my-1,k)+
     .             dt22(3,my  )*u(1,my,k)
         wk2(my)=  dt22(1,my  )*u(2,my-2,k)+dt22(2,my  )*u(2,my-1,k)+
     .             dt22(3,my  )*u(2,my,k)

         call banbks(prem3,my,wk1)
         do j=1,my
            du(1,j,k)= wk1(j)
         enddo

         call banbks(prem3,my,wk2)

         do j=1,my
            du(2,j,k)=wk2(j)
         enddo

      enddo

      end

c**************************************************
c*    SUBRUTINAS PARA TRABAJAR CON n-DIAGONALES	
c*    (Numerical Recipes in FORTRAN)	
c*	banmul, bandec, banbks	
c*      Solo para liso---
c**************************************************
      SUBROUTINE banmul(a,n,m1,m2,np,mp,x,b)
      INTEGER m1,m2,mp,n,np
      REAL*8 a(np,mp),b(n),x(n)
      INTEGER i,j,k
      do 12 i=1,n
        b(i)=0d0
        k=i-m1-1
        do 11 j=max(1,1-k),min(m1+m2+1,n-k)
          b(i)=b(i)+a(i,j)*x(j+k)
11      continue
12    continue
      return
      END


c     Adapted to liso

      SUBROUTINE banbks(a,n,b)
      INTEGER n
      REAL*8 a(n,5),b(n)
      INTEGER i,k

      do k=1,n-2
         b(k+1) = b(k+1)-a(k,4)*b(k)
         b(k+2) = b(k+2)-a(k,5)*b(k)
      enddo
      b(n) = b(n)- a(n-1,4)*b(n-1)

c     back substitution

      b(n) = b(n)*a(n,1)
      b(n-1) = (b(n-1)-a(n-1,2)*b(n))*a(n-1,1)
      do i=n-2,1,-1
         b(i) = (b(i)-a(i,2)*b(1+i)-a(i,3)*b(2+i))*a(i,1)
      enddo
      
      return
      END


      SUBROUTINE bandec(a,n)
      INTEGER n
      REAL*8 a(n,5)
      INTEGER j,k
      
      do j=1,3
         a(1,j)=a(1,j+2)
      enddo
      do j=1,4
         a(2,j)=a(2,j+1)
      enddo


      do k=1,n-2
        a(k,1)   = 1./a(k,1)

        a(k,4)   = a(k+1,1)*a(k,1)

        a(k+1,1) = a(k+1,2)-a(k,4)*a(k,2)
        a(k+1,2) = a(k+1,3)-a(k,4)*a(k,3)
        a(k+1,3) = a(k+1,4)
c       
        a(k,5)   = a(k+2,1)*a(k,1)

        a(k+2,1) = a(k+2,2)-a(k,5)*a(k,2)
        a(k+2,2) = a(k+2,3)-a(k,5)*a(k,3)
        a(k+2,3) = a(k+2,4) 
        a(k+2,4) = a(k+2,5)
      enddo

c     k=n-1

      a(n-1,1) = 1./a(n-1,1)

      a(n-1,4)=a(n,1)*a(n-1,1)

      a(n,1) = a(n,2)-a(n-1,4)*a(n-1,2)
      a(n,2) = a(n,3)-a(n-1,4)*a(n-1,3)
      a(n,3) = a(n,4)

c     the next loop will be used in banbk

      a(n,1)=1./a(n,1)


      END


      SUBROUTINE banbks7(a,n,b)
      INTEGER n
      REAL*8 a(n,7),b(n)
      INTEGER i,k

      do k=1,n-3
         b(k+1) = b(k+1)-a(k,5)*b(k)
         b(k+2) = b(k+2)-a(k,6)*b(k)
         b(k+3) = b(k+3)-a(k,7)*b(k)
      enddo

c     n-2

      b(n-1) = b(n-1)-a(n-2,5)*b(n-2)
      b(n)   = b(n)  -a(n-2,6)*b(n-2)

c     n-1

      b(n) = b(n)    -a(n-1,5)*b(n-1)

c     back substitution

      b(n) = b(n)*a(n,1)
      b(n-1) = (b(n-1)-a(n-1,2)*b(n))*a(n-1,1)
      b(n-2) = (b(n-2)-a(n-2,2)*b(n-1)-a(n-2,3)*b(n))*a(n-2,1)

      do i=n-3,1,-1
         b(i) = (b(i)-a(i,2)*b(1+i)-a(i,3)*b(2+i)-a(i,4)*b(3+i))*a(i,1)
      enddo

      return
      END



      SUBROUTINE bandec7(a,n)
      INTEGER n
      REAL*8 a(n,7)
      INTEGER j,k
c
      do j=1,4
         a(1,j)=a(1,j+3)
      enddo

      do j=1,5
         a(2,j)=a(2,j+2)
      enddo

      do j=1,6
         a(3,j)=a(3,j+1)
      enddo


c     LU

      do k=1,n-3
        a(k,1)   = 1./a(k,1)

        a(k,5)   = a(k+1,1)*a(k,1)
        a(k+1,1) = a(k+1,2)-a(k,5)*a(k,2)
        a(k+1,2) = a(k+1,3)-a(k,5)*a(k,3)
        a(k+1,3) = a(k+1,4)-a(k,5)*a(k,4)
        a(k+1,4) = a(k+1,5)

        a(k,6)   = a(k+2,1)*a(k,1)
        a(k+2,1) = a(k+2,2)-a(k,6)*a(k,2)
        a(k+2,2) = a(k+2,3)-a(k,6)*a(k,3)
        a(k+2,3) = a(k+2,4)-a(k,6)*a(k,4)
        a(k+2,4) = a(k+2,5)
        a(k+2,5) = a(k+2,6)

        a(k,7)   = a(k+3,1)*a(k,1)
        a(k+3,1) = a(k+3,2)-a(k,7)*a(k,2)
        a(k+3,2) = a(k+3,3)-a(k,7)*a(k,3)
        a(k+3,3) = a(k+3,4)-a(k,7)*a(k,4)
        a(k+3,4) = a(k+3,5)
        a(k+3,5) = a(k+3,6)
        a(k+3,6) = a(k+3,7)

      enddo

c     k=n-2
      a(n-2,1) = 1./a(n-2,1)

      a(n-2,5) = a(n-1,1)*a(n-2,1)

      a(n-1,1) = a(n-1,2)-a(n-2,5)*a(n-2,2)
      a(n-1,2) = a(n-1,3)-a(n-2,5)*a(n-2,3)
      a(n-1,3) = a(n-1,4)-a(n-2,5)*a(n-2,4)
      a(n-1,4) = a(n-1,5)

      a(n-2,6) = a(n,1)*a(n-2,1)

      a(n,1)   = a(n,2)-a(n-2,6)*a(n-2,2)
      a(n,2)   = a(n,3)-a(n-2,6)*a(n-2,3)
      a(n,3)   = a(n,4)-a(n-2,6)*a(n-2,4)
      a(n,4)   = a(n,5)
      a(n,5)   = a(n,6)

c     k=n-1
     
      a(n-1,1) = 1./a(n-1,1)
      a(n-1,5) = a(n,1)*a(n-1,1)

      a(n,1)   = a(n,2)-a(n-1,5)*a(n-1,2)
      a(n,2)   = a(n,3)-a(n-1,5)*a(n-1,3)
      a(n,3)   = a(n,4)-a(n-1,5)*a(n-1,4)
      a(n,4)   = a(n,5)

c     the next loop will be used in banbk

      a(n,1)=1./a(n,1)

      return
      END



c**************************************************
c*    SUBRUTINA PARA RESOVER SIST LINEALES
c*    (Numerical Recipes in FORTRAN)	
c**************************************************
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
  


      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1d-20)
      INTEGER i,imax,j,k
      REAL*4 aamax,dum,sum,vv(NMAX)
      d=1d0
      do 12 i=1,n
        aamax=0d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0d0) pause 'singular matrix in ludcmp'
        vv(i)=1d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
