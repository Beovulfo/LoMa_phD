!=======================================================================
!       PAQUETE DE SUBRUTINAS PARA DIFERENCIAS FINITAS COMPACTAS
!       O.F, nov 2002
!       S.H. may 2005
!
!=======================================================================
! --------------------------Subrutinas --------------------------------

      subroutine malla8(n,y,g)

      implicit none

      integer j,n
      real*8 y(0:*),g(0:*)
      character*99 text


      open(13,file='malla2.dat',status='old')
      do j=0,n-1
         read(13,'(a)') text
         read(text,*) y(j),g(j)
      enddo

      close(13)

      end



!**********************************************************************
! pred1m60:
!
! Genera los coeficientes alpha,beta,delta y
!      a0,a1,a2 y a3 para el esquema de
! diferencias finitas compactas con malla uniforme.
!               derivada primera
!
! delta*(du(j+2)+du(j-2))beta*(du(j+2)+du(j-2))+
!      alpha*(du(j+1)+du(j-1))+du(j) =
!    = a1*(u(j+1)-u(j-1)) + a2(u(j+2)-u(j-2)) + a3(u(j+2)-u(j-3))
!
!    En la frontera del dominio:
!    alpha*du(2) + du(1) =
!     = a1*u(1) + a2*u(2) + a3*u(3) + a4*u(4) + a5*u(5)
!
!    alpha*(du(3)+du(1)) + du(2) =
!     = a1*(u(3)-u(1))
!
!*********************************************************************
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


!-----------------------------------------------------------------------
!     centro del dominio
!     IMPONGO LOS COEFICIENTES A PARTIR DE CALCULOS EXTERNOS
! ----------------------------------------------------------------------

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


!-----------------------------------------------------------------------
!
!       fronteras del dominio
!       ESTANDAR:SOLO IMPONGO ECUACIONES DE ORDEN
!       USANDO TODA LA INFORMACION DISPONIBLE:
!
!-----------------------------------------------------------------------


      j=1
      den(j) =  0d0
      ben(j) =  0d0
      aln(j) =  0d0
      alp(j) =  9d0
      bep(j) =  alp(j)
      dep(j) =  1.d0

      a3n(j) =  0d0
      a2n(j) =  0d0
      a1n(j) =  0d0
      a0 (j) =  -11d0/3d0/h
      a1p(j) =  -9d0/h
      a2p(j) =  -a1p(j)
      a3p(j) =  -a0(j)

      j=2
      den(j) =  0d0
      ben(j) =  0d0
      aln(j) =  1d0/16d0
      alp(j) =  9d0/4d0
      bep(j) =  1d0
      dep(j) =  aln(j)

      a3n(j) =  0d0
      a2n(j) =  0d0
      a1n(j) =  -25d0/96d0/h
      a0 (j) =  -5d0/3d0/h
      a1p(j) =  0d0
      a2p(j) =  -a0(j)
      a3p(j) =  -a1n(j)


      j=3
      den(j) =   0d0
      ben(j) =   2.538032824d-2
      aln(j) =   0.347040189d0
      alp(j) =   1d0
      bep(j) =   aln(j)
      dep(j) =   ben(j)

      a3n(j) =   0d0
      a2n(j) =   -9.841716629d-2/h
      a1n(j) =   -0.5830132425d0/h
      a0 (j) =   -0.5037154692d0/h
      a1p(j) =    -a0(j)
      a2p(j) =    -a1n(j)
      a3p(j) =    -a2n(j)

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


!***********************************************************************
! prederiv2:
!
! Genera los coeficientes alpha,beta,a0,a1 y a2 para el esquema de
! diferencias finitas compactas. Derivada segunda
!
! beta*(d2u(j+2)+d2u(j-2))+alpha*(d2u(j+1)+d2u(j-1))+d2u(j) =
!  = a0*u(j)+a1*(u(j+1)+u(j-1)) + a2(u(j+2)+u(j-2))
!
!
!
!***********************************************************************
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

!               Ahora calculo los coeficientes del centro del dominio

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

!               call LSLRG(9,c,9,b,1,sol)
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

!-----------------------------------------------------------------------
!
!                       fronteras del dominio
!-----------------------------------------------------------------------

      j=1
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
!      c(1,5)=1d0
!      c(1,6)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j+1)-x(j))
      c(2,4)=(x(j+2)-x(j))
!      c(2,5)=(x(j+3)-x(j))
!      c(2,6)=(x(j+4)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j+1)-x(j))**2
      c(3,4)=-(x(j+2)-x(j))**2
!      c(3,5)=-(x(j+3)-x(j))**2
!      c(3,6)=-(x(j+4)-x(j))**2
      b(3)=-2d0

      do n=3,3
      c(n+1,1)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,2)=0d0
      c(n+1,3)=-(x(j+1)-x(j))**n
      c(n+1,4)=-(x(j+2)-x(j))**n
!      c(n+1,5)=-(x(j+3)-x(j))**n
!      c(n+1,6)=-(x(j+4)-x(j))**n
      b(n+1)=0d0
      enddo

      call gaussj(c,4,9,b,1,1)

      ben(j)=0d0
      aln(j)=0d0
      alp(j)=b(1)
      bep(j)=0d0
!      a2n(j)=b(6)
!      a1n(j)=b(5)
      a2n(j)=0d0
      a1n(j)=0d0
      a0(j)=b(2)
      a1p(j)=b(3)
      a2p(j)=b(4)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------

      j=my
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
!      c(1,5)=1d0
!      c(1,6)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=(x(j-2)-x(j))
!      c(2,5)=(x(j-3)-x(j))
!      c(2,6)=(x(j-4)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=-(x(j-2)-x(j))**2
!      c(3,5)=-(x(j-3)-x(j))**2
!      c(3,6)=-(x(j-4)-x(j))**2
      b(3)=-2d0

!      do n=3,5
      do n=3,3
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=0d0
      c(n+1,3)=-(x(j-1)-x(j))**n
      c(n+1,4)=-(x(j-2)-x(j))**n
!      c(n+1,5)=-(x(j-3)-x(j))**n
!      c(n+1,6)=-(x(j-4)-x(j))**n
      b(n+1)=0d0
      enddo

      call gaussj(c,4,9,b,1,1)

      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=0d0
      bep(j)=0d0
      a2n(j)=b(4)
      a1n(j)=b(3)
      a0(j)=b(2)
!      a1p(j)=b(5)
!      a2p(j)=b(6)
      a1p(j)=0d0
      a2p(j)=0d0


      end

!************************************************************************
!                                                                       !
!    Subroutine deryr2                                                  !
!                                                                       !
!    computes the first derivative in y of u, with size(u)=(2,my,m)     !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be complex.                  !
!         m: size of the third dimension of u                           !
!    Output:                                                            !
!        du: First derivative of u                                      !
!                                                                       !
!************************************************************************

      subroutine deryr2(u,du,m)
      implicit none
      include "ctes3D"

      integer m,i,j,k
      real*4 u(2,my,m),du(2,my,m)


      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $             trp(my),mss(my)
      save   /fis/

      real*8 prem1,dt12,endder
      common /cfdiff/ prem1(7,my),dt12(7,my),endder(my)
      save   /cfdiff/

      real*8  wk1(my),wk2(my)

!               !!!  calculo el termino independiente
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
     .          dt12(4,3)*u(2,3,k)+dt12(5,3)*u(2,4,k)+
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
!              backsubstitution, storage and MAPPING

         call banbks7(prem1,my,wk2)

         do j=1,my
            du(2,j,k)=wk2(j)*fmap(j)
         enddo

      enddo

      endsubroutine

!************************************************************************
!                                                                       !
!    Subroutine deryr                                                   !
!                                                                       !
!    computes the first derivative in y of u, with size(u)=(my)         !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be real of size my           !
!    Output:                                                            !
!        du: First derivative of u                                      !
!                                                                       !
!************************************************************************

      subroutine deryr(u,du)

      implicit none
      include "ctes3D"
      integer i,j

      real*8 u(my),du(my)
      real*8 fwk(my)

      real*8 prem1,dt12,endder
      common /cfdiff/ prem1(7,my),dt12(7,my),endder(my)
      save   /cfdiff/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/

!-------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 casi-pentadiagonal!!


      fwk(1)=   dt12(4,1)*u(1)+dt12(5,1)*u(2) +
     .          dt12(6,1)*u(3)+dt12(7,1)*u(4)

      fwk(2)=   dt12(3,2)*u(1)+dt12(4,2)*u(2) +
     .          dt12(5,2)*u(3)+dt12(6,2)*u(4) +
     .          dt12(7,2)*u(5)

      fwk(3)=   dt12(2,3)*u(1)+dt12(3,3)*u(2) +
     .          dt12(4,3)*u(3)+dt12(5,3)*u(4) +
     .          dt12(6,3)*u(5)+dt12(7,3)*u(6)

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

!************************************************************************
!                                                                       !
!    Subroutine deryyr                                                 !
!                                                                       !
!    Computes the second derivative in y of u, with size(u)=(my)        !
!    and performs some operations intended for SMR CODE,                !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be real of size my           !
!      rkn1,dalre,dtgamma: Constants of RK                              !
!                                                                       !
!    Output:                                                            !
!        nl: nolinear term of Modes00 equation                          !
!                                                                       !
!************************************************************************

      subroutine deryyr(u,nl,rkn1,dalre,dtgamma)
      implicit none
      include "ctes3D"

      integer i,j

      real*8 u(my),nl(my),dalre,dtgamma,fwk2(my)
      real*8 fwk(my),rkn1

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

!------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 pentadiagonal!!!

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

! ---- primera  derivada ----

      do j=1,my
         u(j) = -rkn1*(u(j)+dalre*fwk(j)+dtgamma*nl(j))
      enddo

      end

!************************************************************************
!                                                                       !
!    Subroutine deryyr2                                                 !
!                                                                       !
!    computes the second derivative in y of u, with size(u)=(2,my,m)    !
!    and creates the rhs of vor and phi equations when the CFL Changes  !
!                                                                       !
!    Input:                                                             !
!         u : Field to derivate. Assumed to be complex                  !
!        nl : Nonlinear term                                            !
!        xpl: x-plane                                                   !
!        rkn1,dalre,dtgamma: Parameters of R_K                          !
!                                                                       !
!    Output:                                                            !
!        du: R.H.S of  vor and phi equations                            !
!                                                                       !
!************************************************************************

      subroutine deryyr2(u,du,nl,m,xpl,rkn1,dalre,dtgamma)

      implicit none
      include "ctes3D"
      integer m,i,j,k,xpl

      real*4 u(2,my,m),du(2,my,m),nl(2,my,m)
      real*8 wk1(my),wk2(my),rkn1
      real*8 dalre,dtgamma,coef

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
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
         coef = 1d0-dalre*(bet2(icx(k-1))+alp2(xpl-1))
         do j=1,my
            du(1,j,k)=-rkn1*(coef*u(1,j,k)+dalre*wk1(j)
     .                 +dtgamma*nl(1,j,k))
         enddo

         call banbks(prem3,my,wk2)

         do j=1,my
            du(2,j,k)=-rkn1*(coef*u(2,j,k)+dalre*wk2(j)
     .                 +dtgamma*nl(2,j,k))
         enddo

      enddo

      end


!************************************************************************
!
!     Subroutine premass
!
!     Compute coefficients for the mean of the second derivative
!
!          If dv/dyy is computed by a.dvyy= dt2.v
!     where  a and dt2 are in common /d2y/ and a is already in a=LU
!          we want to compute m such that
!              m'.v  = u'.dvyy
!     where u is some averaging vector
!
!     It can be shown that  a'.q = u,   and  m= dt2'.q
!
!     the linear system is solved using that a'= U'L'
!
!     This is used to compute the mass correction from the viscous term
!
!************************************************************************

      subroutine premass(u,m)
      implicit none
      include "ctes3D"

      integer i,j,k

      real*8 u(my),m(my)
      real*8 fwk(my)

      real*8 a,dt21,dt22
      common /d2y/ a(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/



!   -------  this is banbks adapted to work on the transpose ---


      fwk(1)= a(1,1)*u(1)
      fwk(2)= a(1,2)*(u(2)-a(2,1)*fwk(1))

      do j=3,my
         fwk(j) = a(1,j)*(u(j)-a(2,j-1)*fwk(j-1)-a(3,j-2)*fwk(j-2))
      enddo

!     back substitution

      fwk(my-1) = fwk(my-1)-a(4,my-1)*fwk(my)
      do j=my-2,1,-1
         fwk(j) = fwk(j)-a(4,j)*fwk(j+1)-a(5,j)*fwk(j+2)
      enddo

!     ---  this multiplies the right hand side

         m(1) = fwk(1)*dt22(3,1)+fwk(2)*dt22(2,2)
     $        + fwk(3)*dt22(1,3)
         m(2) = fwk(1)*dt22(4,1)+fwk(2)*dt22(3,2)
     $        + fwk(3)*dt22(2,3)+fwk(4)*dt22(1,4)

      do j=3,my-2
         m(j) = fwk(j-2)*dt22(5,j-2)+fwk(j-1)*dt22(4,j-1)
     $        + fwk(j  )*dt22(3,j  )+fwk(j+1)*dt22(2,j+1)
     $        + fwk(j+2)*dt22(1,j+2)
      enddo

      m(my-1) = fwk(my-3)*dt22(5,my-3)+fwk(my-2)*dt22(4,my-2)
     $        + fwk(my-1)*dt22(3,my-1)+fwk(my  )*dt22(2,my  )
      m(my  ) = fwk(my-2)*dt22(5,my-2)+fwk(my-1)*dt22(4,my-1)
     $        + fwk(my  )*dt22(3,my  )


      endsubroutine


!************************************************************************
!
!     Subroutine endderiv
!
!     Compute coefficients for the derivative at first point
!
!          If dv/dy is computed by a.dvyy= dt12.v
!     where  a and dt12 are in common /cfdiff/ and a is already in a=LU
!          we want to compute m such that
!              m'.v  = u'.dvyy
!     where u is the Kronecker delta(k,1)
!
!     It can be shown that  a'.q = u,   and  m= dt12'.q
!
!     the linear system is solved using that a'= U'L'
!
!     This is used to compute the mass correction from the viscous term
!
!************************************************************************

      subroutine endderiv()
      implicit none
      include "ctes3D"

      integer i,j,k

      real*8 fwk(my)

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $             trp(my),mss(my)
      save   /fis/

      real*8 a,dt12,endder
      common /cfdiff/ a(7,my),dt12(7,my),endder(my)
      save   /cfdiff/



!   -------  this is banbks7 adapted to work on the transpose ---


      fwk(1)= a(1,1)                    ! ---- this is the delta ------
      fwk(2)= -a(1,2)*a(2,1)*fwk(1)     ! -- everything else has no rhs
      fwk(3)= -a(1,3)*(a(2,2)*fwk(2)+a(3,1)*fwk(1))

      do j=4,my
         fwk(j) = -a(1,j)*(a(2,j-1)*fwk(j-1)+a(3,j-2)*fwk(j-2)
     $                   + a(4,j-3)*fwk(j-3))
      enddo

!     back substitution

      fwk(my-1) = fwk(my-1)-a(5,my-1)*fwk(my)
      fwk(my-2) = fwk(my-2)-a(5,my-2)*fwk(my-1)-a(6,my-2)*fwk(my)
      do j=my-3,1,-1
         fwk(j) = fwk(j)-a(5,j)*fwk(j+1)-a(6,j)*fwk(j+2)
     $                  -a(7,j)*fwk(j+3)
      enddo

!     ---  this multiplies the right hand side

      do j=1,my
         endder(j) = 0d0
         do i=max(1,j-3),min(my,j+3)
            endder(j) = endder(j) +fwk(i)*dt12(j-i+4,i)
         enddo
         endder(j)=endder(j)*fmap(1)
      enddo


      endsubroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    prepares commons and things for using cfdiff and laps later
c                             jj/may/05
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivadas(d11,d12,d21,d22)

      implicit none
      include "ctes3D"

      real*8 d11(my,7),d12(my,7),d21(my,5),d22(my,5)
      integer i,j

c     ------   commons --------------

      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12,endder
      common /cfdiff/ prem1(7,my),dt12(7,my),endder(my)
      save   /cfdiff/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $             trp(my),mss(my)
      save   /fis/

c   ---  start doing something -------------

      call pred1m6o (my,d11(1,1),d11(1,2),d11(1,3),d11(1,4),d11(1,5),
     .                 d11(1,6),d11(1,7),d12(1,1),d12(1,2),d12(1,3),
     .                 d12(1,4),d12(1,5),d12(1,6),d12(1,7))

      call prederiv2(my,d21(1,1),d21(1,2),d21(1,3),d21(1,4),d21(1,5),
     .               d22(1,1),d22(1,2),d22(1,3),d22(1,4),d22(1,5),y(1))

c       Copiamos, transponemos, y calculamos LU

      do j=1,my
         do i=1,7
           prem1(i,j) = d11(j,i)
           dt12(i,j)  = d12(j,i)
         enddo
      enddo

      do j=1,my
         do i=1,5
           prem3(i,j) = d21(j,i)
           dt21(i,j)  = d21(j,i)
           dt22(i,j)  = d22(j,i)
         enddo
      enddo

c     Calculo de traspuestas

      call bandec7(prem1,my)
      call bandec (prem3,my)

      call endderiv()

      end


!***********************************************************************
!    SUBRUTINAS PARA TRABAJAR CON n-DIAGONALES
!    (Numerical Recipes in FORTRAN)
!      bandec, banbks
!      Solo para liso---
!***********************************************************************

!     Adapted to liso

      SUBROUTINE banbks(a,n,b)
      INTEGER n
      REAL*8 a(5,n),b(n)
      INTEGER i,k

      do k=1,n-2
         b(k+1) = b(k+1)-a(4,k)*b(k)
         b(k+2) = b(k+2)-a(5,k)*b(k)
      enddo
      b(n) = b(n)- a(4,n-1)*b(n-1)

!     back substitution

      b(n) = b(n)*a(1,n)
      b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
      do i=n-2,1,-1
         b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i))*a(1,i)
      enddo

      return
      END


      SUBROUTINE bandec(a,n)
      INTEGER n
      REAL*8 a(5,n)
      INTEGER j,k

      do j=1,3
         a(j,1)=a(j+2,1)
      enddo
      do j=1,4
         a(j,2)=a(j+1,2)
      enddo


      do k=1,n-2
        a(1,k)   = 1d0/a(1,k)

        a(4,k)   = a(1,k+1)*a(1,k)

        a(1,k+1) = a(2,k+1)-a(4,k)*a(2,k)
        a(2,k+1) = a(3,k+1)-a(4,k)*a(3,k)
        a(3,k+1) = a(4,k+1)
!
        a(5,k)   = a(1,k+2)*a(1,k)

        a(1,k+2) = a(2,k+2)-a(5,k)*a(2,k)
        a(2,k+2) = a(3,k+2)-a(5,k)*a(3,k)
        a(3,k+2) = a(4,k+2)
        a(4,k+2) = a(5,k+2)
      enddo

!     k=n-1

      a(1,n-1) = 1d0/a(1,n-1)

      a(4,n-1)=a(1,n)*a(1,n-1)

      a(1,n) = a(2,n)-a(4,n-1)*a(2,n-1)
      a(2,n) = a(3,n)-a(4,n-1)*a(3,n-1)
      a(3,n) = a(4,n)

!     the next loop will be used in banbk

      a(1,n)=1d0/a(1,n)


      END


      SUBROUTINE banbks7(a,n,b)
      INTEGER n
      REAL*8 a(7,n),b(n)
      INTEGER i,k

      do k=1,n-3
         b(k+1) = b(k+1)-a(5,k)*b(k)
         b(k+2) = b(k+2)-a(6,k)*b(k)
         b(k+3) = b(k+3)-a(7,k)*b(k)
      enddo

!     n-2

      b(n-1) = b(n-1)-a(5,n-2)*b(n-2)
      b(n)   = b(n)  -a(6,n-2)*b(n-2)

!     n-1

      b(n) = b(n)    -a(5,n-1)*b(n-1)

!     back substitution

      b(n) = b(n)*a(1,n)
      b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
      b(n-2) = (b(n-2)-a(2,n-2)*b(n-1)-a(3,n-2)*b(n))*a(1,n-2)

      do i=n-3,1,-1
         b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i)-a(4,i)*b(3+i))*a(1,i)
      enddo

      return
      END



      SUBROUTINE bandec7(a,n)
      INTEGER n
      REAL*8 a(7,n)
      INTEGER j,k
!
      do j=1,4
         a(j,1)=a(j+3,1)
      enddo

      do j=1,5
         a(j,2)=a(j+2,2)
      enddo

      do j=1,6
         a(j,3)=a(j+1,3)
      enddo


!     LU

      do k=1,n-3
        a(1,k)   = 1d0/a(1,k)

        a(5,k)   = a(1,k+1)*a(1,k)
        a(1,k+1) = a(2,k+1)-a(5,k)*a(2,k)
        a(2,k+1) = a(3,k+1)-a(5,k)*a(3,k)
        a(3,k+1) = a(4,k+1)-a(5,k)*a(4,k)
        a(4,k+1) = a(5,k+1)

        a(6,k)   = a(1,k+2)*a(1,k)
        a(1,k+2) = a(2,k+2)-a(6,k)*a(2,k)
        a(2,k+2) = a(3,k+2)-a(6,k)*a(3,k)
        a(3,k+2) = a(4,k+2)-a(6,k)*a(4,k)
        a(4,k+2) = a(5,k+2)
        a(5,k+2) = a(6,k+2)

        a(7,k)   = a(1,k+3)*a(1,k)
        a(1,k+3) = a(2,k+3)-a(7,k)*a(2,k)
        a(2,k+3) = a(3,k+3)-a(7,k)*a(3,k)
        a(3,k+3) = a(4,k+3)-a(7,k)*a(4,k)
        a(4,k+3) = a(5,k+3)
        a(5,k+3) = a(6,k+3)
        a(6,k+3) = a(7,k+3)

      enddo

!     k=n-2
      a(1,n-2) = 1d0/a(1,n-2)

      a(5,n-2) = a(1,n-1)*a(1,n-2)

      a(1,n-1) = a(2,n-1)-a(5,n-2)*a(2,n-2)
      a(2,n-1) = a(3,n-1)-a(5,n-2)*a(3,n-2)
      a(3,n-1) = a(4,n-1)-a(5,n-2)*a(4,n-2)
      a(4,n-1) = a(5,n-1)

      a(6,n-2) = a(1,n)*a(1,n-2)

      a(1,n)   = a(2,n)-a(6,n-2)*a(2,n-2)
      a(2,n)   = a(3,n)-a(6,n-2)*a(3,n-2)
      a(3,n)   = a(4,n)-a(6,n-2)*a(4,n-2)
      a(4,n)   = a(5,n)
      a(5,n)   = a(6,n)

!     k=n-1

      a(1,n-1) = 1d0/a(1,n-1)
      a(5,n-1) = a(1,n)*a(1,n-1)

      a(1,n)   = a(2,n)-a(5,n-1)*a(2,n-1)
      a(2,n)   = a(3,n)-a(5,n-1)*a(3,n-1)
      a(3,n)   = a(4,n)-a(5,n-1)*a(4,n-1)
      a(4,n)   = a(5,n)

!     the next loop will be used in banbk

      a(1,n)=1d0/a(1,n)

      return
      END



!**************************************************
!*    SUBRUTINA PARA RESOLVER SIST LINEALES
!*    (Numerical Recipes in FORTRAN)
!**************************************************
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
        pivinv=1d0/a(icol,icol)
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

