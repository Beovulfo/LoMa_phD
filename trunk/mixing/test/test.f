      program main

      implicit none

      include "ctes3D"

      real*8     gama(3), alpha(4), beta(3), ibeta(3), xi(3)
      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
      parameter (alpha=(/ 29d0/96d0, -3d0/40d0, 1d0/6d0, 29d0/96d0/))
      parameter (beta =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /))
      parameter (ibeta=(/ 160d0/37d0, 24d0/5d0,   6d0     /))
      parameter (xi   =(/ -17d0/60d0, -5d0/12d0,  0d0     /))


!--------------COMMONS--------------------------------------------------!      
      real*8 prem3,dt21,dt22
      common /d2y/ prem3(5,my),dt21(5,my),dt22(5,my)
      save   /d2y/

      real*8 prem1,dt12
      common /cfdiff/ prem1(7,my),dt12(7,my)
      save /cfdiff/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),trp(my),mss(my)
      save   /fis/


      integer iax,icx
      real*4 alp2,bet2,ralp,rbet
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .              ralp(0:mx1),rbet(0:mz1),iax(mx),icx(0:mz1)
      save  /wave/

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save   /tem/
!------------end commons-----------------------------------------!
!-------------VARIABLES------------------------------------------!
     
      real*8  u00(my),u00wk(my),rf0u(my),rhs(my),rhs2(my),nl(my) 
      real(4),dimension(2,my) :: vorwk,  vor,  phi,   v, dvdy,phiwk
      real*8  const1,const2 
      real*8  rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire
      real*4 pi,zero
      integer myid,istep,irun,rkstep,i,k,j,k1,ierr,nanerror,kk

      real*8, allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
      real*8  two
      real*4  u(2,my,mz),du(2,my,mz)
!---------------------------------------------------------------------!
      Pi=4.*atan(1.)
      zero = 0.0
! CREATE MESH manually, only vector y and vector size my needed
      do i=1,my
         y(i)=-1.+(i-1)*2./(my-1)
      enddo
      do i=1,my
         y(i)=1.02*tanh(y(i)*atanh(1/1.02))
      enddo 
      y(1)=-1
      y(my)=1
      write(*,*) "my", my
      
      allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))
      
      call derivadas(d11,d12,d21,d22)      
            
      deallocate(d11,d12,d21,d22)
! write fmap
!      write(*,*) '-----FMAP---- '      
!      do i=1,my
!         write(*,*) fmap(i)
!      enddo
!      write(*,*) '-----FMAP---- '      

!------------------------------------------------------------------!
!                 TESTING FOR RK00
!------------------------------------------------------------------!
      Re=1;
      Deltat=0.0001;
      ire=1/Re;
!tan     const1=10d0
!      const2=atan(const1)
!      do j=1,my
!         u00(j) = 1./const2*atan(const1*y(j)) 
!      enddo
      do j=1,my
         u00(j) = 5./2.*(y(j)**3.-3./5.*y(j)**5)
      enddo
!      write(*,*) '-----y---- '      
!      do j=1,my
!         write(*,*) y(j)
!      enddo
!      write(*,*) '-----y---- '      

!      write(*,*) '-----u00---- '      
!      do j=1,my
!         write(*,*) u00(j)
!      enddo
!      write(*,*) '-----u00---- '      

!      rf0u=0;
! Some u00 has to be created in order to estimat LHS,
! Then u00 in the last step should be also created in order
! to obtain NL (rf0u) analitically.
! With this we can know RHS so that the u00 in the next step
! is the one we imposed at the beginning
! CAUTION: Check BC for given u00 !!
      rkstep=2; 
      write(*,*) "rkstep", rkstep
         rkn1    = Re*ibeta(rkstep)/Deltat 
         dalre   = Deltat*alpha(rkstep)*ire
         dtgamma = Deltat*gama(rkstep)
         dalbe   = 1d0+alpha(rkstep+1)*ibeta(rkstep)
         dtri    = Deltat*alpha(rkstep+1)*ire
         dtxi    = Deltat*xi(rkstep)
   
       write(*,*) 'rkn1',rkn1,'dalre',dalre,'dtgamma',dtgamma,
     .            'dalbe', dalbe, 'dtri',dtri,'dtxi',dtxi
!for rkstep.ne.1 
      do j=1,my
!rf0u no linear term
         rf0u(j)=100.*sin(Pi*y(j))
!u00wk not used as input in rkstep.ne.1
         u00wk(j) = u00(j) 
!u00'' first
         rhs(j)=15./2.*(2.*y(j)-4.*y(j)**3)
!         rhs(j)=const1/const2*(-2.*const1*const1*y(j))/
!     .   (1+const1*const1*y(j)*y(j))
!now build RHS
         rhs(j)=rhs(j)-rkn1*u00(j)
!finally obtain rhs
         rhs(j)=-rhs(j)/rkn1-dtgamma*rf0u(j)

      enddo


!      write(*,*) "before calling rk00"
!      do j=1,my
!        write(*,*) u00wk(j)
!      enddo
!       write(*,*) "rhs(1),rhs(2)",rhs(1),rhs(2)
                       
!        call rk00(u00wk,rhs,rf0u,rkn1,
!     &       dalre,dtgamma,dalbe,dtri,dtxi,rkstep,1)
!u00wk should be the same as u00
!      write(*,*) "after calling rk00"
!      write(*,*) "difference between u00wk-u00"
!      do j=1,my
!        write(*,*) u00wk(j)-u00(j)
!      enddo
!       write(*,*) "rhs(1),rhs(2)",rhs(1),rhs(2)
!----------------------------------------------------
!-----------VISC (COUETTE BC)
      do i=1,2
      do j=1,my
!forcing for omega mode 1,1
          vorwk(i,j) = 2-(2+rkn1)*(y(j)**2-1)
! vor is just an output...
          vor(i,j)   = zero
! phi is just an output...
          phi(i,j)   = zero
!NL for omega
          v(i,j) =   sin(2*Pi*y(j))
!NL for vor
          dvdy(i,j)= sin(2*Pi*y(j))
!forcing for phi
          phiwk(i,j) = -24.*y(j)**2.+32.-(2.+rkn1)*
     .                 (-2*y(j)**4.+16.*y(j)**2.-6.)
          do k=1,mz
             u(i,j,k)=sin(2*Pi*y(j))
          enddo
        enddo
      enddo
      

      two = 2d0;

!      call deryr2(u,du,mz)
!test for first derivative
!      write(*,*) 'du from deryr2'
!      write(*,*) du(1,:,1)

!      write(*,*) 'analytical first derivative'
!      write(*,*) 2*Pi*cos(2*Pi*y(:))-du(1,:,1)
! FIRST DERIVATIVE WITH DERYR2 WORKS!!
! but trying to do second derivative using 2 times deryr2 fails in both
! extremes :(
!      call deryr2(du,du,mz)
!now u keeps the second derivative of u!!
!      call deryyr2(u,u,u,mz,0d0,0d0,0d0,0d0)
!      write(*,*) 'checking deryyr2'
!      write(*,*) du(1,:,1)
!      do j=1,my
!      write(*,*) du(1,j,1)-(-4*Pi**2*sin(2*Pi*y(j)))
!      enddo
!deryyr2 gives really good results calculating second derivative
      
!---------------------------------------------------
!-----------VISC (ROBIN mixing layer BC's)



!--------------------------------------------------
!      call visc(vorwk(1,1),vor(1,1), phi(1,1),v(1,1),dvdy(1,1),
!     .     phiwk(1,1),two+rkn1,two,dalbe,dtri,dtxi,1)                     
!----------------------------------------------------------------!
!      write(*,*) "vor after visc compared with analitical"
!      do j=1,my
!the expected result
!        v(1,j) = y(j)**2.-1. 
!        write(*,*) vorwk(1,j)-v(1,j)
!      enddo
!------------------------------------------------------------
!      write(*,*) "velocity v compared with analitical"
!      do j=1,my
!!the expected result
!        v(2,j) = (y(j)**2.-1.)**2.
!        write(*,*) v(1,j), v(2,j)
!      enddo

!SHOW HY
      write(*,*) hy(:)


      end program

     
