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
!     ----------------------- Commons  ---------------------------------

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
      use matrices
!     Input  :: f
!     Output :: rhslap

      implicit none
      include "ctes3D"

      integer i,j
      real(8), dimension(my):: f,rhsd1
        

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
      use fis
      use matrices
      use wave
      use worksp

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


      zero = 0d0

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
      use matrices
      use fis
      use worksp
      use wave
 
      implicit none
      include "ctes3D" 
!    ----------------------- IN & OUT ----------------- 
    
      real(8),dimension(my),intent(inout) :: f, u
      real(8),    intent(in)              :: rk,dalbe,dtri
      integer(4), intent(in)              :: ind


!--------- ----------- BC functions----------------------------
      real(8) :: C_phi_top1, C_vvv_top1,C_dyv_top1,BC_phi_top1
      real(8) :: C_phi_top2, C_vvv_top2,C_dyv_top2,BC_phi_top2
      real(8) :: C_phi_bot1, C_vvv_bot1,C_dyv_bot1,BC_phi_bot1
      real(8) :: C_phi_bot2, C_vvv_bot2,C_dyv_bot2,BC_phi_bot2
      real(8) :: C_ome_top,  C_dyo_top,BC_ome_top
      real(8) :: C_ome_bot,  C_dyo_bot,BC_ome_bot
!     ----------------------- work --------------------     

      integer(4)              ::  i,j
      real(8),dimension(my)   ::  df,fwkh1,fwkh2,dfwk,dfwkh1,dfwkh2
      real(8)                 ::  coef(2),AA(2,2)

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
       
 
