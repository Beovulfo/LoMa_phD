! ----------------------------------------------------------------------!
!                                                                       !
!    Hace las cuentas en y para planos kz=cte                           !
!                                                                       !
! ----------------------------------------------------------------------!

      subroutine uvwyzx(vor,phi,u,v,w,dudy,dvdy,nvecy,i)

      implicit none
      include "ctes3D"
      
      ! -------------------- Commons ------------------------------------!

      integer      icx
      real*4       alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     &              icx(0:mz1)
      save  /wave/      
      
      ! ----------------------- constants -------------------------------!
      integer,intent(in)::i
     
      complex*8,intent(in),dimension(my,0:mz1,0:nvecy-1):: vor,phi
      complex*8,intent(out),dimension(my,0:mz1,0:nvecy-1)::
     &          u,v,w,dudy,dvdy

        
      !------------------------Aux arrays    ----------------------------!
      integer k,j,k1,ind,ii,ip,nvecy
      complex*8 hgg (my),hvv (my),xa,xb,dyvor(my),o1(my),o3(my)
      real*8 rK

      ! ------------------------ ! Computes u,v,w ! -- ------------------!
       
      do ind=0,nvecy-1
         ip = i+ind
         do k=0,mz1
            k1 = icx(k)
            rK = bet2(k1)+alp2(ip-1)

            !call Lapvdv(phi(1,k,ind),hgg,hvv,rK)
            call Lapvdvhom(phi(1,k,ind),hgg,hvv,rK)
            call deryr2(vor(1,k,ind),dyvor,1)

            do j=1,my
               u   (j,k,ind) = phi(j,k,ind)
               v   (j,k,ind) = hgg(j)
               dvdy(j,k,ind) = hvv(j)
               w   (j,k,ind) = vor(j,k,ind)
               dudy(j,k,ind) = dyvor(j)
            enddo
         enddo
      enddo

      endsubroutine uvwyzx
      
      
!!   !***********************************************************************!
!!   !                                                                       !
!!   !	OPERACIONES EN LINEAS EN X                                           !
!!   !                                                                       !
!!   !***********************************************************************!
!!   
!!         subroutine presrhs (phi,v,vor,dvordy,dvdy,pr,ps,wk,
!!        %          wk1,wk2,wk3,wk1r,wk2r,wk3r,u00,w00,du00,dw00,j,ind)    
!!         implicit none
!!         include "ctes3D"
!!   
!!         ! -------------------- Commons ------------------------------------!
!!   
!!         integer icx
!!         real*4 alp2,bet2
!!         complex*8    xalp, xbet
!!         common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
!!        &              icx(0:mz1)
!!         save /wave/      
!!   
!!         ! ------------------------ I&O ----------------------------------!
!!   
!!         integer,intent(in)::j
!!         integer ind
!!         
!!         complex*8 phi (0:mx1,0:mz1),     v(0:mx1,0:mz1),
!!        %          vor (0:mx1,0:mz1),dvordy(0:mx1,0:mz1),
!!        %          dvdy(0:mx1,0:mz1),
!!        &          pr  (0:mx1,0:mz1),ps    (0:mx1,0:mz1) 
!!   
!!         complex*8 wk1 (0:mx1,0:mz1),wk2 (0:mx1,0:mz1),wk3 (0:mx1,0:mz1)
!!   
!!         real*4    wk1r(mgalx+2 ,mgalz),wk2r(mgalx+2 ,mgalz),
!!        %          wk3r(mgalx+2 ,mgalz),wk(0:mx1,0:mz1)
!!   
!!         real*4    u00(my),w00(my),du00(my),dw00(my)
!!   
!!         ! ----------------------- constants -------------------------------!
!!         
!!         integer i,k,kk,nlen,buc,fj
!!         complex*8 xa,xb,zero
!!         real*4    temp
!!   
!!         !------------------------------------------------------------------!
!!         
!!         zero=cmplx(0,0)
!!         
!!         do k=0,mz1
!!            do i=0,mx1
!!               wk1(i,k) = xalp(i)*v(i,k)
!!            enddo
!!         enddo
!!   
!!         do k=0,mz1
!!            do i=0,mx1
!!               wk2(i,k) = xbet(k)*v(i,k)
!!            enddo
!!         enddo
!!   
!!         ! 2. Transformo a Fourier
!!      
!!         call fourxz(wk1,wk1,1,1)    ! BUG corregido
!!         call fourxz(wk2,wk2,1,1)    !
!!      
!!         ! 3. RHS OF SLOW PRESS -------------! En Fourier 2 y 8.
!!         ! Libre 3
!!      
!!         ! 3.1 , dudy and v
!!          
!!         do k=0,mz1
!!            xb = xbet(k)
!!            do i=0,mx1
!!               xa = xalp(i)
!!               temp = wk(i,k)
!!               wk3(i,k) = xa*v(i,k)-(dvordy(i,k)*xb -phi(i,k)*xa)*temp  ! 4
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk3,wk3,1,1)
!!      
!!         !Calculo primer sumando y sobreescribo en 1
!!      
!!         do k=1,mgalz
!!            do i=1,mgalx
!!               wk1r(i,k) = 2*wk1r(i,k)*wk3r(i,k)   ! 2 and 4
!!            enddo
!!         enddo
!!      
!!         ! Calculo 6
!!      
!!         do k=0,mz1
!!            xb = xbet(k)
!!            do i=0,mx1
!!               xa = xalp(i)
!!               temp= wk(i,k)
!!               wk3(i,k) = xb*v(i,k)+(dvordy(i,k)*xa+phi(i,k)*xb)*temp  ! 4
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk3,wk3,1,1)
!!      
!!         do k=1,mgalz
!!            do i=1,mgalx
!!               wk1r(i,k) = wk1r(i,k) + 2*wk2r(i,k)*wk3r(i,k)  ! 6 and 8
!!            enddo
!!         enddo
!!      
!!         ! Libres wk2 y wk3 ! calculo 3 y 7
!!         ! dwdx dzdu
!!      
!!         do k=0,mz1
!!            xb = xbet(k)
!!            do i=0,mx1
!!               xa = xalp(i)
!!               temp = wk(i,k)
!!               wk2(i,k) = xa*temp*(dvdy(i,k)*xb+vor(i,k)*xa)
!!               wk3(i,k) = xb*temp*(dvdy(i,k)*xa-vor(i,k)*xb)
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk2,wk2,1,1)
!!         call fourxz(wk3,wk3,1,1)
!!      
!!         do k=1,mgalz
!!            do i=1,mgalx
!!                wk1r(i,k) = wk1r(i,k) + 2*wk2r(i,k)*wk3r(i,k)  ! 3 and 7
!!            enddo
!!         enddo
!!      
!!         ! 5 and 9
!!      
!!         do k=0,mz1
!!            do i=0,mx1
!!               wk2(i,k) = dvdy(i,k) !uyzx(nplany,0:mz1,mx1+1)
!!            enddo
!!         enddo
!!      
!!         do k=0,mz1
!!            xb = xbet(k)
!!            do i=0,mx1
!!               xa = xalp(i)
!!               temp = wk(i,k)
!!               wk3(i,k) = xb*temp*(dvdy(i,k)*xb+vor(i,k)*xa)
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk2,wk2,1,1)
!!         call fourxz(wk3,wk3,1,1)
!!      
!!         do k=1,mgalz
!!            do i=1,mgalx
!!               wk1r(i,k) = wk1r(i,k) + wk2r(i,k)*wk2r(i,k)
!!        &                            + wk3r(i,k)*wk3r(i,k)! 5 and 9
!!            enddo
!!         enddo
!!      
!!         ! 'ultimo, 1
!!      
!!         do k=0,mz1
!!            xb = xbet(k)
!!            do i=0,mx1
!!               xa = xalp(i)
!!               temp= wk(i,k) 
!!               wk2(i,k) = xa*temp*(dvdy(i,k)*xa-vor(i,k)*xb)
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk2,wk2,1,1)
!!      
!!         do k=1,mgalz
!!            do i=1,mgalx
!!                 wk1r(i,k) =-(wk1r(i,k) + wk2r(i,k)*wk2r(i,k)) ! 1
!!            enddo
!!         enddo
!!      
!!         call fourxz(wk1,wk1,-1,1)
!!      
!!         do k=0,mz1
!!            do i=1,mx1
!!               ps(i,k) = wk1(i,k)  
!!            enddo
!!         enddo
!!   
!!   
!!       !!!  ! Uncomment the following routine if you want velocities or 
!!       !!!  ! Vorticities in physics space
!!       !!!  ! 
!!   
!!   
!!       !!!  ! fvor(1,1,1) -> wx,wy,wz
!!       !!!  fj=0
!!       !!!  if (ind.le.nspec) then
!!       !!!     write(*,*)'chequeando',ind,j,jspecy(ind)
!!       !!!     if (jspecy(ind)==j) then
!!       !!!        fj=1
!!       !!!        ind=ind+1
!!       !!!     endif
!!       !!!  endif
!!   
!!       !!!  
!!       !!!  write(*,*) 'entrando escrphys',j,u00(j)
!!       !!!  call escrphys(phi,v,vor,dvdy,dvordy,wk,wk1,wk1,
!!       !!! &              u00,w00,du00,dw00,j,fj)
!!       !!!  write(*,*) 'saliendo de escrphys',j,u00(j)
!!         
!!         endsubroutine presrhs
!!         
!!   
!!   !***********************************************************************!
!!   !                                                                       !
!!   !	 PRESION: calculo presion y balance                                  !
!!   !                                                                       !
!!   !***********************************************************************!
!!         
!!         subroutine pressure(rapid,slow,stok,bcpress,nvecy,i)
!!   
!!         implicit none
!!         include "ctes3D"
!!         
!!         ! -------------------- Commons ------------------------------------!
!!         
!!         integer icx
!!         real*4 alp2,bet2
!!         complex*8    xalp, xbet
!!         common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
!!        &              icx(0:mz1)
!!         save /wave/     
!!    
!!         ! ----------------------- constants -------------------------------!
!!             
!!         integer i,k,j,kk,nlen,k1,iplan,jj,nvecy,buc
!!   
!!         complex*8 slow (0:my-1,0:mz1,0:nvecy-1),
!!        &          rapid(0:my-1,0:mz1,0:nvecy-1),
!!        &          stok (0:my-1,0:mz1,0:nvecy-1)
!!   
!!         complex*8 bcpress(2,0:mz1,mx)
!!   
!!         complex*8 bct,bcb,zero
!!         real*8 rk
!!        
!!         !------------------------------------------------------------------!
!!         
!!         zero=cmplx(0.,0.)
!!         
!!         do jj = 0,nvecy-1
!!            iplan = i+jj
!!            if (iplan==1) then ! from k=1
!!            
!!               do j=0,my1
!!                  slow (j,0,jj) = zero
!!                  rapid(j,0,jj) = zero
!!                  stok (j,0,jj) = zero
!!               enddo
!!               do k=1,mz1
!!                  k1 = icx(k)
!!                  rK = bet2(k1)+alp2(iplan-1)
!!   
!!                  bcb = bcpress(1,k,iplan)
!!                  bct = bcpress(2,k,iplan)
!!                  call visc(stok (0,k,jj),bcb,bct,rK,1)
!!                  bcb = zero
!!                  bct = zero
!!                  call visc(slow (0,k,jj),bcb,bct,rK,1)
!!               enddo
!!            else
!!               do k=0,mz1
!!                  k1 = icx(k)
!!                  rK = bet2(k1)+alp2(iplan-1)
!!                  bcb = bcpress(1,k,iplan)
!!                  bct = bcpress(2,k,iplan)             
!!                  call visc(stok (0,k,jj),bcb,bct,rK,1)
!!                  bcb = zero
!!                  bct = zero 
!!                  call visc(slow (0,k,jj),bcb,bct,rK,1)
!!   
!!               enddo
!!            endif      
!!         enddo      
!!   
!!            
!!         endsubroutine pressure
