!***********************************************************************!
!                                                                       !
!               writes a intermediate solution                          !
!    single jjs  4/01/01                                                !
!    plane  lines shc 01/03/05                                          !
!                                                                       !
!                                                                       !
!                                                                       !
!***********************************************************************!

      subroutine escrphys(phi,v,vor,dvdy,dvordy,wk,wk1,wk1r,
     %                    u00,w00,du00,dw00,j,fj)
      ! Esta subrutina escribe en disco un plano yz, un yx y los
      ! marcados en ctes3d para j. 

      implicit none
      include 'ctes3D'
      
      ! -------------- Common ---------!

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      integer iinp,iswap
      character*50 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(7)
      save /ficheros/

      real*4  Re,alp,bet,a0
      real*8  y,fmap
      common /fis/ Re,alp,bet,a0,y(my),fmap(my)
      save   /fis/

      integer icx
      real*4 alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     &              icx(0:mz1)
      save /wave/      
!     --------------------------------------------------------
      integer,intent(in)::j,fj
      integer plyx,plyz !planes for writing
            
      complex*8,dimension(0:mx1,0:mz1)::phi,v,vor,dvordy,dvdy,wk1
!      real*4    u00(my),w00(my),du00(my),dw00(my)
      real*8    u00(my),w00(my),du00(my),dw00(my)
      real*4 wk1r(mgalx+2 ,mgalz),wk(0:mx1,0:mz1)
      integer i,k,kk,nacum,jj
      complex*8 xa,xb
      character*80 fname
      character*5  extp(21)
      ! ------------------

      kk=0; nacum=1
      plyx=48  !plyx RP
      plyz=1  !plyz we anto to write

!      write(*,*)'entrando',j
      extp=(/'.upyz','.upyx','.upxz','.vpyz','.vpyx','.vpxz',
     &       '.wpyz','.wpyx','.wpxz','.oxyz','.oxyx','.oxxz',
     &       '.oyyz','.oyyx','.oyxz','.ozyz','.ozyx','.ozxz',
     &       '.oz2d','.oz3d','.z3xz'/)

      ! Velocidades       

      if(j.eq.1) then
         do k=1,21
            fname= filout(1:index(filout,' ')-1)//extp(k)
            open (1000+k,file=fname,status='unknown',form='unformatted')
            write(1000+k) time,Re,alp,bet,mgalx,my,mgalz,nspec,nacum 
           write(1000+k)(jspecy(jj),jj=1,nspec),(y(jj),fmap(jj),jj=1,my)
         enddo
      endif
!      write(*,*)'despues',j
 
      ! u

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
!            if(j==my) write(*,*)i,k,wk1(i,k) 
            xa = xalp(i)
            wk1(i,k) =  (dvdy(i,k)*xa-vor(i,k)*xb)*wk(i,k)
! calculate it on FOU and...
         enddo
      enddo
!      WRITE(*,*)'primer',j

      wk1r(1,1) = u00(j)
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
! ... change it to phys
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)
!      if (fj==1) write(1000+kk)((dvdy(i,k),i=0,mx1),k=0,mz1)


      ! v

      do k=0,mz1
         do i=0,mx1
            wk1(i,k) = v(i,k)
         enddo
      enddo

      wk1r(1,1) = 0d0
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)

      ! w

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            wk1(i,k) = (dvdy(i,k)*xb+vor(i,k)*xa)*wk(i,k)
         enddo
      enddo
!      WRITE(*,*)'segun',j

      wk1r(1,1) = w00(j)
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)
!      if (fj==1) write(1000+kk)((vor(i,k),i=0,mx1),k=0,mz1)


      ! ox
      
      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            wk1(i,k) = (dvordy(i,k)*xa+phi(i,k)*xb)*wk(i,k)
         enddo
      enddo

!IMPORTANT: if we want to plot perturbations set to zero all
!          zero modes
      wk1r(1,1) = dw00(j)
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)      
!      if (fj==1) write(1000+kk)((dvordy(i,k),i=0,mx1),k=0,mz1)      

      ! oy

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            wk1(i,k) = vor(i,k)
         enddo
      enddo
      wk1r(1,1) = 0d0
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)      



      ! oz

      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            wk1(i,k) = (dvordy(i,k)*xb-phi(i,k)*xa)*wk(i,k)
         enddo
      enddo
      wk1r(1,1) = -du00(j)
      wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
!special plane for oz
      plyx=1  !plyx BP 96 is z=0 approx
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)     
 

! ----------------------------------------------------------------!
!Adding calculation of wz2d and wz3d(aaf)
!    oz2d & oz3d
!-------------------------!
!   oz2d (only with kz=0)
      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            if (k==0) then 
              wk1(i,k) = (dvordy(i,k)*xb-phi(i,k)*xa)*wk(i,k)
            else
              wk1(i,k) = 0.
            endif
         enddo
      enddo
      wk1r(1,1) = -du00(j)
      wk1r(2,1) = 0d0
      call fourxz(wk1,wk1,1,1)
      kk=kk+1
!special plane for oz
      plyx=1  !plyx BP 96 is z=0 approx
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)

!oz3D  (kz not zero),kz=0 zero
      do k=0,mz1
         xb = xbet(k)
         do i=0,mx1
            xa = xalp(i)
            if (k==0) then
              wk1(i,k) = 0.
            else
              wk1(i,k) = (dvordy(i,k)*xb-phi(i,k)*xa)*wk(i,k)
            endif
         enddo
      enddo
      wk1r(1,1) = 0d0
      wk1r(2,1) = 0d0
      call fourxz(wk1,wk1,1,1)
      kk=kk+1
!special plane for oz in plane YX
      plyx=1  !plyx BP 96 is z=0 approx
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
!special plane for oz in plane XZ
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)     


      end  

