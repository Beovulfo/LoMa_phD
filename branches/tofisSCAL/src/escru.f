!***********************************************************************!
!                                                                       !
!               writes a intermediate solution                          !
!    single jjs  4/01/01                                                !
!    plane  lines shc 01/03/05                                          !
!                                                                       !
!                                                                       !
!                                                                       !
!***********************************************************************!

      subroutine escru1phys(v,wk1,wk1r,j,fj)
      ! Esta subrutina escribe en disco un plano yz, un yx y los
      ! marcados en ctes3d para j. 

      implicit none
      include 'ctes3D'
      
      ! -------------- Common ---------!

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      integer iinp,iswap
      character*80 filinp,filout,filswap
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
            
      complex*8,dimension(0:mx1,0:mz1)::v,wk1
!      real*4    u00(my),w00(my),du00(my),dw00(my)
      real*4 wk1r(mgalx+2 ,mgalz)
      integer i,k,kk,nacum,jj
      complex*8 xa,xb
      character*80 fname
      character*5  extp(3)
      ! ------------------

      kk=0; nacum=1
      plyx=48  !plyx RP
      plyz=1  !plyz we anto to write

!      write(*,*)'entrando',j
      extp=(/'.scyz','.scyx','.scxz'/)

      ! Velocidades       

      if(j.eq.1) then
         do k=1,3
            fname= filout(1:index(filout,' ')-1)//extp(k)
            open (1000+k,file=fname,status='unknown',form='unformatted')
            write(1000+k) time,Re,alp,bet,mgalx,my,mgalz,nspec,nacum 
           write(1000+k)(jspecy(jj),jj=1,nspec),(y(jj),fmap(jj),jj=1,my)
         enddo
      endif
 

      ! v

      do k=0,mz1
         do i=0,mx1
            wk1(i,k) = v(i,k)
         enddo
      enddo

      !wk1r(1,1) = 0d0
      !wk1r(2,1) = 0d0

      call fourxz(wk1,wk1,1,1)
      kk=kk+1
      write(1000+kk) (wk1r(plyz,k),k=1,mgalz)
      kk=kk+1
      write(1000+kk) (wk1r(i,plyx),i=1,mgalx)
      kk=kk+1
      if (fj==1) write(1000+kk)((wk1r(i,k),i=1,mgalx),k=1,mgalz)



      end  

