!******************************************************************!
!                                                                  !
!               writes a intermediate solution                     !
!    single jjs  4/01/01                                           !
!    plane  lines shc 01/03/05                                     !
!                                                                  !
!                                                                  !
!                                                                  !
!******************************************************************!
      !subroutine escru(vor,phi,u00,w00,v00,myid)
      subroutine escru(vor,phi,psi,scal,u00,w00,v00,myid)
      use tem
      use fis
      use ficheros
      use point

      implicit none

      include "ctes3D"

c  ------------------------- Variables ----------------------------
     
      integer iproc,i,j,k,leng,jj,kk,ii
      integer myid
      
      real(4) :: vor(0:(2*my-1),0:mz-1,1:nplanes),
     .        phi(0:(2*my-1),0:mz-1,1:nplanes),
     .        psi(0:(2*my-1),0:mz-1,1:nplanes),
     .        scal(0:(2*my-1),0:mz-1,1:nplanes)

      real(8) ::  u00(my),w00(my),v00(my)
      real(8) ::  fac,wk(my,13)      
      real(4) ::  wk00(3*my)
      real(4)     Ree,alpe,bete,a0e

      character*3 ext1
      character*100 fnameima,fnamespe,fnamesta
     
      Ree  = Re
      alpe = alp
      bete = bet
      a0e  = a0
c    ------------------------  Program      -----------------

      leng = 2*my*mz 
           
      fnameima=filout(1:index(filout,' ')-1)
         
      open (iout,file=fnameima,status='unknown',
     .      form='unformatted')

      rewind(iout)
         
      do j=1,my
         wk00(3*j-2) = u00(j)
         wk00(3*j-1) = v00(j)
          wk00(3*j)   = w00(j)
      enddo

      write(iout) time,Ree,alpe,bete,a0e,mx,my,mz,
     .               (y(j),fmap(j),j=1,my),
     .               (wk00(j),j=1,3*my)
 
         
      do i=1,nplanes
            !write(iout) ((vor(j,k,i),phi(j,k,i),j=0,2*my-1),k=0,mz1)
            write(iout) ((vor(j,k,i),phi(j,k,i),psi(j,k,i),
     .                   scal(j,k,i),j=0,2*my-1),k=0,mz1)
      enddo
         
      close(iout)
         
      
    
      end subroutine escru

!=======================================================!
!AAF subroutine to write any size variable
! to a file for matlab
      subroutine write2file(fwk,varname,tamy,tamz,tamx,myid)
      use ficheros
      use fis

      implicit none
      include "ctes3D"

!---------Variables------------------------------------!

      integer fid 
      integer tamy,tamz,tamx,myid
      real(4) :: fwk(2,1:tamy,1:tamz,1:tamx)
      character*100 fnameima,fnamespe,fnamesta
      character*3 varname
      character*3 ext1
      integer i,j,k

!.......Program.......................................!
      
      write(ext1,'(i3.3)') id22 
      fnamesta=filstt(1:index(filstt,' ')-1)//'_'//
     .                               varname//'.'//ext1
      open (fid,file=fnamesta,status='unknown',form='unformatted')
      rewind(fid)
      write(fid) (y(j),fmap(j),j=1,my)
      write(fid) tamy,tamz,tamx
!write real part
      write(fid)(((fwk(1,j,k,i),j=1,tamy),k=1,tamz),i=1,tamx)
!write imaginary part
      write(fid)(((fwk(2,j,k,i),j=1,tamy),k=1,tamz),i=1,tamx)
      close(fid)

      end
 


!=======================================================!
!AAF subroutine to write any size variable
! to a file for matlab
      subroutine write2file1(fwk,varname,tamy,tamz,tamx,myid)
      use ficheros
      use fis

      implicit none

      include "ctes3D"
!-------------------------------------------------------!      

      integer fid 
      integer tamy,tamz,tamx,myid
      real(8) :: fwk(1:tamy,1:tamz,1:tamx)
      character*100 fnameima,fnamespe,fnamesta
      character*3 varname
      character*3 ext1
      integer i,j,k
!-------------------------------------------------------!

      write(ext1,'(i3.3)') id22 
      fnamesta=filstt(1:index(filstt,' ')-1)//'_'//
     .                              varname//'.'//ext1
      open (fid,file=fnamesta,status='unknown',form='unformatted')
      rewind(fid)
      write(fid) (y(j),fmap(j),j=1,my)
      write(fid) tamy,tamz,tamx
!write real part
      write(fid)(((fwk(j,k,i),j=1,tamy),k=1,tamz),i=1,tamx)
!write imaginary part
      close(fid)

      end
 



