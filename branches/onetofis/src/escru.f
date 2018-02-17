
!--------------------------------------------------------------------!
! Subroutine WRITESTATS
! AAF June 2014
! AAF February 2016 for onetofis-getstats
!....................................................................!
      subroutine writestats(myid)
      
      use statis
      use point
      use fis
      use ficheros
      use MPI_GROUPS 
      use tem
      use spectra

      implicit none 
      include "mpif.h"
      include "ctes3D"

!-------------Variables----------------------------!
      integer myid,i,j,ii,jj,k,ierr,istat(MPI_STATUS_SIZE)
      real*4 fac
      character*3 ext1
      real*4 Ree,alpe,bete,a0e
      integer istep,leng,kk,iproc
    
      !extension number of sta file & increment on it
      !write(ext1,'(i3.3)') ista 
      Ree  = Re
      alpe = alp
      bete = bet
      a0e  = a0
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !Compute stats
      !Note:MEAN VALUES do not need to REDUCE because only MASTER PROC 
      ! accumulates the data
      call MPI_ALLREDUCE(ff,wkstats(1,1),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ff2,wkstats(1,2),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ff3,wkstats(1,3),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ff4,wkstats(1,4),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)

      !WRITE STATS
      if (myid.eq.0) then
         !variables accumulated onf phys space need to be averaged
         !by fac
         fac = 1./dble(mgalx*mgalz)
         !Calculate stats vectors for this field
         !mean:
         do j=1,my
              mean(j) = wkstats(j,1)*fac
         enddo 
         !variance
         do j=1,my
              stddev(j) = wkstats(j,2)*fac-mean(j)**2
         enddo
         !skewness
         do j=1,my
              if (stddev(j).eq.0.0) then
                skew(j) = stddev(j)
              else
                skew(j) = (wkstats(j,3)*fac-3*mean(j)*stddev(j)-
     .      mean(j)**3)/stddev(j)**(1.5)               
              endif
         enddo
         !kurt
         do j=1,my
              if (stddev(j).eq.0.0) then
                kurt(j) = stddev(j)
              else
                kurt(j) = (wkstats(j,4)*fac-4*mean(j)*wkstats(j,3)*fac
     .    +6*mean(j)**2*stddev(j)+3*mean(j)**4)/stddev(j)**2
              endif
         enddo

         write(*,*) "Writing stats at  time = ",time 
         open (isn,file=fnamesta,status='unknown',form='unformatted')
         rewind(isn)
         nacum = mgalx*mgalz
         write(isn) my,time,Ree,alpe,bete,a0e,nacum
         write(isn) (y(j),fmap(j),j=1,my)
         write(isn) (mean(j),stddev(j),skew(j),kurt(j),
     .              wkstats(j,1),wkstats(j, 2),
     .              wkstats(j,3),wkstats(j, 4),
     .                                   j=1,my)
         close(isn)
      endif

    
      end subroutine
!-------------------------------------------------------------------!

!-----------------------------------------------------------------------!
      subroutine resetstats()

      use statis
      use point
      implicit none
      include "ctes3D"
      include "mpif.h"
    
      integer i,j,k,jj,kk,ii,myid

      nacum = 0
      if (myid.eq.0) then
         write(*,*) " Reset stats..."
         mean(j) = 0d0
         stddev(j) = 0d0
         skew(j) = 0d0
         kurt(j) = 0d0 !rhoum
      endif
      !clean mean buffers
      do j=1,my
         ff(j) = 0d0 
         ff2(j) = 0d0
         ff3(j) = 0d0
         ff4(j) = 0d0
      enddo
      !cleaning wkstats (just in case)
      do j=1,my 
         wkstats(j,1) = 0.0
         wkstats(j,2) = 0.0
         wkstats(j,3) = 0.0 
         wkstats(j,4) = 0.0 
         wkstats(j,5) = 0.0
         wkstats(j,6) = 0.0
         wkstats(j,7) = 0.0
         wkstats(j,8) = 0.0
         wkstats(j,9) = 0.0 
         wkstats(j,10)= 0.0 
      enddo


      end subroutine
!-------------------------------------------------------------------!
 


!=======================================================!
!AAF subroutine to write any size variable
! to a file for matlab
      subroutine write2file(fwk,varname,tamy,tamz,tamx,myid)
      use ficheros
      use fis

      implicit none
      include "mpif.h"
      include "ctes3D"

!---------Variables------------------------------------!

      integer fid 
      integer tamy,tamz,tamx,myid
      real(4) :: fwk(2,1:tamy,1:tamz,1:tamx)
      character*100 fnameima,fnamespe
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

      include "mpif.h"
      include "ctes3D"
!-------------------------------------------------------!      

      integer fid 
      integer tamy,tamz,tamx,myid
      real(8) :: fwk(1:tamy,1:tamz,1:tamx)
      character*100 fnameima,fnamespe
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
      close(fid)

      end
 
