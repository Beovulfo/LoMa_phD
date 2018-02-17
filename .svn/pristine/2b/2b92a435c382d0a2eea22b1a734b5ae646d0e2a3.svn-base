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
      !use MPI_GROUPS 

      implicit none

      include "mpif.h"
      include "ctes3D"

c  ------------------------- Variables ----------------------------
     
      integer iproc,i,j,k,leng,jj,kk,ii
      integer myid,istat(MPI_STATUS_SIZE),ierr
      
      real(4) :: vor(0:(2*my-1),0:mz-1,pb:pe),
     .        phi(0:(2*my-1),0:mz-1,pb:pe),
     .        psi(0:(2*my-1),0:mz-1,pb:pe),
     .        scal(0:(2*my-1),0:mz-1,pb:pe)

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
      write(ext1,'(i3.3)') id22  
      id22=id22+1         
c    ------------------------  Program      -----------------

        ! save procs wait for comp. ones ! 
        
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
c    ----------------------- Statistics ----------------------         

      ! Save procs receive data from calc proc !
      

      if(myid.eq.0) then 
         call MPI_SEND(u00,my,MPI_DOUBLE_PRECISION,
     .                  numerop,1,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w00,my,MPI_DOUBLE_PRECISION,
     .                  numerop,2,MPI_COMM_WORLD,ierr)
         call MPI_SEND(v00,my,MPI_DOUBLE_PRECISION,
     .                  numerop,3,MPI_COMM_WORLD,ierr)
      elseif(myid.eq.numerop) then 
         call MPI_RECV(u00,my,MPI_DOUBLE_PRECISION,
     .                  0,1,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(w00,my,MPI_DOUBLE_PRECISION,
     .                  0,2,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(v00,my,MPI_DOUBLE_PRECISION,
     .                  0,3,MPI_COMM_WORLD,istat,ierr)
      endif
      
      call mpi_barrier(MPI_COMM_WORLD,ierr)            
      leng = 2*my*mz 
           
      if (myid.lt.numerop) then    ! calc proc
         do i=pb,pe
            call MPI_SEND(vor(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr) 
                 
            call MPI_SEND(phi(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr)

            call MPI_SEND(psi(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr)

            call MPI_SEND(scal(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr)
         enddo
         return
      else
         do i = pb,pe
            call MPI_RECV(vor(0,0,i),leng,MPI_REAL,plcalc(i),
     .                 MPI_any_tag,MPI_COMM_WORLD,istat,ierr)
              
            call MPI_RECV(phi(0,0,i),leng,MPI_REAL,plcalc(i),
     .                 MPI_any_tag,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(psi(0,0,i),leng,MPI_REAL,plcalc(i),
     .                 MPI_any_tag,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(scal(0,0,i),leng,MPI_REAL,plcalc(i),
     .                 MPI_any_tag,MPI_COMM_WORLD,istat,ierr)
         enddo
      endif
      
           
      if(myid.eq.numerop) then 

         !   open files    !

         if(id22.gt.999.or.id22.lt.0) then
             write(*,*) 'number of images out of range'
             stop
         endif
         fnameima=filout(1:index(filout,' ')-1)//'.'//ext1
         
         open (iout,file=fnameima,status='unknown',
     .         form='unformatted')

         rewind(iout)
         
         do j=1,my
            wk00(3*j-2) = u00(j)
            wk00(3*j-1) = v00(j)
            wk00(3*j)   = w00(j)
         enddo

         write(iout) time,Ree,alpe,bete,a0e,mx,my,mz,
     .               (y(j),fmap(j),j=1,my),
     .               (wk00(j),j=1,3*my)
     
         ! master write its fields
         
         do i=pb,pe
            !write(iout) ((vor(j,k,i),phi(j,k,i),j=0,2*my-1),k=0,mz1)
            write(iout) ((vor(j,k,i),phi(j,k,i),psi(j,k,i),
     .                   scal(j,k,i),j=0,2*my-1),k=0,mz1)
         enddo
         
         ! then calls slaves
         
         do iproc=numerop+1,numtot-1
         
            call MPI_RECV(leng,1,MPI_INTEGER,
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(vor(0,0,1),leng,MPI_REAL,
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            call MPI_RECV(phi(0,0,1),leng,MPI_REAL,
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            call MPI_RECV(psi(0,0,1),leng,MPI_REAL,
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            call MPI_RECV(scal(0,0,1),leng,MPI_REAL,
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            write(*,*) iproc
            do i=1,pend(iproc-numerop)-pbeg(iproc-numerop)+1
            !  write(iout)((vor(j,k,i),phi(j,k,i),j=0,2*my-1),k=0,mz1)
               write(iout) ((vor(j,k,i),phi(j,k,i),psi(j,k,i),
     .                       scal(j,k,i),j=0,2*my-1),k=0,mz1)
            enddo
            
         enddo
         close(iout)
         
      else
      
         ! slaves send their fields
         leng = 2*my*mz*mmp
         call MPI_SEND(leng,1,MPI_INTEGER,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)           

         call MPI_SEND(vor(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr) 
     
         call MPI_SEND(phi(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)

         call MPI_SEND(psi(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)

         call MPI_SEND(scal(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)
     
      endif
      
    
      end subroutine escru

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

      include "mpif.h"
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
 

!--------------------------------------------------------------------!
! Subroutine WRITESTATS
! AAF June 2014
!....................................................................!
      subroutine writestats(sp,spwk,myid,istep)
      
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
      character*84 fnamesta,fnamespe
      character*3 ext1
      real*4 Ree,alpe,bete,a0e
      integer istep,leng,kk,iproc
!spectra
      real*4 sp  (0:nz1,1:2*nspec+1,12,pb:pe),
     .       spwk(0:nz1,1:  nspec+1,12,pb:pe)

    
      !extension number of sta file & increment on it
      write(ext1,'(i3.3)') ista 
      ista = ista +1
      Ree  = Re
      alpe = alp
      bete = bet
      a0e  = a0

      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       
      !Compute stats
      !Note:MEAN VALUES do not need to REDUCE because only MASTER PROC 
      ! accumulates the data
      if (myid.lt.numerop) then 
         call MPI_ALLREDUCE(up,wkstats(1,1),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(vp,wkstats(1,2),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(wp,wkstats(1,3),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(uvr,wkstats(1,4),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(uwr,wkstats(1,5),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(vwr,wkstats(1,6),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w1p,wkstats(1,7),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w2p,wkstats(1,8),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w3p,wkstats(1,9),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(Tp,wkstats(1,10),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(ep,wkstats(1,11),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(ruu,wkstats(1,12),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(ruv,wkstats(1,13),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(ruw,wkstats(1,14),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(rvv,wkstats(1,15),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(rvw,wkstats(1,16),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(rww,wkstats(1,17),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(rhom,wkstats(1,18),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
 
       endif

      !WRITE STATS
      if (myid.eq.0) then
         !variables accumulated onf phys space need to be averaged
         !by fac
         fac = 1./dble(mgalx*mgalz)
         write(*,*) "Writing stats at istep = ",istep, 
     .              "nacum = ",nacum
         fac = 1.0/(mgalx*mgalz)
         fnamesta = filstt(1:index(filstt,' ')-1)//'_'//ext1//'.sta'
         open (isn,file=fnamesta,status='unknown',form='unformatted')
         rewind(isn)
         write(isn) my,time,Ree,alpe,bete,a0e,nacum
         write(isn) (y(j),fmap(j),j=1,my)
         write(isn) (um(j),vm(j),wm(j),
     .              rum(j),rvm(j),rwm(j),
     .              Tm(j),
     .              w1m(j),w2m(j),w3m(j),
     .              wkstats(j,1),wkstats(j, 2),
     .              wkstats(j,3),wkstats(j, 4),
     .              wkstats(j,5),wkstats(j, 6),
     .              wkstats(j,7),wkstats(j, 8),
     .              wkstats(j,9),wkstats(j,10),
     .              wkstats(j,11),
     .              wkstats(j,12)*fac,
     .              wkstats(j,13)*fac,wkstats(j,14)*fac,
     .              wkstats(j,15)*fac,wkstats(j,16)*fac,
     .              wkstats(j,17)*fac,wkstats(j,18)*fac,
     .                                             j=1,my)
         close(isn)
      endif
c   ----------------------- Spectra --------------------
           ! All processors computes spectra     !

c  Warning!!! All processors have to compute the same number of planes

      leng = (mx1+1)*(nz1+1)*12
      do i=pb,pe
         ! velocity and vorticity spectra are symmetric !
         do kk=1,3
            do j=2*nspec+1,nspec+1,-1
               jj = 2*nspec+2 - j
               do k = 0,nz1
                  spwk(k,jj,kk,i) =.5*(sp(k,j,kk,i) + sp(k,jj,kk,i))
               enddo
            enddo
         enddo
         do kk=6,12
            do j=2*nspec+1,nspec+1,-1
               jj = 2*nspec+2 - j
               do k = 0,nz1
                  spwk(k,jj,kk,i) =.5*(sp(k,j,kk,i) + sp(k,jj,kk,i))
               enddo
            enddo
         enddo

c------------- velocity cospectra are skew-symmetric
         do kk=4,5
            do j=2*nspec+1,nspec+1,-1
               jj = 2*nspec+2 - j
               do k=0,nz1
                  spwk(k,jj,kk,i) =.5*(-sp(k,j,kk,i) + sp(k,jj,kk,i))
               enddo
            enddo
         enddo
      enddo     ! i

c------- everybody sends data to the master

      if (myid.ne.0.and.myid.lt.numerop) then
         leng = (nspec+1)*(nz1+1)*12*mmp
         call MPI_SEND(spwk(0,1,1,pb),leng,MPI_REAL,
     .        0,myid,MPI_COMM_WORLD,ierr)

         write(*,*) 'envio espectros',myid

      elseif(myid.eq.0) then

c-----  the master first writes its stuff

         fnamespe=filstt(1:index(filstt,' ')-1)//'_'//ext1//'.spe'
         open (ispf,file=fnamespe,status='unknown',
     .         form='unformatted')
         rewind(ispf)
         write(*,*) ispf
         write(ispf) time,Ree, alpe, bete, mx,my,mz,nspec+1,
     .               nacumsp
        write(ispf) (jsptot(j), j=1,nspec+1),(y(j),fmap(j),j=1,my)

         write(*,*) 'escribe espectro, proc ',myid,
     .               'pb= ',pb

         do i=1,mmp
            write(ispf)(spwk(j,1,1,i),j=0,12*(nz1+1)*(nspec+1)-1)
         enddo

c----- then receives everything from everybody

         leng = (nspec+1)*(nz1+1)*12*mmp
         do iproc=1,numerop-1
            call MPI_RECV(spwk(0,1,1,1),leng,MPI_REAL,iproc,
     .                    iproc,MPI_COMM_WORLD,istat,ierr)

c------ and writes it
            do i=1,mmp
               write(ispf) (spwk(j,1,1,i),j=0,12*(nz1+1)*(nspec+1)-1)
            enddo
         enddo
         write(*,*) 'closing files'
         close(ispf)
      endif

      do i=pb,pe
         do ii=1,12
            do j=1,2*nspec+1
               do k=0,nz1
                  sp(k,j,ii,i) = 0.0
               enddo
            enddo
         enddo
      enddo

      nacumsp  = 0


    
      end subroutine
!-------------------------------------------------------------------!

!-----------------------------------------------------------------------!
      subroutine resetstats(sp,myid)

      use statis
      use point
      use spectra
      implicit none
      include "ctes3D"
      include "mpif.h"
    
      integer i,j,k,jj,kk,ii,myid
      real*4 sp  (0:nz1,1:2*nspec+1,12,pb:pe)

      nacum = 0
      nacumsp = 0
      if (myid.eq.0) then
         write(*,*) " Reset stats..."
      endif
      !clean mean buffers
      do j=1,my
         um(j) = 0d0
         vm(j) = 0d0
         wm(j) = 0d0
         rum(j) = 0d0 !rhoum
         rvm(j) = 0d0 !rhovm
         rwm(j) = 0d0 !rhowm
         Tm(j) = 0d0
         w1m(j) = 0d0
         w2m(j) = 0d0
         w3m(j) = 0d0
      enddo
      !Clean ALL STATIS BUFFERS
      do j=1,my
          up(j) = 0d0
          vp(j) = 0d0
          wp(j) = 0d0
         uvr(j) = 0d0
         uwr(j) = 0d0
         vwr(j) = 0d0
         w1p(j) = 0d0
         w2p(j) = 0d0
         w3p(j) = 0d0
         Tp(j) = 0d0
         ep(j) =0d0
         ruu(j) = 0d0
         ruv(j) = 0d0
         ruw(j) = 0d0
         rvv(j) = 0d0
         rvw(j) = 0d0
         rww(j) = 0d0
         rhom(j) = 0d0
      enddo

      do i=pb,pe
         do ii=1,12
            do j=1,2*nspec+1
               do k=0,nz1
                  sp(k,j,ii,i) = 0.0
               enddo
            enddo
         enddo
      enddo


      end subroutine
!-------------------------------------------------------------------!
 





