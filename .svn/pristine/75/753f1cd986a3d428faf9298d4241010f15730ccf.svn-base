!******************************************************************!
!                                                                  !
!               writes a intermediate solution                     !
!    single jjs  4/01/01                                           !
!    plane  lines shc 01/03/05                                     !
!                                                                  !
!                                                                  !
!                                                                  !
!******************************************************************!
      subroutine escru(vor,phi,u00,w00,sp,spwk,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

c ------------------Commons --------------------------------

      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .        ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm     
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/

      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save /mesh/

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/
 
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/

      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save /ficheros/

      integer nacumsp,jsptot
      common/spectra/   nacumsp,jsptot(2*nspec+1)
      save/spectra/
      
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
     
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     
          
c  ------------------------- Variables ----------------------------
     
      integer iproc,i,j,k,leng,jj,kk,ii
      integer myid,istat(MPI_STATUS_SIZE),ierr
      
      real*4  vor(0:2*my-1,0:mz-1,pb:pe),
     .        phi(0:2*my-1,0:mz-1,pb:pe)

      real*4 sp  (0:nz1,1:2*nspec+1,8,pb:pe),
     .       spwk(0:nz1,1:  nspec+1,8,pb:pe)

      
      real*8 u00(my),w00(my)
      real*8 fac,wk(my,13)      
      real*4 wk00(2*my),Ree,alpe,bete,a0e

      character*3 ext1
      character*84 fnameima,fnamespe,fnamesta
     
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
         
            
c        Computes statistics
      if (myid.lt.numerop) then
         call MPI_ALLREDUCE(up,wk(1,1),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(vp,wk(1,2),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(wp,wk(1,3),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w1p,wk(1,4),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w2p,wk(1,5),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(w3p,wk(1,6),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(uvr,wk(1,7),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(uwr,wk(1,8),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(vwr,wk(1,9),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(ep, wk(1,10),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(uuv,wk(1,11),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(wwv,wk(1,12),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
         call MPI_ALLREDUCE(vvv,wk(1,13),my,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_CALC,ierr)
      
      endif

c             !       write statistics       !

      if (myid.eq.0) then 
         fac = 1./dble(mgalx*mgalz)
         fnamesta=filstt(1:index(filstt,' ')-1)//'_'//ext1//'.sta'
         open (isn,file=fnamesta,status='unknown',form='unformatted')
         rewind(isn)
         write(isn) nacum,Wx0a/nacum,Wz0a/nacum
         write(isn) my,time,Ree,alpe,bete,a0e
         write(isn) (y(j),fmap(j),j=1,my)
         write(isn)(um(j) ,vm(j) ,wm(j) ,wk(j,1),wk(j,2),wk(j,3),
     .              w1m(j),w2m(j),w3m(j),wk(j,4),wk(j,5),wk(j,6),
     .              wk(j,7),wk(j,8),wk(j,9),wk(j,10),
     .              wk(j,11)*fac,wk(j,12)*fac,wk(j,13)*fac,j=1,my)
         close(isn)
      endif
 
      do j=1,my
         um(j)  = 0.
         vm(j)  = 0.
         wm(j)  = 0.
         up(j)  = 0.
         vp(j)  = 0.
         wp(j)  = 0.
         uvr(j) = 0.
         uwr(j) = 0.
         vwr(j) = 0.
         w1m(j) = 0.
         w2m(j) = 0.
         w3m(j) = 0.
         w1p(j) = 0.
         w2p(j) = 0.
         w3p(j) = 0.
         ep(j)  = 0.
         uuv(j) = 0.
         wwv(j) = 0.
         vvv(j) = 0.
      enddo

c   ----------------------- Spectra --------------------
     
     
           ! All processors computes spectra     !

c  Warning!!! All processors have to compute the same number of planes

      leng = (mx1+1)*(nz1+1)*8
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
         do kk=6,8
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
         leng = (nspec+1)*(nz1+1)*8*mmp
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

         write(*,*) 'escribe espectro, proc ',iproc,
     .               'pb= ',pb

         do i=1,mmp
            write(ispf)(spwk(j,1,1,i),j=0,8*(nz1+1)*(nspec+1)-1)
         enddo    

c----- then receives everything from everybody

         leng = (nspec+1)*(nz1+1)*8*mmp
         do iproc=1,numerop-1
            call MPI_RECV(spwk(0,1,1,1),leng,MPI_REAL,iproc,
     .                    iproc,MPI_COMM_WORLD,istat,ierr)

c------ and writes it
            do i=1,mmp
               write(ispf) (spwk(j,1,1,i),j=0,8*(nz1+1)*(nspec+1)-1)
            enddo 
         enddo
         write(*,*) 'closing files'
         close(ispf)            
      endif

      do i=pb,pe
         do ii=1,8
            do j=1,2*nspec+1
               do k=0,nz1
                  sp(k,j,ii,i) = 0.0
               enddo
            enddo
         enddo
      enddo     

      nacumsp  = 0
      nacum = 0

      ! Save procs receive data from calc proc !
      
      

      if(myid.eq.0) then 
         call MPI_SEND(u00,my,MPI_DOUBLE_PRECISION,
     .                  numerop,1,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w00,my,MPI_DOUBLE_PRECISION,
     .                  numerop,2,MPI_COMM_WORLD,ierr)
      elseif(myid.eq.numerop) then 
         call MPI_RECV(u00,my,MPI_DOUBLE_PRECISION,
     .                  0,1,MPI_COMM_WORLD,istat,ierr)
     
         call MPI_RECV(w00,my,MPI_DOUBLE_PRECISION,
     .                  0,2,MPI_COMM_WORLD,istat,ierr)
      endif
      
      call mpi_barrier(MPI_COMM_WORLD,ierr)            
      leng = 2*my*mz 
           
      if (myid.lt.numerop) then    ! calc proc
         do i=pb,pe
         
            call MPI_SEND(vor(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr) 
                 
            call MPI_SEND(phi(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr)
         enddo
         return
      else
         do i = pb,pe
                            
            call MPI_RECV(vor(0,0,i),leng,MPI_REAL,plcalc(i),
     .                 MPI_any_tag,MPI_COMM_WORLD,istat,ierr)
              
            call MPI_RECV(phi(0,0,i),leng,MPI_REAL,plcalc(i),
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
            wk00(2*j-1) = u00(j)
            wk00(2*j)   = w00(j)
         enddo
         
         write(iout) time,Ree,alpe,bete,a0e,mx,my,mz,
     .               (y(j),fmap(j),j=1,my),
     .               (wk00(j),j=1,2*my)
     
         ! master write its fields
         
         do i=pb,pe
            write(iout) ((vor(j,k,i),phi(j,k,i),j=0,2*my-1),k=0,mz1)
         enddo
         
         ! then calls slaves
         
         do iproc=numerop+1,numtot-1
         
            call MPI_RECV(leng,1,MPI_INTEGER,
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(vor(0,0,1),leng,MPI_REAL,                   
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            call MPI_RECV(phi(0,0,1),leng,MPI_REAL,                   
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)            
            write(*,*) iproc
            do i=1,pend(iproc-numerop)-pbeg(iproc-numerop)+1
               write(iout)((vor(j,k,i),phi(j,k,i),j=0,2*my-1),k=0,mz1)
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
     
      endif
      
         
         

    
      end subroutine escru

!******************************************************************!
!                                                                  !
!               writes a intermediate solution                     !
!   adapted by AAF from escru                                      !
!                                                                  !
!                                                                  !
!******************************************************************!
      subroutine escruscal(scal,iscal,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

c ------------------Commons --------------------------------

      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .        ep,uuv,wwv,vvv,wkst,Wx0a,Wz0a,dm     
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),wkst(my),
     .                  Wx0a,Wz0a,dm,
     .                  istati,ntimes,nacum,nstart
      save /statis/

      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save /mesh/

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/
 
      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/

      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                              filinp,filout,filstt
      save /ficheros/

      integer nacumsp,jsptot
      common/spectra/   nacumsp,jsptot(2*nspec+1)
      save/spectra/
      
      integer lbeg,lend,pbeg,pend,lb,le,pb,pe,mmp,mml,procs,plcalc,plsav
      common /point /lbeg(0:numerop-1),lend(0:numerop-1),
     .               pbeg(0:numerop-1),pend(0:numerop-1),
     .               plcalc(nplanes),plsav(nplanes),
     .               pb,pe,lb,le,mmp,mml,procs
      save /point/ 
      
     
      integer MPI_GROUP_WORLD
      integer MPI_GROUP_CALC,MPI_COMM_CALC
      integer MPI_GROUP_SAVE,MPI_COMM_SAVE
      common /MPI_GROUPS/ MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save /MPI_GROUPS/
     
          
c  ------------------------- Variables ----------------------------
     
      integer iproc,i,j,k,leng,jj,kk,ii,iscal
      integer myid,istat(MPI_STATUS_SIZE),ierr
      
      real*4  scal(0:2*my-1,0:mz-1,pb:pe)

      
      real*8 fac,wk(my,13)      

      character*3 ext1,scalstr
      character*84 fnameima
     
      write(ext1,'(i3.3)') id22-1 
      write(scalstr,'(i3.3)') iscal
c    ------------------------  Program      -----------------

        ! save procs wait for comp. ones ! 
        
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
c    ----------------------- Statistics ----------------------         
!aaf  Nothing to do about it...        
!
c   ----------------------- Spectra --------------------
!aaf  Nothing to do about this    
!     

c------- everybody sends data to the master

c-----  the master first writes its stuff

c----- then receives everything from everybody

c------ and writes it


      ! Save procs receive data from calc proc !
      
      
      leng = 2*my*mz 
           
      if (myid.lt.numerop) then    ! calc proc
         do i=pb,pe
         
            call MPI_SEND(scal(0,0,i),leng,MPI_REAL,plsav(i),
     .                    plsav(i),MPI_COMM_WORLD,ierr) 
         enddo
         return
      else
         do i = pb,pe
                            
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
         fnameima=filout(1:index(filout,' ')-1)//'sc'//
     .                               scalstr//'.'//ext1
         
         open (iout,file=fnameima,status='unknown',
     .         form='unformatted')

         rewind(iout)
         
         
!         write(iout) time,Ree,alpe,bete,a0e,mx,my,mz,
!     .               (y(j),fmap(j),j=1,my)
     
         ! master write its fields
         
         do i=pb,pe
            write(iout) ((scal(j,k,i),j=0,2*my-1),k=0,mz1)
         enddo
         
         ! then calls slaves
         
         do iproc=numerop+1,numtot-1
         
            call MPI_RECV(leng,1,MPI_INTEGER,
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(scal(0,0,1),leng,MPI_REAL,                   
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            write(*,*) iproc
            do i=1,pend(iproc-numerop)-pbeg(iproc-numerop)+1
               write(iout)((scal(j,k,i),j=0,2*my-1),k=0,mz1)
            enddo
            
         enddo
         close(iout)
         
      else
      
         ! slaves send their fields
         leng = 2*my*mz*mmp
         call MPI_SEND(leng,1,MPI_INTEGER,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)           

         call MPI_SEND(scal(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr) 
     
     
      endif
      
         
         
      end subroutine

 


