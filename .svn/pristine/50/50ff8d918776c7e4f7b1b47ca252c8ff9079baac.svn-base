!******************************************************************!
!                                                                  !
!               writes a intermediate solution                     !
!    single jjs  4/01/01                                           !
!    plane  lines shc 01/03/05                                     !
!                                                                  !
!                                                                  !
!                                                                  !
!******************************************************************!
      subroutine escru(vor,phi,u00,w00,v00,myid)
!      subroutine escru(vor,phi,u00,w00,sp,spwk,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"

c ------------------Commons --------------------------------

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
      
      real*4  vor(0:(2*my-1),0:mz-1,pb:pe),
     .        phi(0:(2*my-1),0:mz-1,pb:pe)

!      real*4 sp  (0:nz1,1:2*nspec+1,8,pb:pe),
!     .       spwk(0:nz1,1:  nspec+1,8,pb:pe)

      
      real*8 u00(my),w00(my),v00(my)
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
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            call MPI_RECV(phi(0,0,1),leng,MPI_REAL,
     .            iproc,iproc,MPI_COMM_WORLD,istat,ierr)            
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
      
      real*4 wk00(2*my),Ree,alpe,bete,a0e

      
      real*8 fac,wk(my,13)      

      character*3 ext1,scalstr
      character*84 fnameima
      integer ioutscal
!-------------------------------------------------------------!
      Ree  = Re
      alpe = alp
      bete = bet
      a0e  = a0
      ioutscal = 66
     
      write(ext1,'(i3.3)') id22 
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
      
!      call mpi_barrier(MPI_COMM_WORLD,ierr) 
         do j=1,my
            wk00(2*j-1) = 0.0
            wk00(2*j)   = 0.0
         enddo
 
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
         
         open (ioutscal,file=fnameima,status='unknown',
     .         form='unformatted')

         rewind(ioutscal)
         
         
         ! master write its fields
         !First write the heading
         write(ioutscal) time,Ree,alpe,bete,a0e,mx,my,mz,
     .               (y(j),fmap(j),j=1,my),
     .               (wk00(j),j=1,2*my)
         
         do i=pb,pe
            write(ioutscal) ((scal(j,k,i),j=0,2*my-1),k=0,mz1)
         enddo
         
         ! then calls slaves
         
         do iproc=numerop+1,numtot-1
         
            call MPI_RECV(leng,1,MPI_INTEGER,
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)

            call MPI_RECV(scal(0,0,1),leng,MPI_REAL,                   
     .              iproc,iproc,MPI_COMM_WORLD,istat,ierr)            

            write(*,*) iproc
            do i=1,pend(iproc-numerop)-pbeg(iproc-numerop)+1
               write(ioutscal)((scal(j,k,i),j=0,2*my-1),k=0,mz1)
            enddo
            
         enddo
         close(ioutscal)
         
      else
      
         ! slaves send their fields
         leng = 2*my*mz*mmp
         call MPI_SEND(leng,1,MPI_INTEGER,
     .                 numerop,myid,MPI_COMM_WORLD,ierr)           

         call MPI_SEND(scal(0,0,pb),leng,MPI_REAL,
     .                 numerop,myid,MPI_COMM_WORLD,ierr) 
     
     
      endif
      
         
         
      end subroutine escruscal


!=======================================================!
!AAF subroutine to write any size variable
! to a file for matlab
      subroutine write2file(fwk,varname,tamy,tamz,tamx,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"
!-------------------------------------------------------!      
      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save /ficheros/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/


      integer fid 
      integer tamy,tamz,tamx,myid
      real*4 fwk(2,1:tamy,1:tamz,1:tamx)
      character*84 fnameima,fnamespe,fnamesta
      character*3 varname
      character*3 ext1
      integer i,j,k

      
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

      endsubroutine
 


!=======================================================!
!AAF subroutine to write any size variable
! to a file for matlab
      subroutine write2file1(fwk,varname,tamy,tamz,tamx,myid)

      implicit none

      include "mpif.h"
      include "ctes3D"
!-------------------------------------------------------!      
      integer iinp,iout,id22,isn,ispf
      character*70 filinp,filout,filstt,filscainp
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                  filinp,filout,filstt,filscainp
      save /ficheros/

      real*8  Re,alp,bet,a0
      real*8  y,hy,fmap,trp,mss
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),
     $            trp(0:my-1),mss(0:my-1)
      save   /fis/


      integer fid 
      integer tamy,tamz,tamx,myid
      real*8 fwk(1:tamy,1:tamz,1:tamx)
      character*84 fnameima,fnamespe,fnamesta
      character*3 varname
      character*3 ext1
      integer i,j,k


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

      endsubroutine
 


