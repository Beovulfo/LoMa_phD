      subroutine chpl2ln(recta,plano,wa,myid)
      use MPI_GROUPS
      use point
      use MPI_datatype
      use timers

      implicit none
      include "mpif.h"
      include "ctes3D"

c     ----------------------- Variables -------------------------  
     
      integer istat(MPI_STATUS_SIZE),ierr,myid,i,j,ii
      integer iproc,pproc,nsetotr,nsetotrid,numprocs,block
      
      complex*8  plano(my*mgalz,pb:pe)
      complex*8  recta(mx1+1,lb:le)
      complex*8  wa(lb:le,1:mx1+1)

c     ----------------------- Programa --------------------------       

      if (myid.eq.0) then
         commtimer = commtimer-MPI_WTIME()
      endif

      do i=pb,pe
         do j=lb,le
            wa(j,i) = plano(j,i)
         enddo
      enddo
      nsetotr=2*(le-lb+1)
      do pproc=1,pnodes-1
         iproc=ieor(myid,pproc)          
         if (iproc<numerop) then
            nsetotrid = 2*(lend(iproc)-lbeg(iproc)+1)
            do ii=pb,pe
               call MPI_SENDRECV(plano(lbeg(iproc),ii),nsetotrid,
     .           MPI_REAL,iproc,0,
     .           wa(lb,pbeg(iproc)+ii-pb),nsetotr,MPI_REAL,iproc,0,
     .           MPI_COMM_WORLD,istat,ierr)
            enddo
         endif
      enddo
      
      if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
      endif

      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif

      block = 64
      do ii = 1,mx1+1,block
         do j=lb,le     
            do i=ii,ii+block-1
               recta(i,j) = wa(j,i)
            enddo
         enddo
      enddo

      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
      endif


      end subroutine chpl2ln

c     ---------------------------------------------------------
c     ---------------------------------------------------------
c     ---------------------------------------------------------

      subroutine chln2pl(recta,plano,wa,myid)
      use point
      use MPI_datatype
      use timers
      use MPI_GROUPS

      implicit none

      include "mpif.h"
      include "ctes3D"

c     ----------------------- Variables -------------------------         

      integer istat(MPI_STATUS_SIZE),ierr,myid,jj,ii
      integer iproc,pproc,nsetotr,i,j,k,mmz2,block,nsetotrid
      
      complex*8  plano(my*mgalz,pb:pe)
      complex*8  recta(mx1+1,lb:le)
      complex*8  wa(lb:le,1:mx1+1)

      
c     ----------------------- Programa -------------------------- 
      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif

      block = 64
      do ii = 1,mx1+1,block
         do j=lb,le     
            do i=ii,ii+block-1
               wa(j,i) = recta(i,j)
            enddo
         enddo
      enddo
      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
      endif 
           
      if (myid.eq.0) then
         commtimer = commtimer-MPI_WTIME()
      endif
      
      nsetotr=2*(le-lb+1)
      do pproc=1,pnodes-1
         iproc=ieor(myid,pproc)          
         if (iproc<numerop) then
            nsetotrid = 2*(lend(iproc)-lbeg(iproc)+1)
            do ii=pb,pe
               call MPI_SENDRECV(wa(lb,pbeg(iproc)+ii-pb),nsetotr
     .           ,MPI_REAL,iproc,0
     .           ,plano(lbeg(iproc),ii),nsetotrid,MPI_REAL,iproc,0
     .           ,MPI_COMM_WORLD,istat,ierr)
            enddo
         endif
      enddo
            
      do i=pb,pe
         do j=lb,le
            plano(j,i) = wa(j,i)
         enddo
      enddo

      if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
      endif


      end


      subroutine pointers_calc(pb,pe,lb,le,procs,myid)
c =============================================================
c     jjs,     aug/2001  (bug in old version)
c     shc,     aug/2005      plane-pencil
c =============================================================

      implicit none

      include "ctes3D"
   
      integer pb(numerop),pe(numerop),lb(numerop),le(numerop)
      integer n,n1,n2,i,procs,myid                 
c    nplanes planes (mx1+1)
!    nplanes = mx/2 from ctes3D
      n1=nplanes/numerop
      n2=nplanes-numerop*n1  
      pb(1)=1
      do n=1,n2
         pe(n)  = pb(n)+n1
         pb(n+1)= pe(n)+1
      enddo
      do n=n2+1,numerop-1
         pe(n)  =pb(n)+n1-1
         pb(n+1)= pe(n)+1
      enddo
      pe(numerop) = pb(numerop)+n1-1
      
      ! lb will be used as wk in the next loop
      
      n1=nplanes/numerosa
      n2=nplanes-numerosa*n1
      lb(1)=1
      do n=1,n2
         le(n)  = lb(n)+n1
         lb(n+1)= le(n)+1
      enddo
      do n=n2+1,numerosa-1
         le(n)  =lb(n)+n1-1
         lb(n+1)= le(n)+1
      enddo
      le(numerosa) = lb(numerosa)+n1-1
      
      do i=1,numerosa
          if (pb(myid+1).gt.lb(numerosa)-1) procs = numerosa+numerop-1
      enddo 
          
      

c     nlines lines (mgalz*my)

      n1=nlines/numerop
      n2=nlines-numerop*n1

      lb(1)=1

      do n=1,n2
         le(n)  = lb(n)+n1
         lb(n+1)= le(n)+1
      enddo
      do n=n2+1,numerop-1
         le(n)  =lb(n)+n1-1
         lb(n+1)= le(n)+1
      enddo
      le(numerop)=lb(numerop)+n1-1

      end
      
            
      subroutine pointers_save(pb,pe)
c =============================================================
c     jjs,     aug/2001  (bug in old version)
c     shc,     aug/2005      plane-pencil
c =============================================================

      implicit none

      include "ctes3D"

      integer pb(numerosa),pe(numerosa)
      integer n,n1,n2                 

c    nplanes planes (mx1+1)

      n1=nplanes/numerosa
      n2=nplanes-numerosa*n1
      pb(1)=1
      do n=1,n2
         pe(n)  = pb(n)+n1
         pb(n+1)= pe(n)+1
      enddo
      do n=n2+1,numerosa-1
         pe(n)  =pb(n)+n1-1
         pb(n+1)= pe(n)+1
      enddo
      pe(numerosa) = pb(numerosa)+n1-1
      
      end
