      subroutine chpl2ln(recta,plano,wa,myid)
      use timers
      use point

      implicit none
      include "ctes3D"

c     ----------------------- Variables -------------------------  
     
      integer myid,i,j,ii
      integer block
      
      complex*8  plano(my*mgalz,pb:pe)
      complex*8  recta(mx1+1,lb:le)
      complex*8  wa(lb:le,1:mx1+1)

c     ----------------------- Programa --------------------------       

      do i=pb,pe
         do j=lb,le
            wa(j,i) = plano(j,i)
         enddo
      enddo
      
      block = 64
      do ii = 1,mx1+1,block
         do j=lb,le     
            do i=ii,ii+block-1
               recta(i,j) = wa(j,i)
            enddo
         enddo
      enddo

      end subroutine chpl2ln

c     ---------------------------------------------------------
c     ---------------------------------------------------------
c     ---------------------------------------------------------

      subroutine chln2pl(recta,plano,wa,myid)
      use point
      use timers

      implicit none

      include "ctes3D"

c     ----------------------- Variables -------------------------         

      integer myid,jj,ii
      integer i,j,k,mmz2,block
      
      complex*8  plano(my*mgalz,pb:pe)
      complex*8  recta(mx1+1,lb:le)
      complex*8  wa(lb:le,1:mx1+1)

      
c     ----------------------- Programa -------------------------- 

      block = 64
      do ii = 1,mx1+1,block
         do j=lb,le     
            do i=ii,ii+block-1
               wa(j,i) = recta(i,j)
            enddo
         enddo
      enddo
      
      do i=pb,pe
         do j=lb,le
            plano(j,i) = wa(j,i)
         enddo
      enddo

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
