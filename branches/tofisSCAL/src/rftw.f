c ==============================================================
c              initialize rft_w
c     uses fftw_3.0, supposed to work as the old rftsingle 
c              tested only in aeolos 
c     compile with flags: 
c     
c     ifort -O3 -tpp7 -xW -align dcommons -c rftw3.f \
c         -lfftw3f -lm -I/usr/local/include
c
c ==============================================================
      subroutine rfti(n)
      implicit none
      include "fftw3.f"

      integer*4  n,nmax
      parameter (nmax=8192)

      real*4     fdum,bdum,dnf,dnb
      integer*4  nf,nb
      integer*8  planf,planb
      common /rfttofour/ planf,nf,dnf,fdum(nmax+2)
      common /rfttophys/ planb,nb,dnb,bdum(nmax+2)
      save /rfttofour/, /rfttophys/

      nb=n
      nf=n
      dnf=1d0/n
      dnb=1d0
      
      call sfftw_plan_dft_r2c_1d(planf,nf,fdum,fdum,FFTW_MEASURE )
      call sfftw_plan_dft_c2r_1d(planb,nb,bdum,bdum,FFTW_MEASURE )

      end


c ================================================================
c             the real  rftsingle 
c ================================================================
      subroutine rft(c,isa,m,iopt)
      implicit none
      integer isa,m,iopt,i,j
      real*4 c(isa,*)

      real*4     fdum,bdum,dnf,dnb
      integer*4  nf,nb
      integer*8  planf,planb
      common /rfttofour/ planf,nf,dnf,fdum(1)  
      common /rfttophys/ planb,nb,dnb,bdum(1)
      save /rfttofour/, /rfttophys/

      if (iopt<0) then

         do j=1,m
            do i=1,nf 
               fdum(i)=c(i,j)
            enddo
            call sfftw_execute(planf)
            do i=1,nf+2 
               c(i,j)=dnf*fdum(i)
            enddo
         enddo

      else   

         do j=1,m
            do i=1,nb+2 
               bdum(i)=c(i,j)
            enddo
            call sfftw_execute(planb)
            do i=1,nb 
               c(i,j)=bdum(i)
            enddo            
            c(nb+1,j)=0.
            c(nb+2,j)=0.
         enddo

      endif

      end
