c ==============================================================
c              initialize cft_w
c     uses fftw_3.0, supposed to work as the old rftsingle 
c              tested only in aeolos 
c     compile with flags: 
c     
c     ifort -O3 -tpp7 -xW -align dcommons -c rftw3.f \
c         -lfftw3f -lm -I/usr/local/include
c
c ==============================================================
      subroutine cfti(n)
      implicit none
      include "fftw3.f"

      integer*4  n,nmax,nth,i
      parameter (nmax=4096*2)

      real*4     dnf,dnb 
      complex*8  fdum,bdum
      integer*4  nf,nb
      integer*8  planf,planb
      common /cfttofour/ planf,nf,dnf,fdum(nmax+2)
      common /cfttophys/ planb,nb,dnb,bdum(nmax+2)
      save /cfttofour/, /cfttophys/
      

      nb=n
      nf=n
      dnf=1d0/n
      dnb=1d0
      
      call sfftw_plan_dft_1d(planf,nf,fdum(1),fdum(1),
     .                                FFTW_FORWARD,FFTW_MEASURE )
      call sfftw_plan_dft_1d(planb,nb,bdum(1),bdum(1),
     .                                FFTW_BACKWARD,FFTW_MEASURE )
      
      end


c ================================================================ c
c             the real  rftsingle                                  c
c ================================================================ c

      subroutine cft(c,ise,isa,m,iopt)
      implicit none
      integer isa,m,iopt,i,j,ise,ise1,isa1,j0,j00,nth,nmax
      complex*8 c(*)

      parameter (nmax=4096*2)
      
      real*4     dnf,dnb 
      complex*8  fdum,bdum
      integer*4  nf,nb
      integer*8  planf,planb
      common /cfttofour/ planf,nf,dnf,fdum(nmax)
      common /cfttophys/ planb,nb,dnb,bdum(nmax)
      save   /cfttofour/, /cfttophys/

      isa1 = isa/2
      ise1 = ise/2
 
      j00  = 1
      
      if (iopt<0) then

         do j=1,m
            j0 = j00
            do i=1,nf 
               fdum(i)=c(j0)
               j0 = j0 + ise1
            enddo

            call sfftw_execute(planf)

            j0 = j00
            do i=1,nf
               c(j0)=dnf*fdum(i)
               j0 = j0 + ise1
            enddo
            j00 = j00 + isa1
         enddo

      else   

         do j=1,m
            j0 = j00
            do i=1,nb 
               bdum(i)=c(j0)
               j0 = j0 + ise1
            enddo

            call sfftw_execute(planb)
            j0 = j00
            do i=1,nb 
               c(j0)=bdum(i)
               j0 = j0 + ise1
            enddo   
            j00 = j00 + isa1
         enddo

      endif

      end
