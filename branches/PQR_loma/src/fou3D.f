c/********************************************************************/
c/*   FAST TRANSFORMS INTERFACE PACKAGE                              */
c/*..................................................................*/
c/*  MODULE FOR MPI SP2                                              */
c/*..................................................................*
c/*                                                                  */
c/********************************************************************/
        subroutine fourxz(fou,phys,iopt,mmy)
        implicit none
        include"ctes3D"
        integer iopt,mmy
        real*4 fou(*),phys((mgalx+2)*mgalz,*)
        integer j,ipy
c/********************************************************************/
c/*                                                                  */
c/*    makes a 3-dimensional real to complex fourier transform       */
c/*      from F-F-Physical   (array fou)                             */
c/*  to/ from Phys-Phys-Phys (array phys)                            */
c/*                                                                  */
c/*       iopt >=0    ===>  inverse transforms( fou ---> fis)        */
c/*       iopt < 0    ===>  direct  transforms( fis ---> fou)        */
c/*                                                                  */
c/*       does the expanding and compacting itself                   */
c/*                                                                  */
c/*  NOTE:                                                           */
c/*    fis ---> fou : supposes high frecuency modes are not          */
c/*                   needed (dealiasing) and throws them away       */
c/*                                                                  */
c/********************************************************************/
c                                /* physical to fourier transforms  */
        if (iopt.lt.0) then
           call rft  (phys,mgalx+2,mmy*mgalz,-1)
           do 20 j=1,mmy
              call cft(phys(1,j),mgalx+2,2,mx1+1,-1)
 20        continue
           call pack0(phys,fou,mmy)
c                                /* fourier to physical transforms  */
         else
           call fill00(fou,phys,mmy)
           do 30 j=1,mmy
              call cft(phys(1,j),mgalx+2,2,mx1+1,1)
 30        continue
           call rft(phys,mgalx+2,mmy*mgalz,1)
         endif
         end


c/********************************************************************/
c/*                                                                  */
c/*         fills high armonics of a 3D fourier representation       */
c/*         (f-f-T) with  zeroes to get extra points enough to       */
c/*         perform dealiasing.                                      */
c/*                                                                  */
c/*     input:                                                       */
c/*       fou: variable at fourier space                             */
c/*                                                                  */
c/*    output:                                                       */
c/*      dfou: variable at fourier space ready to dealiasing         */
c/*                                                                  */
c/********************************************************************/
      subroutine fill00(fou,foud,mmy)
      implicit none
      include "ctes3D"
      integer mmy
      complex*8 fou(0:mx1,0:mz1,mmy),foud(0:mgx,0:mgalz1,mmy),zero
      integer i,j,k,kk
      zero = cmplx(0.,0.)
cDIR$ NOVECTOR
      do 10 j=mmy,1,-1
cDIR$ NOVECTOR
         do 70 k=mgalz1,nz2,-1
            kk = k-mgalz1+mz1
cDIR$ NOVECTOR
            do 80 i=mx1,0,-1
               foud(i,k,j)=fou(i,kk,j)
 80         continue
 70      continue

cDIR$ NOVECTOR
         do 20 k=nz1,0,-1
cDIR$ NOVECTOR
            do 30 i=mx1,0,-1
               foud(i,k,j)=fou(i,k,j)
 30         continue
 20      continue

cDIR$ NOVECTOR
         do 50 k=nz2-1,nz1+1,-1
cDIR$ NOVECTOR
            do 60 i=mx1,0,-1
               foud(i,k,j)=zero
 60         continue
 50      continue

cDIR$ NOVECTOR
         do 25 k=mgalz1,0,-1
cDIR$ NOVECTOR
            do 35 i=mgx,mx1+1,-1
               foud(i,k,j)=zero
 35         continue
 25      continue

 10   continue

      end



c/********************************************************************/
c/*                                                                  */
c/*         fills high armonics of a 3D fourier representation       */
c/*         (f-f-T) with  zeroes to get extra points enough to       */
c/*         perform dealiasing.                                      */
c/*                                                                  */
c/*     input:                                                       */
c/*      dfou: variable at fourier space after dealiasing with       */
c/*            zeroes at high frecuencies                            */
c/*                                                                  */
c/*    output:                                                       */
c/*       fou: packed variable at fourier space (if fou=fdou         */
c/*            fou overwrites fdou)                                  */
c/*                                                                  */
c/********************************************************************/
      subroutine pack0(foud,fou,mmy)
      implicit none
      include "ctes3D"
      integer mmy
      complex*8 fou(0:mx1,0:mz1,mmy),foud(0:mgx,0:mgalz1,mmy)
      integer i,j,k,kk
      do 10 j=1,mmy

         do 20 k=0,nz1
            do 30 i=0,mx1
               fou(i,k,j)=foud(i,k,j)
 30         continue
 20      continue

         do 70 k=nz2,mgalz1
            kk = k-mgalz1+mz1
            do 80 i=0,mx1
               fou(i,kk,j)=foud(i,k,j)
 80         continue
 70      continue

 10   continue

      end


c/********************************************************************/
c/*                                                                  */
c/*      computes the derivative of u respect chi and gives it       */
c/*      back en du (for an array of the form array(mx,mz,my)        */
c/*      in f-f-t space)                                             */
c/*                                                                  */
c/*     input:                                                       */
c/*         u: array with fourier coefficients                       */
c/*       chi: indicates the direction of the derivative             */
c/*            (options chi='x', or 'z')                             */
c/*    output:                                                       */
c/*        du: derivative coefficients (if du=u then du over-        */
c/*            write u)                                              */
c/*                                                                  */
c/********************************************************************/
      subroutine deriv(u,du,chi,mmy)
      implicit none
      include "ctes3D"
      integer mmy
      complex*8 u(0:mx1,0:mz1,mmy),du(0:mx1,0:mz1,mmy)
      character chi
      complex*8    dk

      integer iax,icx
      real*4 alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     >              iax(mx),icx(0:mz1)
      save /wave/

      complex*8 ci,zero
      integer i,j,k

      ci=cmplx(0.,1.)
      zero = cmplx (0e0,0e0)

      if(chi.eq.'x') then

         do 10 j=1,mmy
            do 20 k=0,mz1
               do 30 i=0,mx1
                  du(i,k,j)=xalp(i)*u(i,k,j)
 30            continue
 20         continue
 10      continue

      elseif(chi.eq.'z') then

         do 40 j=1,mmy
            do 50 k=0,mz1
               dk=xbet(k)
               do 60 i=0,mx1
                  du(i,k,j)=dk*u(i,k,j)
 60            continue
 50         continue
 40      continue

      else
          write(*,*) 'error en deriv'
      endif
      end
