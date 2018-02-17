      subroutine rfti(mgal)
c*********************************************************************/
c*     initialises coefficients for real to complex transforms       */
c*                                                                   */
c*      transforms to be of size:        mgal                        */
c*                                                                   */
c*      assumes:                                                     */
c*           mgal <= nmaxt                                           */
c*                                                                   */
c*********************************************************************/

      parameter(nmaxt=4096)
      implicit real*4(a-h,o-z)

      common /four1/ wsaf1(2*nmaxt+18)
      save /four1/

      wsaf1(1)=mgal
      wsaf1(2)=1e0/mgal
      call rffti(mgal,wsaf1(3))

      end


      subroutine rft(c,isa,m,iopt)
c*********************************************************************/
c*   given  m    ==> real vectors cm(1:2*n+2), stored in c as
c*          isa  ==> stride between initial elements of vectors
c*                   isa >= 2*n+2
c*          n    ==> length of transform, passed in common aux (rfti)
c*
c*   considering them as complex coefficients
c*      cc(k)= cm(2*k+1)+i*cm(2*k+2), (k=0,n)
c*   computes fourier expansion
c*      c(j) = cc(0) + cc(n) *(-1)**(j-1) +
c*         sum from j=1 to n-1 of cc(k)*exp(i*(j-1)*k*2*pi/2*n ) +
c*             complex conjg. of sum
c*   (or viceversa)
c*
c*   note that the imaginary part of cc(0) and cc(n) are assumed to
c*        be zero on the inverse transform, and that
c*        imag(cc(0)) =  cm(2)
c*        imag(cc(n)) =  cm(2n+2) ==> are set to zero on direct trans.
c*
c*   note also that there are 2*n+2 numbers on fourier, but only
c*                            2*n   numbers on physical
c*
c*         iopt >= 0  ==> inverse transform (x fourier --> c phys.)
c*         iopt <  0  ==> direct  transform (x physic. --> c four.)
c*             jimenez/    may-90
c*********************************************************************/
      implicit real*4 (a-h,o-z)
      dimension c(isa,*)

      common /four1/ aux(3)
      save /four1/

      nn = aux(1)
      dn1= aux(2)


c*********************************************************************/
c*            forward transform  (to fourier)                        */
c*********************************************************************/
      if (iopt.lt.0) then
         do 106 j=1,m
106         call rfftf(nn,c(1,j),aux(3))

         do 105 j=1,m
            do 100 i=nn,2,-1
100            c(i+1,j)=dn1*c(i,j)
            c(   1,j) = dn1*c(1,j)
            c(nn+2,j)=0.
            c(   2,j)=0.
105      continue
c*********************************************************************/
c*            inverse transform  (to physical)                       */
c*********************************************************************/
      else
         do 103 j=1,m
            do 103 i=2,nn
103            c(i,j)=c(i+1,j)

         do 107 j=1,m
107         call rfftb(nn,c(1,j),aux(3))
      endif

      end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c*  from here down is rcfti.packet.ncar adapted to double precision
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine rffti (n,wsave)
c**********************************************************************
c
c    subroutine rffti(n,wsave)
c
c**********************************************************************
c
c    subroutine rffti initializes the array wsave which is used in
c    both rfftf and rfftb. the prime factorization of n together with
c    a tabulation of the trigonometric functions are computed and
c    stored in wsave.
c
c    input parameter
c
c    n       the length of the sequence to be transformed.
c
c    output parameter
c
c    wsave   a work array which must be dimensioned at least 2*n+15.
c            the same work array can be used for both rfftf and rfftb
c            as long as n remains unchanged. different wsave arrays
c            are required for different values of n. the contents of
c            wsave must not be changed between calls of rfftf or rfftb.
c
c
c
c**********************************************************************

      implicit real*4 (a-h,o-z)
      dimension       wsave(*)
      if (n .eq. 1) return
      call rffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end




      subroutine rffti1 (n,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension       wa(*)      ,ifac(*)    ,ntryh(4)
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
      save

      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      tpi = 2e0*atan2(0e0,-1e0)
      argh = tpi/float(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = float(ld)*argh
            fi = 0e0
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1e0
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      return
      end


      subroutine rfftf (n,r,wsave)
c******************************************************************
c
c    subroutine rfftf(n,r,wsave)
c
c       modified from fftpack.ncar by j.jimenez   dec/87
c
c       n must be even and  only multiples allowed are 2,3,4
c
c******************************************************************
c
c    subroutine rfftf computes the fourier coefficients of a real
c    perodic sequence (fourier analysis). the transform is defined
c    below at output parameter r.
c
c    input parameters
c
c    n       the length of the array r to be transformed.  the method
c            is most efficient when n is a product of small primes.
c            n may change so long as different work arrays are provided
c
c    r       a real array of length n which contains the sequence
c            to be transformed
c
c    wsave   a work array which must be dimensioned at least 2*n+15.
c            in the program that calls rfftf. the wsave array must be
c            initialized by calling subroutine rffti(n,wsave) and a
c            different wsave array must be used for each different
c            value of n. this initialization does not have to be
c            repeated so long as n remains unchanged thus subsequent
c            transforms can be obtained faster than the first.
c            the same wsave array can be used by rfftf and rfftb.
c
c
c    output parameters
c
c    r       r(1) = the sum from i=1 to i=n of r(i)
c
c            if n is even set l =n/2   , if n is odd set l = (n+1)/2
c
c              then for k = 2,...,l
c
c                 r(2*k-2) = the sum from i = 1 to i = n of
c
c                      r(i)*cos((k-1)*(i-1)*2*pi/n)
c
c                 r(2*k-1) = the sum from i = 1 to i = n of
c
c                     -r(i)*sin((k-1)*(i-1)*2*pi/n)
c
c            if n is even
c
c                 r(n) = the sum from i = 1 to i = n of
c
c                      (-1)**(i-1)*r(i)
c
c    *****  note
c                 this transform is unnormalized since a call of rfftf
c                 followed by a call of rfftb will multiply the input
c                 sequence by n.
c
c    wsave   contains results which must not be destroyed between
c            calls of rfftf or rfftb.
c
c******************************************************************
      implicit real*4 (a-h,o-z)
      dimension       r(*)       ,wsave(*)
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end




      subroutine rfftf1 (n,c,ch,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)

      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
c        write(6,*) k1,ip,l1,l2,ido
         if (ip .eq. 4) then
            ix2 = iw+ido
            ix3 = ix2+ido
            if (na .eq. 0) then
               call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            else
               call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
            endif

         else if (ip.eq.2) then
            if (na .eq. 0) then
               call radf2 (ido,l1,c,ch,wa(iw))
            else
               call radf2 (ido,l1,ch,c,wa(iw))
            endif

         else if (ip.eq.3) then
            ix2 = iw+ido
            if (na .eq. 0) then
               call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
            else
               call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
            endif

         else
            write(6,*) 'error in factors rfftf',ip
            stop

         endif

         l2 = l1

111   continue

      if (na .eq. 0) then
         do 112 i=1,n
112         c(i) = ch(i)
      endif

      end




      subroutine radf2 (ido,l1,cc,ch,wa1)
      implicit real*4 (a-h,o-z)
      dimension       ch(ido,2,l1)           ,cc(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
      end




      subroutine radf3 (ido,l1,cc,ch,wa1,wa2)
      implicit real*4 (a-h,o-z)
      dimension       ch(ido,3,l1)           ,cc(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      save taur,taui

      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
      return
      end




      subroutine radf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,l1,4)           ,ch(ido,4,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data hsqt2 /.7071067811865475/
      save hsqt2

      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
      end



      subroutine rfftb (n,r,wsave)
c******************************************************************
c
c    subroutine rfftb(n,r,wsave)
c
c    n must be even and have only 2 & 3 as multiples
c
c           adapted by j.jimenez, dec/87.
c
c
c******************************************************************
c
c    subroutine rfftb computes the real perodic sequence from its
c    fourier coefficients (fourier synthesis). the transform is defined
c    below at output parameter r.
c
c    input parameters
c
c    n       the length of the array r to be transformed.  the method
c            is most efficient when n is a product of small primes.
c            n may change so long as different work arrays are provided
c
c    r       a real array of length n which contains the sequence
c            to be transformed
c
c    wsave   a work array which must be dimensioned at least 2*n+15.
c            in the program that calls rfftb. the wsave array must be
c            initialized by calling subroutine rffti(n,wsave) and a
c            different wsave array must be used for each different
c            value of n. this initialization does not have to be
c            repeated so long as n remains unchanged thus subsequent
c            transforms can be obtained faster than the first.
c            the same wsave array can be used by rfftf and rfftb.
c
c
c    output parameters
c
c    r       for n even and for i = 1,...,n
c
c                 r(i) = r(1)+(-1)**(i-1)*r(n)
c
c                      plus the sum from k=2 to k=n/2 of
c
c                       2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                      -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c
c     *****  note
c                 this transform is unnormalized since a call of rfftf
c                 followed by a call of rfftb will multiply the input
c                 sequence by n.
c
c**********************************************************************

      implicit real*4 (a-h,o-z)
      dimension r(*), wsave(*)
      if (n .eq. 1) return
      call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end




      subroutine rfftb1 (n,c,ch,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension  ch(*),c(*),wa(*),ifac(*)
      nf = ifac(2)
      na = 1
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         na = 1-na

         if (ip .eq. 4) then
            ix2 = iw+ido
            ix3 = ix2+ido
            if (na .eq. 0) then
               call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            else
               call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
            endif

         else if (ip .eq. 2) then
            if (na .eq. 0) then
               call radb2 (ido,l1,c,ch,wa(iw))
            else
               call radb2 (ido,l1,ch,c,wa(iw))
            endif

         else if (ip .eq. 3) then
            ix2 = iw+ido
            if (na .eq. 0) then
               call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
            else
               call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
            endif

         else
            write(6,*) 'error in factors rfftf',ip
            stop

         endif

         l1 = l2
         iw = iw+(ip-1)*ido
  116 continue

      if (na .eq. 0) then
         do 117 i=1,n
117         c(i) = ch(i)
      endif

      end




      subroutine radb2 (ido,l1,cc,ch,wa1)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
      end




      subroutine radb3 (ido,l1,cc,ch,wa1,wa2)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      save taur,taui

      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
      return
      end




      subroutine radb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data sqrt2 /1.414213562373095/
      save sqrt2

      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
      end
