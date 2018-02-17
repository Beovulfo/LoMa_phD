      subroutine cft(c,ise,isa,m,iopt)
c*********************************************************************/
c*   given  m    ==> complex vectors cm(0:n-1), stored in c as
c*          ise  ==> stride between elements in each vectors
c*          isa  ==> stride between initial elements of vectors
c*
c*   computes fourier expansion                                      */
c*
c*   cm(j) = sum from j=0 to n-1 of cm(k)*exp(i*j*k*2*pi/n )
c*
c*   (or viceversa, normalising factor is taken care of)
c*         iopt >= 0  ==> inverse transform (x fourier --> c phys.)  */
c*         iopt <  0  ==> direct  transform (x physic. --> c four.)  */
c*             jimenez/    may 90                                    */
c*********************************************************************/
      implicit real*4 (a-h,o-z)
      dimension c(*)

c******* coefficients and savearea for fourier   expansions **********/
      common /fouc1/ aux(3)
      save /fouc1/

      nn = aux(1)
      dn1= aux(2)
      n2 = 2*nn

      j00 = 1

      do 105 k=1,m
c                            *********  copy to contiguous ******
         j0 = j00
         do 110 j=1,n2,2
            aux(2+j) = c(j0)
            aux(3+j) = c(j0+1)
110         j0 = j0+ise

c-------   forward transform  (to fourier)  --------------
         if (iopt.lt.0) then

            call cfftf(nn,aux(3),aux(4+n2))
c                            *********  copy to original ********
            j0 = j00
            do 100 j=1,n2,2
               c(j0)   = dn1*aux(2+j)
               c(j0+1) = dn1*aux(3+j)
100            j0 = j0+ise

c-------   inverse transform  (to physical) --------------
         else

            call cfftb(nn,aux(3),aux(4+n2))
c                            *********  copy to original ********
            j0 = j00
            do 102 j=1,n2,2
               c(j0)   = aux(2+j)
               c(j0+1) = aux(3+j)
102            j0 = j0+ise

         endif

         j00 = j00+isa

105   continue

      end

      subroutine cfti(ngal)
c*********************************************************************/
c*     initialises coefficients for complex fourier derivatives      */
c*       and transforms.                                             */
c*                                                                   */
c*      transforms to be of size:  complex ngal                      */
c*      assumes:                                                     */
c*           ngal <= nmaxt                                           */
c*                                                                   */
c*********************************************************************/

      parameter (nmaxt=4096)

c*********************************************************************/
      complex*8 fder
      real*4 wsac1,xx

c******* coefficients and savearea for fourier   expansions **********/
      common /fouc1/ wsac1(10*nmaxt+19)
      save /fouc1/


c***********************  transforms (ngal) **************************/
      wsac1(1)=ngal
      wsac1(2)=1e0/wsac1(1)
      call cffti(ngal,wsac1(2*ngal+4))


      end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c*      from here down is cffti.packet.ncar adapted to real*4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine cffti (n,wsave)
c***********************************************************************
c
c     subroutine cffti(n,wsave)
c
c***********************************************************************
c
c     subroutine cffti initializes the array wsave which is used in
c     both cfftf and cfftb. the prime factorization of n together with
c     a tabulation of the trigonometric functions are computed and
c     stored in wsave.
c
c     input parameter
c
c     n       the length of the sequence to be transformed
c
c     output parameter
c
c     wsave   a work array which must be dimensioned at least 4*n+15
c             the same work array can be used for both cfftf and cfftb
c             as long as n remains unchanged. different wsave arrays
c             are required for different values of n. the contents of
c             wsave must not be changed between calls of cfftf or cfftb.
c
c***********************************************************************

      implicit real*4 (a-h,o-z)
      dimension       wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cffti1 (n,wsave(iw1),wsave(iw2))
      return
      end




      subroutine cffti1 (n,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension       wa(*)      ,ifac(*)    ,ntryh(4)
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
      save ntryh

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
      tpi = 8.0*atan(1.0)
      arg1 = tpi/float(n)
      dc = cos(arg1)
      ds = sin(arg1)
      wa(1) = dc
      wa(2) = ds
      nt = n+n
      do 108 i=4,nt,2
         wa(i-1) = dc*wa(i-3)-ds*wa(i-2)
         wa(i) = ds*wa(i-3)+dc*wa(i-2)
  108 continue
      return
      end




      subroutine cfftb (n,c,wsave)
c***********************************************************************
c
c     subroutine cfftb(n,c,wsave)
c
c***********************************************************************
c
c     subroutine cfftb computes the backward complex discrete fourier
c     transform (the fourier synthesis). equivalently , cfftb computes
c     a complex periodic sequence from its fourier coefficients.
c     the transform is defined below at output parameter c.
c
c     a call of cfftf followed by a call of cfftb will multiply the
c     sequence by n.
c
c     the array wsave which is used by subroutine cfftb must be
c     initialized by calling subroutine cffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes.
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls cfftb. the wsave array must be
c             initialized by calling subroutine cffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by cfftf and cfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(i*(j-1)*(k-1)*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine cfftf or cfftb
c
c***********************************************************************
      implicit real*4 (a-h,o-z)
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end




      subroutine cfftb1 (n,c,ch,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
      nf = ifac(2)
      l1 = 1
      do 105 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 101
         ix2 = l1+l1
         ix3 = ix2+l1
         call passb4 (idot,l1,idl1,ix2,ix3,c,c,c,ch,ch,wa,wa,wa)
         go to 104
  101    if (ip .ne. 2) go to 102
         call passb2 (idot,l1,idl1,c,c,c,ch,ch,wa)
         go to 104
  102    if (ip .ne. 3) go to 103
         ix2 = l1+l1
         call passb3 (idot,l1,idl1,ix2,c,c,c,ch,ch,wa,wa)
         go to 104
  103    call passb (idot,ip,l1,idl1,c,c,c,ch,ch,wa)
  104    l1 = l2
  105 continue
      return
      end




      subroutine passb2 (ido,l1,idl1,cc,c1,c2,ch,ch2,wa1)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,2,l1)           ,c1(ido,l1,2)           ,
     1                c2(idl1,2) ,ch(ido,l1,2)           ,ch2(idl1,2),
     2                wa1(l1,1)
      idot = ido/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ch(i,k,2) = cc(i,1,k)-cc(i,2,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ch(i,k,2) = cc(i,1,k)-cc(i,2,k)
  104    continue
  105 continue
  106 do 107 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  107 continue
      do 108 k=1,l1
         c1(1,k,2) = ch(1,k,2)
         c1(2,k,2) = ch(2,k,2)
  108 continue
      if (ido .eq. 2) return
      if (idot .lt. l1) go to 111
      do 110 k=1,l1
         do 109 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
  109    continue
  110 continue
      return
  111 do 113 i=4,ido,2
         do 112 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
  112    continue
  113 continue
      return
      end



      subroutine passb3 (ido,l1,idl1,ix2,cc,c1,c2,ch,ch2,wa1,wa2)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,3,l1)           ,c1(ido,l1,3)           ,
     1                c2(idl1,3) ,ch(ido,l1,3)           ,ch2(idl1,3),
     2                wa1(l1,1)  ,wa2(ix2,1)
      data taur,taui /-.5,.866025403784439/
      save taur,taui

      idot = ido/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
            ch(i,k,2) = cc(i,2,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)-cc(i,3,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
            ch(i,k,2) = cc(i,2,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)-cc(i,3,k)
  104    continue
  105 continue
c
  106 do 107 ik=1,idl1
         c2(ik,1) = ch2(ik,1)+ch2(ik,2)
         c2(ik,2) = ch2(ik,1)+taur*ch2(ik,2)
         c2(ik,3) = taui*ch2(ik,3)
  107 continue
      do 108 ik=2,idl1,2
         ch2(ik-1,2) = c2(ik-1,2)-c2(ik,3)
         ch2(ik-1,3) = c2(ik-1,2)+c2(ik,3)
  108 continue
      do 109 ik=2,idl1,2
         ch2(ik,2) = c2(ik,2)+c2(ik-1,3)
         ch2(ik,3) = c2(ik,2)-c2(ik-1,3)
  109 continue
      do 111 j=2,3
         do 110 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  110    continue
  111 continue
      if (ido .eq. 2) return
      if (idot-1 .lt. l1) go to 114
      do 113 k=1,l1
         do 112 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)+
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)-
     1                    wa2(ix2,i-2)*ch(i,k,3)
  112    continue
  113 continue
      return
  114 do 116 i=4,ido,2
         do 115 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)+
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)-
     1                    wa2(ix2,i-2)*ch(i,k,3)
  115    continue
  116 continue
      return
      end




      subroutine passb4 (ido,l1,idl1,ix2,ix3,cc,c1,c2,ch,ch2,wa1,wa2,
     1                   wa3)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,4,l1)           ,c1(ido,l1,4)           ,
     1                c2(idl1,4) ,ch(ido,l1,4)           ,ch2(idl1,4),
     2                wa1(l1,1)  ,wa2(ix2,1) ,wa3(ix3,1)
      idot = ido/2
c
      if (ido .lt. l1) go to 106
      do 103 k=1,l1
         do 101 i=2,ido,2
            ch(i-1,k,4) = cc(i,4,k)-cc(i,2,k)
  101    continue
         do 102 i=2,ido,2
            ch(i,k,4) = cc(i-1,2,k)-cc(i-1,4,k)
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,2) = cc(i,1,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)+cc(i,4,k)
            ch(i,k,1) = cc(i,1,k)-cc(i,3,k)
  104    continue
  105 continue
      go to 111
  106 do 108 i=2,ido,2
         do 107 k=1,l1
            ch(i-1,k,4) = cc(i,4,k)-cc(i,2,k)
            ch(i,k,4) = cc(i-1,2,k)-cc(i-1,4,k)
  107    continue
  108 continue
      do 110 i=1,ido
         do 109 k=1,l1
            ch(i,k,2) = cc(i,1,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)+cc(i,4,k)
            ch(i,k,1) = cc(i,1,k)-cc(i,3,k)
  109    continue
  110 continue
  111 do 112 ik=1,idl1
         c2(ik,1) = ch2(ik,2)+ch2(ik,3)
  112 continue
      do 113 ik=1,idl1
         ch2(ik,3) = ch2(ik,2)-ch2(ik,3)
  113 continue
      do 114 ik=1,idl1
         ch2(ik,2) = ch2(ik,1)+ch2(ik,4)
  114 continue
      do 115 ik=1,idl1
         ch2(ik,4) = ch2(ik,1)-ch2(ik,4)
  115 continue
      do 117 j=2,4
         do 116 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  116    continue
  117 continue
      if (ido .eq. 2) return
      if (idot .lt. l1) go to 120
      do 119 k=1,l1
         do 118 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)+
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)-
     1                    wa2(ix2,i-2)*ch(i,k,3)
            c1(i,k,4) = wa3(ix3-1,i-2)*ch(i,k,4)+
     1                  wa3(ix3,i-2)*ch(i-1,k,4)
            c1(i-1,k,4) = wa3(ix3-1,i-2)*ch(i-1,k,4)-
     1                    wa3(ix3,i-2)*ch(i,k,4)
  118    continue
  119 continue
      return
  120 do 122 i=4,ido,2
         do 121 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)+wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)-
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)+
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)-
     1                    wa2(ix2,i-2)*ch(i,k,3)
            c1(i,k,4) = wa3(ix3-1,i-2)*ch(i,k,4)+
     1                  wa3(ix3,i-2)*ch(i-1,k,4)
            c1(i-1,k,4) = wa3(ix3-1,i-2)*ch(i-1,k,4)-
     1                    wa3(ix3,i-2)*ch(i,k,4)
  121    continue
  122 continue
      return
      end




      subroutine passb (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit real*4 (a-h,o-z)
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      l1t = l1+l1
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 do 113 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  113 continue
      idj = 0
      do 115 j=2,ipph
         jc = ipp2-j
         idj = idj+idl1
         do 114 ik=1,idl1
            c2(ik,j) = ch2(ik,1)+wa(idj-1)*ch2(ik,2)
            c2(ik,jc) = wa(idj)*ch2(ik,ip)
  114    continue
  115 continue
      do 117 j=2,ipph
         do 116 ik=1,idl1
            c2(ik,1) = c2(ik,1)+ch2(ik,j)
  116    continue
  117 continue
c
      idl = 0
      do 120 l=2,ipph
         lc = ipp2-l
         idl = idl+idl1
         idlj = idl
         do 119 j=3,ipph
            jc = ipp2-j
            idlj = idlj+idl
            if (idlj .gt. nt) idlj = idlj-nt
            war = wa(idlj-1)
            wai = wa(idlj)
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
c
      do 122 j=2,ipph
         jc = ipp2-j
         do 121 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 ik=2,idl1,2
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  123    continue
  124 continue
c
      do 126 j=2,ip
         do 125 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  125    continue
  126 continue
      if (ido .eq. 2) return
      idj = 0
      if (idot .gt. l1) go to 130
      do 129 j=2,ip
         idj = idj+l1t
         idij = 0
         do 128 i=4,ido,2
            idij = idij+idj
            do 127 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  127       continue
  128    continue
  129 continue
      return
  130 do 134 j=2,ip
         idj = idj+l1t
         do 133 k=1,l1
            idij = 0
            do 131 i=4,ido,2
               idij = idij+idj
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
  131       continue
            idij = 0
            do 132 i=4,ido,2
               idij = idij+idj
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  132       continue
  133    continue
  134 continue
      return
      end




      subroutine cfftf (n,c,wsave)
c***********************************************************************
c
c     subroutine cfftf(n,c,wsave)
c
c***********************************************************************
c
c     subroutine cfftf computes the forward complex discrete fourier
c     transform (the fourier analysis). equivalently , cfftf computes
c     the fourier coefficients of a complex periodic sequence.
c     the transform is defined below at output parameter c.
c
c     the transform is not normalized. to obtain a normalized transform
c     the output must be divided by n. otherwise a call of cfftf
c     followed by a call of cfftb will multiply the sequence by n.
c
c     the array wsave which is used by subroutine cfftf must be
c     initialized by calling subroutine cffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes. n
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls cfftf. the wsave array must be
c             initialized by calling subroutine cffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by cfftf and cfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(-i*j*k*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine cfftf or cfftb
c
c***********************************************************************
      implicit real*4 (a-h,o-z)
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end




      subroutine cfftf1 (n,c,ch,wa,ifac)
      implicit real*4 (a-h,o-z)
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
      nf = ifac(2)
      l1 = 1
      do 105 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 101
         ix2 = l1+l1
         ix3 = ix2+l1
         call passf4 (idot,l1,idl1,ix2,ix3,c,c,c,ch,ch,wa,wa,wa)
         go to 104
  101    if (ip .ne. 2) go to 102
         call passf2 (idot,l1,idl1,c,c,c,ch,ch,wa)
         go to 104
  102    if (ip .ne. 3) go to 103
         ix2 = l1+l1
         call passf3 (idot,l1,idl1,ix2,c,c,c,ch,ch,wa,wa)
         go to 104
  103    call passf (idot,ip,l1,idl1,c,c,c,ch,ch,wa)
  104    l1 = l2
  105 continue
      return
      end




      subroutine passf2 (ido,l1,idl1,cc,c1,c2,ch,ch2,wa1)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,2,l1)           ,c1(ido,l1,2)           ,
     1                c2(idl1,2) ,ch(ido,l1,2)           ,ch2(idl1,2),
     2                wa1(l1,1)
      idot = ido/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ch(i,k,2) = cc(i,1,k)-cc(i,2,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ch(i,k,2) = cc(i,1,k)-cc(i,2,k)
  104    continue
  105 continue
  106 do 107 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  107 continue
      do 108 k=1,l1
         c1(1,k,2) = ch(1,k,2)
         c1(2,k,2) = ch(2,k,2)
  108 continue
      if (ido .eq. 2) return
      if (idot .lt. l1) go to 111
      do 110 k=1,l1
         do 109 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
  109    continue
  110 continue
      return
  111 do 113 i=4,ido,2
         do 112 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
  112    continue
  113 continue
      return
      end




      subroutine passf3 (ido,l1,idl1,ix2,cc,c1,c2,ch,ch2,wa1,wa2)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,3,l1)           ,c1(ido,l1,3)           ,
     1                c2(idl1,3) ,ch(ido,l1,3)           ,ch2(idl1,3),
     2                wa1(l1,1)  ,wa2(ix2,1)
      data taur,taui /-.5,-.866025403784439/
      save taur,taui

      idot = ido/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
            ch(i,k,2) = cc(i,2,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)-cc(i,3,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
            ch(i,k,2) = cc(i,2,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)-cc(i,3,k)
  104    continue
  105 continue
c
  106 do 107 ik=1,idl1
         c2(ik,1) = ch2(ik,1)+ch2(ik,2)
         c2(ik,2) = ch2(ik,1)+taur*ch2(ik,2)
         c2(ik,3) = taui*ch2(ik,3)
  107 continue
      do 108 ik=2,idl1,2
         ch2(ik-1,2) = c2(ik-1,2)-c2(ik,3)
         ch2(ik-1,3) = c2(ik-1,2)+c2(ik,3)
  108 continue
      do 109 ik=2,idl1,2
         ch2(ik,2) = c2(ik,2)+c2(ik-1,3)
         ch2(ik,3) = c2(ik,2)-c2(ik-1,3)
  109 continue
      do 111 j=2,3
         do 110 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  110    continue
  111 continue
      if (ido .eq. 2) return
      if (idot-1 .lt. l1) go to 114
      do 113 k=1,l1
         do 112 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)-
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)+
     1                    wa2(ix2,i-2)*ch(i,k,3)
  112    continue
  113 continue
      return
  114 do 116 i=4,ido,2
         do 115 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)-
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)+
     1                    wa2(ix2,i-2)*ch(i,k,3)
  115    continue
  116 continue
      return
      end




      subroutine passf4 (ido,l1,idl1,ix2,ix3,cc,c1,c2,ch,ch2,wa1,wa2,
     1                   wa3)
      implicit real*4 (a-h,o-z)
      dimension       cc(ido,4,l1)           ,c1(ido,l1,4)           ,
     1                c2(idl1,4) ,ch(ido,l1,4)           ,ch2(idl1,4),
     2                wa1(l1,1)  ,wa2(ix2,1) ,wa3(ix3,1)
      idot = ido/2
c
      if (ido .lt. l1) go to 106
      do 103 k=1,l1
         do 101 i=2,ido,2
            ch(i-1,k,4) = cc(i,2,k)-cc(i,4,k)
  101    continue
         do 102 i=2,ido,2
            ch(i,k,4) = cc(i-1,4,k)-cc(i-1,2,k)
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,2) = cc(i,1,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)+cc(i,4,k)
            ch(i,k,1) = cc(i,1,k)-cc(i,3,k)
  104    continue
  105 continue
      go to 111
  106 do 108 i=2,ido,2
         do 107 k=1,l1
            ch(i,k,4) = cc(i-1,4,k)-cc(i-1,2,k)
            ch(i-1,k,4) = cc(i,2,k)-cc(i,4,k)
  107    continue
  108 continue
      do 110 i=1,ido
         do 109 k=1,l1
            ch(i,k,2) = cc(i,1,k)+cc(i,3,k)
            ch(i,k,3) = cc(i,2,k)+cc(i,4,k)
            ch(i,k,1) = cc(i,1,k)-cc(i,3,k)
  109    continue
  110 continue
  111 do 112 ik=1,idl1
         c2(ik,1) = ch2(ik,2)+ch2(ik,3)
  112 continue
      do 113 ik=1,idl1
         ch2(ik,3) = ch2(ik,2)-ch2(ik,3)
  113 continue
      do 114 ik=1,idl1
         ch2(ik,2) = ch2(ik,1)+ch2(ik,4)
  114 continue
      do 115 ik=1,idl1
         ch2(ik,4) = ch2(ik,1)-ch2(ik,4)
  115 continue
      do 117 j=2,4
         do 116 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  116    continue
  117 continue
      if (ido .eq. 2) return
      if (idot .lt. l1) go to 120
      do 119 k=1,l1
         do 118 i=4,ido,2
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)-
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)+
     1                    wa2(ix2,i-2)*ch(i,k,3)
            c1(i,k,4) = wa3(ix3-1,i-2)*ch(i,k,4)-
     1                  wa3(ix3,i-2)*ch(i-1,k,4)
            c1(i-1,k,4) = wa3(ix3-1,i-2)*ch(i-1,k,4)+
     1                    wa3(ix3,i-2)*ch(i,k,4)
  118    continue
  119 continue
      return
  120 do 122 i=4,ido,2
         do 121 k=1,l1
            c1(i,k,2) = wa1(l1-1,i-2)*ch(i,k,2)-wa1(l1,i-2)*ch(i-1,k,2)
            c1(i-1,k,2) = wa1(l1-1,i-2)*ch(i-1,k,2)+
     1                    wa1(l1,i-2)*ch(i,k,2)
            c1(i,k,3) = wa2(ix2-1,i-2)*ch(i,k,3)-
     1                  wa2(ix2,i-2)*ch(i-1,k,3)
            c1(i-1,k,3) = wa2(ix2-1,i-2)*ch(i-1,k,3)+
     1                    wa2(ix2,i-2)*ch(i,k,3)
            c1(i,k,4) = wa3(ix3-1,i-2)*ch(i,k,4)-
     1                  wa3(ix3,i-2)*ch(i-1,k,4)
            c1(i-1,k,4) = wa3(ix3-1,i-2)*ch(i-1,k,4)+
     1                    wa3(ix3,i-2)*ch(i,k,4)
  121    continue
  122 continue
      return
      end




      subroutine passf (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit real*4 (a-h,o-z)
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      l1t = l1+l1
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 do 113 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  113 continue
      idj = 0
      do 115 j=2,ipph
         jc = ipp2-j
         idj = idj+idl1
         do 114 ik=1,idl1
            c2(ik,j) = ch2(ik,1)+wa(idj-1)*ch2(ik,2)
            c2(ik,jc) = -wa(idj)*ch2(ik,ip)
  114    continue
  115 continue
      do 117 j=2,ipph
         do 116 ik=1,idl1
            c2(ik,1) = c2(ik,1)+ch2(ik,j)
  116    continue
  117 continue
c
      idl = 0
      do 120 l=2,ipph
         lc = ipp2-l
         idl = idl+idl1
         idlj = idl
         do 119 j=3,ipph
            jc = ipp2-j
            idlj = idlj+idl
            if (idlj .gt. nt) idlj = idlj-nt
            war = wa(idlj-1)
            wai = wa(idlj)
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
c
      do 122 j=2,ipph
         jc = ipp2-j
         do 121 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 ik=2,idl1,2
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  123    continue
  124 continue
c
      do 126 j=2,ip
         do 125 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  125    continue
  126 continue
      if (ido .eq. 2) return
      idj = 0
      if (idot .gt. l1) go to 130
      do 129 j=2,ip
         idj = idj+l1t
         idij = 0
         do 128 i=4,ido,2
            idij = idij+idj
            do 127 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  127       continue
  128    continue
  129 continue
      return
  130 do 134 j=2,ip
         idj = idj+l1t
         do 133 k=1,l1
            idij = 0
            do 131 i=4,ido,2
               idij = idij+idj
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
  131       continue
            idij = 0
            do 132 i=4,ido,2
               idij = idij+idj
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  132       continue
  133    continue
  134 continue
      return
      end
