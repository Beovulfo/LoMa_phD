C  cft precision simple
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
      implicit double precision (a-h,o-z)
      dimension c(*)

c******* coefficients and savearea for fourier   expansions **********/
      common /fouc1/ aux(3)

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

            call zfftf(nn,aux(3),aux(4+n2))
            !call cfftf(nn,aux(3),aux(4+n2))
c                            *********  copy to original ********
            j0 = j00
            do 100 j=1,n2,2
               c(j0)   = dn1*aux(2+j)
               c(j0+1) = dn1*aux(3+j)
100            j0 = j0+ise

c-------   inverse transform  (to physical) --------------
         else

            call zfftb(nn,aux(3),aux(4+n2))
            !call cfftb(nn,aux(3),aux(4+n2))
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

      parameter (nmaxt=2048*2)

c*********************************************************************/
      double precision wsac1,xx

c******* coefficients and savearea for fourier   expansions **********/
      common /fouc1/ wsac1(10*nmaxt+19)


c***********************  transforms (ngal) **************************/
      wsac1(1)=ngal
      wsac1(2)=1d0/wsac1(1)
      !call cffti(ngal,wsac1(2*ngal+4))
      call zffti(ngal,wsac1(2*ngal+4))


      end


