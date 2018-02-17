C rft precision simple
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

      parameter(nmaxt=2048*8)
      implicit real*8(a-h,o-z)

      common /four1/ wsaf1(2*nmaxt+18)

      wsaf1(1)=mgal
      wsaf1(2)=1d0/mgal
      call dffti(mgal,wsaf1(3))
      !call rffti(mgal,wsaf1(3))

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
      implicit real*8 (a-h,o-z)
      dimension c(isa,*)

      common /four1/ aux(3)

      nn = aux(1)
      dn1= aux(2)


c*********************************************************************/
c*            forward transform  (to fourier)                        */
c*********************************************************************/
      if (iopt.lt.0) then
         do 106 j=1,m
106         call dfftf(nn,c(1,j),aux(3))
!106         call rfftf(nn,c(1,j),aux(3))

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
	   
         do j=1,m
            do i=2,nn
                 c(i,j)=c(i+1,j)
	        enddo
	   enddo

         do 107 j=1,m
107         call dfftb(nn,c(1,j),aux(3))
!107         call rfftb(nn,c(1,j),aux(3))
      endif

      end


