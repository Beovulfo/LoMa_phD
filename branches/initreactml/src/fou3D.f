!=======================================================================!
!                                                                       !
!                 FAST TRANSFORMS INTERFACE PACKAGE                     !
!                                                                       !
!=======================================================================!

!***********************************************************************!
!                                                                       !
!    Subrotine Fourx                                                    !
!                                                                       !
!    makes a 1-dimensional real to complex fourier transform            !
!      from F   -Phys-F   (array fou)                                   !
!  to/ from Phys-Phys-F   (array phys)                                  !
!                                                                       !
!       iopt >=0    ===>  inverse transforms( fou ---> fis)             !
!       iopt < 0    ===>  direct  transforms( fis ---> fou)             !
!                                                                       !
!       does the expanding and compacting itself                        !
!                                                                       !
!  Remarks:                                                             !
!    fis ---> fou : supposes high frecuency modes are not               !
!                   needed (dealiasing) and throws them away            !
!    out in phys in real*8                                              !
!                                                                       !
!                                                                       !
!***********************************************************************!


        subroutine fourx(phys,phys8,iopt)

        implicit none
        include"ctes3D"

        integer iopt,j

        real(4) :: phys(mgalx+2)
        real(8) :: phys8(mgalx+2)

! ---------------------- Program ---------------------------------------  

        if (iopt.lt.0) then !  physical to fourier transforms 
        
          call rft(phys,mgalx+2,1,-1)          
          
        else                !  fourier to physical transforms !

          do j=mx+1,mgalx+2
             phys(j)=0.
          enddo

          call rft(phys,mgalx+2,1,1)

          do j=1,mgalx+2
             phys8(j)=phys(j)
          enddo

        endif


        end


!***********************************************************************!
!                                                                       !
!    Subroutine fourz                                                   !
!                                                                       !
!    makes a 1-dimensional complex to complex fourier transform         ! 
!                                                                       !
!       iopt >=0    ===>  inverse transforms( fou ---> fis)             !
!       iopt < 0    ===>  direct  transforms( fis ---> fou)             !
!                                                                       !
!       does the expanding and compacting itself                        !
!                                                                       !
!  NOTE:                                                                !
!    fis ---> fou : supposes high frecuency modes are not               !
!                   needed (dealiasing) and throws them away            !
!                                                                       !
!***********************************************************************!

      subroutine fourz(phys,iopt)

      implicit none
      include"ctes3D"

      integer iopt
      real(4) :: phys(2*my*mgalz)

      if (iopt.lt.0) then  !   physical to fourier transforms  
         call cft(phys,2,2*mgalz,my,-1)
      else                !    fourier to physical transforms  
         call cft(phys,2,2*mgalz,my,1)
      endif


       end

      
!     ---------------------------------------------------------
!     ------------ LOCAL TRANSPOSES        --------------------
!     ---------------------------------------------------------


      subroutine localyz2zy(fieldyz,fieldzy,chwk)
      use point
      
      implicit none
      include "ctes3D"
      
      complex*8:: fieldyz(0:my-1,0:mgalz-1,pb:pe) ! 0:mz1
      complex*8:: fieldzy(0:mgalz-1,0:my-1 ,pb:pe)
      complex*8:: chwk(0:mz-1,0:my-1)
      
      complex*8 zero 
      integer i,j,k,kk,jj,blocking  
          
      zero =(0.,0.)
      kk = mz1 - mgalz1
      
      blocking = 32
      
      do i=pb,pe
      
         do jj=0,mz1,blocking
            do j=0,my1
               do k=jj,min(jj+blocking-1,mz1)
                  chwk(k,j) = fieldyz(j,k,i)
               enddo
            enddo
         enddo               
         
         do j=0,my1
            do k=0,nz1
               fieldzy(k,j,i) = chwk(k,j)
            enddo
            do k=nz2,mgalz1
               fieldzy(k,j,i) = chwk(k+kk,j)
            enddo
            do k=nz1+1,nz2-1
               fieldzy(k,j,i) = zero
            enddo
         enddo

      enddo
       
      end subroutine localyz2zy
      
!     ---------------------------------------------------------

      subroutine localzy2yz(fieldyz,fieldzy,chwk)
      use point
      
      implicit none
      include "ctes3D"
      
      complex*8 :: fieldyz(0:my-1,0:mgalz-1,pb:pe)   ! 
      complex*8 :: fieldzy(0:mgalz-1,0:my-1,pb:pe)   !
      complex*8 :: chwk(0:mz-1,0:my-1)
      
      integer i,j,k,kk,jj,blocking

      kk = mz1-mgalz1 
      blocking = 32
      
      do i=pb,pe
         
         do j=0,my1
            do k=0,nz1
               chwk(k,j) = fieldzy(k,j,i)
            enddo
         enddo
         
         do j=0,my1
            do k=nz2,mgalz1
               chwk(k+kk,j) = fieldzy(k,j,i)
            enddo
         enddo
         
         do jj=0,my1,blocking
            do k=0,mz1
               do j=jj,min(jj+blocking-1,my1)
                  fieldyz(j,k,i) = chwk(k,j)
               enddo
            enddo
         enddo

    
      enddo
      
      end subroutine localzy2yz
