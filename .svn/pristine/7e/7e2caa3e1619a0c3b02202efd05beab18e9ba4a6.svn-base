         program main
         use matrices
         use fis
         use ficheros
         implicit none
         include 'ctes3D'
         !.............................!
         integer::i,j,k,opt,z1,z2
         real(8):: fun(my),dfun(my),d2fun(my),funb(my)
         real(8):: dfunb(my),d2funb(my),d3fun(my),d3funb(my)
         complex(8),dimension(mgalz) :: fouc
         real(8),dimension(mgalx) :: four
         complex(8),dimension(mgalz) :: phys
         complex(8),dimension(mgalx) :: physx
         complex(8) :: box
         real(8),allocatable::d11(:,:),d12(:,:),d21(:,:),d22(:,:)
         real(8) :: ky,fref(my)
         integer :: nn
         real(8) :: e1,e2,swin

         !Reading mesh
         iinp=2
         filinp='mesh.txt'
!        READ MESH
         write(*,*) "Reading mesh and mean profiles.."
         open(iinp,file=filinp,status='old')
         read(iinp,*) (y(j),j=1,my)
         rewind(iinp)
         close(iinp)

         
         write(*,*) "Let's check some FFT routines"
         call cfti(mgalz) !complex Fourier Transform init
         call rfti(mgalx) !real Fourier Transform init
         write (*,*) 'precalculo matrices ...'
      
         allocate(d11(my,7),d12(my,7),d21(my,5),d22(my,5))
      
         call derivadas(d11,d12,d21,d22)      
            
         deallocate(d11,d12,d21,d22)
      
         write (*,*) 'FIN precalculo matrices ...'
         swin=8.0
         !Creating fun=sin(y)
         do nn=1,320,20
           ky=nn*2*3.1415/380
           do j=1,my
              fun(j) = exp(-(y(j)/swin)**2)*sin(ky*y(j))
              fref(j) = exp(-(y(j)/swin)**2)*(-2/swin**2*sin(ky*y(j))
     .              + 4*y(j)**2/swin**4*sin(ky*y(j))
     .              - 2*y(j)/swin**2*ky*cos(ky*y(j))
     .              - ky**2*sin(ky*y(j))-2*y(j)/swin**2*ky*cos(ky*y(j)))
           enddo
           call deryr(fun,dfun) 
           call deryr(dfun,d2fun) 

           call deryyr(fun,d2funb) 
           e1 = 0.0
           e2 = 0.0
           do j=1,my
              e1 = e1+abs(d2fun(j)-fref(j))
              e2 = e2+abs(d2funb(j)-fref(j))
           enddo
           write(*,*) "nn=",nn
           write(*,*) "e1=",e1
           write(*,*) "e2=",e2
        enddo

         
  
        
         !phys to fourier transfrom (-1) 
         !------------------------
         !opt = -1
         !four X 
         !four Z 

         !PHYS to CFFT 
         !
         z1 = mgalz/3
         z2 = 2*mgalz/3
         !rectangle function
         do k=1,mgalz
            if ((k.eq.z1).or.(k.eq.z2)) then
                 fouc(k) = 0.5
            elseif ((k.gt.z1).and.(k.lt.z2)) then
                 fouc(k) = 1.0
            else
                fouc(k)=0d0
            endif
         enddo
         !copy
         do k=1,mgalz
            phys(k) = fouc(k)
         enddo

         !Check complex fourier transform
         !-----------------------------------
         !transform to complexfourier
         opt = -1
         call fourz(fouc,opt)

         !Write to file
         open(unit=1,file='cfttest_double.txt')
         write(1,5)(real(phys(k)),real(fouc(k)),
     .                 aimag(fouc(k)),k=1,mgalz)
         close(1)
         !Go back to phys
         opt = 1
         call fourz(fouc,opt)

         !Check that is the same
         box=0.0
         do k=1,mgalz
            box = box+(fouc(k)-phys(k))
         enddo
         write(*,*) "Total error between phys-fou-phys=",box

         !transform to realfourier
         z1 = mgalx/3
         z2 = 2*mgalx/3
         !rectangle function
         do i=1,mgalx
            if ((i.eq.z1).or.(i.eq.z2)) then
                 four(k) = 0.5
            elseif ((i.gt.z1).and.(i.lt.z2)) then
                 four(i) = 1.0
            else
                four(i)=0d0
            endif
         enddo
         !copy
         do i=1,mgalx
            physx(i) = four(i)
         enddo

 
         opt = -1
         call fourx(four,opt)

         !Write to file
         open(unit=2,file='rffttest_double.txt')
         write(2,6)(real(physx(i)),four(i),i=1,mgalx)
         close(2)
         !Go back to phys
         opt = 1
         call fourx(four,opt)

         !Check that is the same
         box=0.0
         do i=1,mgalx
            box = box+(four(i)-physx(i))
         enddo
         write(*,*) "Total error between phys-fou-phys=",box
5        format (5E14.7, 5E14.7, 5E14.7)
6        format (5E14.7, 5E14.7)



          
         end program main




!---!----!----------------------------------------------------

      subroutine fourz(phys,iopt)

      implicit none
      include 'ctes3D'

      integer iopt
      real(8) :: phys(2*my*mgalz)

      if (iopt.lt.0) then  !   physical to fourier transforms  
         call cft(phys,2,2*mgalz,my,-1)
      else                !    fourier to physical transforms  
         call cft(phys,2,2*mgalz,my,1)
      endif


       end
 

        subroutine fourx(phys,iopt)

        implicit none
        include 'ctes3D'

        integer iopt,j
        real(8) :: phys(mgalx+2)

! ---------------------- Program ---------------------------------------  

        if (iopt.lt.0) then !  physical to fourier transforms 
        
          call rft(phys,mgalx+2,1,-1)          
          
        else                !  fourier to physical transforms !

          do j=mx+1,mgalx+2
             phys(j)=0.
          enddo

          call rft(phys,mgalx+2,1,1)

        endif


        end



