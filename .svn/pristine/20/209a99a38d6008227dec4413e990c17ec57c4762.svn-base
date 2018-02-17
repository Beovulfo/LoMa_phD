! ----------------------------------------------------------------------!
!                                                                       !
!    Hace las cuentas en y para planos kz=cte                           !
!                                                                       !
! ----------------------------------------------------------------------!

      subroutine uvwyzx(vor,phi,psi,scal,mfz,u,v,w,dudy,dvdy,mfzbis,
     .    dwdy, dscal,dmfz,wk,nvecy,i)
      use wave
      use combustion
      use fis

      implicit none
      include "ctes3D"
      
      ! ----------------------- constants -------------------------------!
      integer,intent(in)::i
     
      complex*8,intent(in),dimension(my,0:mz1,0:nvecy-1):: vor,phi,psi,
     .                                                     scal,mfz
      complex*8,intent(out),dimension(my,0:mz1,0:nvecy-1)::
     &          u,v,w,dudy,dvdy,dwdy,dscal,dmfz,mfzbis
        
      !------------------------Aux arrays    ----------------------------!
      integer k,j,k1,ind,ii,ip,nvecy
      complex*8 hgg (my),hvv (my),dyvor(my),o1(my),o3(my),
     .          dpsidy(my),dTdy(my),dZdy(my)
      real*8 rK
      complex*8 xa,xb
      complex*8 tempa,tempb 
      real*4 wk(0:mx1,0:mz1)
      !----------------------------------------------------------------------


      ! ------------------------ ! Computes ru,rv,rw ! -- ------------------!
       
      do ind=0,nvecy-1
         ip = i+ind
         xa = xalp(ip-1)
         do k=0,mz1
            xb = xbet(k)
            k1 = icx(k)
            rK = bet2(k1)+alp2(ip-1)

            !1) Solve Lap(phi)=my
            call Lapvdv(phi(1,k,ind),hgg,hvv,rK) !hgg=my, hvv=d(my)/dy
            !2) Calculating d(PSI)/dy
            call deryr2(psi(1,k,ind),dpsidy,1)
            !3) Calculating d(scal)/dy
                !First need to compute T
            !do j=1,my
            !    scal(j,k,ind)=Temper(scal(j,k,ind),mfz(j,k,ind))
            !enddo
            call deryr2(scal(1,k,ind),dTdy,1)
            call deryr2(mfz(1,k,ind),dZdy,1)
            !Compute rhou,rhov,rhow
            !xbet(k)/rK, xalp(ip-1)/rK-> to solve mx,mz
            do j=1,my
                !calculate mx = (ikx*dmydy-ikz*vor)/rk
                !rhou = mx(j,k,ind) + dOmega/dx
                u  (j,k,ind) = (xa*hvv(j)-
     .          xb*vor(j,k,ind))*wk(ip-1,k)
     .          + xa*psi(j,k,ind) 
                !mz =(ikz*dmdy+ikx*vor)/rk
                !rhow = mz + dOmega/dz
                w  (j,k,ind) = (xb*hvv(j)+
     .          xa*vor(j,k,ind))*wk(ip-1,k)
     .          + xb*psi(j,k,ind) 
                !rhov = my + dOmega/dy
                v  (j,k,ind) = hgg(j) + dpsidy(j)   !v keeps rhov
                !dT/dy
                dscal   (j,k,ind) = dTdy(j) !dscal = dT/dy
                dmfz   (j,k,ind) = dZdy(j) !dscal = dT/dy
            enddo
         enddo
       enddo
       !Loop again for drhoU/dy, drhoW/dy
       do ind=0,nvecy-1
          ip = i+ind
          do k=0,mz1
          !Calculating derivatives of rhoU,rhoV
             call deryr2(u(1,k,ind),o1,1)
             call deryr2(w(1,k,ind),o3,1)
             do j=1,my
                !rhov = my + d(psi)/dy
                dudy(j,k,ind) = o1(j)    !drhou/dy 
                dwdy(j,k,ind) = o3(j)    !drhow/dy 
                dvdy(j,k,ind) = scal(j,k,ind) !T 
                mfzbis(j,k,ind) = mfz(j,k,ind) !T 
             enddo
         enddo
      enddo

      endsubroutine uvwyzx
      
      
