!=======================================================================!
!                                                                       !
!                    swap  Routines Package i                           !
! ---  Escribimos trozos de lineas                                      !
! ---  enrique 2003                                                     !
! ---  modificado a mono-variable:  oscar 2003                          !
! ---  Modificado a planos yz    : SHC dec 2005                         !
! ---  getfil Direct Access Routine: Reads pl-pl file in l-p style      !
!                                                                       !
!                                                                       !
!=======================================================================!
!***********************************************************************!
!                                                                       !
!     Escribe ficheros de swap en yz                                    !
!     Opuesto de getswapyz, leido por getswapyz2xz                      !
!                                                                       !
!     SHC dec 2005                                                      !
!     Entradas: var. -> Variable a guardar                              !
!               aux  -> auxiliar, comparte con auxr                     !
!               iplan -> i inicial a guardar                            !
!               nvecy -> Numero de bloques nvecy                        !
!               plx   -> Numero de planos x                             !
!               ivar  -> Variables                                      !
!                                                                       !
!***********************************************************************!
      subroutine swap(var,aux,auxr,iplan,plx,ivar)

      implicit none
      include 'ctes3D'
      
      integer iinp,iswap
      character*50 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(15)
      save /ficheros/
      
      ! ------------------------ variables -----------------------------!
            
      integer j,k,ivar,iplan,jj,irec,plx,wk,nregt
      complex*8  var(my,mz,plx),aux(mz,plx,my)
      real*4 auxr(2*mz,plx,my)
      
      character*70 fnameima
      character*3  ext
      ! ------------------------ Program -------------------------------!

      nregt = (mx1+1)/plx
      irec = (iplan-1)/plx+1
      iswap = 12

      if (iplan==1) then
         fnameima=filswap(1:index(filswap,' ')-1)//extensiones(ivar)
         open(unit=iswap+ivar,file=fnameima,status='unknown',
     &        form='unformatted',access='direct',recl=plx*mz*2)
      endif
         

      do j=1,my
         do jj=1,plx
            do k=1,mz
               aux(k,jj,j) = var(j,k,jj)
            enddo
         enddo
      enddo

      do j=1,my
         write(iswap+ivar,rec=irec+(j-1)*nregt)
     .        ((auxr(k,jj,j),k=1,2*mz),jj=1,plx)
      enddo

      if ((iplan-1+plx)==(mx1+1)) then
!          write(*,*)'he entrado'
         close(iswap+ivar)
      endif

      endsubroutine swap

!***********************************************************************!
!                                                                       !
!     Lee ficheros de swap en yz                                        !
!                                                                       !
!     Entradas: var. -> Variable a guardar                              !
!               aux  -> auxiliar, comparte con auxr                     !
!               iplan -> i inicial a guardar                            !
!               nvecy -> Numero de bloques nvecy                        !
!               plx   -> Numero de planos x                             !
!               ivar  -> Variables                                      !
!     SHC dec 2005                                                      !
!                                                                       !
!***********************************************************************!
      subroutine getswapyz(var,aux,auxr,iplan,plx,ivar)

      implicit none
      include 'ctes3D'
      
      integer iinp,iswap
      character*50 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(15)
      save /ficheros/
      
      ! ------------------------ variables -----------------------------!
            
      integer j,k,ivar,iplan,jj,irec,plx,wk
      complex*8  var(my,mz,plx),aux(mz,plx,my)
      real*4 auxr(2*mz,plx,my)
      
      character*70 fnameima
      character*3  ext 
      ! ------------------------ Program -------------------------------!
      
      iswap = 12
      irec = (iplan-1)/plx+1

      do j=1,my
         write(ext,'(i3.3)') j
         fnameima=filswap(1:index(filswap,' ')-1)//ext//
     &   '/swap.'//extensiones(ivar)

         open(unit=iswap,file=fnameima,status='unknown',
     &   form='unformatted',access='direct',recl=plx*mz*2) ! Problema !!!
         read(iswap,rec=irec)((auxr(k,jj,j),k=1,2*mz),jj=1,plx)
         close(iswap)
         
      enddo

      do jj=1,plx
         do k=1,mz
            do j=1,my
               var(j,k,jj) = aux(k,jj,j)
            enddo
         enddo
      enddo


      endsubroutine getswapyz

      
!***********************************************************************!
!                                                                       !
!     Lee ficheros de swap en x (lineas de tamanyo mx)                  !
!     Entradas: var. -> Variable a guardar                              !
!               aux  -> auxiliar, comparte con auxr                     !
!               iplan -> i inicial a guardar                            !
!               nvecy -> Numero de bloques nvecy                        !
!               plx   -> Numero de planos x                             !
!               ind   -> Variables                                      !
!     SHC dec 2005                                                      !
!                                                                       !
!***********************************************************************!

      subroutine getswapyz2xz(var,aux,auxc,j,ivar)

      implicit none
      include 'ctes3D'
      
      integer iinp,iswap
      character*50 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(15)
      save /ficheros/
      
      ! ------------------------ variables -----------------------------!

      integer,intent(in)::j
      integer i,ivar,jj,k,jp,wk
      real*4 aux(2*mz,mx1+1)
      complex*8 var(mx1+1,mz),auxc(mz,mx1+1)
      character*52 fnameima
      character*3 ext
      
      ! ------------------------ Program -------------------------------!

      iswap = 12
      if (j==1) then
         fnameima=filswap(1:index(filswap,' ')-1)//extensiones(ivar)
         open(unit=iswap+ivar,file=fnameima,status='unknown',
     &           form='unformatted',access='direct',recl=mx*mz)
      endif

      read(iswap+ivar,rec=j)((aux(k,i),k=1,2*mz),i=1,mx1+1)

      do k=1,mz
         do i=1,mx1+1
            var(i,k) = auxc(k,i)
         enddo
      enddo

      if (j==my) then
         close(iswap+ivar)
      endif
    
      
      
      endsubroutine getswapyz2xz

      
!***********************************************************************!
!                                                                       !
!     Escribe ficheros de swap en yz desde lineas de tamanyo mx         !
!                                                                       !
!     SHC dec 2005                                                      !
!                                                                       !
!***********************************************************************!
      subroutine swapxz2yz(var,aux,auxc,j,ivar)
      implicit none
      include 'ctes3D'
      
      integer iinp,iswap
      character*50 filinp,filout,filswap
      character*4  extensiones
      common /ficheros/ iinp,iswap,filinp,filout,filswap,extensiones(15)
      save /ficheros/
      
      ! ------------------------ variables -----------------------------!
      
      integer,intent(in)::j
      integer i,ivar,k,jp,jj,wk
      real*4 aux(2*mz,mx1+1)
      complex*8 var(mx1+1,mz),auxc(mz,mx1+1)
      character*52 fnameima     
      character*3 ext
      
      ! ------------------------ Program -------------------------------!

      iswap = 12
      do k=1,mz
         do i=1,mx1+1
            auxc(k,i) = var(i,k)
         enddo
      enddo

       write(ext,'(i3.3)') j

      fnameima=filswap(1:index(filswap,' ')-1)//ext//'/swap.'//
     %        extensiones(ivar)
      open(unit=iswap,file=fnameima,status='unknown',
     &   form='unformatted',access='direct',recl=mx*mz)

      write(iswap,rec=1)((aux(k,i),k=1,2*mz),i=1,mx1+1)
      
      endsubroutine swapxz2yz
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!      subroutine check_nan(u,n,flag)
!      implicit none
!
!      real*4 u(n)
!      integer j,i,n,flag
!
!      do j=1,n
!         if (isnan(u(j)))write(*,*) 'encontre nan'
!         if (isnan(u(j))) exit
!      enddo
!         
!      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!      subroutine check_nan8(u,n,flag)
!      implicit none
!
!      real*8 u(n)
!      integer j,i,n,flag
!
!      do j=1,n
!         if (u(j).ne.u(j))write(*,*) 'encontre nan' 
!      enddo
!         
!      end
!

