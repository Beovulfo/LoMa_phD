c-------------------------------------------------------------------
c diag AAF
      module diag
      implicit none
      real(4) :: ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL,energy(9)
      real(8) :: dm
      end module diag
c-------------------------------------------------------------------
c ficheros AAF
      module ficheros
      implicit none
      integer :: iinp,iout,id22,isn,ispf,ista
      character(100) :: filinp,filout,filstt
      end module ficheros
c-------------------------------------------------------------------
c fis AAF
      module fis
      implicit none
      private
      include 'ctes3D'
      real(8), public ::  Re,alp,bet,a0,y(my),hy(my),
     .                    fmap(my),trp(my),mss(my),rhocrit,
     .                    Pr,peclet,sigma
      end module fis

cc rkcoef AAF
      module rkcoef
      implicit none
      real(8)  :: gama(3),alpha(3),beta(3),ibeta(3),xi(3),p0,p1,b3
      real(8)  :: dtbeta
      real(8) ::  rk,rkn1,dalbe,dtri,dtxi,dtgamma,dalre,ire
      !parameter (p0 = -0.2d0, p1=0.4d0)
      parameter (p0 = 5d0/66d0, p1=-1.0d0)
      parameter (xi   =(/      0d0, -17d0/60d0, -5d0/12d0 /)) !aaf 
      parameter (gama= (/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
      parameter (b3 = (p1*(gama(1)*gama(3)+gama(2)*gama(3)+
     .   xi(2)*gama(3)+gama(1)
     .  *xi(3))-gama(1)*gama(2)*gama(3)-gama(1)*gama(2)*xi(3))/
     .   (p1-gama(1)*gama(2)))
       parameter (beta =(/gama(1)-p0, p1, b3/))
       parameter (alpha=(/p0,gama(2)+xi(2)-p1,gama(3)+xi(3)-beta(3)/))
      end module rkcoef      
c-------------------------------------------------------------------
c-------------------------------------------------------------------
!MPI_GROUPS aaf
      module MPI_GROUPS
      implicit none
      integer :: MPI_GROUP_WORLD,MPI_GROUP_CALC,MPI_COMM_CALC,
     .                    MPI_GROUP_SAVE,MPI_COMM_SAVE
      save
      end module MPI_GROUPS
!c-------------------------------------------------------------------
c MPI_datatype aaf
      module MPI_datatype 
      implicit none
      private
      include 'ctes3D'
      integer, public :: myslice(0:numerop-1),myplane(0:numerop-1),
     .                   plbypl(0:numerop-1),req(0:2*numerop-1)
      end module MPI_datatype
!c-------------------------------------------------------------------
c point AAF
      module point
      implicit none
      private
      include 'ctes3D'
      integer, public:: lbeg(0:numerop-1),lend(0:numerop-1),
     .                  pbeg(0:numerop-1),pend(0:numerop-1),
     .                  plcalc(nplanes),plsav(nplanes),
     .                  pb,pe,lb,le,mmp,mml,procs
      end module point
c-------------------------------------------------------------------
      module spectra
      implicit none
      private
      include 'ctes3D'
      integer, public:: nacumsp,jsptot(nspec+1)
      end module spectra
!c-------------------------------------------------------------------
      module statis
      implicit none
      private
      include 'ctes3D'
      real(8), public:: um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  rhom(my),drhom(my),rho2(my),
     .                  rum(my),rvm(my),rwm(my),Tm(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  Tp(my),ep(my),ruu(my),rvv(my),
     .                  rww(my),ruv(my),ruw(my),
     .                  rhop(my),num(my),chi(my),
     .                  rvw(my),them(my),thep(my),theup(my),
     .                  mum(my),Zm(my),Hm(my),Zp(my),Hp(my)
      real(8), public:: wkstats(my,25) !working var
      !real(8), public:: wkstats(my,18) !working var
      integer, public:: istati,ntimes,nacum,nstart
      end module statis
!c-------------------------------------------------------------------
c tem AAF no private
      module tem
      implicit none
      real(4) :: Deltat,CFL,time,dtr
      end module tem
c-------------------------------------------------------------------
c timacc AAF no private
      module timacc
      implicit none
      integer :: nimag,nstep,nhist,ihist,icfl,ncfl
      end module timacc
c-------------------------------------------------------------------
c timers AAF
      module timers
      real(8) :: commtimer,transtimer,totaltimer,
     .           ctm,ctm2,ttm,writetimer
      end module timers
c-------------------------------------------------------------------
c cnan AAF 
      module cnan
      implicit none
      integer,public:: nanerror,nanproblem
      end module cnan
c-------------------------------------------------------------------
c wave AAF
      module wave
      implicit none
      private 
      include 'ctes3D'
      integer, public :: iax(mx),icx(0:mz1)
      real(4), public :: alp2(0:mx1),bet2(0:mz1)
      complex(8), public :: xalp(0:mx1),xbet(0:mz1)
      end module wave
c-------------------------------------------------------------------
c matrices AAF
      module matrices
      implicit none
      private
      include 'ctes3D'
      real(8), public :: prem1(7,my),dt12(7,my),endder(my),
     .                   prem3(5,my),dt21(5,my),dt22(5,my),
     .                   dt11(7,my)                
      end module matrices
c-------------------------------------------------------------------
c wkhvect AAF 
      module wkhvect
      implicit none
      private
      include 'ctes3D'
      real(4),public :: up1wk(0:mgalx+1),up2wk(0:mgalx+1),
     .                  up3wk(0:mgalx+1),ten11wk(0:mgalx+1),
     .                  ten22wk(0:mgalx+1),ten33wk(0:mgalx+1),
     .                  ten12wk(0:mgalx+1),ten13wk(0:mgalx+1),
     .                  ten23wk(0:mgalx+1),rhstwk(0:mgalx+1),
     .                  tau11wk(0:mgalx+1),tau22wk(0:mgalx+1),
     .                  tau33wk(0:mgalx+1),tau12wk(0:mgalx+1),
     .                  tau13wk(0:mgalx+1),tau23wk(0:mgalx+1),
     .                  rhszwk(0:mgalx+1),lapzwk(0:mgalx+1),
     .                  drhowk(0:mgalx+1),tnextwk(0:mgalx+1)
      real(8),public :: up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .                  up3wk8(0:mgalx+1),
     .                  ten11wk8(0:mgalx+1),ten22wk8(0:mgalx+1),
     .                  ten33wk8(0:mgalx+1),
     .                  ten12wk8(0:mgalx+1),ten13wk8(0:mgalx+1),
     .                  ten23wk8(0:mgalx+1),rhstwk8(0:mgalx+1),
     .                  tmp1wk8(0:mgalx+1),tmp2wk8(0:mgalx+1),
     .                  tmp3wk8(0:mgalx+1),tnwk8(0:mgalx+1),
     .                  znwk8(0:mgalx+1),dtwk8(0:mgalx+1),
     .                  rhszwk8(0:mgalx+1),tnextwk8(0:mgalx+1),
     .                  drhowk8(0:mgalx+1),lapzwk8(0:mgalx+1)
      real(8),public :: hache(0:mgalx-1),zeta(0:mgalx-1)
      real(8),public :: Dh(0:mgalx-1),Dzb(0:mgalx-1)
      real(8),public :: DTemp(0:mgalx-1)
      end module wkhvect
c-------------------------------------------------------------------
c worksp AAF
      module worksp
      implicit none
      private
      include 'ctes3D'
      real(8),public :: wk1(5,my),fwk1(my),fwk2(my),fwk(my),
     .                dphi1(my),dphi2(my),
     .                phipr(my),phipi(my),vpr(my),vpi(my),v1(my),v2(my),
     .                dvpr(my),dvpi(my),phi1(my),dv1(my),dv2(my)
      end module worksp
c-------------------------------------------------------------------
c fis AAF
      !Now Z is actually Zbar
      module combustion
      implicit none
      !private
      !include 'ctes3D'
      real(8), public ::  gam,Zs,Zsbar,betha,iLm,Lf,S,T0f,Hf
      real(8), public ::  diffcoef,diff1,diff2,diff3

!     FUNCTIONS 
      contains
        function heav(x)
        real(8), intent(in)::x
        real(8):: heav
        heav = (tanh(betha*x)+1.0)/2
        end function
        
        function dheav(x)
        real(8), intent(in)::x
        real(8):: dheav
        dheav = betha*(1.0-tanh(betha*x)**2)/2
        end function



        function Temper(H,Z)
        real(8),intent(in):: H,Z
        real(8)::Temper
        real(8):: T1,T2,A
        A = gam/((1.0-Zs)*Zsbar)
        T1 = A*Z
        T2 = gam/(1.0-Zs)
        !T1 =  A*Z
        !T2 =  gam/(1.0-Zs)
        !Temper = 1.0 + H + T2*(1+(Z/Zsbar-1)*heav(Zsbar-Z))
        Temper= 1.0 + H + T1*heav(Zsbar-Z)+T2*heav(Z-Zsbar)
        end function

        function dTdZ(Z)
        real(8),intent(in):: Z
        real(8)::dTdZ
        real(8)::A,T1,T2
        A = gam/((1.0-Zs)*Zsbar)
        T1 = A*Z
        T2 = gam/(1.0-Zs)
        dTdZ = -T1*dheav(Zsbar-Z) + A*heav(Zsbar-Z) +
     .          T2*dheav(Z-Zsbar)

        !T1 = gam/((1-Zs)*Zsbar)
        !T2 = 0
        !dTdZ = T1*heav(Zsbar-Z)
        end function


        function dZbdZ(Z)
        real(8),intent(in):: Z
        real(8)::dZbdZ
        real(8):: Z1,Z2
        Z1 = Zsbar/Zs
        Z2 = (1.0-Zsbar)/(1.0-Zs)
        dZbdZ = Z1*heav(Zsbar-Z)+Z2*heav(Z-Zsbar)
        end function


      end module combustion
c-------------------------------------------------------------------
!
