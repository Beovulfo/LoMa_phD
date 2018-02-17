c-------------------------------------------------------------------
c inparams
      module inparams
      implicit none
      private
      real(8), public :: euutarget,evvtarget,ewwtarget,dm0,u0
     .                   ,swin,sden,DT
      integer, public::  nymax, nhmin,nhmax
      end module inparams
c-------------------------------------------------------------------
c diag AAF
      module diag
      implicit none
      real(4) :: ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
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
     .                    fmap(my),trp(my),mss(my),A1,A2
      end module fis
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
!c spectra AAF -> this module has to be reworked
!      module spectra
!      implicit none
!      private
!      include 'ctes3D'
!      integer, public:: nacumsp,jsp(my),
!     .                  jspiproc(nspec),
!     .                  jspbeg(0:numerop-1), jspend(0:numerop-1),
!     .                  jspb,jspe,jspbb
!      end module spectra
!c-------------------------------------------------------------------
!c statis AAF -> this module has to be reworked
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
     .                  Tp(my),ep(my)
      real(8), public:: wkstats(my,11) !working var
      integer, public:: istati,ntimes,nacum,nstart
      end module statis

!c-------------------------------------------------------------------
      module base
      implicit none
      private
      include 'ctes3D'
      real(8),public :: Um0(my), V0(my), H0(my), Z0(my), rho0(my)
      real(8),public :: upp(my), vpp(my),rhop(my),hp(my),zp(my)
      end module base
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
     .           ctm,ctm2,ttm
      end module timers
c-------------------------------------------------------------------
c cnan AAF 
      module cnan
      implicit none
      integer nanerror,nanproblem
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
     .                  ten23wk(0:mgalx+1),rhstwk(0:mgalx+1)
      real(8),public :: up1wk8(0:mgalx+1),up2wk8(0:mgalx+1),
     .                  up3wk8(0:mgalx+1),
     .                  ten11wk8(0:mgalx+1),ten22wk8(0:mgalx+1),
     .                  ten33wk8(0:mgalx+1),
     .                  ten12wk8(0:mgalx+1),ten13wk8(0:mgalx+1),
     .                  ten23wk8(0:mgalx+1),rhstwk8(0:mgalx+1),
     .                  tmp1wk8(0:mgalx+1),tmp2wk8(0:mgalx+1),
     .                  tmp3wk8(0:mgalx+1)
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

!c-------------------------------------------------------------------
!      module mass 
!      implicit none
!      private
!      include 'ctes3D'
!      integer, public :: imass
!      real(4), public :: mub,put
!      real(4), public :: trp(my)
!      real(8), public :: trp2(my)
!      end module mass
!c-------------------------------------------------------------------
c      module mesh 
c      implicit none
c      real(4) gamma
c      integer imesh
c      end module mesh 

!c-------------------------------------------------------------------
!      module strat
!      implicit none
!      integer istra
!      real(4) Pr,Fr,iFr2
!      end module strat
!!c-------------------------------------------------------------------
!      module override
!      implicit none
!      integer ioverride
!      real(4) overridemub,overrideput
!      end module override

c-------------------------------------------------------------------
c-------------------------------------------------------------------
      !Now Z is actually Zbar
      module combustion
      implicit none
      !private
      !include 'ctes3D'
      real(8), public ::  gam,Zs,Zsbar,betha,iLm,Lf,S,T0f,Hf
      real(8), public ::  diffcoef,diff1,diff2,diff3

!     FUNCTIONS 
      contains
        function dTdZ(Z)
        real(8),intent(in):: Z
        real(8)::dTdZ
        dTdZ = gam/((1d0-Zs)*Zsbar)*0.5d0*(dtanh(betha*(Zsbar-Z))+1d0)
        !dTdZ = gam/((1d0-Zs)*Zs)*0.5d0*(dtanh(betha*(Zs-Z))+1d0)
        end function

!        function Temper(H,Z)
!        real(8),intent(in):: H,Z
!        real(8)::Temper
!        Temper = 1 + H +0.5d0*(gam/((1-Zs)*Zsbar))*(Z +1/betha*(-dlog
!     .   (dcosh( betha*(Z-Zsbar)))+dlog(dcosh(betha*Zsbar))))
!!        Temper = 1 + H +0.5d0*(gam/((1-Zs)*Zs))*(Z +1/betha*(-dlog
!!     .   (dcosh( betha*(Z-Zs)))+dlog(dcosh(betha*Zs))))
!        end function
        function Temper(H,Z)
        real(8),intent(in):: H,Z
        real(8)::Temper
        Temper = 1 + H + gam/((1-Zs)*Zsbar)*Z
        if (Z.ge.Zsbar) then
           Temper = 1 + H + gam/(1-Zs)
        endif
        end function


        function Zb(Z)
        real(8),intent(in):: Z
        real(8)::Zb
        Zb = Z*Zsbar/Zs 
        if (Z.ge.Zs) then
           Zb = Zsbar + (1-Zsbar)/(1-Zs)*(Z-Zs)
        endif
        end function


!        function Zb(Z)
!        real(8),intent(in):: Z
!        real(8)::Zb
!        Zb = Z*0.5*((1.0-Zsbar)/(1.0-Zs) +  Zsbar/Zs) 
!     .    +((1.0-Zsbar)/(1.0-Zs) -Zsbar/Zs)/(2.0*betha)
!     .    *(dlog(dcosh(betha*(Z-Zs)))
!     .     -dlog(dcosh(betha*Zs))) 
!!        dZbdZ = Zsbar/Zs+((1d0-Zsbar)/(1d0-Zs)-Zsbar/Zs)*
!!     .             0.5d0*(dtanh(betha*(Z-Zs))+1d0)
!        end function


      end module combustion
c-------------------------------------------------------------------
!

!
