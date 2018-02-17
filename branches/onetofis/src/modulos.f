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
      character(100) :: filinp,filout,filstt,fnamesta
      end module ficheros
c-------------------------------------------------------------------
c fis AAF
      module fis
      implicit none
      private
      include 'ctes3D'
      real(8), public ::  Re,alp,bet,a0,y(my),hy(my),
     .                    fmap(my),trp(my),mss(my),rhocrit
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
c


c--- point AAF
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
      real(8), public:: ff(my),ff2(my),ff3(my),ff4(my)
      real(8), public:: mean(my),stddev(my),skew(my),kurt(my)
      real(8), public:: wkstats(my,10) !working var
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


