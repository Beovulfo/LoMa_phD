      program Binary_Plane_to_Plot
       !This program writes a binary Tecplot or VTK file in single precision
       !from a plane of data.
       !Modification History
       !-----------------------------------------------------------------------
       !------------------------------------------------------------------------
       implicit none
       !REAL TYPES
       integer,parameter :: r8 = 8 
       integer,parameter :: r4 = 4 
       !Tem module
       real(4) :: Deltat,CFL,time,dtr
       !fis module
       real(4) :: Re,alp,bet,a0
       real(8), allocatable ::  y(:),hy(:),
     .             fmap(:),trp(:),mss(:)
       real(4), allocatable ::  xvec(:),zvec(:)
       real(4) ::  lengthx, lengthz
       character*3 :: nfill
       character*4 :: extp
       character*2 :: plnij

      
       !PLANE OUTPUT
       integer                             :: nstep
       integer                             :: np1, np2
       logical                             :: single
       real(r8),allocatable,dimension(:,:) :: DPplane
       real(r4),allocatable,dimension(:,:) :: SPplane
       real(r4),allocatable,dimension(:,:,:) :: SP3D
       real(r4),allocatable,dimension(:,:) :: SPplanet
       integer                             :: nstart, nend, nskip
       integer                             :: dir, index1,iu, iv, iw 
       logical                             :: List,oldnew
       real(r8)                            :: cloc,eloc
       real(4)                             :: cL,eL
       !GRID
       real(r8),allocatable,dimension(:)    :: gc1, gc2, ge1, ge2

       !TECPLOT
       !real(r4),allocatable,dimension(:,:) :: tecvar
       !real(r4),allocatable,dimension(:,:) :: g1tec,g2tec
       !logical                             :: tecplot_output
       !character(len=150)                  :: ss1, ss2, ss3
      
       !VTK
       real(r4),allocatable,dimension(:,:) :: vtkvar
       real(r4),allocatable,dimension(:)   :: g1vtk,g2vtk
       logical                             :: paraview_output
     
       !FILENAMES
       character(len=100)                   :: InFileName, OutFileName 
       character(len=100)                   :: dataDIR, outputDIR
       character(len=100)                   :: basename
       !OTHER STRINGS
       character(len=25)                   :: ss
       character(len=25)                   :: TITLE
       character(len=100)                  :: stemp
       integer                             :: l1,l2,l3
       logical                             :: slogical
       character(len=25),allocatable,dimension(:) :: Sname, Gname
       integer,allocatable,dimension(:)    :: group
      
       !LOOPING VARIABLES
       integer                             :: i,j,k,n,jj
      
      
       !STATUS VARIABLES
       integer                            :: s1
       
       !DEBUG
      ! logical,parameter                          :: debug=.false.
       logical,parameter                  :: debug=.true.
       !ctes3D reading parameters
       integer :: mgalx,mgalz,my,mx,mz,nplanes,nspec,plyx,plyz
       integer,allocatable :: jspecy(:)
         
      
       !write(6,'(a)') 'Data Directory:'
       !read(5,*) dataDir
       ! if ( debug ) write(6,*) dataDir
      
       !write(6,'(a)') 'Output Directory:'
       !read(5,*) outputDIR
       ! if ( debug ) write(6,*) outputDir
      
       !write(6,'(a)') 'Tecplot:'
       !read(5,*) tecplot_output
       ! if ( debug ) write(6,*) tecplot_output
      
!       write(6,'(a)') 'Paraview:'
!       read(5,*) paraview_output
!        if ( debug ) write(6,*) paraview_output
      
       write(6,'(a)') 'Base Name:'
       read(5,*) BaseName
        if ( debug ) write(6,*) BaseName

       write(6,'(a)') 'Extension:'
       read(5,*) extp
        if ( debug ) write(6,*) extp
        plnij = extp(3:4)
        if ( debug ) write(6,*) plnij
      
       write(6,'(a)') 'Start:'
       read(5,*) nstart
        if ( debug ) write(6,*) nstart
      
       write(6,'(a)') 'End:'
       read(5,*) nend
        if ( debug ) write(6,*) nend
      
       write(6,'(a)') 'Stride:'
       read(5,*) nskip
        if ( debug ) write(6,*) nskip
      
       plyz = 1 !always 1 if not changing 
       write(*,*) "EOF"
       do n=nstart,nend,nskip
        write(nfill,'(I3.3)') n
        write(*,*) 'image number:',nfill
        write(InFileName,'(a,a,a)')
     . trim(BaseName)//"_",nfill,"."//trim(extp)
        if ( debug ) write(*,*) InFileName
        open(unit=500,file=InFileName,status='old',
     .        form='unformatted',iostat=s1,action='read')
       read(500) time,Re,alp,bet,mgalx,my,mgalz,nspec,plyx
       if ( debug ) write(*,*) 'header 1 read'
       allocate (jspecy(1:nspec))
       allocate (y(1:my))
       allocate (fmap(1:my))
       allocate (xvec(1:mgalx))
       allocate (zvec(1:mgalz))
       read(500) (jspecy(jj),jj=1,nspec),(y(jj),fmap(jj),jj=1,my)
       if ( debug ) write(*,*) 'header 2 read'
       
       if (plnij.eq.'yx') then
         allocate( SPplane(1:my,1:mgalx))
         allocate( SPplanet(1:mgalx,1:my))
         do jj=1,my
            read(500) (SPplane(jj,i),i=1,mgalx)
         enddo
         !tranpose
         do i = 1,mgalx
           do jj=1,my
              SPplanet(i,jj) = SPplane(jj,i)
           enddo
         enddo

       elseif (plnij.eq.'xz') then
          allocate( SPplane(1:mgalx,1:mgalz))
          allocate( SP3D(1:mgalx,1:nspec,1:mgalz))
          do jj=1,nspec
             read(500) ((SPplane(i,k),i=1,mgalx),k=1,mgalz)
             do i=1,mgalx
                do k=1,mgalz
                  SP3D(i,jj,k) = SPplane(i,k)
                enddo
              enddo
             !save in 3D
          enddo
     
       endif
          
       if ( debug ) write(*,*) 'pln ',plnij,' read'

       !*********************************************
       !***************WRITE VTK FILE****************
       !*********************************************
       write(OutFileName,'(a)')
     . trim(InFileName)//".vtk" 
 
       open(unit=13,file=OutFileName,access='stream',form='unformatted',status='new',
     &          convert='big_endian',iostat=s1)
       if ( debug ) write(*,*) OutFileName
      
        !HEADER: note termination with char(10)
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(InFileName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)
        if ( debug ) write(*,*) 'VTK header written'
     
         if (plnij.eq.'yx') then
            k = 1
            write(ss,fmt='(A10,3I5)') "DIMENSIONS",mgalx,my,1
            write(13) ss//char(10)
            !X-grid
            write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",mgalx," float"
            write(13) char(10)//ss//char(10)
            do j = 1, mgalx
             xvec(j) = 2*3.141592/alp*(-0.5 + (j-1)/(mgalx-1))
             write(13) xvec(j)
            enddo
            !Y-grid
            write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",my," float"
            write(13) char(10)//ss//char(10)
            do j = 1, my
             write(13) real(y(j))
            enddo
            !Z-grid
            write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
            write(13) char(10)//ss//char(10)
            zvec(plyz) = -2*3.141592/bet
            write(13) zvec(plyz)

           !Field
           write(ss,fmt='(A10,I15)') "POINT_DATA",my*mgalx
           write(13) char(10)//ss//char(10)
           write(13) "SCALARS Fraction float 1"//char(10)
           write(13) "LOOKUP_TABLE default"//char(10)
           write(13) SPplane(:,:)
           if ( debug ) write(*,*) 'VTK file written'
         elseif (plnij.eq.'xz') then
            k = 1
            write(ss,fmt='(A10,3I5)') "DIMENSIONS",mgalx,nspec,mgalz
            write(13) ss//char(10)
            !X-grid
            write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",mgalx," float"
            write(13) char(10)//ss//char(10)
            do j = 1, mgalx
             xvec(j) = 2*3.141592/alp*(-0.5 + (j-1)/(mgalx-1))
             write(13) xvec(j)
            enddo
            !Y-grid
            write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",nspec," float"
            write(13) char(10)//ss//char(10)
            do jj = 1, nspec
             write(13) real(y(jspecy(jj)))
            enddo
            !Z-grid
            write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",mgalz," float"
            write(13) char(10)//ss//char(10)
            do j = 1, mgalz
             zvec(j) = 2*3.141592/bet*(-0.5 + (j-1)/(mgalz-1))
             write(13) zvec(j)
            enddo

           !Field
           write(ss,fmt='(A10,I15)') "POINT_DATA",nspec*mgalx*mgalz
           write(13) char(10)//ss//char(10)
           write(13) "SCALARS Field float 1"//char(10)
           write(13) "LOOKUP_TABLE default"//char(10)
           write(13) SP3D
           if ( debug ) write(*,*) 'VTK file written'
 


         endif
         
         !Close VTK File
        close(13)
      
       if (allocated(SPplane) ) deallocate(SPplane)
       if (allocated(jspecy) ) deallocate(jspecy)
       if (allocated(y) ) deallocate(y)
       if (allocated(fmap) ) deallocate(fmap)
       if (allocated(xvec) ) deallocate(xvec,zvec)
       if (allocated(SPplanet) ) deallocate(SPplanet)
       if (allocated(SP3D) ) deallocate(SP3D)
      
       enddo !LOOP OVER FILES
      
      stop
      end program
      
      
      !  l1=1
      !  l3=len(header)
      ! do i=1,nstats
      !   stemp=header(l1:l3)
      !   l2=scan(stemp,' ',slogical)
      !   Names(i)=stemp(1:l2)
      !   write(6,'(i3,3x,a)') i, trim(Names(i)) 
      !   l1=l1+l2
      !  enddo
      
      
      
