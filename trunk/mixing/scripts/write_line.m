%****************************************************************************%
% Matlab script to create starting file for mixingLISO			     %
% Antonio Almagro							     %
% June 2013								     %
%****************************************************************************%
clear all;
close all;
clc
%Relevant variables
%Relevant velocity this time would me umax
umax=1;
Re=10000;
% we need to scale the output
factor=0.1;
%Variables (for record 1)
time=0;
Ree= Re;%Re_max because our profile is for u_max=1.
alpe=1;%tamaño de caja periodica en X (alpe=2PI/Lx)
bete=2;%tamaño de caja periodiza en z (bete=2PI/Lz)
aOe=0;%velocidad relativa de observacion
%mxe and mze must be some multiple of some kind.. mgalx 384 works
mye=256*2
mgalx=384
mgalz=mgalx

mxe=mgalx*2/3
mze=mxe
my=mye
%Dimension half of the height of the channel (delta)
h=1; %do we need to change this deffinition?

%FMAP(j)
fmap=ones(mye,1); %fmap is indeed recalculated by program
%FMAP is changed by code
%-------------------------------------------------%
%Creating vector y ---GRID--------
%Using user matlab function to create the grid
%Check matlab script "create_mesh1D" for more info
y=create_mesh1D('hype',1.02,mye);

%-----------MODE 1 ------------------%
%kx, kz
nx1=25;%any disturbance with a large enough wavel enght
       %wich causes the shear layer to be locally disp.
       %in the y-direction, will trigger such an instability
nz1=0;

kx1=nx1*alpe
kz1=nz1*bete

%Calculation for each desired mode
[v,phi,vor,ee1,vecty]=ossq_for_toni(kx1,kz1,Re,mye-2);
phimod1=phi(:,1);%first mode for phi
vormod1=vor(:,1);%first mode for vor
%Must be turned cause ossq writes from 1 to -1
% and our code reads from -1 to 1.
%vecty does not include extremes!
%--------------MODE 2---------------
%kx, kz
nx2=2;%in this case nx2 must be different than nx1
nz2=1;

kx2=nx2*alpe
kz2=nz2*bete
%Calculation for each desired mode
[v,phi,vor,ee2,vecty]=ossq_for_toni(kx2,kz2,Re,mye-2);
phimod2=phi(:,1);%first mode for phi
vormod2=vor(:,1);%first mode for vor


%interpolating
phi1=interp1(vecty,phimod1,y,'spline','extrap');
vor1=interp1(vecty,vormod1,y,'spline','extrap');
phi2=interp1(vecty,phimod2,y,'spline','extrap');
vor2=interp1(vecty,vormod2,y,'spline','extrap');
%scaling
phi1=phi1*factor;
vor1=vor1*factor;
phi2=phi2*factor;
vor2=vor2*factor;
%wk1==> line for kx=1, kz=0
% Order: Re(vor),Im(vor),Re(phi),Im(phi)
%Same way as assing put data:
wk1=zeros(2,2*mye,mze);
linez1 = nz1 + 1
planex1 = nx1+1
 %kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
% but matlab goes from 1 to mze
	for j=1:mye
		wk1(1,2*j-1,linez1) =real(vor1(j));
		wk1(1,2*j  ,linez1) =imag(vor1(j));
       	        wk1(2,2*j-1,linez1) =real(phi1(j));
 		wk1(2,2*j,  linez1) =imag(phi1(j));
	end

wk2=zeros(2,2*mye,mze);
linez2 = nz2 + 1
planex2 = nx2+1
 %kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
% but matlab goes from 1 to mze
	for j=1:mye
		wk2(1,2*j-1,linez2) =real(vor2(j));
		wk2(1,2*j  ,linez2) =imag(vor2(j));
       	        wk2(2,2*j-1,linez2) =real(phi2(j));
 		wk2(2,2*j,  linez2) =imag(phi2(j));
	end
%check size of wk1 (ntotr)	
%size(wk1(:))
ntotr=4*mye*mze %FOR EACH PLANE
%In this case kz=0 is the first line for kz
%Creating a complete plane for this input mode:
wkplane1=zeros(4*mye,mze);
wkplane1=wk1(:);
wkplane2=zeros(4*mye,mze);
wkplane2=wk2(:);
%wkplane(:,linez)=wk1(:);

%All other lines/planes should be zero
wkplane0=0*wk1(:); %empty
%wkplane2=zeros(4*mye,mze);

%---------------------------------------------------------%
%==================================================%
%ZERO MODES (Laminar)
uOO=y;
ubulk=trapz(y,uOO)/2
%Checking profile u00(y)
figure(1)
plot(y,uOO,'r*')
title('u00 vs y')
%w00 = 0;
wOO=uOO*0;


%============================================================%
%===================WRITING TO FILE==========================%
%save record 1 values
finp=fopen('../data/finp.dat','w','l'); %open/create file tab_var.dat
%Record 1 - counting bytes
rec1_real4=5+2*mye;%time is real 4
%Ree,alpe,bete,a0e + y(j),fmap(j) (1:my)+wk1(j)(1:2my)
rec1_real8=mye*2;
%integer mxe mye mze
rec1_integer=3;
%total record 1 size
rec1_size=4*rec1_real4+8*rec1_real8+4*rec1_integer;
%----------------------Record 1----------------------------%
%Buffer size    
fwrite(finp, rec1_size,'uint');
%record 1 data 
    fwrite(finp,time,'float');
    fwrite(finp,Ree,'float');
    fwrite(finp,alpe,'float');
    fwrite(finp,bete,'float');
    fwrite(finp,aOe,'float');
    fwrite(finp, mxe,'int');
    fwrite(finp, mye,'int');
    fwrite(finp, mze,'int');
%write y(j) and fmap(j)
for (j=1:mye)
    fwrite(finp,y(j),'double');
    fwrite(finp,fmap(j),'double');
end
%wk1==>u00, w00, odd/even
for (j=1:mye)
	fwrite(finp,uOO(j),'float');
	fwrite(finp,wOO(j),'float');
end
%write buffer size again
fwrite(finp, rec1_size,'uint');
%------end of record 1---------------------------------------%

%----------RECORD 2 ---FOR EACH PLANES ALL VALUES------------%
%Record 2 - vor,phi %NEED REVISION

rec2_size=ntotr*4; %SAME SIZE FOR EACH PLANES
%Buffer size
%for each plane write wk2 if not planex
%First plane is i=0, but we start from one!
%WARNING WITH THE ELSEIF!!! not working as wished :P

%CAUTION!!!! 0 x !!!!
if 1 % ACHTUNG this is a stupid test
   for i=1:mxe
	if (i==planex1)	%write first mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,0.*wkplane1(:),'float');
		fwrite(finp,rec2_size,'int');
        elseif (i==planex2) %write second mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,0.*wkplane2(:),'float');
		fwrite(finp,rec2_size,'int');
	else
		fwrite(finp,rec2_size,'int');
		fwrite(finp,0.*wkplane0(:),'float');
		fwrite(finp,rec2_size,'int');
	end
end
 
else
for i=1:mxe
	if (i==planex1)	%write first mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane1(:),'float');
		fwrite(finp,rec2_size,'int');
        elseif (i==planex2) %write second mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane2(:),'float');
		fwrite(finp,rec2_size,'int');
	else
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane0(:),'float');
		fwrite(finp,rec2_size,'int');
	end
end
end
	
    fclose(finp);
%============================================================%
%===========END OF WRITING TO A FILE=========================%
%============================================================%


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
%Calculations for Re_tau
%dudy_wall=3; %Parabolic u(y) profile!!!
%nu=1/Re;
%u_tau=sqrt(nu*dudy_wall);
%Re_tau=u_tau/nu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================================================%
%CHECKING
%FIRST MODE
figure(2)
plot(y,real(phi1),'r*')
hold on
plot(y,imag(phi1),'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 1 PHI')

figure(3)
plot(y,real(vor1),'r*')
hold on
plot(y,imag(vor1),'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 1 VOR')

%SECOND MODE
figure(5)
plot(y,real(phi2),'r*')
hold on
plot(y,imag(phi2),'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 2 PHI')

figure(6)
plot(y,real(vor2),'r*')
hold on
plot(y,imag(vor2),'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 2 VOR')
%Checking grid, Wall units!
%figure(7)
%plot(diff(y)*Re_tau)
%title('Grid in wall units (Re_tau)')
%axis 'equal'
%xlabel('Grid y')
%ylabel('Dy*')
%===============================================================================%
%Select type of initial condition
%%TURBULENT PROFILE%%%%
%U00 == torroja
%W00 == zeros
%%	nghost=0;
%	dmax=2; %dimenion del canal 2x1
%	turbprof=importdata('Re180.dat');%import data from file torroja
%	%xprof are the coordinates y/h
%	xprof=turbprof(:,1); 
%	%yprof are the mean velocities
%	yprof=turbprof(:,3);
%	%Start the mirror,we have only half of the profile
%	lxprof=length(xprof);
%	xprofmirror=zeros(lxprof-1,1);
%	yprofmirror=zeros(lxprof-1,1);
%	%create the other half!
%	for (i=(lxprof-1):-1:1)
%		xprofmirror(i)=dmax-xprof(lxprof-i);
%		yprofmirror(i)=yprof(lxprof-i);
%	end
%	xprof=[xprof;xprofmirror];
%	yprof=[yprof;yprofmirror];
%	%xi_int keeps the coordinates as taken from our run and sums 1 to make them monotically positive
%	nghost=0;
%	xi_int=1+y(1:mye);
	%interpolate the velocity in our points (from the 49 points taken from file Re180.dat to the nytab points).
%	u_interp=interp1q(xprof,yprof,xi_int');
	%Save the u_interp into the u matrix
 %      uOO=u_interp;
%============================================================%
% COmpute wavenumbers like LISO
%as defined in ctes3D (change for matlab indexes (1:n))
%NEED REVISION!!!
%nz=mze/2;
%nz1=nz;
%mgalx=mxe;
%mx=2*mgalx/3;
%
%mx1=mxe/2;
%mz1=mze;
%%First half
%for k=1:nz1
%	xbet(k)=complex(0,bete*(k-1));
%	rbet(k)=bete*(k-1);
%	icx(k)=k-1;
%end
%
%for k=nz1+1:mz1
%	xbet(k)=complex(0,-bete*(mz1-k));
%	rbet(k)=-bete*(mz1-k);
%end
%
%for k=1:nz1
%	icx(mze-(k-1))=k-1;
%end
%
%for i=1:mx1
%	iax(2*i-1)=i-1;%revisar
%	iax(2*i)=i-1;%revisar
%	xalp(i)=complex(0,alpe*(i-1));
%	ralp(i)=alpe*(i-1);
%	alp2(i)=-xalp(i)^2;
%end
%
%for j=1:mz1
%	bet2(j)=-xbet(j)^2; %revisar
%end


