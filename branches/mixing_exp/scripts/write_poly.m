%************************************************************************************%
% Matlab script to create starting file 					     %
% Antonio Almagro 								     %
% May 2013									     %
% Script for LISO 								     %
%************************************************************************************%
clear all;
close all;
clc
%Relevant variables
%El perfil del modod 0 esta normalizado con 
% Ub=1, por tanto el Ree que le voy a meter al programa
% sera el Reb que quiero
Reb=10000;
%Reb=3800;
Ub=1;
%Variables (for record 1)
time=0;
Ree= Reb;
alpe=1;%tamaño de caja periodica en X (alpe=2PI/Lx)
bete=2;%tamaño de caja periodiza en z (bete=2PI/Lz)
aOe=0;%velocidad relativa de observacion
%mxe and mze must be some multiple of some kind.. mgalx 384 workss
mye=256
mgalx=384
mgalz=mgalx

mxe=mgalx*2/3
mze=mxe
my=mye

%Creating vector y (perpendicular to walls)
 %duct type -> chebysev
h=1;
%Desfase de theta para que no salgan muy junts
%los puntos cercanos a la pared
%for (j=1:mye)
%	theta=(mye-j)*pi/(mye-1);
%	y_cheb(j)=h*cos(theta);
%end
%y=y_cheb;

%FMAP(j)
fmap=ones(mye,1); %all ones not importante yet
%-------------------------------------------------%
%-----------MODE 1 for kx=1,kz=0------------------%
%kx, kz
kx=1;
kz=0;
[v,phi,vor,ee,vecty]=ossq_for_toni(kx,kz,Reb,mye-2);
phimod1=phi(:,1);%first mode for phi
vormod1=vor(:,1);%first mode for vor
%Must be turned cause ossq writes from 1 to -1
% and our code reads from -1 to 1.
%ii=[mye:-1:1]; %inverse order
%vecty=vecty(ii)
%vecty does not include extremes!

%Adding 0 --- 0
%Creating vector y
delta=3; %strecthing factor
eps=linspace(0,2,mye);
y=(1+tanh(delta*(eps/2-1/2))/tanh(delta/2))-1;
%Y can be manually set
%phimod1=phimod1(ii);
%vormod1=vormod1(ii);

%interpolating
phi1=interp1(vecty,phimod1,y,'spline','extrap');
vor1=interp1(vecty,vormod1,y,'spline','extrap');
%SCALING
%the BASE flow is ubulk=1; 
% we need to scale the output
factor=1E-4;
phi1=phi1*factor;
vor1=vor1*factor;
%INTERPOLATE FOR 
%!!!!!!!!!FUTURE WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% Interpolate between vecty given by ossq and   %
% the desired "y".                              %
% FOR SIMPLICITY now we work with y=vecty       %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%CHECKING
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

%wk1==> line for kx=1, kz=0
% Order: Re(vor),Im(vor),Re(phi),Im(phi)
%TEST same way as assing put data!
wk1=zeros(2,2*mye,mze);
k=1
	for j=1:mye
		wk1(1,2*j-1,k) =real(vor1(j));
		wk1(1,2*j  ,k) =imag(vor1(j));
        wk1(2,2*j-1,k) =real(phi1(j));
 		wk1(2,2*j,  k) =imag(phi1(j));
	end
	
size(wk1(:))

%wk1=[]; 
%for k=1:mze
%for j=1:mye
%	wk1(4*(j-1)+1:4*j)=[real(vormod1(j)),imag(vormod1(j)),real(phimod1(j)),imag(phimod1(j))];
%        wk1 = [wk1, real(vor1(j)),imag(vor1(j)),real(phi1(j)),imag(phi1(j))]; 
%        wk1 = [wk1, real(vor1),imag(vor1),real(phi1),imag(phi1)]; 
%end
%wk1 = [real(vormod1(:));imag(vormod1(:));real(phimod1(:)); imag(phimod1(:))]; 
%wk1=wk1; 
%wk1 = wk1(:); 



ntotr=4*mye*mze %FOR EACH PLANE
%In this case kz=0 is the first line for kz
linez = round(kz/bete) + 1
planex = round(kx/alpe)
%Creating a complete plane for this input mode:
wkplane=zeros(4*mye,mze);
wkplane=wk1(:);
%wkplane(:,linez)=wk1(:);

%All other lines/planes should be zero
wkplane2=0*wk1(:); %empty
%wkplane2=zeros(4*mye,mze);

%---------------------------------------------------------%
%==================================================%
%ZERO MODES
%LAMINAR PROFILE (Re_tau=180), Re_b=2800.
%uOO=3/2*Ub*(1-y.^2);
%uOO=y;
uOO=5/2*(y.^3-3/5*y.^5);
ubulk=trapz(y,uOO)/2
%Checking profile u00(y)
figure(1)
plot(y,uOO,'r*')
title('u00 vs y')
%w00 = 0;
wOO=uOO*0;
%============================================================%
%Re_tau=180;
%figure(4)
%plot(diff(y)*Re_tau)
%title('diff(y)*Re_tau')


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

for i=0:mxe-1
	if (i==planex)	
		fwrite(finp,rec2_size,'int');
		fwrite(finp,0*wkplane(:),'float');
		fwrite(finp,rec2_size,'int');
	else
		fwrite(finp,rec2_size,'int');
		fwrite(finp,0*wkplane2(:),'float');
		fwrite(finp,rec2_size,'int');
	end
end
	
%fwrite(finp,rec2_size,'uint');
%Write in wk1 the plane we are creating in the order
%according to LISO getfil().
%fwrite(finp,rec2_size,'uint');

%------------END OF RECORD 2---------------------------------%
    %shift
    fclose(finp);
%============================================================%
%===========END OF WRITING TO A FILE=========================%
%============================================================%


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

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
%Calculations for Re_tau
dudy_wall=3;
nu=1/Reb;
u_tau=sqrt(nu*dudy_wall);
Re_tau=u_tau/nu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

