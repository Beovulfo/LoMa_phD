%****************************************************************************%
% Matlab script to create starting file for mixingLISO			     %
% Antonio Almagro							     %
% June 2013								     %
%****************************************************************************%
clear all;
close all;
clc
%Relevant variables
%The reference speed will be u2-u1=2*u2=u0
u0=1;
Re=500;
%Re0=U*dw0/nu
ymax=5; %half heigth of the box

d=1; %shear layer thickness

% we need to scale the disturbance output

%Variables (for record 1)
time=0;
Ree= Re;%Re_max because our profile is for u_max=1.
%alpe=1;%tamaño de caja periodica en X (alpe=2PI/Lx)
%As stated in paper 3d-evolution of a plane mixing layer (moser):
%Streamwise and spanwise fundamental wavelengths chosen as the most
%unstable wavelengths from linear theory.
lambdax=2.32*pi;
lambdaz=0.6*lambdax;
N=1;
M=1;
Lx=lambdax*N;
Lz=lambdaz*M;
alpe=2*pi/Lx %tamaño de caja periodica en X (alpe=2PI/Lx)
bete=2*pi/Lz %tamaño de caja periodiza en z (bete=2PI/Lz)
aOe=0;%velocidad relativa de observacion
%mxe and mze must be some multiple of some kind.. mgalx 384 works
mye=256
mgalx=384
mgalz=mgalx

mxe=mgalx*2/3
mze=mxe
my=mye

Volum=Lx*Lz*ymax*2;
%Scealing to do integrated energy = 1
factor=1*Volum;%this gives energy=1
%factor=0.1*Volum/Re;
%Dimension half of the height of the computational box
%h=5;

%FMAP(j)
fmap=ones(mye,1); %fmap is indeed recalculated by program
%FMAP is changed by code
%-------------------------------------------------%
%Creating vector y ---GRID--------
%Using user matlab function to create the grid
%Check matlab script "create_mesh1D" for more info
y=create_mesh1D(ymax,'mixl',0.4,mye);
%y=create_mesh1D(length,'hype',5,mye);

%============================================================%
%ZERO MODES
%Wissink - On unconditional conservation of kinetic energy
%by finite-difference discretizations of the linear
%and non-linear convection equation
%Error function
U=u0;%1/2*(U_top-U_bot)
dw0=1; %vorticiy thickness, set to 1
z1=sqrt(pi).*y/dw0;
dz1dy=sqrt(pi)/dw0;
u00=U*erf(z1);
u00p=2/sqrt(pi)*dz1dy*U*exp(-z1.^2);%first derivative
u00pp=u00p.*(-2*z1*dz1dy); %second derivative


%Checking profile u00(y)
figure(1)
hold on
plot(y,u00,'b-')
plot(y,u00p,'g-')
plot(y,u00pp,'r-')

title('u00 vs y')
hold off
%----------------------------__ACCHTUNG NO MEAN FLOWWWW!!!!!-----%
w00=u00*0;
%==================================================%


%-----------one MODE for Disturbance ------------------%
%kx, kz
nx1=1;%first mode
nz1=0;

kx1=nx1*alpe
kz1=nz1*bete



%Calculation for each desired mode
%alp=kx1;bet=kz1;Reb=Re;N=mye;yy=y;uu=u00;Ly=heigth;uup=u00p;uup2=u00pp;
[phi,vor,ee1,vecty]=ossq_for_mixing(kx1,kz1,Re,mye,y,u00,u00p,u00pp,ymax);
u00=u00*0;
%We want a eigenvector with vor=0 and phi.neq.0
phimod1=phi(:,1);%first mode for phi
vormod1=vor(:,1);%first mode for vor
%clean the noise created on extremes
ndel=10;
cleanext=[zeros(1,ndel) ones(1,mye-2*ndel) zeros(1,ndel)];
phimod1=phimod1.*cleanext';
vormod1=vormod1.*cleanext';
%SCALING (integrated energy unity??)

%  phimod1=phimod1/;%scaling to 1
%  vormod1=vormod1/;%scaling to 1


%kx, kz
nx2=0;%in this case nx2 must be different than nx1
nz2=1;

kx2=nx2*alpe
kz2=nz2*bete

%interpolating /REAL and IMAGINARY and scaling
phi1real=interp1(vecty,real(phimod1),y,'spline','extrap')*factor;
phi1imag=interp1(vecty,imag(phimod1),y,'spline','extrap')*factor;
vor1real=interp1(vecty,real(vormod1),y,'spline','extrap')*factor;
vor1imag=interp1(vecty,imag(vormod1),y,'spline','extrap')*factor;


%2nd mode is (0,1) , phi=(0,-kz2*ome_x) --> Moses
c=1*Volum/Re;%needs to be calculated with energy
phi2real=phi1real*0;
phi2imag=-c*kz2*exp(-pi*y.^2)
vor2real=vor1real*0;
vor2imag=vor1imag*0;
%mode (0,-1) %we can use just -phi2imag and skip a variable
phi3imag=-phi2imag;


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
		wk1(1,2*j-1,linez1) =vor1real(j);
		wk1(1,2*j  ,linez1) =vor1imag(j);
%----------------------------------------------------!
                wk1(2,2*j-1,linez1) =phi1real(j);
 		wk1(2,2*j,  linez1) =phi1imag(j);
	end

wk2=zeros(2,2*mye,mze);
linez2 = nz2 + 1
planex2 = nx2+1
%kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
%% but matlab goes from 1 to mze
	for j=1:mye
		wk2(1,2*j-1,linez2) =vor2real(j);
		wk2(1,2*j  ,linez2) =vor2imag(j);
%----------------------------------------------------%
       	        wk2(2,2*j-1,linez2) =phi2real(j);
 		wk2(2,2*j,  linez2) =phi2imag(j);
	end
%add the mode (0,-1)---------------------------------------VIP
%mod_index=mze;
% De momento no SEGURO SI ESTE MODO ESTA BIEN METIDO!!
	for j=1:mye
		wk2(1,2*j-1,mze) =vor2real(j);
		wk2(1,2*j  ,mze) =vor2imag(j);
%----------------------------------------------------%
       	        wk2(2,2*j-1,mze) =phi2real(j);
 		wk2(2,2*j,  mze) =-phi2imag(j);
	end
%check size of wk1 (ntotr)	
%check size of wk1 (ntotr)	
%size(wk1(:))
ntotr=4*mye*mze %FOR EACH PLANE
%In this case kz=0 is the first line for kz
%Creating a complete plane for this input mode:
wkplane1=zeros(4*mye,mze);
wkplane1=wk1(:);
wkplane2=zeros(4*mye,mze);
wkplane2=wk2(:);

%All other lines/planes should be zero
wkplane0=0*wk1(:); %empty


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
	fwrite(finp,u00(j),'float');
	fwrite(finp,w00(j),'float');
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
% Only working if planex1 is different from planex2
%check
if (planex1==planex2)
   disp('CAUTION!!! planex1 and planex2 MUST BE DIFFERENT with this code')
else
   disp('Good job selecting the modes, they have different planex!')
end

%Write to file finp.dat
for i=1:mxe
	if (i==planex1)	%write first mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane1(:),'float');
		fwrite(finp,rec2_size,'int');
%Only 1 mod in order to calculate energy!
%        elseif (i==planex2) %write second mode
%		fwrite(finp,rec2_size,'int');
%		fwrite(finp,wkplane2(:),'float');
%		fwrite(finp,rec2_size,'int');
	else
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane0(:),'float');
		fwrite(finp,rec2_size,'int');
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
plot(y,phi1real,'r*')
hold on
plot(y,phi1imag,'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 1 PHI')

figure(3)
plot(y,vor1real,'r*')
hold on
plot(y,vor1imag,'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 1 VOR')

%CHECKING
%SECOND MODE
figure(4)
hold on
plot(y,phi2imag,'b*')
plot(y,phi3imag,'r*')
%plot(y,phi2real,'r*')
xlabel('y')
ylabel('Re/Im')
title('IMAG part of 2nd PHI mode:mode(0,1) and mode (0,-1)')

figure(5)
%plot(y,vor2real,'r*')
plot(y,vor2imag,'r*')
hold on
plot(y,vor2imag,'b*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 2 VOR')
