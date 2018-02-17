%****************************************************************************%
% Matlab script to create starting file for LOMA			     %
% Antonio Almagro							     %
% June 2013								     %
%****************************************************************************%
clear all;
close all;
clc
theta=pi/4
%Relevant variables
%The reference speed will be u2-u1=2*u2=u0
u0=1;
Re=500;
%Re0=U*dw0/nu
%ACHTUNG MESH NOW FROM -1 TO 1
ymax=7; %half heigth of the box
%init time
time=0;
Ree= Re;
%factor de onda, define tamaño de caja (alpe=2PI/Lx)
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
%mxe and mze must be some multiple of 64*3/2
mye=256
mgalx=192
mgalz=mgalx
%physical mesh size
Dx=Lx/(mgalx-1)
Dz=Lz/(mgalz-1)
%convert to number of Fourier modes
mxe=mgalx*2/3
mze=mgalz*2/3-1
my=mye
%Size of Comp. Box
Volum=Lx*Lz*ymax*2;
%--------------
%FMAP(j)
fmap=ones(mye,1); %fmap is indeed recalculated by program
%FMAP is changed by code
%-------------------------------------------------%
%Creating vector y ---GRID--------
%Using user matlab function to create the grid
%Check matlab script "create_mesh1D" for more info
y=create_mesh1D(ymax,'mixl',0.38,mye); 
%figure
%plot(diff(y),'r.')
%hold on
y=smooth(y,5)'; %make the mesh smoother!!
%plot(diff(y),'b.')
%y=create_mesh1D(length,'hype',5,mye);
%Calculate hy vector
hy(1)=y(2)-y(1);
for (j=2:mye-1)
    hy(j)=0.5*(y(j+1)-y(j-1));
end
hy(mye)=y(mye)-y(mye-1);
hy=hy;
%size of Dy
Dy=min(hy)
%============================================================%
%ZERO MODES
%by Rogers&Moser
%Error function
U=u0;%1/2*(U_top-U_bot)
dw0=1; %vorticiy thickness, set to 1
z1=sqrt(pi).*y/dw0;
dz1dy=sqrt(pi)/dw0;
u00=U*erf(z1);
u00p=2/sqrt(pi)*dz1dy*U*exp(-z1.^2);%first derivative
u00pp=u00p.*(-2*z1*dz1dy); %second derivative


%Checking profile u00(y)
%figure(1)
%hold on
%plot(y,u00,'b-')
%plot(y,u00p,'g-')
%plot(y,u00pp,'r-')
dm0=0.25*sum((1-u00.^2).*hy)

%title('u00 vs y')
%hold off
%w zero mode is zero.
w00=u00*0;


%===========================================================%
%            Calculation of initial perturbances
%--------------------MODE 1-------------------------------
%MODE (1,+-1) --> leads to the 3d disturbance
%kx, kz
nx1=1;%first mode
nz1=1;
kx1=nx1*alpe;
kz1=nz1*bete;
%Solve the Orr-Sommerfeld Equation
[phi,vor,ee1,vecty,v,dvdy]=ossq_for_mixing(kx1,kz1,Re,mye,y,u00,u00p,u00pp,ymax);
% Need to clean the solution to zero instead of noise!
% we need an eigenvector
mod=1;%THIS MODE IS WHAT MOSES says:
%positive phi in y=0, symmetric for real(phi) and asymmetric for imag(phi).
%We founded manually a complex number that makes the real part of phi symmetrical 
% so that the vortex will be in the middle of our box

%-pi/2;%THIS THETA CANT BE DIFFERENT

phimod1= exp(sqrt(-1)*theta)*phi(:,mod);%first mode for phi=0 and vor neq zero
vormod1=  exp(sqrt(-1)*theta)*vor(:,mod);        
v1     =   exp(sqrt(-1)*theta)* v(:,mod);
dvdy1  = exp(sqrt(-1)*theta)*dvdy(:,mod);

phimod1=smooth(phimod1,5)
vormod1=smooth(vormod1,5)        
v1     =smooth( v1  ,5)   
dvdy1  =smooth(dvdy1,5)

%Now clean noise, only for ymax=5 and u00=erf!!
%indleft=find(vecty<=-4);indright=find(vecty>=4);
%phimod1(indleft)=complex(0,0);
%phimod1(indright)=complex(0,0);
%vormod1(indleft)=complex(0,0);
%vormod1(indright)=complex(0,0);;
%v1(indleft)=complex(0,0);;
%v1(indright)=complex(0,0);;
%dvdy1(indleft)=complex(0,0);;
%dvdy1(indright)=complex(0,0);;

phi1 = interp1(vecty,phimod1,y,'spline','extrap');
vor1 = interp1(vecty,vormod1,y,'spline','extrap');
dvdy1= interp1(vecty,  dvdy1,y,'spline','extrap');
v1   = interp1(vecty,     v1,y,'spline','extrap');

%this calculations are only used for energy scaling, so it is not necessary
%to apply the phase change yet.

%from here we want to compute u and w solving an alegebra eq system
%1) From vor_y definition:
% i kz u - i kx w = vor_y
%2) From Continuity equation:
% i kx u + dvdy + ikz w = 0
% kx1,kz1=0 simplifies system
w1    = dvdy1*sqrt(-1)/(kz1+kx1^2/kz1);
u1    = w1*kx1/kz1;
%%%vorz1 =-phimod1*sqrt(-1)/kx1;

% Energy calculation of mode 1,1
A11=(2*2*sum(u1.*conj(u1).*hy+v1.*conj(v1).*hy+w1.*conj(w1).*hy))^0.5;
factor1=1/A11; %factor1 is the factor we have to multiply all in order to
                 %normalize energy to 1;

%normalizing phi1 and applying the phase change:
phi1=phi1*factor1;
A03D=0.0612;

phi1=phi1*A03D;

phi1real =real(phi1);
phi1imag =imag(phi1);

%Now clean noise, only for ymax=5 and u00=erf!!
%indleft=find(y<=-2.4);indright=find(y>=2.4);
%phi1real(indleft)=0;
%phi1real(indright)=0;

%for imag part
%indleft=find(y<=-2.4);indright=find(y>=2.4);
%phi1imag(indleft)=0;
%phi1imag(indright)=0;


vor1real =0*real(vor1);%ZERO
vor1imag =0*imag(vor1); %ZERO

% Order: Re(vor),Im(vor),Re(phi),Im(phi)
%Same way as assing put data:
wk1=zeros(2,2*mye,mze);
linez1 = nz1 + 1
planex1 = nx1+1
 %kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
% but matlab goes from 1 to mze
%MODE 0,1
	for j=1:mye
		wk1(1,2*j-1,linez1) =vor1real(j);
		wk1(1,2*j  ,linez1) =vor1imag(j);
%----------------------------------------------------!
                wk1(2,2*j-1,linez1) =phi1real(j);
 		wk1(2,2*j,  linez1) =phi1imag(j);
	end
%MODE 0,-1
	for j=1:mye
		wk1(1,2*j-1,mze+2-linez1) =vor1real(j);
		wk1(1,2*j  ,mze+2-linez1) =vor1imag(j);
%----------------------------------------------------!
                wk1(2,2*j-1,mze+2-linez1) =phi1real(j);
 		wk1(2,2*j,  mze+2-linez1) =phi1imag(j);
	end

%check size of wk1 (ntotr)	
%check size of wk1 (ntotr)	
%size(wk1(:))
ntotr=4*mye*mze %FOR EACH PLANE
%In this case kz=0 is the first line for kz
%Creating a complete plane for this input mode:
wkplane1=zeros(4*mye,mze);
wkplane1=wk1(:);
%wkplane2=zeros(4*mye,mze);
%wkplane2=wk2(:);

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

%nplanes only counts for half of the planes
nplanes=mxe/2;

%Write to file finp.dat
for i=1:nplanes
	if (i==planex1)	%write first mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane1(:),'float');
		fwrite(finp,rec2_size,'int');
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
figure(1)
plot(y,phi1real,'b*')
hold on
plot(y,phi1imag,'r*')
xlabel('y')
ylabel('Re/Im')
title('REAL+IMAG of mode 1 PHI')

%figure(3)
%plot(y,vor1real,'b*')
%hold on
%plot(y,vor1imag,'r*')
%xlabel('y')
%ylabel('Re/Im')
%title('REAL+IMAG of mode 1 VOR')

%CHECKING
%SECOND MODE
%figure(4)
%hold on
%plot(y,phi2real,'b*')
%plot(y,phi2imag,'r*')
%xlabel('y')
%ylabel('Re/Im')
%title('REAL+IMAG of mode 2 PHI')
%title('IMAG part of 2nd PHI mode:mode(0,1) and mode (0,-1)')

%figure(5)
%plot(y,vor2real,'r*')
%plot(y,vor2imag,'r*')
%hold on
%plot(y,vor2imag,'b*')
%xlabel('y')
%ylabel('Re/Im')
%title('REAL+IMAG of mode 2 VOR')



