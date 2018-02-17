%****************************************************************************%
% Matlab script to create starting file for LOMA			     %
% Antonio Almagro							     %
% June 2013								     %
% July 2014 --> modified for Compressible LoMA  			     %
%****************************************************************************%
clear all;
close all;
clc

%-----------------------------------------------------------------%
%Relevant variables
% ----------------------------------------------------------------%
%The reference speed will be u2-u1=2*u2=u0
u0=1;
Re=500;
mye=513;mgalx=192;mgalz=192;
%Re0=U*dw0/nu
ymax=10; %half heigth of the box
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


%-------------------------------------------------%
%Creating vector y ---GRID--------
%Using user matlab function to create the grid
%Check matlab script "create_mesh1D" for more info
%y=create_mesh1D(ymax,'mixl',0.38,mye); 
y=create_mesh1D(ymax,'line',0.38,mye); 
%y=create_mesh1D(ymax,'manl',0.38,mye); 
figure
plot(diff(y),'r.')
%hold on
%y=smooth(y,5)'; %make the mesh smoother!!
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
%00 MODES
%by Rogers&Moser
%Error function
%============================================================%
%00 MODES
U=u0;%1/2*(U_top-U_bot)
dw0=1; %vorticiy thickness, set to 1
z1=sqrt(pi).*y/dw0;
dz1dy=sqrt(pi)/dw0;
%Now u00 is rho*u00

%TEMPERATURE 
%-----------------------------------------------------------------%
% Defining Temperature Stuff (SCAL)                               %
% ----------------------------------------------------------------%
Tmean = 1;
DT = 0; %INCOMPRESSIBLE CASE AGAIN
DT= DT/2;
T00=Tmean+DT*erf(sqrt(pi).*y);

figure(10)
plot(y,T00)
title('T00')
%-----------------------------------------------------------------%

up00 = U.*erf(z1);
figure(11)
plot(y,up00)
title('up00')

r00=1./T00;
figure(12)
drho = max(r00)-min(r00) 
plot(y,r00)
title('rho00')

u00=r00.*U.*erf(z1);
%u00=U.*erf(z1); %Incompressible
figure(13)
plot(y,u00)
title('rhou00')
%Incompressible:
%-------------
%u00=U.*erf(z1);
u00p=2/sqrt(pi)*dz1dy*U*exp(-z1.^2);%first derivative
u00pp=u00p.*(-2*z1*dz1dy); %second derivative

dm0=0.25*sum((1-up00.^2).*hy)

%Other 00 profiles
v00=u00*0;
w00=u00*0;

%===========================================================%
%            Calculation of initial perturbances
%===========================================================%

%--------------------MODE 1-------------------------------

%MODE (1,0) --> leads to the Kelvin-Helmholtz rollup
%kx, kz
nx1=1;%first mode
nz1=0;
kx1=nx1*alpe;
kz1=nz1*bete;
%Solve the Orr-Sommerfeld Equation
[phi,vor,ee1,vecty,v,dvdy]=ossq_for_mixing(kx1,kz1,Re,mye,y,u00,u00p,u00pp,ymax);
% Need to clean the solution to zero instead of noise!

%this is used in order to get the stability eigenfunction for the vorticity.
%we want to introduce vorticity z (wz=A10*Real(f(y))
% we need an eigenvector
mod=1;%THIS MODE IS WHAT MOSES says:
%mod=29;%THIS MODE IS WHAT MOSES says:
%positive phi in y=0, symmetric for real(phi) and asymmetric for imag(phi).
%We founded manually a complex number that makes the real part of phi symmetrical 
% so that the vortex will be in the middle of our box
theta=0.412-pi/2;

phimod1= exp(sqrt(-1)*theta)*phi(:,mod);%first mode for phi=0 and vor neq zero
vormod1=  exp(sqrt(-1)*theta)*vor(:,mod);        
v1     =   exp(sqrt(-1)*theta)* v(:,mod);
dvdy1  = exp(sqrt(-1)*theta)*dvdy(:,mod);

%phimod1=smooth(phimod1,5)
%vormod1=smooth(vormod1,5)        
%v1     =smooth( v1  ,5)   
%dvdy1  =smooth(dvdy1,5)

%Now clean noise, only for ymax=5 and u00=erf!!
indleft=find(vecty<=-4);indright=find(vecty>=4);
phimod1(indleft)=complex(0,0);
phimod1(indright)=complex(0,0);
vormod1(indleft)=complex(0,0);
vormod1(indright)=complex(0,0);;
v1(indleft)=complex(0,0);;
v1(indright)=complex(0,0);;
dvdy1(indleft)=complex(0,0);;
dvdy1(indright)=complex(0,0);;

phi1 = interp1(vecty,phimod1,y,'spline','extrap');
vor1 = interp1(vecty,vormod1,y,'spline','extrap');
dvdy1= interp1(vecty,  dvdy1,y,'spline','extrap');
v1   = interp1(vecty,     v1,y,'spline','extrap');

%this calculations are only used for energy scaling, so it is not necessary
%to apply the phase change yet.

%from here we want to compute u and w solving an alegebra eq system
%1) From vor_y definition:
% i kz u - i kz w = vor_y
%2) From Continuity equation:
% i kx u + dvdy + ikz w = 0
% kx1,kz1=0 simplifies system
w1    = sqrt(-1)*vor1/kx1;
u1    = sqrt(-1)*dvdy1/kx1;
%vorz1 =-phimod1*sqrt(-1)/kx1;

% Energy calculation of mode 1
A10=(2*sum(u1.*conj(u1).*hy+v1.*conj(v1).*hy+w1.*conj(w1).*hy))^0.5;
factor1=1/A10; %factor1 is the factor we have to multiply all in order to
                 %normalize energy to 1;

%normalizing phi1 and applying the phase change:
phi1=phi1*factor1;


%from Rogers&Moser phi calculated from ome_z should be
% phi=f(y)*i*kx*A01 (i*f(y) is the eigenvector of omega
phi1=phi1*0.1;

phi1real =real(phi1);
phi1imag =imag(phi1);

vor1real =0*real(vor1);%ZERO
vor1imag =0*imag(vor1); %ZERO



%---------------------MODE 2----------------------------
%kx, kz
nx2=0;%in this case nx2 must be different than nx1
nz2=1;
kx2=nx2*alpe
kz2=nz2*bete

%As read from Moser, second mode is used to give 3D
% ome_x=A*exp(-pi*y^2), we need to introduce a phi able to
% aproach the BC -> phi

%2nd mode is (0,1) , phi=-kz2*A01*c*exp(-pi*y^2) --> Moses
%needs to be calculated with energy
%
rhs1                 = -exp(-pi*y.^2);%(0, 1) phi
%rhs2                 = -exp(-pi*y.^2);%(0,-1) phi
[vecty2 ,v2 ,dvdy2 ] = lapsolver(kx2, kz2,mye,y,ymax,rhs1);
%[vecty2n,v2n,dvdy2n] = lapsolver(kx2,-kz2,mye,y,ymax,rhs2);

%interpolate
v2    = interp1(vecty2,   v2,y,'spline','extrap');
dvdy2 = interp1(vecty2,dvdy2,y,'spline','extrap');
u2    = zeros(1,mye);%u2 is zero
w2    = sqrt(-1)*dvdy2/kz2;
%mode -kz2
%v2n    = interp1(vecty2n,   v2n,y,'spline','extrap');
%dvdy2n = interp1(vecty2n,dvdy2n,y,'spline','extrap');
%u2n    = zeros(1,mye);%u2 is zero
%w2n    = -sqrt(-1)*dvdy2/kz2;
A01   = (2*sum(v2.*conj(v2).*hy+w2.*conj(w2).*hy))^0.5;%ml4
c=1/A01; %c achieves integrated energy equals 1.
%Calculate c:
%Now we need to compute A3d in order to scale properly
A3D0=0.0166;%A3D0=A01 (LOROLL)
%A3D0=0.0;%A3D0=A01 (HIROLL)
%A3D0=0.3327;%A3D0=A01 (HIROLL)
A01=A3D0; %in this case

phi2real =  A01*c*rhs1;
phi2imag =  0*A01*c*rhs1;%ml1

vor2real =  0*c*vor1real;
vor2imag =  0*c*vor1imag;
%mode (0,-1)
phi3real =  A01*c*rhs1;
phi3imag =  phi2imag;
vor3real =  vor2real;
vor3imag =  vor2imag;
%====================================================================%
%wk00==> line for kx=0, kz=0
% 		MODE 00 PLANE (needed for SCAL only)
%		--------------
% Order: Re(vor),Im(vor),Re(phi),Im(phi),Re(psi),Im(Psi),Re(scal),Im(scal)
%Same way as assing put data:
wk00=zeros(4,2*mye,mze);
 %kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
% but matlab goes from 1 to mze
	for j=1:mye
                %SCAL 00 MODE
%----------------------------------------------------%
       	        wk00(4,2*j-1, 1) = T00(j);%real SCAL
 		wk00(4,2*j  , 1) = 0.0;
	end



%wk1==> line for kx=1, kz=0
% Order: Re(vor),Im(vor),Re(phi),Im(phi),Re(psi),Im(Psi),Re(scal),Im(scal)
%Same way as assing put data:
wk1=zeros(4,2*mye,mze);
linez1 = nz1 + 1
planex1 = nx1+1
 %kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
% but matlab goes from 1 to mze
	for j=1:mye

                %VOR
%----------------------------------------------------%
		wk1(1,2*j-1,linez1) =vor1real(j);
		wk1(1,2*j  ,linez1) =vor1imag(j);

                %PHI
%----------------------------------------------------!
                wk1(2,2*j-1,linez1) =phi1real(j);
 		wk1(2,2*j,  linez1) =phi1imag(j);

                %PSI
%----------------------------------------------------%
       	        wk1(3,2*j-1, linez1) = 0.0;
 		wk1(3,2*j  , linez1) = 0.0;

                %SCAL
%----------------------------------------------------%
       	        wk1(4,2*j-1, linez1) = 0.0;
 		wk1(4,2*j  , linez1) = 0.0;
	end

%planex 2 = 1 --> THIS IS THE FIRST PLANE
%VIP :: THIS PLANE INCLUDES 00 MODE OF SCAL!!!
wk2=zeros(4,2*mye,mze);
linez2 = nz2 + 1
planex2 = nx2+1
%kz=0 means k=1 (first plane) 
%Z lines goes from 0 (first) to mze-1 (last) [FORTRAN], 
%% but matlab goes from 1 to mze
	for j=1:mye
                %VOR
%----------------------------------------------------%
		wk2(1,2*j-1, linez2) = vor2real(j);
		wk2(1,2*j  , linez2) = vor2imag(j);

                %PHI
%----------------------------------------------------%
       	        wk2(2,2*j-1, linez2) = phi2real(j);
 		wk2(2,2*j  , linez2) = phi2imag(j);

                %PSI
%----------------------------------------------------%
       	        %wk2(3,2*j-1, linez2) = 0.0;
 		%wk2(3,2*j  , linez2) = 0.0;

                %SCAL: 00 MODE
%----------------------------------------------------%
       	        wk2(4,2*j-1, 1) = T00(j);
 		wk2(4,2*j  , 1) = 0.001;

	end

%add the mode (0,-1)---------------------------------------VIP
%mod_index=mze;
	for j=1:mye

                %VOR
%----------------------------------------------------%
		wk2(1,2*j-1, mze) = vor3real(j);
		wk2(1,2*j  , mze) = vor3imag(j);

                %PHI
%----------------------------------------------------%
       	        wk2(2,2*j-1, mze) = phi3real(j);
 		wk2(2,2*j  , mze) = phi3imag(j); %ml3

                %PSI
%----------------------------------------------------%
       	        wk2(3,2*j-1, mze) = 0.0;
 		wk2(3,2*j  , mze) = 0.0;

                %SCAL
%----------------------------------------------------%
       	        wk2(4,2*j-1, mze) = 0.0;
 		wk2(4,2*j  , mze) = 0.0;
	end

%check size of wk1 (ntotr)	
%check size of wk1 (ntotr)	
%size(wk1(:))
ntotr=8*mye*mze %FOR EACH PLANE
%In this case kz=0 is the first line for kz
%Creating a complete plane for this input mode:
wkplane0=zeros(8*mye,mze); %ikx=0 (FIRST plane)
wkplane0=wk00(:);

wkplane1=zeros(8*mye,mze);
wkplane1=wk1(:);

wkplane2=zeros(8*mye,mze);
wkplane2=wk2(:);

%All other lines/planes should be zero
zero=0*wk1(:); %empty plane


%============================================================%
%===================WRITING TO FILE==========================%
%save record 1 values
finp=fopen('../data/finp.dat','w','l'); %open/create file tab_var.dat
%Record 1 - counting bytes
rec1_real4=5+3*mye;%time is real 4
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
%wk1==>u00, v00, w00
for (j=1:mye)
	fwrite(finp,u00(j),'float');
	fwrite(finp,v00(j),'float');
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

%nplanes only counts for half of the planes
nplanes=mxe/2;

%Write to file finp.dat
%--------------------------------------

for i=1:nplanes
	if (i==planex1)	%write first mode
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane1(:),'float');
		%fwrite(finp,zero(:),'float');
		fwrite(finp,rec2_size,'int');
        elseif (i==planex2) %write second mode 
                %FIRST PLANE including 00 mode of SCAL
		fwrite(finp,rec2_size,'int');
		fwrite(finp,wkplane2(:),'float');
		%fwrite(finp,wkplane0(:),'float');%For 1D case we still need SCAL
		fwrite(finp,rec2_size,'int');
	else
		fwrite(finp,rec2_size,'int');
		fwrite(finp,zero(:),'float');
		fwrite(finp,rec2_size,'int');
	end
end
	
fclose(finp);
%============================================================%
%===========END OF WRITING TO A FILE=========================%
%============================================================%



