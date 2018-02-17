%****************************************************************************%
% Matlab script to create starting file for LOMA			     %
% Antonio Almagro							     %
% June 2013								     %
%****************************************************************************%
clear all;
close all;
clc
%IC noise
%C=0.001; %for phi
%C2=0.001; %C2 for omega y
pert=0.1;%perturbation between modes x
C=0.001; %Controls magnitud for phi
Cy=0.1; %Controls variation introduced as function of y
%pert=0.001;%Controls variations with respect modes 
%C2=0.000002;%C2 for omega y
C2=0.;%C2 for omega y
Nmaxy=2;%maximum mode number @ y
%k0=0.174; %k0=kmax/48
k0=48;
%Relevant variables
%The reference speed will be u2-u1=2*u2=u0
%u0=0.5;
%The reference speed will be DU=1
DU=1
u0=DU/2;
%Re=1778;%Rem=800 as R&M
%Re=1370;%Rew=1370 as R&M
%Re0=U*dw0/nu
%ACHTUNG MESH NOW FROM -1 TO 1
%init time
time=0;
%-------------------------------------------------%
%Caso MLPANTANO01
%----------------
Ree=160;%Rew=1370 as R&M
alpe=0.0182 %Pantano 345dm0 
bete=0.0731 %Pantano 86 dm0 
%ymax=120; %half heigth of the box
ymax=172/2; %half heigth of the box
%mxe and mze must be some multiple of 64*3/2
%mye=513 mlpantano13
mye=1025 %mlpantano13y2
mgalx=768
%mgalz=128*3/2 %mlpantano13
mgalz=288  
%alpe=2*pi/Lx %tamaño de caja periodica en X (alpe=2PI/Lx)
%bete=2*pi/Lz %tamaño de caja periodiza en z (bete=2PI/Lz)
%Re0=U*dw0/nu
%-------------------------------------------------%
aOe=0;%velocidad relativa de observacion
%physical mesh size
Lx=2*pi/alpe;
Lz=2*pi/bete;
Dx=Lx/(mgalx-1)
Dz=Lz/(mgalz-1)
%convert to number of Fourier modes
mxe=mgalx*2/3
mze=mgalz*2/3-1
my=mye
%Size of Comp. Box
Volum=Lx*Lz*ymax*2;
%kx kz vectors
nxvector=(0:mxe/2-1);
nzvector=[(0:(mze-1)/2) (-(mze-1)/2:1:-1)];
kxvector=nxvector*alpe;
kzvector=nzvector*bete;
%--------------
%FMAP(j)
fmap=ones(mye,1); %fmap is indeed recalculated by program
%FMAP is changed by code
%-------------------------------------------------%
%Creating vector y ---GRID--------
%Using user matlab function to create the grid
%Check matlab script "create_mesh1D" for more info
%y=create_mesh1D(ymax,'mixl', 0.3,mye); 
y=create_mesh1D(ymax,'file',1 ,mye); 
%y=create_mesh1D(ymax,'line',1.0,mye); 
%y=smooth(y,5)'; %make the mesh smoother!!
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
U=u0;%1/2*(U_top-U_bot)
%dw0=1; %vorticiy thickness, set to 1
%Error function
%z1=sqrt(pi).*y/dw0;
%dz1dy=sqrt(pi)/dw0;
%u00=U*erf(z1);
%u00p=2/sqrt(pi)*dz1dy*U*exp(-z1.^2);%first derivative
%u00pp=u00p.*(-2*z1*dz1dy); %second derivative
%----
%tanh
dm0=1.0;
u00=u0*tanh(-y./(2*dm0));
%dm0 set to 1.
dm0=sum((0.25-(u00./DU).^2).*hy)

%Checking profile u00(y)
%figure(1)
%hold on
%plot(y,u00,'b-')
%plot(y,u00p,'g-')
%plot(y,u00pp,'r-')
dm0=sum((1/4-(u00/DU).^2).*hy)

%title('u00 vs y')
%hold off
%w zero mode is zero.
w00=u00*0;

%===========================================================%
%            Calculation of initial perturbances
% Order: Re(vor),Im(vor),Re(phi),Im(phi)
%============================================================%
%===================WRITING TO FILE==========================%
%save record 1 values
finp=fopen('../data/finpINCturb.dat','w','l'); %open/create file tab_var.dat
%finp=fopen('../data/finp2.dat','w','l'); %open/create file tab_var.dat
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

ntotr=4*mye*mze %FOR EACH PLANE
rec2_size=ntotr*4; %SAME SIZE FOR EACH PLANES
%Buffer size
%nplanes only counts for half of the planes
nplanes=mxe/2;
wkplane1=zeros(4*mye,mze);
wk1=zeros(2,2*mye,mze);
indreal=(2*(1:mye)-1);
%CASES---------------------------------------------------------------
%--------------------------------------------------------------------
%MLPANTANO01
%%C=0.0289*(0.01/0.0145)^0.5
%%k0=8;%Really important factor to control magnitud of turbulence
%--------------------------
%MLPANTANO02
%%C=0.001*(0.01/5E-7)^0.5; %kfac*phi
%%C2=0.;
%%k0=2.557^4; %k0=kmax

%prueba
%Really important factor to control magnitud of turbulence
%C= 0.0240*(0.015/0.0005415)^(0.5);
%k0=9;%Really important factor to control magnitud of turbulence
%--------------------------------------------------------------------

%Write to file finp.dat
%figure(1)
%hold on
%modecut=4;
contador=0;
modecut=2;
wk1=zeros(2,2*mye,mze);
%figure;hold on;
for i=1:nplanes
   if (nxvector(i)<=modecut)
      wk1=zeros(2,2*mye,mze);
   else
      for k=1:mze
         if (abs(nzvector(k))<1) %eliminate 2D modes
       		wk1(1,indreal,k)  =zeros(my,1);
       		wk1(1,indreal+1,k)=zeros(my,1);
       		wk1(2,indreal,k)  =zeros(my,1);
       		wk1(2,indreal+1,k)=zeros(my,1);
      	 else
         %fill the wk1
      		contador=contador+1;
      		kfac=sqrt(nxvector(i).^2+nzvector(k).^2)./k0;
      		kvector(contador)=(kxvector(i)^2+kzvector(k).^2)^0.5;
%      kwave=-2*kfac^2;%pantano
      		kwave=-2*kfac^2;
      %		kwave=-2*kfac^2;
        %adding waves on "y2
        %initializing
              funcy1=ones(mye,1);
              funcy2=funcy1;
              funcy3=funcy1;
              funcy4=funcy1;
              for (l=1:Nmaxy)
                 funcy1=funcy1+Cy*(2*rand(1)-1)*cos(l*y'/(2*pi)+2*rand(1)*pi);
                 funcy2=funcy2+Cy*(2*rand(1)-1)*cos(l*y'/(2*pi)+2*rand(1)*pi);
                 funcy3=funcy3+Cy*(2*rand(1)-1)*cos(l*y'/(2*pi)+2*rand(1)*pi);
                 funcy4=funcy4+Cy*(2*rand(1)-1)*cos(l*y'/(2*pi)+2*rand(1)*pi);
              end
        %exponential decay to restrict noise within vorticity thickness
                sigma=2.934*dm0;
      		funcy1=funcy1.*exp(-0.5*(y'./sigma).^2);
      		funcy2=funcy2.*exp(-0.5*(y'./sigma).^2);
      		funcy3=funcy3.*exp(-0.5*(y'./sigma).^2);
      		funcy4=funcy4.*exp(-0.5*(y'./sigma).^2);
%      		funcy1=(pert*(-1+2*rand(mye,1)')+1).*exp(-(y./(2.*dm0)).^2);
%      		funcy2=(pert*(-1+2*rand(mye,1)')+1).*exp(-(y./(2.*dm0)).^2);
%      		funcy3=(pert*(-1+2*rand(mye,1)')+1).*exp(-(y./(2.*dm0)).^2);
%      		funcy4=(pert*(-1+2*rand(mye,1)')+1).*exp(-(y./(2.*dm0)).^2);
%                %VOR y
%      		wk1(1,indreal,k)  =C2*(-1+2*rand(1)).*funcy3.*exp(kwave);
%      		wk1(1,indreal+1,k)=C2*(-1+2*rand(1)).*funcy4.*exp(kwave);
%                %PHI
%      		wk1(2,indreal,k)  =C*(-1+2*rand(1)).*funcy1.*exp(kwave);
%      		wk1(2,indreal+1,k)=C*(-1+2*rand(1)).*funcy2.*exp(kwave);
         %       %VOR y C2=0 neglected.
      	%	wk1(1,indreal,k)  =C2*(2*rand(1)-1)*funcy3.*exp(kwave);
      %		wk1(1,indreal+1,k)=C2*(2*rand(1)-1)*funcy4.*exp(kwave);
                %VOR y white noise C2=0 neglected.
      		wk1(1,indreal,k)  =C2*(2*rand(1)-1)*funcy3.*exp(kwave);
      		wk1(1,indreal+1,k)=C2*(2*rand(1)-1)*funcy4.*exp(kwave);
                %PHI
%      		wk1(2,indreal,k)  =C*(2*rand(1)-1)*funcy1.*exp(kwave);
%      		wk1(2,indreal+1,k)=C*(2*rand(1)-1)*funcy2.*exp(kwave);
                wk1(2,indreal,k)  =C*(kvector(contador))^2*(2*rand(1)-1)*funcy1.*exp(kwave);
%.*exp(kwave);
        	wk1(2,indreal+1,k)=C*(kvector(contador))^2*(2*rand(1)-1)*funcy2.*exp(kwave);
%;.*exp(kwave);
%                %PHI (whitenoise)
%      		wk1(2,indreal,k)  =C*(kvector(contador))^2*(2*rand(1)-1)*funcy1;
%      		wk1(2,indreal+1,k)=C*(kvector(contador))^2*(2*rand(1)-1)*funcy2;
%mlpantano16      		wk1(2,indreal,k)  =C*kfac^4*(2*rand(1)-1)*funcy1.*exp(kwave);
%      		wk1(2,indreal+1,k)=C*kfac^4*(2*rand(1)-1)*funcy2.*exp(kwave);
%mlpantano15	wk1(2,indreal,k)  =C*kfac^2*(2*rand(1)-1)*funcy1.*exp(kwave);
%      	  	wk1(2,indreal+1,k)=C*kfac^2*(2*rand(1)-1)*funcy2.*exp(kwave);
      	 end
      end
%----------------------------------------------------!
%Creating a complete plane for this input mode:
   end
%Writing
     
   % plot(y,wk1(2,indreal,41));

    fwrite(finp,rec2_size,'int');
    fwrite(finp,(pert*(2*rand(1)-1)+1)*wk1(:),'float');
    fwrite(finp,rec2_size,'int');
end
	
    fclose(finp);
%============================================================%
%===========END OF WRITING TO A FILE=========================%
%============================================================%
