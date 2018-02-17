%****************************************************************************%
% Matlab script to create starting file for LOMA			     %
% Antonio Almagro							     %
% June 2013								     %
%****************************************************************************%
clear all;
close all;
clc
%Relevant variables
%%The reference speed will be u2-u1=2*u2=u0
%u0=0.5;
%The reference speed will be DU
DU=1;
u0=DU/2;%u0 is umax
time=0;
%-------------------------------------------------%
%Caso MLPANTANO01
%----------------
Ree=160;%Rew=1370 as R&M
alpe=0.0182 %Pantano 345dm0
bete=0.0731/2 %Pantano 86 dm0
ymax=172/2; %half heigth of the box
%mxe and mze must be some multiple of 64*3/2
mye=513;
mgalx=768
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
y=create_mesh1D(ymax,'line',1.0,mye); 
y=smooth(y,5)'; %make the mesh smoother!!
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
u00=u0*tanh(y./(2*dm0));
%dm0 set to 1.
dm0=sum((0.25-(u00./DU).^2).*hy)

%w zero mode is zero.
w00=u00*0;
u=zeros(mgalx,mgalz,mye);
v=zeros(mgalx,mgalz,mye);
w=zeros(mgalx,mgalz,mye);

%test
%kk=zeros(256,mgalz);
%x=linspace(0,2*pi,256);
%create random fields
c0=max(u00(:))/10;%desviation for velocity
for k=1:mgalz
   for i=1:mgalx
%      u(i,k,:)=c0.*(-1+2.*rand(mye,1)').*exp(-(y./(2.*dm0)).^2);
      v(i,k,:)=c0.*(-1+2.*rand(mye,1)').*exp(-(y./(2.*dm0)).^2);
%      w(i,k,:)=c0.*(-1+2.*rand(mye,1)').*exp(-(y./(2.*dm0)).^2);
   end
   %%end
end
%uf1=fft2(u,mxe,mze);
vf1=fft2(v,mxe,mze);
%wf1=fft2(w,mxe,mze);


%===========================================================%
%            Calculation of initial perturbances
% Order: Re(vor),Im(vor),Re(phi),Im(phi)
%============================================================%
%===================WRITING TO FILE==========================%
%save record 1 values
finp=fopen('../data/finp.dat','w','l'); %open/create file tab_var.dat
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
wk1=zeros(2,2*mye,mze);
indreal=(2*(1:mye)-1);
%CASES---------------------------------------------------------------
%--------------------------------------------------------------------

dvdy=zeros(mye,1);
Dy=min(hy(:));
for i=1:nplanes
   for k=1:mze
       dvdy(2:end-1)=(vf1(i,k,3:end)-2.*vf1(i,k,2:end-1)+vf1(i,k,1:end-2))./(Dy.^2);
%fill the wk1 corresponding to real phi only
%!acc rding to LES paper      kwave=-2*sqrt(kxvector(i)^2+kzvector(k)^2)/k0;
%according Pantano is different
      kfac=(sqrt(kxvector(i)^2+kzvector(k)^2))^2;
%adding vory too
      wk1(1,indreal,k)  =zeros(mye,1);
      wk1(1,indreal+1,k)=zeros(mye,1);
%Adding phi
      wk1(2,indreal  ,k)=-kfac.*squeeze(real(vf1(i,k,:)))+real(dvdy);
      wk1(2,indreal+1,k)=-kfac.*squeeze(imag(vf1(i,k,:)))+imag(dvdy);
     %plot(y,wk1(2,indreal,k))
     %pause(0.01)
%      end
   end
%----------------------------------------------------!
 %  if (i<=modecut)
 %   wkplane=zeros(4*mye,mze);
 %  else
%Creating a complete plane for this input mode:
%Writing
    fwrite(finp,rec2_size,'int');
    fwrite(finp,wk1(:),'float');
    fwrite(finp,rec2_size,'int');
end
	
    fclose(finp);
%============================================================%
%===========END OF WRITING TO A FILE=========================%
%============================================================%


%
