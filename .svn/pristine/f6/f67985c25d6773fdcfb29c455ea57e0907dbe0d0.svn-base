%This program writes once SCALAR FIELD
%CAUTION!! MAKE SURE to run FIRST the writing of finp
%AAlmagro
%=================================================================
%============================================================%
%===================WRITING TO FILE==========================%
%save record 1 values
finp=fopen('../data/finpsca.dat','w','l'); %open/create file tab_var.dat
%Record 1 - counting bytes

rec1_real4=5+2*mye;%time is real 4
%Ree,alpe,bete,a0e + y(j),fmap(j) (1:my)+wk1(j)(1:2my)
rec1_real8=mye*2;
%integer mxe mye mze
rec1_integer=3;
%total record 1 size
rec1_size=4*rec1_real4+4*rec1_integer;
%------end of record 1---------------------------------------%

%WRITE SCALARS SEPARETELY
ntotrSCA=2*mye*mze;%size of a plane for each scalar
iscal=1;
wkSCA=zeros(2*mye,mze);
%wkSCA=rand(2*mye,mze);
%try random number scalar field
recSCA_size=ntotrSCA*4;
C0=0.5;
%Creating a plane for scalar, separating real part and imaginary part
%     for k=1:mze
%	for j=1:mye
%		wkSCA(2*j-1, k) = 0;%REAL PART
%		wkSCA(2*j  , k) = 0; %IMAG PART
%	end
%     end

%----------------------------------------------------%
%Now only debugging, all scalars look the same initially
% furthermore, all planes x are the same
  cvarwkplane(:,iscal)=wkSCA(:);
  %MODE ZERO,ZERO and ZERO,k
  i=1;
  k=1;
  %draw profile of mode zero
  dm0=1; %momentum thickness, set to 1
  %z1=sqrt(pi).*y/dw0;
%ROGERS-MOSER CASE KELVIN-HELMHOLTZ

U=1;%1/2*(U_top-U_bot)
dw0=1; %vorticiy thickness, set to 1
z1=sqrt(pi).*y/dw0;
l=100;%not important?
%dz1dy=sqrt(pi)/dw0;
C00=0.5*(l+erf(z1));

%OPTION1
%  C00=smooth(cos(0.5*y).*(unitstep(y+2.5)-unitstep(y-2.5)),5);
%OPTION 2
%  C00=C0*(1+tanh(-y./(2*dm0)));
%  C00=C00.*(unitstep(y+50)-unitstep(y-50));%BC set to zero this way

  for j=1:mye
     wkSCA00(2*j-1,k)=C00(j);
     wkSCA00(2*j  ,k)=0.0;
  end
  %defining plane of mode 00
  cvar00plane(:,iscal)=wkSCA(:);
  cvar00plane(1:2*my,iscal)=wkSCA00(:);

  size(cvar00plane)
  %write modes 0,0:mz-1
       fwrite(finp,recSCA_size,'int');
       fwrite(finp,cvar00plane(:,iscal),'float');
       fwrite(finp,recSCA_size,'int');
  
  %write the rest of the modes
  for i=2:nplanes
       fwrite(finp,recSCA_size,'int');
       fwrite(finp,cvarwkplane(:,iscal),'float');
       fwrite(finp,recSCA_size,'int');
  end

fclose(finp)          

figure
plot(y,C00,'r*')



