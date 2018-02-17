%=============================================================
%This program writes one SCALAR FIELD
%aalmagro
%=================================================================
%create/open file format: finpsc###.dat, where ### is the iscal number.
finp=fopen('../data/finpsca.dat','w','l'); %open/create file tab_var.dat
mye=513;mgalx=192;mgalz=192;my=mye;
Ree=500;alpe=0.8621;bete=1.4368;

fmap=ones(mye,1);

ymax=10;
y=create_mesh1D(ymax,'line',0.38,mye); 

time=0;a0e=0;
u00=zeros(1,mye);
w00=zeros(1,mye);

mxe=mgalz*2/3;mze=mgalz*2/3-1;
nplanes=mxe/2;
%==========WRITE HEADER==========
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
%fwrite(finp, rec1_size,'uint');
%record 1 data 
%    fwrite(finp,time,'float');
%    fwrite(finp,Ree,'float');
%    fwrite(finp,alpe,'float');
%    fwrite(finp,bete,'float');
%    fwrite(finp,a0e,'float');
%    fwrite(finp, mxe,'int');
%    fwrite(finp, mye,'int');
%   fwrite(finp, mze,'int');
%write y(j) and fmap(j)
%for (j=1:mye)
%    fwrite(finp,y(j),'double');
%    fwrite(finp,fmap(j),'double');
%end
%wk1==>u00, w00, odd/even
%for (j=1:mye)
%	fwrite(finp,u00(j),'float');
%	fwrite(finp,w00(j),'float');
%end
%write buffer size again
%fwrite(finp, rec1_size,'uint');
%------end of record 1---------------------------------------%



%======END OF HEADER====================%



ntotrSCA=2*mye*mze;%size of a plane for each scalar
wkSCA=zeros(2*mye,mze);
cvarwkplane=zeros(2*mye,mze);
cvar00plane=zeros(2*mye,mze);
recSCA_size=ntotrSCA*4;
%ROGERS-MOSER CASE KELVIN-HELMHOLTZ
%MODE 00 SCALAr
U=1;%1/2*(U_top-U_bot)
dw0=1; %vorticiy thickness, set to 1
z1=sqrt(pi).*y/dw0;
l=1;%not important,RM use 0.5
%dz1dy=sqrt(pi)/dw0;
%C00=0.5*(l+erf(z1));%iscal=001
C00=l+0.01*erf(z1);%iscal=002
%defining mode 00
  for j=1:mye
     wkSCA00(2*j-1,1)=C00(j);
     wkSCA00(2*j  ,1)=0.0;
  end
  %defining plane of mode 00
  cvar00plane(1:2*my)=wkSCA00(:);%Mode (0,0)
  %====writing modes (0,0:mz-1)
  fwrite(finp,recSCA_size,'int');
  fwrite(finp,cvar00plane(:),'float');
  fwrite(finp,recSCA_size,'int');
  
  %====writing the rest of the modes
  for i=2:nplanes
       fwrite(finp,recSCA_size,'int');
       fwrite(finp,cvarwkplane(:),'float');
       fwrite(finp,recSCA_size,'int');
  end

fclose(finp)          

figure
plot(y,C00,'r*')



