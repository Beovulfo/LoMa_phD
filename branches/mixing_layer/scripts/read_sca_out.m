%****************************************************************************%
% only for LISOuc3m			                                     %
% Antonio Almagro							     %
% September 2013								     %
% Script for LISO							     %
%****************************************************************************%
%clear all;
%clc
%===================READING FROM FILE==========================%
%filename0='/data2/toni/mldiffSCAL01sc001.';
%filename0='/data2/toni/mlturbSCAL02sc001.';
%filename0='/data2/toni/mlRM2scal/hiroll01sc002.';
filename0='/data2/toni/test02sc002.';
%filename0='/data2/toni/mlRMSCAL01sc001.';
ibeg=1
iend=1
figure(1)
hold on
for ii=ibeg:iend
  if ii<10
     string2=strcat('00',num2str(ii));
  elseif ii<100
     string2=strcat('0',num2str(ii));
  else
     string2=num2str(ii);
  end
  filename=strcat(filename0,string2)

  finp=fopen(filename,'r','l');
%----------------------Record 1----------------------------%
%%Buffer size    
dummy1=fread(finp, 1,'int');
%%record 1 data 
   time= fread(finp,1,'float');
   Re= fread(finp,1,'float');
   alp= fread(finp,1,'float');
   bet= fread(finp,1,'float');
   a0= fread(finp,1,'float');
   mx= fread(finp,1,'int');
   my= fread(finp,1,'int');
   mz= fread(finp,1,'int');
   y=zeros(my,1);
   fmap=y;

%wk1==>u00, w00, odd/even
for (j=1:my)
     y(j) = fread(finp,1,'float64');
     fmap(j) = fread(finp,1,'float64');  
end
%modes cero
for (j=1:my)
        u00(j)=fread(finp,1,'float');
        w00(j)=fread(finp,1,'float');
end
%
dummy2=fread(finp, 1,'int');
%display([dummy1,dummy2])
nplanes=mx/2;
%dummy3=fread(finp,1,'int')
%keyboard
ntotr=2*my*mz;

 for (i=1:nplanes)
     dummy3=fread(finp,1,'int');
     wk1(i,:)=fread(finp,ntotr,'float');
     dummy4=fread(finp,1,'int');
 end

 fclose(finp)

 scalar=zeros(1,my);
 for (j=1:my)
    scalar(j)=wk1(1,2*j-1,1);
 end


plot(y,scalar,'r-');

pause(0.5)
end

mgalx=mx*3/2;
mgalz=(mz+1)*3/2;
%start transformation
scal_phys=zeros(mgalx,my,mgalz);
%1)Sum real part and imaginary part in order to create the complex number
scal_complex=wk1(:,2*(1:my)-1,:)+sqrt(-1)*wk1(:,2*(1:my),:);
%2)For each XZ plane we are going to Transfrom from FOU2FIS
for j=1:my
   tempFF=squeeze(scal_complex(:,j,:));
   %2.1) transform ifft for Z
   tempFP=ifft(tempFF,mgalz,2);
   %2.2) transform ifft for x
   tempPP=ifft(tempFP,mgalx,1,'symmetric');
   scal_phys(:,j,:)=tempPP*mgalz*mgalx;
end

%Calculate 1/scal_phys
	kk=1./scal_phys;
kk_fou=zeros(mx,my,mz);
for j=1:my
        tempPP2=squeeze(kk(:,j,:));
        kk_fou1=fft2(tempPP2,mx,mz)./(mx*mz);
        kk_fou(:,j,:)=kk_fou1(:,:);
end

scal_phys2=zeros(mgalx,my,mgalz);
%1)Sum real part and imaginary part in order to create the complex number
%2)For each XZ plane we are going to Transfrom from FOU2FIS
clear tempFF,tempFP,tempPP;
for j=1:my
   tempFF=squeeze(kk_fou(:,j,:));
   %2.1) transform ifft for Z
   tempFP=ifft(tempFF,mgalz,2);
   %2.2) transform ifft for x
   tempPP=ifft(tempFP,mgalx,1,'symmetric');
   scal_phys2(:,j,:)=tempPP*mgalz*mgalx;
end

checking=max(scal_phys2-kk);

