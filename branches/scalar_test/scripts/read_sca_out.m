%****************************************************************************%
% only for LISOuc3m			                                     %
% Antonio Almagro							     %
% September 2013								     %
% Script for LISO							     %
%****************************************************************************%
%clear all;
%clc
%===================READING FROM FILE==========================%
%
my=513 %mlpantano13y2
ymax=10; %half heigth of the box
%ymax=80; %half heigth of the box
%ymax=10; %half heigth of the box
%mgalx=768;
%mgalz=288;
mgalx=192;
mgalz=192;
y=create_mesh1D(ymax,'manl',1 ,my); 
%mgalx=192;
%mgalz=192;
%y=create_mesh1D(ymax,'line',0.38 ,my); 
mx=mgalx*2/3;
mz=mgalz*2/3-1

%
%filename0='/data2/toni/mldiffSCAL01sc001.';
%filename0='/data2/toni/mlturbSCAL02sc001.';
filename0='/data2/toni/mlRMhiroll02sc001.';
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
%dummy1=fread(finp, 1,'int');
%%record 1 data 
%   time= fread(finp,1,'float');
%   Re= fread(finp,1,'float');
%   alp= fread(finp,1,'float');
%   bet= fread(finp,1,'float');
%   a0= fread(finp,1,'float');
%   mx= fread(finp,1,'int');
%   my= fread(finp,1,'int');
%   mz= fread(finp,1,'int');
%   y=zeros(my,1);
%   fmap=y;

%wk1==>u00, w00, odd/even
%for (j=1:my)
%     y(j) = fread(finp,1,'float64');
%     fmap(j) = fread(finp,1,'float64');  
%end
%dummy2=fread(finp, 1,'int');
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

plot(y,scalar,'b-');

pause(0.5)
end
