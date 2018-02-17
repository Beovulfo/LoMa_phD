%****************************************************************************%
% Matlab script to read .sta file for LISOuc3m		     %
% Antonio Almagro							     %
% May 2013								     %
% 
% Script for LISO							     %
%****************************************************************************%
%clear all;
%close all;
%clc
%============================================================%
%===================READING FROM FILE==========================%
%
function [y,fmap,fwkreal,fwkimag]=leekk(filename,my)
%finp=fopen('/data/toni/mixing/modesx2.055','r','l'); %open/create file tab_var.dat
%filename=input('filename:\n')
%
finp=fopen(filename,'r','l');
%frewind(finp);
%----------------------Record 1----------------------------%
%Buffer size    


dummy1=fread(finp, 1,'int')
   y=zeros(my,1);
   fmap=y;
%wk1==>u00, w00, odd/even
   for (j=1:my)
       y(j) = fread(finp,1,'float64');
       fmap(j) = fread(finp,1,'float64');  
   end
dummy2=fread(finp, 1,'int')

%RECORD 2
dummy3=fread(finp, 1,'int');

    tamy=fread(finp,1,'int');
    tamz=fread(finp,1,'int');
    tamx=fread(finp,1,'int');

fwkreal=zeros(tamy,tamz,tamx);
fwkimag=zeros(tamy,tamz,tamx);
dummy4=fread(finp,1,'int')

%========================================%
dummy5=fread(finp,1,'int')

for (i=1:tamx)
  for (k=1:tamz)
     for (j=1:tamy)
       fwkreal(j,k,i) = fread(finp,1,'float');
     end
  end
end
dummy6=fread(finp, 1,'int')
dummy5=fread(finp,1,'int')

for (i=1:tamx)
  for (k=1:tamz)
     for (j=1:tamy)
       fwkimag(j,k,i) = fread(finp,1,'float');
     end
  end
end
dummy6=fread(finp, 1,'int')


fclose(finp)
