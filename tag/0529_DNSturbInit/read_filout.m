function [y,u00,w00,time]=read_filout(filename)
%****************************************************************************%
% Matlab script to read re-start file header only for LISOuc3m			     %
% Antonio Almagro							     %
% May 2013								     %
% Ready for input 2 different modes
% Script for LISO							     %
%****************************************************************************%
%clear all;
%close all;
%clc
%============================================================%
%===================READING FROM FILE==========================%
%
%finp=fopen('/data/toni/mixing/modesx2.055','r','l'); %open/create file tab_var.dat
%
finp=fopen(filename,'r','l');
%----------------------Record 1----------------------------%
%Buffer size    
dummy1=fread(finp, 1,'int');
%record 1 data 
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
   u00=y;
   w00=y;

%wk1==>u00, w00, odd/even
for (j=1:my)
     y(j) = fread(finp,1,'float64');
     fmap(j) = fread(finp,1,'float64');  
end
for (j=1:my)
	u00(j)=fread(finp,1,'float');
	w00(j)=fread(finp,1,'float');
end
dummy2=fread(finp, 1,'int');
%display([dummy1,dummy2])
fclose(finp)