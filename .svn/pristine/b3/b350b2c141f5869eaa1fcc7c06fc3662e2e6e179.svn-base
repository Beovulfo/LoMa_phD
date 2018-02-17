%****************************************************************************%
% Matlab script to read velocity field                           	     %
% Antonio Almagro							     %
% May 2013								     %
% Script for LISO							     %
%****************************************************************************%
%clear all;
%close all;
%clc
%============================================================%
%===================READING FROM FILE==========================%
%
%filename=ml8.001
function [time,x,y,z,posy,wk1r]=readfieldxz(filename)

finp=fopen(filename,'r','l');
%----------------------Record 1----------------------------%
%Buffer size    
dummy1=fread(finp, 1,'int');
%record 1 data 
   time= fread(finp,1,'float');
   Re= fread(finp,1,'float');
   alp= fread(finp,1,'float');
   bet= fread(finp,1,'float');
%   a0= fread(finp,1,'float');
   mgalx= fread(finp,1,'int')
   my= fread(finp,1,'int')
   mgalz= fread(finp,1,'int')
   nspec=fread(finp,1,'int');
   nacum=fread(finp,1,'int');

dummy2=fread(finp, 1,'int');
   y=zeros(my,1);
   fmap=y;
   u00=y;
   w00=y;

dummy3=fread(finp, 1,'int');
%wk1==>u00, w00, odd/even
for (jj=1:nspec)
    jspec(jj)=fread(finp,1,'int');
end

for (j=1:my)
     y(j) = fread(finp,1,'float64');
     fmap(j) = fread(finp,1,'float64');  
end
%for (j=1:my)
%	u00(j)=fread(finp,1,'float');
%	w00(j)=fread(finp,1,'float');
%end
dummy4=fread(finp, 1,'int');
wk1r=zeros(nspec,mgalx,mgalz);
% Here strarts the real wriging
%-------------------------------------
for (jj=1:nspec)
dummy5=fread(finp,1,'int');
dummy=fread(finp,[mgalx,mgalz],'float');
wk1r(jj,:,:)=dummy(:,:);
%for (k=1:mgalz)
%    for (i=1:mgalx)
%      			wk1r(jj,i,k)=fread(finp,1,'float');
%    		end
%	end

dummy6=fread(finp,1,'int');
end
%display([dummy1,dummy2])
fclose(finp)
Lx=2*pi/alp;
Lz=2*pi/bet;
x=(-Lx/2:Lx/(mgalx-1):Lx/2);
z=(-Lz/2:Lz/(mgalz-1):Lz/2);
%Ly=(max(y)-min(y));
posy=y(jspec);
