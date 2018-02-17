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
function [y,time,temp,Re,alp,bet]=lee_sta(filename)
%finp=fopen('/data/toni/mixing/modesx2.055','r','l'); %open/create file tab_var.dat
%filename=input('filename:\n')
%
finp=fopen(filename,'r','l');
%frewind(finp);
%----------------------Record 1----------------------------%
%Buffer size    


dummy1=fread(finp, 1,'int');
%record 1 data 
 nacum= fread(finp,1,'int');
 Wx0a1= fread(finp,1,'float64');% Wx0a/nacum
 Wz0a1= fread(finp,1,'float64');% Wz0a/nacum

dummy2=fread(finp, 1,'int');
%end of record 1

%RECORD 2
dummy3=fread(finp, 1,'int');
   my= fread(finp,1,'int');
   time=fread(finp,1,'float');
   Re=fread(finp,1,'float');
   alp= fread(finp,1,'float');
   bet= fread(finp,1,'float');
   a0= fread(finp,1,'float');
dummy4=fread(finp, 1,'int');


   y=zeros(my,1);
   fmap=y;
   um=y;
   vm=um;
   wm=um;
   w1m=um;
   w2m=um;
   w3m=um;

   wk=zeros(my,13);

%RECORD3
dummy5=fread(finp, 1,'int');
%wk1==>u00, w00, odd/even
for (j=1:my)
     y(j) = fread(finp,1,'float64');
     fmap(j) = fread(finp,1,'float64');  
end
dummy6=fread(finp, 1,'int');


dummy7=fread(finp, 1,'int');
for (j=1:my)
     um(j) = fread(finp,1,'float64');
     vm(j) = fread(finp,1,'float64');
     wm(j) = fread(finp,1,'float64');
     wk(j,1)=fread(finp,1,'float64');
     wk(j,2)=fread(finp,1,'float64');
     wk(j,3)=fread(finp,1,'float64');
     
     w1m(j)=fread(finp,1,'float64');
     w2m(j)=fread(finp,1,'float64');
     w3m(j)=fread(finp,1,'float64');
     wk(j,4)=fread(finp,1,'float64');
     wk(j,5)=fread(finp,1,'float64');
     wk(j,6)=fread(finp,1,'float64');
     wk(j,7)=fread(finp,1,'float64');
     wk(j,8)=fread(finp,1,'float64');
     wk(j,9)=fread(finp,1,'float64');
     wk(j,10)=fread(finp,1,'float64');
     wk(j,11)=fread(finp,1,'float64');
     wk(j,12)=fread(finp,1,'float64');
     wk(j,13)=fread(finp,1,'float64');

end
dummy8=fread(finp, 1,'int');
%display([dummy1,dummy2])
fclose(finp)

U  = um/nacum;
V  = vm/nacum;
W  = wm/nacum;
uu = wk(:,1)/nacum;
vv = wk(:,2)/nacum;
ww = wk(:,3)/nacum;
uv = wk(:,7)/nacum;

wxm= w1m/nacum;
wym= w2m/nacum;
wzm= w3m/nacum;

wx = wk(:,4)/nacum;
wy = wk(:,5)/nacum;
wz = wk(:,6)/nacum;
%
urms = sqrt(uu-U.^2);
vrms = sqrt(vv-V.^2);
wrms = sqrt(ww-W.^2);

wxrms = sqrt(wx-wxm.^2);
wyrms = sqrt(wy-wym.^2);
wzrms = sqrt(wz-wzm.^2);

%
uv = uv - U.*V;
%

ep=wk(:,10)/nacum;
%
temp=[U,V,W,urms,vrms,wrms,uv,wxm,wym,wzm,wxrms,wyrms,wzrms,ep];

%Interesting to calculate all in wall units to draw
% U+=f(y+)!
%yplus=y/
