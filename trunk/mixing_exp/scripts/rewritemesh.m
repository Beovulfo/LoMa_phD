%****************************************************************************%
% Matlab script to read re-start file header and write again after changes   % 
% only for LISOuc3m			                                     %
% Antonio Almagro							     %
% September 2013								     %
% Script for LISO							     %
%****************************************************************************%
clear all;
close all;
clc
%===================READING FROM FILE==========================%
%
%
filename='../data/finpmlpantano09.dat'
finp=fopen(filename,'r','l');
rfactor=1;
%new mesh
my2=513*2;
ymax=172/2; %half heigth of the box
%y2=create_mesh1D(ymax,'mixl', 0.3,my2); 
y2=create_mesh1D(ymax,'line', 0.3,my2); 
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
   fmap2=zeros(my2,1);
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
nplanes=mx/2;
%dummy3=fread(finp,1,'int')
%keyboard
ntotr=4*my*mz;

% for (j=1:2*my)
%   for (k=1:(2*mz-1))
%                vor(j , k) = fread(finp,1,'float64');
%                vor(j , k) = fread(finp,1,'float64');
%%----------------------------------------------------%
%                phi(j , k) = fread(finp,1,'float64');
%                phi(j , k) = fread(finp,1,'float64');
%   end
% end
for (i=1:nplanes)
dummy3=fread(finp,1,'int')
 wk1(i,:)=fread(finp,ntotr,'float');
dummy4=fread(finp,1,'int')
end
dummy5=fread(finp,1,'int')


fclose(finp)


%===================WRITING TO FILE==========================%
%CAUTION DOING TIME=0
time=0.
%save record 1 values
finp=fopen('../data/finpmlpantano09y2.dat','w','l'); %open/create file tab_var.dat
my=my2
%Record 1 - counting bytes
rec1_real4=5+2*my;%time is real 4
%Ree,alpe,bete,a0e + y(j),fmap(j) (1:my)+wk1(j)(1:2my)
rec1_real8=my*2;
%integer mxe mye mze
rec1_integer=3;
%total record 1 size
rec1_size=4*rec1_real4+8*rec1_real8+4*rec1_integer;
%----------------------Record 1----------------------------%
%Buffer size    
fwrite(finp, rec1_size,'uint');
%record 1 data 
    fwrite(finp,time,'float');
    fwrite(finp,Re,'float');
    fwrite(finp,alp,'float');
    fwrite(finp,bet,'float');
    fwrite(finp,a0,'float');
    fwrite(finp, mx,'int');
    fwrite(finp, my,'int');
    fwrite(finp, mz,'int');
%write y(j) and fmap(j)
for (j=1:my)
    fwrite(finp,y2(j),'double');
    fwrite(finp,fmap2(j),'double');
end
%wk1==>u00, w00, odd/even
%interpolate u00 and clean w00
u002=interp1(y,u00,y2,'spline','extrap');
for (j=1:my)
	fwrite(finp,u002(j),'float');
	fwrite(finp,0*u002(j),'float');
end
%write buffer size again
fwrite(finp, rec1_size,'uint');
%------end of record 1---------------------------------------%

%----------RECORD 2 ---FOR EACH PLANES ALL VALUES------------%
%ntotr=4*my*mz %FOR EACH PLANE
rec2_size=ntotr*4; %SAME SIZE FOR EACH PLANES
indreal=(2*(1:my)-1);
for i=1:nplanes
    %reshape for each plane i
    %wk2=squeeze(reshape(wk1(i,:),2,2*my,mz));
    for k=1:mz
%scale and save vor 
      wk2(1,indreal  ,k)=interp1(y,wk1(1,indreal  ,k),y2,'linear',0);
      wk2(1,indreal+1,k)=interp1(y,wk1(1,indreal+1,k),y2,'linear',0);
%write PHI
      wk2(2,indreal  ,k)=interp1(y,wk1(2,indreal  ,k),y2,'linear',0);
      wk2(2,indreal+1,k)=interp1(y,wk1(2,indreal+1,k),y2,'linear',0);
    end
%----------------------------------------------------!
%Writing
    fwrite(finp,rec2_size,'int');
    fwrite(finp,wk2(:),'float');
    fwrite(finp,rec2_size,'int');
end
	
    fclose(finp);
%
