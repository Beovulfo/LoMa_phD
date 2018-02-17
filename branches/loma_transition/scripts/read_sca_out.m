%****************************************************************************%
% only for LISOuc3m			                                     %
% Antonio Almagro							     %
% September 2013								     %
% Script for LISO							     %
%****************************************************************************%
%clear all;
%clc
%===================READING FROM FILE==========================%
function [y,mode00,modenxnz]=read_sca_out(filename0,iframe,nx,nz)
%filename0='/data2/toni/mldiffSCAL01sc001.';G
%filename0='/data2/toni/mlturbSCAL02sc001.';
%filename0='/data2/toni/testcomp02RK2sc001.';
%filename0='/data2/toni/mlRMSCAL01sc001.';
%ibeg=1
%iend=10
my=513;mgalx=192;mgalz=192;
Ree=500;alpe=0.8621;bete=1.4368;
mye=my;

fmap=ones(mye,1);

ymax=10;
y=create_mesh1D(ymax,'manl',0.38,my); 

time=0;a0e=0;
u00=zeros(1,mye);
w00=zeros(1,mye);


mx=mgalz*2/3;mz=mgalz*2/3-1;
nplanes=mx/2;
%
  if iframe<10
     string2=strcat('00',num2str(iframe));
  elseif iframe<100
     string2=strcat('0',num2str(iframe));
  else
     string2=num2str(iframe);
  end
  filename=strcat(filename0,string2)

  finp=fopen(filename,'r','l');
%Read header
%----------------------Record 1----------------------------%
%%Buffer size    
dummy1=fread(finp, 1,'int');
%%record 1 data %
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
     if (i==(nx+1))
   	 mode=squeeze(reshape(wk1(i,:),2*my,mz));
     end
 end
 I=1:2:2*my
 modenxnz=mode(I,nz+1);
 mode00  =wk1(1,I,1);

 fclose(finp)

