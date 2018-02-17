%------------------------------------%
% M file ready to read velocities and plot them
% at the same plots for comparition 
%file='/data2/toni/inc03cte_01_';
%file='/data2/toni/inc03nocte_01_';
%file='/data2/toni/inc02cte_01_';
%file='/data2/toni/inc02nocte_01_';
%file='/data2/toni/inc03ctenewdy_01_';
file='/data2/toni/inc03noctenewdy_01_';
mode='v11';
my=257;
%---------------------------------%
figure(1);hold on;%for velocity
figure(2);hold on;%for dv/dy
figure(3);hold on;%for second derivatives
figure(4);hold on;%for third derivatives

% VELOCITY
filetail='000v11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(1);plot(y,v,'r',y,vimag,'b')
title('Velocity, real part (red),imaginary part(blue)')

%FIRST DERIVATIVE
filetail='00dv11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(2);plot(y,v,'r',y,vimag,'b')
title('d(v)/dy, real part (red),imaginary part(blue)')

%SECOND DERIVATIVES
filetail='dvdv11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(3);plot(y,v,'r*',y,vimag,'b*')
title('second derivatives, real part (red),imaginary part(blue)')

vaux=v;
vimagaux=vimag;

filetail='0d2v11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(3);plot(y,v,'r',y,vimag,'b')

%Estimating correlation
correlation2=corrcoef(v(:),vaux(:))
correlation2imag=corrcoef(vimag(:),vimagaux(:))

%Third derivatives
filetail='0d3v11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(4);plot(y,v,'r*',y,vimag,'b*')
title('Third derivatives, real part (red),imaginary part(blue)')

vaux=v;
vimagaux=vimag;

filetail='d2dv11.kk';
filename=[file,filetail];
[y,fmap,v,vimag]=leekk(filename,my);
figure(4);plot(y,v,'r',y,vimag,'b')

%Estimating correlation
correlation3=corrcoef(v(:),vaux(:))
correlation3imag=corrcoef(vimag(:),vimagaux(:))

correlation2
correlation2imag

