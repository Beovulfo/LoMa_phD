%[y,tau,data,dm,epsilon]=readstaseries2(string0,ibeg,iend,my,mgalx,mgalz,naverages)
%clear all
%close all
function [y,tau,data,dm,epsilon]=readstaseries2(string0,ibeg,iend,my,mgalx,mgalz,naverages)
%File 01 limits and name
%string0='ml05turb';
%string0='ml04turb';
%string0='ml09turb';
% string0='mlpantano06';
% %string0='mlnozeromodes';
% %string0='ml01loroll';
% %Re=1370;
% %Re=675;
% %Re=1370;
% ibeg=1;
% iend=43;
% my=513;
% mgalx=768;
% mgalz=288;
imid=ceil(my/2)

%itau=[50 51 52  53 54 55];
%itau=[38 39 40];
%itau=[38 39 40];
%itau=[22:28];
%itau=[30 31 32];
%itau=[45 46 47];
itau=iend-naverages+1:iend;
%itau=52:57;
%itau=44:56;
%itau=[3]
string1=strcat(string0,'_01_');

for ii=1:iend
    ind=ii;
    if ii<10
        string2=strcat('00',num2str(ii));
    elseif ii<100
        string2=strcat('0',num2str(ii));       
    else 
        string2=num2str(ii);
    end
    filename=strcat(string1,string2,'.sta');
    %check if 01 doesn't exist to increase to 02
    if exist(filename,'file')==0
        string1=strcat(string0,'_02_');
        filename=strcat(string1,string2,'.sta');
    end
    
    if exist(filename,'file')==0
        string1=strcat(string0,'_03_');
        filename=strcat(string1,string2,'.sta');
    end
    
    [y,time,temp,Re,alp,bet]=lee_sta(filename);
    if ind == 1;
        my =size(temp,1);
        nvar = size(temp,2);
        data=zeros(my,nvar,iend-ibeg);
    end
    data(:,:,ind) =temp;
    timev(ind)=time;
    clear temp
end

% stringcf1=strcat(string0,'_01.cf');
% fid = fopen(stringcf1);
% file = textscan(fid, '%s', 'delimiter', '\n', ...
%                 'whitespace', '');
% fclose(fid);
nu=1/Re;

%textscan returns a 1-by-1 cell array, file,
%that contains a 31-by-1 cell array:
% %Show # of lines in the file:
% lines1 = length(file{1});
% [t,ener,dm]=read_cf(stringcf1);
% if lines1<(iend-ibeg+2)
%     dm1=dm(1:lines1);
% else
%     dm1=dm(1:iend+1);
% end
% 
% if lines1<(iend-ibeg+2)
%     stringcf2=strcat(string0,'_02.cf');
%     fid = fopen(stringcf2);
%     file = textscan(fid, '%s', 'delimiter', '\n', ...
%                 'whitespace', '');
%     fclose(fid);
%     lines2 = length(file{1});
%     %si con 2 archivos no tengo suficiente
%     if (lines2+lines1-1<iend-ibeg+2)
%         %read cf3
%         %the num lines i need from cf2
%         last2=lines2;
%         last3=iend-ibeg+4-lines1-lines2;    
%         
%     else
%          last2=iend-ibeg+3-lines1; 
%         
%     end
%     %read cf 02
%         [t,ener,dm]=read_cf(stringcf2);
%         dm1=[dm1;dm(2:last2)];
%      if (lines2+lines1-1<iend-ibeg+2)  
%          
%         stringcf3=strcat(string0,'_03.cf');
%         [t,ener,dm]=read_cf(stringcf3);
%         dm1=[dm1;dm(2:last3)];
%         
%      end
% %               
% end

%Epsilon (pseudo)
%figure
%hold on
%xlim([-5 5])
%need reynolds for scale

epsilon1=zeros(1,iend-ibeg+1);
mye=length(y);
hy(1)=y(2)-y(1);
for (j=2:mye-1)
    hy(j)=0.5*(y(j+1)-y(j-1));
end
hy(mye)=y(mye)-y(mye-1);
Dy=min(hy(:));

dm1=sum(hy*(0.25-(squeeze(data(:,1,:))).^2),1);

tau1=timev/dm1(1);
%ep1=zeros(my,(iend-ibeg+1));
size(dm1)
ep1=nu*squeeze(data(:,end,:));%This has u^3 dimension
epsilon1=squeeze(ep1'*hy'); %this has u^3/l dimension
%be careful with this definitioon, it has length dimension
%epsilon1=sum(epsilon1,2)
%the velocity as it is with mean and pert
uu1=squeeze(data(:,1,1:(iend-ibeg+1)));
%urms1=squeeze(data(:,4,1:(iend-ibeg+1)));
%Trying with V
urms1=squeeze(data(:,4,1:(iend-ibeg+1)));
vrms1=squeeze(data(:,5,1:(iend-ibeg+1)));
wrms1=squeeze(data(:,6,1:(iend-ibeg+1)));
wm=squeeze(data(:,3,1:(iend-ibeg+1)));


%kturb=0.5*(squeeze(mean(urms1.^2,1))+squeeze(mean(vrms1.^2,1))+squeeze(mean(wrms1.^2,1)));
%TKE=mean(squeeze(data(:,end-1,1:(iend-ibeg+1))),1);
TKE=0.5*(urms1(imid,:).^2+vrms1(imid,:).^2+wrms1(imid,:).^2);
%save(string1,'y','timev','data')
%TEST
%epsilon1=ep1(imid,:)*Dy./dm1';
%--------------JOIN FILES-------------------------------%
epsilon=[epsilon1];% epsilon3]
tau    =[tau1];% tau3]
dm     =[dm1];%;dm3];
ep=[ep1]; %ep2 %ep3];
uu=[uu1]; %uu2 uu3];
urms=[urms1];% urms2 urms3];
%-----------------PLOTS----------------------------------%


figure
subplot(2,2,1)
plot(tau,500*epsilon,'k*')
hold on
plot(tau,dm./dm(1),'b-')
title('Momentum thickness and eps evolution')
xlabel('t*')
ylabel('500*eps,dm/dm0')
p=polyfit(tau(end-naverages+1:end),dm(end-naverages+1:end),1)


subplot(2,2,2)
hold on
for n=1:length(itau)
    plot(y/dm(itau(n)),urms(:,itau(n)).^2)
%     plot(y/dm(itau(n)),abs(wm(:,itau(n))))
  
end
title('Self-similarity of urms*^2 vs y*')
xlabel('y*')
ylabel('urms*^2')
xlim([-5 5])

subplot(2,2,3)
hold on;
for n=1:length(itau)
    plot(y/dm(itau(n)),uu(:,itau(n)))
end
title('Self-similarity of u mean* vs y*')
xlabel('y*')
ylabel('um*')
xlim([-5 5])


subplot(2,2,4)
hold on;
xlim([-5 5])
for n=1:length(itau)
    %plot(y/dm(itau(n)),ep(:,itau(n)).*dm(itau(n)))
     plot(y/dm(itau(n)),vrms1(:,itau(n)).^2)
end
title('vrms^2 vs y* self-similarity')

figure(200)
hold on
Dx=2*pi/(alp*(mgalx-1));
Dz=2*pi/(bet*(mgalz-1));
Dy=min(diff(y));
plot(tau,(nu^3./(ep(imid,:))).^0.25)
plot(tau,min(diff(y)),'k*')
plot(tau,Dx,'r*',tau,Dz,'gsquare')
title('Evolution of kolmogorov scale')
xlabel('t*')
ylabel('eta')

figure(666)
plot(y,hy'./(nu^3./(squeeze(ep(:,end)))).^0.25,'ko')
xlabel('$y$','Interpreter','latex','FontSize',14)
ylabel('$\frac{\Delta y}{\eta}(y)$','Interpreter','latex','FontSize',14)



%urmsy0=squeeze(urms(257,:));%urms in the middle
%lambda=(epsilon./(15*nu.*(urmsy0))).^(-0.5);

Relambda=2.*TKE.*sqrt(5./(nu.*ep(imid,:)));%NOT SURE
%Relambda=2.*TKE.*sqrt(5./(nu.*epsilon1'));%NOT SURE
figure
plot(tau,Relambda)
title('Reynolds lambda evolution');
xlabel('tau'); ylabel('Re-lambda');


Redm=Re.*dm;
figure
plot(tau(1:end-1),diff(dm)./diff(tau),'bo')
%plot(tau,Redm,'b*')
title('Slope of dm vs t')
%title('Reynolds number evolution based on momentum thickness');
xlabel('tau'); ylabel('Ddm/Dtau');

figure
plot(tau,TKE,'r*')
xlabel('tau')
ylabel('TKE');


%========================================================================%
%------------------------------------------------------------------------%
%PROCEDURE TO GET AVERAGE OF sqrt(Rij) variables
% close all

%naverages=6;
%istart=52;
istart=iend-naverages;
Dw=5;
figure; hold on;
xlabel('$y/\delta_w(\tau)$','Interpreter','latex','FontSize',13),ylabel('$\frac{\sqrt{R_{12}}}{\Delta u}$','Interpreter','latex','FontSize',14)
R12_1990=load('~/Documents/MASTER/MasterThesis/TFM/R12_1990.dat');
R12_1994=load('~/Documents/MASTER/MasterThesis/TFM/R12_1994.dat');
R12_1971=load('~/Documents/MASTER/MasterThesis/TFM/R12_1971.dat');
for i=1:naverages
    %Dwdm=Dy/max(diff(data(imid-1:imid+1,1,istart+i-1)));
    Dwdm=Dw*dm(istart+i);
    y2=y./Dwdm;
    psi=[-1:0.05:1];
    datascaled=squeeze(interp1(y2,sqrt(abs(data(:,7,istart+i-1))),psi,'spline'));
    %plot(psi,datascaled)
    datascaled2(:,i)=datascaled(:);
end
datamean=squeeze(mean(datascaled2,2));
R12peak=max(datamean(:))
plot(psi,datamean,'k-');
plot(R12_1990(:,1),R12_1990(:,2),'bo')
plot(R12_1994(:,1),R12_1994(:,2),'r--')
plot(R12_1971(:,1),R12_1971(:,2),'gsquare')



figure; hold on;
xlabel('$y/\delta_w(\tau)$','Interpreter','latex','FontSize',13),ylabel('$\frac{\sqrt{R_{11}}}{\Delta u}$','Interpreter','latex','FontSize',14)
R11_1990=load('~/Documents/MASTER/MasterThesis/TFM/R11_1990.dat');
R11_1994=load('~/Documents/MASTER/MasterThesis/TFM/R11_1994.dat');
R11_1971=load('~/Documents/MASTER/MasterThesis/TFM/R11_1971.dat');
for i=1:naverages
    %Dwdm=Dy/max(diff(data(imid-1:imid+1,1,istart+i-1)));
    Dwdm=Dw*dm(istart+i);
    y2=y./Dwdm;
    psi=[-1:0.05:1];
    datascaled=squeeze(interp1(y2,abs(data(:,4,istart+i-1)),psi,'spline'));
    %plot(psi,datascaled)
    datascaled2(:,i)=datascaled(:);
end
datamean=squeeze(mean(datascaled2,2));
R11peak=max(datamean(:))
plot(psi,datamean,'k-');
plot(R11_1990(:,1),R11_1990(:,2),'bo')
plot(R11_1994(:,1),R11_1994(:,2),'r--')
plot(R11_1971(:,1),R11_1971(:,2),'gsquare')



figure; hold on;
xlabel('$y/\delta_w(\tau)$','Interpreter','latex','FontSize',13),ylabel('$\frac{\sqrt{R_{22}}}{\Delta u}$','Interpreter','latex','FontSize',14)
R22_1990=load('~/Documents/MASTER/MasterThesis/TFM/R22_1990.dat');
R22_1994=load('~/Documents/MASTER/MasterThesis/TFM/R22_1994.dat');
R22_1971=load('~/Documents/MASTER/MasterThesis/TFM/R22_1971.dat');
for i=1:naverages
    %Dwdm=Dy/max(diff(data(imid-1:imid+1,1,istart+i-1)));
    Dwdm=Dw*dm(istart+i);
    y2=y./Dwdm;
    psi=[-1:0.05:1];
    datascaled=squeeze(interp1(y2,abs(data(:,5,istart+i-1)),psi,'spline'));
    %plot(psi,datascaled)
    datascaled2(:,i)=datascaled(:);
end
datamean=squeeze(mean(datascaled2,2));
R22peak=max(datamean(:))
plot(psi,datamean,'k-');
plot(R22_1990(:,1),R22_1990(:,2),'bo')
plot(R22_1994(:,1),R22_1994(:,2),'r--')
plot(R22_1971(:,1),R22_1971(:,2),'gsquare')

figure; hold on;
xlabel('$y/\delta_w(\tau)$','Interpreter','latex','FontSize',13),ylabel('$\frac{\sqrt{R_{33}}}{\Delta u}$','Interpreter','latex','FontSize',14)
R33_1990=load('~/Documents/MASTER/MasterThesis/TFM/R33_1990.dat');
R33_1994=load('~/Documents/MASTER/MasterThesis/TFM/R33_1994.dat');
R33_1971=load('~/Documents/MASTER/MasterThesis/TFM/R33_1971.dat');
for i=1:naverages
    %Dwdm=Dy/max(diff(data(imid-1:imid+1,1,istart+i-1)));
    Dwdm=Dw*dm(istart+i);
    y2=y./Dwdm;
    psi=[-1:0.05:1];
    datascaled=squeeze(interp1(y2,abs(data(:,6,istart+i-1)),psi,'spline'));
    %plot(psi,datascaled)
    datascaled2(:,i)=datascaled(:);
end
datamean=squeeze(mean(datascaled2,2));
R33peak=max(datamean(:))
plot(psi,datamean,'k-');
plot(R33_1990(:,1),R33_1990(:,2),'bo')
plot(R33_1994(:,1),R33_1994(:,2),'r--')
plot(R33_1971(:,1),R33_1971(:,2),'gsquare')

Re

figure; hold on;
xlabel('$y/\delta_w(\tau)$','Interpreter','latex','FontSize',13),ylabel('$\bar{u}_1$','Interpreter','latex','FontSize',14)
uself=load('~/Documents/MASTER/MasterThesis/TFM/uselfpantano.dat');
for i=1:naverages
    %Dwdm=Dy/max(diff(data(imid-1:imid+1,1,istart+i-1)))
    Dwdm=Dw*dm(istart+i);
    y2=y./Dwdm;
    psi=[-1:0.05:1];
    datascaled=squeeze(interp1(y2,data(:,1,istart+i-1),psi,'spline'));
    %plot(psi,datascaled)
    datascaled2(:,i)=datascaled(:);
end
datamean=squeeze(mean(datascaled2,2));
plot(uself(:,1),uself(:,2),'ro')
plot(psi,-datamean,'k-');

%Dwdm=Dy/max(diff(data(:,1,istart+i-1)));

figure(666)
plot(y,hy'./(nu^3./(squeeze(ep(:,end)))).^0.25,'k-')
xlabel('$y$','Interpreter','latex','FontSize',14)
ylabel('$\frac{\Delta}{\eta}(y)$','Interpreter','latex','FontSize',14)
hold on
plot(y,Dx./(nu^3./(squeeze(ep(:,end)))).^0.25,'r--')
plot(y,Dz'./(nu^3./(squeeze(ep(:,end)))).^0.25,'b-.')



