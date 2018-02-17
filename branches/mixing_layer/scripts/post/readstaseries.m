clear all
close all

%File 01 limits and name
string0='mlRM64';
Re=500;
nu=1/Re;
ibeg=1;
iend =20;

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
    [y,time,temp]=lee_sta(filename);
    if ind == 1;
        my =size(temp,1);
        nvar = size(temp,2);
        data=zeros(my,nvar,iend-ibeg);
    end
    data(:,:,ind) =temp;
    timev(ind)=time;
    clear temp
end

stringcf1=strcat(string0,'_01.cf');
[t1,ener1,dm1]=read_cf(stringcf1);
dm1=dm1(1:iend-ibeg+2);
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


for i=1:(iend-ibeg+1)
    %epsilon1(i)=nu*sum(data(:,end,i).*hy');
    epsilon1(i)=nu*sum(data(:,end,i).*hy')*dm1(i+1);
    %plot(y/dm1(i),fact*data(:,end,i)*dm1(i))
    
end

tau1=timev*2/dm1(1);
%ep1=zeros(my,(iend-ibeg+1));
ep1=squeeze(data(:,end,1:(iend-ibeg+1)));
%the velocity as it is with mean and pert
uu1=squeeze(data(:,1,1:(iend-ibeg+1)));
urms1=squeeze(data(:,4,1:(iend-ibeg+1)));
%save(string1,'y','timev','data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-------------FILE 02--------------------------------------------%
%
% clear timev data nvar time y
% ibeg=21;
% iend =39;
% 
% 
% string1=strcat(string0,'_02_');
% 
% 
% for ii=ibeg:iend-1
%     ind=ii-ibeg+1;
%     
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end
%     filename=strcat(string1,string2,'.sta');
%     [y,time,temp]=lee_sta(filename);
%     if ind == 1;
%         my =size(temp,1);
%         nvar = size(temp,2);
%         data=zeros(my,nvar,iend-ibeg);
%     end
%     data(:,:,ind) =temp;
%     timev(ind)=time;
%     clear temp
% end
% %
% stringcf2=strcat(string0,'_02.cf');
% [t2,ener2,dm2]=read_cf(stringcf2);
% %neglect first row because is repeated in .cf file
% dm2=dm2(2:iend-ibeg+1);
% epsilon2=zeros(1,iend-ibeg);
% 
% % figure
% % hold on
% % xlim([-5 5])
% for i=1:(iend-ibeg)
% %    epsilon2(i)=nu*sum(data(:,end,i).*hy');
%     epsilon2(i)=nu*sum(data(:,end,i).*hy')*dm2(i);
% end
% 
% tau2=timev*2/dm1(1);
% ep2=squeeze(data(:,end,1:(iend-ibeg)));
% uu2=squeeze(data(:,1,1:(iend-ibeg)));
% urms2=squeeze(data(:,4,1:(iend-ibeg)));
% 
%  
%--------------JOIN FILES-------------------------------%
epsilon=[epsilon1]
tau    =[tau1]
dm     =[dm1];
% figure
% plot(tau,500*epsilon/8)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%..........:FILE 03.................
%save(string1,'y','timev','data')
% clear timev data nvar time y
% ibeg=1;
% iend =10;
% 
% 
% string1=strcat(string0,'_03_');
% 
% 
% for ii=ibeg:iend-1
%     ind=ii-ibeg+1;
%     
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end
%     filename=strcat(string1,string2,'.sta');
%     [y,time,temp]=lee_sta(filename);
%     if ind == 1;
%         my =size(temp,1);
%         nvar = size(temp,2);
%         data=zeros(my,nvar,iend-ibeg);
%     end
%     data(:,:,ind) =temp;
%     timev(ind)=time;
%     clear temp
% end
% %
% stringcf3=strcat(string0,'_03.cf');
% [t3,ener3,dm3]=read_cf(stringcf3);
% %neglect first row because is repeated in .cf file
% dm3=dm3(2:iend-ibeg+1);
% epsilon3=zeros(1,iend-ibeg);
% %hold on
% %xlim([-5 5])
% for i=1:(iend-ibeg)
%    % epsilon3(i)=nu*sum(data(:,end,i).*hy');
%   epsilon3(i)=nu*sum(data(:,end,i).*hy')*dm3(i);
%  
%     
% end
% 
% tau3=timev*2/dm1(1);
% ep3=squeeze(data(:,end,1:(iend-ibeg)));
% uu3=squeeze(data(:,1,1:(iend-ibeg)));
% urms3=squeeze(data(:,4,1:(iend-ibeg)));
%  %plot(y/dm2(i+1),nu/8*data(:,end,i)*dm2(i+1))

%--------------JOIN FILES-------------------------------%
epsilon=[epsilon]
tau    =[tau]
dm     =[dm];
ep=[ep1];
uu=[uu1];
urms=[urms1];
%-----------------PLOTS----------------------------------%


figure
subplot(2,2,1)
plot(tau,500*epsilon/8,'k--')
hold on
plot(tau,dm(2:end)/dm1(1),'b-')


subplot(2,2,2)
hold on;
itau=[1 5 10 15];
for n=1:length(itau)
    plot(y/dm(itau(n)),uu(:,itau(n))/2)
end

subplot(2,2,3)
hold on
for n=1:length(itau)
    plot(y/dm(itau(n)),(urms(:,itau(n))/2).^2)
end

subplot(2,2,4)
hold on
Dx=2*pi/(0.2247*(768-1));
Dz=2*pi/(0.8977*(288-1));
plot(tau,(nu^3./(epsilon.*dm(2:end)')).^0.25)
plot(tau,min(diff(y)),'b--')
plot(tau,Dx,'r-',tau,Dz,'g-')


