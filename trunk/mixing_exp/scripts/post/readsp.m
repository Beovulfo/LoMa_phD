% [yy,kx,kz,SPA]=readsp(fil,8)
function [ysp,kx,kz,SPA,Ener1,Ener2,Ener3,Ener1z,Ener2z,Ener3z]=readsp(fil,nvar)
% [yy,kx,kz,SPA]=readsp(fil,nvar)
fid=fopen(fil,'r','l');
dummy1=fread(fid,1,'int'); 
  t=fread(fid,1,'float'); 
  Re=fread(fid,1,'float'); 
  alp=fread(fid,1,'float'); 
  bet=fread(fid,1,'float');
  nx = fread(fid,1,'int')./2; 
  ny = fread(fid,1,'int');
  nz = (fread(fid,1,'int')+1)./2; 
  np = fread(fid,1,'int');
  nacum = fread(fid,1,'int'); 
  %
  kx=[0:nx-1]*alp; 
  kz=[0:nz-1]*bet; 
dummy2=fread(fid,2,'int'); 
  jsp = fread(fid,np,'int'); 
  kk = fread(fid,[2 ny],'double'); 
  y = kk(1,:); ysp = y(jsp); 
%  yy = sin((2*(yy-1)./(ny-1)-1)*pi./2) + 1; 
fread(fid,1,'int');
%
SPA=zeros(nx,nz,np,nvar);
for i=1:nx 
    i0 = fread(fid,1,'int');
    for j=1:nvar
       SPA(i,:,:,j)=fread(fid,[nz np],'float'); 
    end
    i1 =fread(fid,1,'int'); 
    if i0~=i1, pause; end
end
SPA = squeeze(SPA)./nacum;
fread(fid,1,'int'),
fclose(fid);


%nxvector=(0:mx/2);
%nzvector=[(0:(mz-1)/2) (-(mz-1)/2:1:-1)];
%kxvector=nxvector*alp;
%kzvector=nzvector*bet;
cont=0;
Ener1=zeros(length(kx),1);
Ener2=zeros(length(kx),1);
Ener3=zeros(length(kx),1);

Ener1z=zeros(length(kz),1);
Ener2z=zeros(length(kz),1);
Ener3z=zeros(length(kz),1);

% Ener1(1)=0; %mode zero
for(i=1:length(kx))
    Ener1(i)=squeeze(sum(SPA(i,:,end,1),2))./kx(2);
    Ener2(i)=squeeze(sum(SPA(i,:,end,2),2))./kx(2);
    Ener3(i)=squeeze(sum(SPA(i,:,end,3),2))./kx(2);
end
for(i=1:length(kz))
    Ener1z(i)=squeeze(sum(SPA(:,i,end,1),1))./kz(2);
    Ener2z(i)=squeeze(sum(SPA(:,i,end,2),1))./kz(2);
    Ener3z(i)=squeeze(sum(SPA(:,i,end,3),1))./kz(2);
end



%PLOT Energy SPECTRUM
% %---------------------
% [ysp,kx,kz,SPA,Ener1,Ener2,Ener3,Ener1z,Ener2z,Ener3z]=readsp(fil,8); 
% itime=1;figure;loglog(kx(2:end).*dm(itime),0.5*(Ener3(2:end)+Ener2(2:end)+Ener1(2:end))./dm(itime),'b-');
% hold on; loglog(kz(2:end).*dm(itime),0.5*(Ener3z(2:end)+Ener2z(2:end)+Ener1z(2:end))./dm(itime),'g--');
% loglog(kz(2:end).*dm(itime),0.002*(kz(2:end).*dm(itime)).^(-5/3),'r-');axis tight;
% xlabel('$k_x \delta_m$,$k_z \delta_m$','Interpreter','latex','FontSize',13),ylabel('$E(k)/(\Delta U^2 \delta_m)$','Interpreter','latex','FontSize',14)
% ylim([10^(-6) 0.1]);xlim([0.1 8])
% legend('Espectro kx','Espectro kz','-5/3')


% itime=20;figure;loglog(kx(2:end).*dm(itime),0.5*(Ener3(2:end)+Ener2(2:end)+Ener1(2:end))./dm(itime),'b-');
% hold on; loglog(kz(2:end).*dm(itime),0.5*(Ener3z(2:end)+Ener2z(2:end)+Ener1z(2:end))./dm(itime),'g--');
% loglog(kz(2:end).*dm(itime),0.02*(kz(2:end).*dm(itime)).^(-5/3),'r-');axis tight;
% xlabel('$k_x \delta_m$,$k_z \delta_m$','Interpreter','latex','FontSize',13),ylabel('$E(k)/(\Delta U^2 \delta_m)$','Interpreter','latex','FontSize',14)
% ylim([10^(-6) 0.1]);xlim([0.1 8])
% legend('Espectro kx','Espectro kz','-5/3')


