close all; clear all;
format long
%fname = 'S10LF03';
%fname ='S01LF20g6';
%fname ='S01LF03g6';
fname ='S01LF10g6';
stab = 0
dw = 1.0
%STABILITY
%load(['../../stability/stab_', fname, '.mat'])
load(['SSphys2_', fname,'.mat'])
Zbase = Zs;
Hbase = Hs;
if (stab==1)
    load(['../../stability/stab2_', fname, '.mat'])
    [~, pos1] = max(wi);
    alp = kv(pos1)
    dyv=diff(y);
    dy = dyv(2);
    Ny = length(y);
    load_FD
    vp = V(:,pos1);
    up = D1*vp/(-1i*alp);
%Eigenvalue
    c = (wr(pos1)+1i*wi(pos1))./alp
    drhos = D1*rhos; dZs = D1*Zs;dHs = D1*Hs;   

   % figure(45);
   % plot(y,rhos,y,drhos)
   % ylabel('\rho_s,\rho_s')

    omez = 1i*alp*vp+D2*vp/(1i*alp);

%up = real(up); 
%vp = real(vp);
%A2 = 0.05; A1 = 0.05;
    fact = 1.0/max(abs(up)); %Factor for scaling max(abs(u))=1
    up = up*fact; vp = vp*fact;
    upr = real(up); vpr = real(vp);
    rhop = -1./(1i*alp*(U-c)).*vp.*drhos;
    zp = -1./(1i*alp*(U-c)).*vp.*dZs;
    hp = -1./(1i*alp*(U-c)).*vp.*dHs;

   % figure(46);
    %plot(y,real(rhos)+real(rhop))
    %ylabel('\rho_s +\rho_p')
    rhopr = real(rhop);
    zpr = real(zp); hpr=real(hp);
%A10 = (2.0*sum((up.*conj(up))+(vp.*conj(vp))))^0.5;
%A100 = 0.1;
    ysta = y;
end

%BASE FLOW - MEAN PROFILES
load(['SSphys2_', fname,'.mat'])
Zbase = Zs;
Hbase = Hs;
eval(fname); Zs = 1./(1.0+S);
bet =500;
Heav = @(Z) 0.5*(1+tanh(bet*Z));
%T    = @(Z,H) 1.0 + H + gamma/(Zs*(1-Zs))*Z.*Heav(Zs-Z) + gamma/(1-Zs).*Heav(Z-Zs);
K =gamma./(Zs*(1-Zs));
T    = @(Z,H) 1.0 + H + 0.5.*K.*Z + K/(2.0*bet)*(-log(cosh(bet.*(Z-Zs)))+log(cosh(bet.*Zs)));

rhos  = 1.0./T(Zbase,Hbase);
%figure(32); plot(y,1./rhos,'.');xlim([-50,50]);
%figure(34); plot(y,Hbase,y,Zbase); xlim([-50,50]);
%figure(36);plot(y,U,y,Vs);xlim([-50,50]);




%Create new mesh
%choose my
my = 301;Ly= 10;
alpha = 1.0;beta=0.6; gg=0.1;
%alpha = 2.5;beta=0.6; gg=0.07;
create_mesh; %ynew
min(diff(y))
Unew = interp1(y,U,ynew,'linear','extrap')';
Vnew = interp1(y,Vs,ynew,'linear','extrap')';
Hnew = interp1(y,Hbase,ynew,'linear','extrap')';
Znew = interp1(y,Zbase,ynew,'linear','extrap')';
rhonew = interp1(y,rhos,ynew,'linear','extrap')';

%figure;plot(ynew/dw,Unew,'.')
dlmwrite(['mean_small_', fname, '.txt'],[ynew' Unew Vnew Hnew Znew rhonew],'delimiter','\t','precision',14)
%dlmwrite(['mean_', fname, '.txt'],[ynew' Unew Vnew Hnew Znew rhonew],'delimiter','\t','precision',14)

if (stab==1)
%keyboard
upnew = interp1(ysta*dw,upr,ynew,'linear',0.0)';
vpnew = interp1(ysta*dw,vpr,ynew,'linear',0.0)';
rhopnew = interp1(ysta*dw,rhopr,ynew,'linear',0.0)';
zpnew = interp1(ysta*dw,zpr,ynew,'linear',0.0)';
hpnew = interp1(ysta*dw,hpr,ynew,'linear',0.0)';
%figure(44)
%plot(ynew,rhopnew,'g-');hold on; plot(ynew,hpnew,ynew,zpnew)
%figure(33)
%plot(ynew,upnew,'.',ynew,vpnew,'.')
%xlim([-30,30])
dlmwrite(['pert_', fname, '.txt'],[ynew' upnew vpnew rhopnew hpnew zpnew],'delimiter','\t','precision',14)
end

%exit;
