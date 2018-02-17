fresults = [fpath fname '.mat'];
%Define vector y and Ny
R=15;
Ny =3501;
%Ny = 2000;
dy0 = 2*R/(Ny-1)
dy = dy0
y = -R:dy:R;
Ny = length(y)

base=load(fresults,'u0','Z0','H0','T0','Zs','sigma','gamma','dom','LF','S','domlen');
Zs = base.Zs;
LF = base.LF; S=base.S;
sigma = base.sigma;
gamma = base.gamma;
u0=base.u0;
T0=base.T0;
Z0=base.Z0;
H0=base.H0;
dom=base.domlen
chebfunpref.setDefaults('splitting',true,'maxLength',9000)
tprima = diff(T0);
%demanding igual entraintmen Vinf = -Vsup
%Vinf = sum(tprima.*x.*0.5)*0.5
N = chebop([-dom,dom]); N.op=@(x,v) [diff(v)+tprima.*x.*0.5]; 
N.lbc=@(v) [v];
solV = N\[0];
V0 = solV;
intV0 = V0(end);
Vinf = -intV0*0.5;

%N = chebop([-dom,dom]); N.op=@(x,v) [diff(v)+tprima.*x.*0.5]; 
N.lbc=@(v) [v-Vinf];
solV = N\[0];
V0 = solV;

%dom = base.dom;
%set of dw=1 gives time:
dw=1.0;
%define eta
Re=1000;
%dom=15;
%Np=2000;eta = [linspace(-dom,-10,500) linspace(-10,10,10000) linspace(10,dom,500)];
Np = 8001; eta=linspace(-dom,dom,Np);
%du=diff(u0);
Ueta = evalcheb(eta,u0);
%dudeta = evalcheb(eta(2:end-1),du);
dudeta = diff(Ueta)./diff(eta);
Teta = evalcheb(eta,T0);
Zeta = evalcheb(eta,Z0);
Veta = evalcheb(eta,V0);
%Veta = evalcheb(eta,V0);
Heta = evalcheb(eta,H0);
dumax=max(1.0./Teta(2:end).*dudeta)
t=(dw*dumax*0.5)^2*Re
delta = (1.0/Re*t)^0.5;
xi = delta.*eta;
%keyboard
eta0 = find(eta>=0.0); eta0=eta0(1)
%Z1 = interp1(eta,Zeta,xi,'pchip');
%T1 = interp1(eta,Teta,xi,'pchip');
%U1 = interp1(eta,Ueta,xi,'pchip');
posZs=find(Zeta>Zs);xif=xi(posZs(end));
xif
etaf = posZs(end);
V0 = Veta(eta0);

ymin = 2*V0*delta - trapz(eta(1:eta0), Teta(1:eta0))*delta
ymax = 2*V0*delta + trapz(eta(eta0:end), Teta(eta0:end))*delta
yf=(2*V0 + trapz(eta(eta0:etaf), Teta(eta0:etaf)))*delta
ynew = 2*V0*delta + cumtrapz(eta,Teta)*delta-trapz(eta(1:eta0),Teta(1:eta0))*delta;
%yv = cumtrapz(xi,Teta);%WARNING NOT SURE!!
Ly=ymax-ymin;
%Find pos Z=Zs
%ynew = ynew-yf;
%yv = yv - yf; %Now y=0 is the position of the flame
method ='linear';
U = interp1(ynew,Ueta,y,method,'extrap')';
Ts = interp1(ynew,Teta,y,method,'extrap')';
Hs = interp1(ynew,Heta,y,method,'extrap')';
Zs = interp1(ynew,Zeta,y,method,'extrap')';

Vs = (1.0/(Re*t)).^0.5.*interp1(ynew,Veta,y,method,'extrap')';
%Vs = interp1(yv,Teta,y,method,'extrap')';
rhos = 1.0./Ts;

load_FD;
dUdy = D1*U;
%figure(1); plot(y,1./rhos,'b.-');hold on; plot(y,U,'r.-')

save(['SSphys2_' fname '.mat'],'y','U','rhos', ...
                     'yf','Vs','Hs','Zs', 'dUdy','LF','S');



