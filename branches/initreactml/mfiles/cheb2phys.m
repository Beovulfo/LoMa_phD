fresults = [fpath fname '.mat'];
%Define vector y and Ny
R=15;
Ny =3001;
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
xi = (1.0/Re*t)^0.5.*eta;
%Z1 = interp1(eta,Zeta,xi,'pchip');
%T1 = interp1(eta,Teta,xi,'pchip');
%U1 = interp1(eta,Ueta,xi,'pchip');
posZs=find(Zeta>Zs);xif=xi(posZs(end));
xif
%eta = eta - etaf;
yf=trapz(xi(posZs),Teta(posZs))
yv = cumtrapz(xi,Teta);%WARNING NOT SURE!!
Ly=trapz(xi,Teta);
%Find pos Z=Zs
yv = yv - yf; %Now y=0 is the position of the flame
method ='linear';
U = interp1(yv,Ueta,y,method,'extrap')';
Ts = interp1(yv,Teta,y,method,'extrap')';
Hs = interp1(yv,Heta,y,method,'extrap')';
Zs = interp1(yv,Zeta,y,method,'extrap')';
Vs = (1.0/(Re*t)).^0.5.*interp1(yv,Veta,y,method,'extrap')';
%Vs = interp1(yv,Teta,y,method,'extrap')';
rhos = 1.0./Ts;

load_FD;
dUdy = D1*U;

save(['SSphys_' fname '.mat'],'y','U','rhos', ...
                     'yf','Vs','Hs','Zs', 'dUdy','LF','S');



