clear all
close all
clc 

my=129;
y(1:my)=-1.+((1:my)-1)*2./(my-1);
y=1.02*tanh(y*atanh(1/1.02))
y(1)=-1;
y(my)=1;
k=10;
c=atan(k); 
%tangente
%u00=1/c*atan(k*y);
%k2=5;
%c2=atan(k2);
%y polynomial
u00=5/2*(y.^3-3/5*y.^5);

plot(y,u00)
%vecty2=y.*y
%u00p=k/c./(1+k^2.*(vecty2));
u00p=5/2*(3.*y.^2-3*y.^4);
u00pp=15/2*(2*y-4*y.^3);

hold on
plot(y,u00p,'r-')
%u00pp=k/c*(-2*k^2.*y)./(1+k^2.*vecty2);
%u00prevpp=k2/c2*(-2*k2^2.*y)./(1+k2^2.*vecty2);
plot(y,u00pp,'g-')
%RK constants
Re=1;
ire=1/Re;
Deltat=0.0001;
rkstep=2;

gamma=[8/15,5/12,3/4];
alpha=[29/96,-3/40,1/6,29/96];
beta=[37/160,5/24,1/6];
ibeta=[160/37,24/5,6];
xi=[-17/60,-5/12,0];

rkn1    = Re*ibeta(rkstep)/Deltat
dalre   = Deltat*alpha(rkstep)*ire
dalbe   = 1+alpha(rkstep+1)*ibeta(rkstep)
dtgamma = Deltat*gamma(rkstep)
dtri    = Deltat*alpha(rkstep+1)*ire
dtxi    = Deltat*xi(rkstep)


RHS=u00pp-rkn1*u00 
rf0u=100*sin(pi*y)
rhs=-RHS/rkn1-dtgamma*rf0u
