dom = [-domlen,domlen];

A    = Zbars/Zs;
B    = (1-Zbars)/(1-Zs);
Heav = @(Z) 0.5*(1+tanh(bet*Z));

Zbar = @(Z) Zbars/Zs*Z.*Heav(Zs - Z) + ((1-Zbars)/(1-Zs)*(Z-Zs) + Zbars).*Heav(Z-Zs);
T    = @(Z,H) 1.0 + H + gamma/(Zs*(1-Zs))*Z.*Heav(Zs-Z) + gamma/(1-Zs).*Heav(Z-Zs);
%
%T    = @(Z,H) 1.0 + H + 0.5*K*Z + K/(2.0*bet)*(-log(cosh(bet*(Z-Zs)))+log(cosh(bet*Zs)));

%T    = @(Z,H) 1.0 + H + 0.5*K*Z + K/(2.0*bet)*(-log(cosh(bet*(Z-Zs)))) ...
%        + K/(2.0*bet)*(bet*Zs-log(2)); %This is by Higuera and achieve BC,losing a bit of Tad
%Zbar =@(Z) Zbars/Zs*Z.*Heav(Zs - Z) + ((1-Zbars)/(1-Zs)*(Z-Zs) + Zbars).*Heav(Z-Zs);
%Zbar = @(Z) 0.5*(A+B)*Z +0.5/bet*(B-A)*(log(cosh(bet*(Z-Zs)))-log(cosh(bet*Zs)));

%Zbar = @(Z) 0.5*(A+B)*Z +0.5/bet*(B-A)*(log(cosh(bet*(Z-Zs)))-(bet*Zs-log(2)));
L    = chebop(dom);       
L.op = @(x, u, Z, H) [diff(u).*(x/2)+              diff(T(Z,H).^(sigma-1.0).*diff(u));...
                      diff(Z).*(x/2)+ 1.0/(Lm*Pr).*diff(T(Z,H).^(sigma-1.0).*diff(Zbar(Z))); ...
 		              diff(H).*(x/2)+      1.0/Pr.*diff(T(Z,H).^(sigma-1.0).*diff(H))];  
L.lbc = @(u,Z,H) [u+1.0;     Z-1.0      ;H-HF  ];
L.rbc = @(u,Z,H) [u-1.0;      Z         ;H     ];
%L.init =[chebfun(@(x) erf(x),dom); chebfun(0,dom) ; chebfun(@(x) erf(x),dom)  ; chebfun(1,dom)];
if (fguess==1)
   load(fstartmat,'u0','Z0','H0');
   %load('hig_s05_q40_b1E4.mat','u0','Z0','H0');
   L.init =[u0;      Z0    ;H0];
end
%u0 = chebfun(@(x) tanh(x))
%Z0 = chebfun(@(x) 1+0.5*tanh(x))
%H0 = chebfun(@(x) 0)
%L.init =[u0 ;      Z0    ;H0; ; ];
zerov = [0;0;0];
sol = L\zerov;
u0=sol(1);Z0=sol(2);H0=sol(3);
T0=T(Z0,H0);
save(fnameresult);
%save('hig_q40_s05_S9_Pr1_b20.mat');


%x = chebfun('x',[-domlen,domlen],'splitting',flag,'splitdegree',512,'splitMaxLength',8000);%do we need to say this?
%u = chebfun([-domlen,domlen],'splitting',flag,'splitdegree',512,'splitMaxLength',8000)
%Z = chebfun([-domlen,domlen],'splitting',flag,'splitdegree',512,'splitMaxLength',8000)
%L.op = @(x, u, Z, H) [diff(u).*(x/2)+              diff(TP2(Z,H,Zs,gamma,bet,K)^(sigma-1.0)*diff(u));...
 %		      diff(Z).*(x/2)+(1.0/(Lm*Pr))*diff(TP2(Z,H,Zs,gamma,bet,K)^(sigma-1.0)*diff(ZB2(Z,Zs,Zbars,bet2))); ...
% 		      diff(H).*(x/2)+     (1.0/Pr)*diff(TP2(Z,H,Zs,gamma,bet,K)^(sigma-1.0)*diff(H))];   
