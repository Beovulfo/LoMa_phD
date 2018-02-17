%Define FD matrix and all discretizated terms
%Define Derivative matrix for U and rhos.
%for (i = 1:lenght(y))
D1  = tridiag(0,-0.5/dy,0.5/dy,Ny);
D1(1,1)    =-1.5/dy; D1(1,2)      =2/dy; D1(1,3)      =-0.5/dy;
D1(end,end)=1.5/dy;  D1(end,end-1)=-2.0/dy; D1(end,end-2)= 0.5/dy;

D2 =  tridiag(-2/dy^2,1/dy^2,1/dy^2,Ny);
D2(1,1)    =1/dy^2; D2(1,2)      =-2/dy^2;D2(1,3)      =1/dy^2;
D2(end,end)=1/dy^2; D2(end,end-1)=-2/dy^2;D2(end,end-2)=1/dy^2;

D1v = tridiag(0,-0.5/dy,0.5/dy,Ny);
D1v(1,1)    =-1.5/dy; D1v(1,2)      =2/dy; D1v(1,3)      =-0.5/dy;
D1v(end,end)=1.5/dy;  D1v(end,end-1)=-2.0/dy; D1v(end,end-2)= 0.5/dy;

D2v = tridiag(-2/dy^2,1/dy^2,1/dy^2,Ny);
D2v(1,1)  =1/dy^2; D2v(1,2)    =-2/dy^2; D2v(1,3)    =1/dy^2;
D2v(Ny,Ny)=1/dy^2; D2v(Ny,Ny-1)=-2/dy^2; D2v(Ny,Ny-2)=1/dy^2;

Urhos = U.*rhos;
UrhosD2  = tridiagvec(Urhos(:).*(-2/dy^2),Urhos(2:end).*(1/dy^2), ...
                  Urhos(1:end-1).*(1/dy^2));
Urhosp = U.*(D1*rhos);
UrhospD1  = tridiagvec(Urhosp(:).*(0.0),Urhosp(2:end).*(-0.5/dy), ...
                  Urhosp(1:end-1).*(0.5/dy));
rhosUpp = rhos.*(D2*U);
rhospUp = (D1*rhos).*(D1*U);
rhosD2  = tridiagvec(rhos(:).*(-2/dy^2),rhos(2:end).*(1/dy^2), ...
                  rhos(1:end-1).*(1/dy^2));
rhosp = D1*rhos;
rhospD1 = tridiagvec(rhosp(:).*(0.0),rhosp(2:end).*(-0.5/dy), ...
                  rhosp(1:end-1).*(0.5/dy));



