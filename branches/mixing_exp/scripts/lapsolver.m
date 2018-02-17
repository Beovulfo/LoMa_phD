function [x1,v,dvdy]=lapsolver(kx,kz,N,yy,Ly,rhs)
%this functions solve laplacian for v knowing a rhs
k2=kx^2+kz^2;

[x0,DM] = chebdif(N+2,2);                   % Compute second derivative
D2 = DM(2:N+1,2:N+1,2);                    % Enforce Dirichlet BCs
D1 = DM(2:N+1,2:N+1,1);                    % Enforce Dirichlet BCs
                                     
%[x,D4] = cheb4c(N+2);                      % Compute fourth derivative
%I = eye(size(D4));                         % Identity matrix
%Z = zeros(size(D4));                       % Identity matrix    

[x,D3, D4] = cheb4c_mod(N+2);                      % Compute fourth derivative
I = eye(size(D4));                         % Identity matrix
Z = zeros(size(D4));                       % Identity matrix    
%
fac1 = 1/Ly; % factor to re-scale equations
x1 = x/fac1;
D1 = D1*fac1;
D2 = D2*fac1^2;
D4 = D4*fac1^4;


rhs       = interp1(yy,rhs,x1,'spline');
figure
plot(x1,rhs)

%Left operator => Laplacian operatior
%

O2  = D2-k2*I;
O2i = inv(O2);

v=O2i*rhs;
dvdy=D1*v;
%A_11 = O2i*((D4-2*k2*D2+(k2^2)*I)/Re+1i*alp*diag(usecond)-1i*alp*diag(u)*(D2-k2*I)); 

%A_21 = -1i*bet*diag(ufirst);

%A_22 = -1i*alp*diag(u) + (D2-k2*I)/Re;


% 
% computing eigenvalues and eigenvectors and order them
%[v,e]=eig(A);
%ee=squeeze(diag(e));
%[tmp,is]=sort(real(ee),'descend');
%ee=ee(is);
%v=v(:,is);
%if real(ee(1))>0
%    display('Unstable')
%end
%ome=v(N+1:2*N,:);
%phi=D2*v(1:N,:)-k2*v(1:N,:);
%dvdy=D1*v(1:N,:);
%v=v(1:N,:);
%
%ceros=zeros(1,2*N);
%ome = [ceros;ome;ceros]; %Considering BC's phi=0,ome=0 in both extremes
%phi = [ceros;phi;ceros];
