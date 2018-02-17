function [v,phi,ome,ee,x]=ossq_for_toni(alp,bet,Reb,N)

%  The script file orrsom.m computes the eigenvalues of the Orr-Sommerfeld
%  equation using NxN Chebyshev differentiation matrices.

% S.C. Reddy, J.A.C. Weideman 1998.  Code modified to display output
% by JACW, May 2003.

%N=29;Re=2000;alp=1;
%alp=1;
%bet=0;
k2=alp^2+bet^2;
Re = 3*Reb/2; % Re based on centerline velocity
%Re = 3*Reb/2; % Re based on centerline velocity
%Re=10000;

%N = input(' Order of the differentiation matrix: N = ? ');
%R = input(' Reynolds number: R = ? ');
%i = sqrt(-1);

[x,DM] = chebdif(N+2,2);                   % Compute second derivative
D2 = DM(2:N+1,2:N+1,2);                    % Enforce Dirichlet BCs
D1 = DM(2:N+1,2:N+1,1);                    % Enforce Dirichlet BCs
                                     
%[x,D4] = cheb4c(N+2);                      % Compute fourth derivative
%I = eye(size(D4));                         % Identity matrix
%Z = zeros(size(D4));                       % Identity matrix    

[x,D3, D4] = cheb4c_mod(N+2);                      % Compute fourth derivative
I = eye(size(D4));                         % Identity matrix
Z = zeros(size(D4));                       % Identity matrix    

u       = 1-x.^2;
ufirst  = -2*x;
usecond = 0*x -2; 

%
%
% Set up A and Ac
%

O2  = D2-k2*I;
O2i = inv(O2);

A_11 = O2i*((D4-2*k2*D2+(k2^2)*I)/Re+1i*alp*diag(usecond)-1i*alp*diag(u)*(D2-k2*I)); 

A_21 = -1i*bet*diag(ufirst);

A_22 = -1i*alp*diag(u) + (D2-k2*I)/Re;

A  = [A_11 Z ; A_21 A_22 ];

% 
% computing eigenvalues and eigenvectors and order them
[v,e]=eig(A);
ee=squeeze(diag(e));
[tmp,is]=sort(real(ee),'descend');
ee=ee(is);
v=v(:,is);
if real(ee(1))>0
    display('Unstable')
end
ome=v(N+1:2*N,:);
phi=D2*v(1:N,:)-k2*v(1:N,:);

