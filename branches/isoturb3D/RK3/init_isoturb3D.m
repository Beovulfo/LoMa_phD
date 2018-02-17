%-----------------------------------------------------
% Homegenous isotropic turbulence                     %
%                                                     %
% Generate initial velocity field, following          %
% Mansour & Wray, 1993                                %
%                                                     %
% Input parameters:                                   %
% nx,ny,nz: number of modes                           %
% kref:   peak wavenumber for the i.c.                %
% sigma:  spectrum ~ k.^sigma for k<kref              %
% q:      ref. value of the kinetic energy            %
%                                                     %
% oflores Oct-2014                                    %
% aalmagro mar-2015 converted to funcion %
%-----------------------------------------------------%

%--- list of buffers (32) ---% 
% kx = 
% ky = 
% kz = 
% k2 = 
% k2o= 
% tem = 
% mx_hat=
% my_hat=
% mz_hat=
% ru    =
% rv    = 
% rw    = 

%-----------domain-----------%
function [] = init_isoturb3D(basename,N,Lx,Ly,Lz,ii)
%N = 96
%ii    = 0; 
%Lx    = pi;     % x length
%Lx    = 2*pi;     % x length
%Ly    = pi;     % y length
%Lz    = pi;     % z length

nx    = N*3/2;  % x-modes
ny    = N*3/2;  % y-modes
nz    = N*3/2;  % z-modes
%----------spe params--------%
sigma = 2;
kref  = 6;
%----------------------------%
dx    = Lx/nx; %(nx-1);
dy    = Ly/ny; %(ny-1);
dz    = Lz/nz; %(nz-1);
ndim  = [ny nx nz];

iu=sqrt(-1);  
%---------display------------%
disp(['x length Lx:          ' ,num2str(Lx)])
disp(['y length Ly:          ' ,num2str(Ly)])
disp(['y length Lz:          ' ,num2str(Lz)])
disp(['number of nx modes:   ' ,int2str(nx)])
disp(['number of ny modes:   ' ,int2str(ny)])
disp(['number of nz modes:   ' ,int2str(nz)])


%-----initialize Fourier------%
[kx ky kz] = meshgrid( [0:nx/2-1 -nx/2:-1], [0:ny/2-1 -ny/2:-1], [0:nz/2-1 -nz/2:-1] );

% Cutting frequencies 2/3 rule: from n/3 to n/2-1
dealias = ( abs(kx)<2/3*(nx/2) & abs(ky)<2/3*(ny/2) & abs(kz)<2/3*(nz/2) );
%dealias(1,1,1) = 0;

kx = 2*pi*kx/Lx;
ky = 2*pi*ky/Ly;
kz = 2*pi*kz/Lz;
k2 = kx.^2 + ky.^2 + kz.^2;
k2o = k2; k2o(1,1,1)=1; 

%------initial conditions----%
ii==0,
mx_hat = zeros(size(kx)); 
my_hat = zeros(size(kx)); 
mz_hat = zeros(size(kx)); 


tvec=0; trand=0; trest=0; 
%for i=2:length(kx(:)); % all modes except (0,0,0). 
listmodes = find(dealias(:)==1); 
for i=listmodes(2:end)'; % all non-dealiased modes except (0,0,0). 

    tic; 
    % definition of e1 and e2. Orthonormal vectors, perpendicular to k
    e0 = [kx(i) ky(i) kz(i)]; e0 = e0./sqrt(e0(1).^2 + e0(2).^2 + e0(3).^2); 
    e1 = cross(e0,[1 0 0]); 
    if norm(e1)<1e-3; e1 = cross(e0,[0 1 0]); end 
    e2 = cross(e0,e1); 
    tvec=tvec+toc; tic; 

    % random numbers: 
    theta1 = rand().*2*pi; 
    theta2 = rand().*2*pi; 
    theta3 = rand().*2*pi; 
    trand=trand+toc; tic; 

    % spe: 
    spe = k2(i).^(sigma/2).*exp(-sigma./2.*(k2(i)./kref.^2)); 

    if ~mod(i,100000); 
        %loglog(sqrt(k2(i)),spe,'o'); hold on;
        %pause(0.1); 
        
        disp([ num2str(tvec+trand+trest) 's, ' num2str(100*i./length(kx(:))) '%']); 
    end 

    % alpha 
    alpha = spe./4./pi./k2(i).*exp(iu.*theta1).*cos(theta3); 
    beta  = spe./4./pi./k2(i).*exp(iu.*theta2).*sin(theta3); 

    % velocities: 
    mx_hat(i) =  e1(1).*alpha + e2(1).*beta; 
    my_hat(i) =  e1(2).*alpha + e2(2).*beta; 
    mz_hat(i) =  e1(3).*alpha + e2(3).*beta; 

    trest=trest+toc; 
end 

% ENFORCE A REAL FIELD IN THE FFTN OF mx_hat, my_hat and mz_hat: 
mx_hat = fftn(real(ifftn(mx_hat,ndim)),ndim).*dealias; 
my_hat = fftn(real(ifftn(my_hat,ndim)),ndim).*dealias; 
mz_hat = fftn(real(ifftn(mz_hat,ndim)),ndim).*dealias; 

% The proyection into a real space spoils the null divergence of the vector field. 
% Enforce it again: 
psi_hat = (-iu.*kx.*mx_hat - iu.*ky.*my_hat - iu.*kz.*mz_hat)./(-k2o); 
mx_hat = mx_hat - iu.*kx.*psi_hat; 
my_hat = my_hat - iu.*ky.*psi_hat; 
mz_hat = mz_hat - iu.*kz.*psi_hat; 

% RESCALE VELOCITY SO THAT 0.5*(u^2 + v^2 + w^2) = 1 !!!
e = sum(mx_hat(:).*conj(mx_hat(:)) + ... 
        my_hat(:).*conj(my_hat(:)) + ...   
        mz_hat(:).*conj(mz_hat(:)) )./length(kx(:)).^2; 

mx_hat = mx_hat./sqrt(e); 
my_hat = my_hat./sqrt(e); 
mz_hat = mz_hat./sqrt(e); 

% COMPUTE ru, rv and rw AND SAVE
ru = real(ifftn(mx_hat,ndim)); 
rv = real(ifftn(my_hat,ndim)); 
rw = real(ifftn(mz_hat,ndim)); 
tem= ones(size(ru)); 
fname = [basename '.' int2str(0) '.mat'] 
save(fname,'ru','rv','rw','tem','Lx','Ly','Lz','kref','sigma');

