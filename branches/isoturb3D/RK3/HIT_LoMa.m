%-----------------------------------------------------
% Homegenous isotropic turbulence                     %
%                                                     %
% Solve 3D incompresible N-S equations                %
% using the Low-Mach formulation                      %
%                                                     %
% dphi_hat/dt = RHS_phi_hat                           %
% dome_hat/dt = RHS_ome_hat                           %
% dtem_hat/dt = RHS_tem_hat                           %
% drho_hat/dt = - k^2 psi_hat                         %
%                                                     %
% k^2   = kx^2 + ky^2 + kz^2                          %
%                                                     %
% RHSphi_hat = 
% RHSome_hat = 
% RHStem_hat = 
%                                                     %
% Spatial discretization: Fourier                     %
% Temporal discretization: Explicit/Implicit Euler    %
%                                                     %
% ()_hat = Fourier coefficient                        %
%                                                     %
% adrian  Dic-2012  (incompressible, never worked)    %
% oflores Dic-2014  compressible Explicit/Implicit Eu.%
% aalmagro Feb-2015 compressible RK3                  %
%-----------------------------------------------------%

%--- list of buffers (32) ---% 
% kx =                % phi_hat =         % RHSphi_hat=
% ky =                % ome_hat =         % RHSome_hat=
% kz =                % psi_hat =         % RHStem_hat=
% k2 =                % tem_hat =         % RHSpsi_hat=
% k2o=                % tem   =           % t11_hat =  
% dealias=            % mx_hat=           % t12_hat =  
% xx =                % my_hat=           % t13_hat =  
% yy =                % mz_hat=           % t22_hat =  
% zz =                % ru    =           % t23_hat =  
% kkk =               % rv    =           % t33_hat =   
% qqq =               % rw    =           
% alpha =             %gamma  =
% beta  =             %delta  =

                                          

%----------------------------%
%RK3 coefficients
%this family of RK3 parameters have been created for lomacte
%-----------------------------% 
%p0 = 0.0 %alpha1
%p1 = -2.0 %beta2
p0 = 0.1 %alpha1
p1 = 0.28 %beta2
xi = [ 0.0, -17.0/60.0 , -5.0/12.0]
gama  = [ 8.0/15.0, 5.0/12.0 , 3.0/4.0] 
b3 = (p1*(gama(1)*gama(3)+gama(2)*gama(3)+ xi(2)*gama(3)+gama(1)*xi(3))-gama(1)*gama(2)*gama(3)+gama(1)*gama(2)*xi(3))/(p1-gama(1)*gama(2))
beta =[gama(1)-p0, p1, b3]
alpha = [ p0,gama(2)+xi(2)-p1,gama(3)+xi(3)-beta(3)]
%-----------domain-----------%
ii    = 0; % 0 -> start form Taylor-Green. else, continue from file. 
taylor2d=0; 
taylor3d=0;   
nsave = 500;
nstats= 10; 
nmax  = 500000;   % max steps
Re    = 500;     % inverse of viscosity
Pr    = 0.7;       % ratio between diffusivity and viscosity
r     = 8;       % density ratio
T     = 6;      % total time
CFL   = 0.25;     % CFL
Lx    = pi;    % x length
%Lx    = 2*pi;    % x length
Ly    = pi;    % y length
Lz    = pi;    % z length
nx    = 96*3/2;   % x-modes
ny    = 96*3/2;  % y-modes
nz    = 96*3/2;  % z-modes
basename = 'results/s80pi96Re500RKold'
%----------------------------%
dx    = Lx/nx; %(nx-1);
dy    = Ly/ny; %(ny-1);
dz    = Lz/nz; %(nz-1);
ndim  = [ny nx nz]; % kx,ky,kz are generated with a meshgrid !!!

iu=sqrt(-1);  
%---------display------------%
disp(['Reynolds number:      ' ,num2str(Re)])
disp(['Density ratio:        ' ,num2str(r)])
disp(['x length Lx:          ' ,num2str(Lx)])
disp(['y length Ly:          ' ,num2str(Ly)])
disp(['y length Lz:          ' ,num2str(Lz)])
disp(['CFL number:           ' ,num2str(CFL)])
disp(['number of nx modes:   ' ,int2str(nx)])
disp(['          ny modes:   ' ,int2str(ny)])
disp(['          nz modes:   ' ,int2str(nz)])
disp(['Resolution dx :       ' ,num2str(dx)])
disp(['           dy :       ' ,num2str(dy)])
disp(['           dz :       ' ,num2str(dz)])
disp(['total time simulated: ' ,num2str(T)])


%----------------------------------------------%
%-----initialize Fourier------%
%----------------------------------------------%
[kx ky kz] = meshgrid( [0:nx/2-1 -nx/2:-1], [0:ny/2-1 -ny/2:-1], [0:nz/2-1 -nz/2:-1] );

% Cutting frequencies 2/3 rule: from n/3 to n/2-1
dealias = ( abs(kx)<2/3*(nx/2) & abs(ky)<2/3*(ny/2) & abs(kz)<2/3*(nz/2) );
%dealias(1,1,1) = 0;

kx = 2*pi*kx/Lx;
ky = 2*pi*ky/Ly;
kz = 2*pi*kz/Lz;
k2 = kx.^2 + ky.^2 + kz.^2;

% singular point pressure (0,0,0)-mode undetermined
k2o        = k2;
k2o(1,1,1) = 1;

%----------------------------------------------%
%------initial conditions----%
%----------------------------------------------%
if ii==0,
    t    = 0;
    dt   = 0; % DON'T KNOW dt YET !!!
    j    = 0;

%   Taylor 2D initial conditions (test case): 
    if taylor2d==1; 
    
       A=1;   aa=1; 
       B=-1;  bb=1; 
   %   continuity implies they have to satisfy Aa+Bb+Cc=0; 
       [xx yy zz] = meshgrid([dx:dx:Lx],[dy:dy:Ly],[dz:dz:Lx]); 
       ru = A.*sin(aa.*xx).*cos(bb.*yy); 
       rv = B.*cos(aa.*xx).*sin(bb.*yy); 
       rw = zeros(size(ru)); 
       tem= ones(size(ru)); 

    elseif taylor3d==1; 

       A=1;   aa=1; 
       B=-1;  bb=1; 
       C=0;   cc=1; 
   %   continuity implies they have to satisfy Aa+Bb+Cc=0; 
       [xx yy zz] = meshgrid([dx:dx:Lx],[dy:dy:Ly],[dz:dz:Lx]); 
       ru = A.*cos(aa.*xx).*sin(bb.*yy).*sin(cc.*zz); 
       rv = B.*sin(aa.*xx).*cos(bb.*yy).*sin(cc.*zz); 
       rw = C.*sin(aa.*xx).*sin(bb.*yy).*cos(cc.*zz); 
       tem= ones(size(ru)); 
       %clear xx yy zz

    else 

       [xx yy zz] = meshgrid([dx:dx:Lx],[dy:dy:Ly],[dz:dz:Lx]); 
       fname = [basename '.' int2str(ii) '.mat'];
       load(fname); % reads ru, rv, rw, and tem
       disp(['loading: ',fname])
       
       % prescribe temperature field as 1 + t0*cos(2*pi*z/Lz)
       % where t0 = (r-1)./(r+1); 

       tem = 1.0 + (r-1)./(r+1).*cos(2*pi*zz./Lz); 
       disp(['min. temp:',num2str(min(tem(:)))])
       disp(['max. temp:',num2str(max(tem(:)))])

       % tem = tem-mean(tem(:))
    end 

else
    j     = 0; % nsave*ii;
    fname = [basename '.' int2str(ii) '.mat'];
    load(fname); % reads ru, rv, rw, and tem
    disp(['loading: ',fname])
end 


%----------------------------------------------%
%-----------viscous CFL-------------%   (see end of file form more details)
%----------------------------------------------%
lam_v = 2./Re.*(1./dx.^2 + 1./dy.^2 + 1./dz.^2)./min(1.0./tem(:)); 
%min density modifies visc. time



% Compute phi, ome, and psi from ru,rv,rw and fft everything: 
mx_hat  = fftn(ru ,ndim); mx_hat  = mx_hat.*dealias; 
my_hat  = fftn(rv ,ndim); my_hat  = my_hat.*dealias; 
mz_hat  = fftn(rw ,ndim); mz_hat  = mz_hat.*dealias; 
tem_hat = fftn(tem,ndim); tem_hat = tem_hat.*dealias; 

mx00_hat = squeeze(mx_hat(:,1,1)); 
mz00_hat = squeeze(mz_hat(:,1,1)); 

mx00_hatwk = mx00_hat; %initialization
mz00_hatwk = mz00_hat; %initialization

psi_hat = (iu.*kx.*mx_hat + iu.*ky.*my_hat + iu.*kz.*mz_hat)./(-k2o); 
phi_hat = -k2.*(my_hat - iu.*ky.*psi_hat); 
ome_hat = iu.*kz.*mx_hat - iu.*kx.*mz_hat; 
RHStem_hat = psi_hat;%just some initialization is needed for first substep 
RHSphi_hat = psi_hat;%just some initialization is needed for first substep 
RHSome_hat = psi_hat;%just some initialization is needed for first substep 

% Compute the divergence of the velocity (sanity check)
kkk =(iu.*kx.*fftn(ru.*tem,ndim) + ... 
      iu.*ky.*fftn(rv.*tem,ndim) + ...
      iu.*kz.*fftn(rw.*tem,ndim) ).*dealias; % Div(v)
lastdiv = [min(abs(kkk(:))), max(abs(kkk(:)))]; 

%----------------------------------------------%
%-----------prepares stats and spe files-------%
%----------------------------------------------%
[ksp,sp] = compuspe(kx,ky,kz,mx_hat); 

header0 = '%% Re,  Lx,  Ly,  Lz,  nx,  ny,  nz \n'
header1 =['%%' num2str(Re) ' ' num2str(Lx) ' ' num2str(Ly) ' ' num2str(Lz) ' ' num2str(nx) ' ' num2str(ny) ' ' num2str(nz) '\n']; 

fid=fopen([basename '.stats'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 

fid=fopen([basename '.spe_u'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 
caca=[-1, ksp]; save([basename '.spe_u'],'-ascii','caca','-append'); 

fid=fopen([basename '.spe_v'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 
caca=[-1, ksp]; save([basename '.spe_v'],'-ascii','caca','-append'); 

fid=fopen([basename '.spe_w'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 
caca=[-1, ksp]; save([basename '.spe_w'],'-ascii','caca','-append'); 

fid=fopen([basename '.spe_T'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 
caca=[-1, ksp]; save([basename '.spe_T'],'-ascii','caca','-append'); 

fid=fopen([basename '.spe_D'],'w'); fprintf(fid,header0); fprintf(fid,header1); fclose(fid); 
caca=[-1, ksp]; save([basename '.spe_D'],'-ascii','caca','-append'); 

disp('!------------------------------!')
disp('      starting computation      ')
disp('!------------------------------!')

%----------------------------------------------%
%------temporal loop---------%
%----------------------------------------------%
tic
while t<=T
    
    j = j+1;
    if j>nmax,break,end
        %start substep loop
        for rkstep=1:3
	   	    % compute my: phi_hat = -k^2.*my_hat
	    my_hat  = phi_hat./(-k2o); 
	    err_my = norm(squeeze(my_hat(:,1,1))); % div(m)=0 requires this to be zero
	
	    % mx, mz: 
	    %   ome_hat = i*kz*mx_hat - i*kx*mz_hat    
	    %  -dmy_hat = i*kx*mx_hat + i*kz*mz_hat    
	    mx_hat = ( iu.*kz.*ome_hat + kx.*ky.*my_hat)./(-kx.^2 - kz.^2); 
	    mz_hat = (-iu.*kx.*ome_hat + kz.*ky.*my_hat)./(-kx.^2 - kz.^2); 
	
	    % need to put the (kx,kz)=(0,0) modes back into mx,my and mz: 
	    mx_hat(:,1,1)=mx00_hat; 
	    my_hat(:,1,1) = 0; 
	    mz_hat(:,1,1)=mz00_hat; 


	    % add grad(psi): 
	    mx_hat = mx_hat+iu.*kx.*psi_hat;
	    my_hat = my_hat+iu.*ky.*psi_hat;
	    mz_hat = mz_hat+iu.*kz.*psi_hat;
	
	    % ru,rv,rw and tem to Phys to compute tij_hat and RHStem: 
	    ru     = real( ifftn(mx_hat,ndim) );
	    rv     = real( ifftn(my_hat,ndim) );
	    rw     = real( ifftn(mz_hat,ndim) );
	    tem    = real( ifftn(tem_hat,ndim) );

	
	    % BEGIN write to disk
	    %----------------------------------------------%
	    if (mod(j,nsave)==0 || j==1) && rkstep == 1,
		disp('------------------------------')
	        toc
	        ii = ii+1;
	        disp(['total time simulated: ',num2str(t)] )
	        disp(['dt:   '                ,num2str(dt)])
	        disp(['step: '                ,num2str(j)])
	        
	        save([basename,'.',int2str(ii),'.mat'],'ru','rv','rw','tem','Lx','Ly','Lz','Re','t');
	        tic
	    end
	    %----------------------------------------------%
	    % END write to disk
	    %----------------------------------------------%
	
	
	    % Compute tij = r*u_i*u_j
	    %
	    kkk = ru.*ru.*tem; 
	    t11_hat=fftn(kkk,ndim); t11_hat = t11_hat.*dealias; 
	    %
	    kkk = ru.*rv.*tem; 
	    t12_hat=fftn(kkk,ndim); t12_hat = t12_hat.*dealias; 
	    %
	    kkk = ru.*rw.*tem; 
	    t13_hat=fftn(kkk,ndim); t13_hat = t13_hat.*dealias; 
	    %
	    kkk = rv.*rv.*tem; 
	    t22_hat=fftn(kkk,ndim); t22_hat = t22_hat.*dealias; 
	    %
	    kkk = rv.*rw.*tem; 
	    t23_hat=fftn(kkk,ndim); t23_hat = t23_hat.*dealias; 
	    %
	    kkk = rw.*rw.*tem; 
	    t33_hat=fftn(kkk,ndim); t33_hat = t33_hat.*dealias; 
	    %
	
	    % compute u, v, w and save them uin ru, rv, rw: 
	    ru = ru.*tem; 
	    rv = rv.*tem; 
	    rw = rw.*tem; 
	
            
	    % compute time step (see end of file form more details)
	    if sum(isnan(ru(:)))~=0; disp('nans in ru !!'); return; end 
	    umax = max( abs(ru(:)) );
	    vmax = max( abs(rv(:)) );
	    wmax = max( abs(rw(:)) );
	    lam_c = umax./dx + vmax./dy + wmax./dz; 
	    dt   = CFL./(lam_c+lam_v); 
	
	    % write u,v and w in PHYS into mx, my, mz: 
	    mx_hat=fftn(ru,ndim); mx_hat = mx_hat.*dealias; 
	    my_hat=fftn(rv,ndim); my_hat = my_hat.*dealias; 
	    mz_hat=fftn(rw,ndim); mz_hat = mz_hat.*dealias; 
            % Now: mx_hat = u, my_hat = v, mz_hat = w
	    %----------------------------------------------%
	    % BEGIN evaluate stats: 
	    %----------------------------------------------%
	    if (mod(j,nstats)==0 || j==1) && rkstep==1, 
	
	%      spectra
	       [ksp,sp] = compuspe(kx,ky,kz,mx_hat); 
	       caca=[t, sp]; save([basename '.spe_u'],'-ascii','caca','-append'); 
	
	       [ksp,sp] = compuspe(kx,ky,kz,my_hat); 
	       caca=[t, sp]; save([basename '.spe_v'],'-ascii','caca','-append'); 
	
	       [ksp,sp] = compuspe(kx,ky,kz,mz_hat); 
	       caca=[t, sp]; save([basename '.spe_w'],'-ascii','caca','-append'); 
	       
	       [ksp,sp] = compuspe(kx,ky,kz,tem_hat); 
	       caca=[t, sp]; save([basename '.spe_T'],'-ascii','caca','-append'); 
               
               rhom = mean(1.0./tem(:));
	
	%      pseudo-dissipation 
	       %indeed this eps = mean(rho Eps)
	       eps = 1/Re.*sum(k2(:).*(mx_hat(:).*conj(mx_hat(:)) + ... 
	                               my_hat(:).*conj(my_hat(:)) + ... 
	                               mz_hat(:).*conj(mz_hat(:))))./length(kx(:)).^2;        
               eps = eps./rhom;
	%      volume averages
	       meanu = mx_hat(1,1,1)./length(kx(:)); 
	       meanv = my_hat(1,1,1)./length(kx(:)); 
	       meanw = mz_hat(1,1,1)./length(kx(:)); 
	       meant = tem_hat(1,1,1)./length(kx(:)); 
	       eneru = sum(mx_hat(2:end).*conj(mx_hat(2:end)))./length(kx(:)).^2; 
	       enerv = sum(my_hat(2:end).*conj(my_hat(2:end)))./length(kx(:)).^2; 
	       enerw = sum(mz_hat(2:end).*conj(mz_hat(2:end)))./length(kx(:)).^2; 
	       enert = sum(tem_hat(2:end).*conj(tem_hat(2:end)))./length(kx(:)).^2; 
	
	%      Kolmogorov lengthscale: 
               %nu = mu/rho = 1/(rho*Re)
	       eta = ((1./(rhom*Re)).^3./eps).^(1./4); 
	       %eta = (1./Re.^3./eps).^(1./4); 
	
	%      Taylor microscale: 
	       lam = sqrt(15./Re.*(eneru + enerv + enerw)./3./eps); 
	       relam = Re.*sqrt((eneru + enerv + enerw)./3).*lam; 
	
	%      save and print:
	       caca=[t, meanu, meanv, meanw, meant, eneru, enerv, enerw,enert, eps, relam, eta, err_my];  
	       save([basename '.stats'],'-ascii','caca','-append'); 
	
	       disp([num2str(t    ),',',num2str(dt),' , u00, t00:', ... 
	             num2str(mean([meanu meanv meanw])),' , ', ... 
	             num2str(meant),' , u^2, T^2:', ... 
	             num2str(mean([eneru enerv enerw])),' , ', ... 
	             num2str(enert),' , eps:', ... 
	             num2str(eps  ),' , Re_lam:', ... 
	             num2str(relam),' , meanrho:', ... 
	             num2str(rhom),' , eta:', ... 
	             num2str(eta  ),' , err:', ... 
	             num2str(err_my),',',num2str(lastdiv(1)),',',num2str(lastdiv(2))]); 
	
	      if taylor2d==1, 
	      figure(1); contour(xx(:,:,1),yy(:,:,1),ru(:,:,1),[-0.75:0.25:0.75],'k'); hold on 
	                 contour(xx(:,:,1),yy(:,:,1),exp(-2.*t./Re).*sin(xx(:,:,1)).*cos(yy(:,:,1)),[-0.75:0.25:0.75],'r--'); 
	                 hold off
	      figure(2); contour(xx(:,:,1),yy(:,:,1),rv(:,:,1),[-0.75:0.25:0.75],'k'); hold on 
	                 contour(xx(:,:,1),yy(:,:,1),-exp(-2.*t./Re).*cos(xx(:,:,1)).*sin(yy(:,:,1)),[-0.75:0.25:0.75],'r--'); 
	                 hold off
	
	      pause
	      end 
	
	      if taylor3d==1, 
	
	
	      ruTaylor = tay3(A,aa,Re,t,xx,yy,zz); 
	      figure(1); contour(xx(:,:,7),yy(:,:,7),ru(:,:,7),[-0.75:0.25:0.75],'k'); hold on 
	                 contour(xx(:,:,7),yy(:,:,7),ruTaylor(:,:,7),[-0.75:0.25:0.75],'r--'); 
	                 hold off
	
	      pause
	      end 
	
	    end 
	    %----------------------------------------------%
	    % END evaluate stats: 
	    %----------------------------------------------%



	
	    % compute RHStem: 
	    kkk = iu.*kx.*tem_hat; 
	    kkk = real( ifftn(kkk,ndim) );
	    qqq = -ru.*kkk; 
	
	    kkk = iu.*ky.*tem_hat; 
	    kkk = real( ifftn(kkk,ndim) );
	    qqq = qqq - rv.*kkk; 
	
	    kkk = iu.*kz.*tem_hat; 
	    kkk = real( ifftn(kkk,ndim) );
	    qqq = qqq - rw.*kkk; % qqq is saved for later !!! 


	
	    RHStem_hatp = RHStem_hat; %previous RHStem 
	    RHStem_hat=fftn(qqq,ndim); RHStem_hat = RHStem_hat.*dealias; % convective
	    RHStem_hat = RHStem_hat - k2.*tem_hat./(Re*Pr);              % diffusion
	
	    % Compute Nx, Ny, Nz (RHS of the momentum eq: advec+viscous): 
	    % OVERWRITE m#_hat !!! 
	
	    kkk = iu.*kx.*mx_hat + iu.*ky.*my_hat + iu.*kz.*mz_hat; % Div(v)
	    lastdiv = [min(abs(kkk(:))), max(abs(kkk(:)))]; 
	 
	    mx_hat = -iu.*kx.*t11_hat - iu.*ky.*t12_hat - iu.*kz.*t13_hat + (-k2.*mx_hat + iu.*kx.*kkk./3)./Re; 
	    my_hat = -iu.*kx.*t12_hat - iu.*ky.*t22_hat - iu.*kz.*t23_hat + (-k2.*my_hat + iu.*ky.*kkk./3)./Re; 
	    mz_hat = -iu.*kx.*t13_hat - iu.*ky.*t23_hat - iu.*kz.*t33_hat + (-k2.*mz_hat + iu.*kz.*kkk./3)./Re; 

	    %evolving 00 modes with xi part before overwriting mx00_hatwk
	    mx00_hat = mx00_hat + xi(rkstep)*dt.*mx00_hatwk;
	    mz00_hat = mz00_hat + xi(rkstep)*dt.*mz00_hatwk;

	    %add last RHS part and then overwrite buffer
	    phi_hat  =  phi_hat + xi(rkstep)*dt.*RHSphi_hat;
	    ome_hat  =  ome_hat + xi(rkstep)*dt.*RHSome_hat;
	
	    % Compute RHS of ome and phi
	    RHSome_hat = iu.*kz.*mx_hat - iu.*kx.*mz_hat; 
	    RHSphi_hat = -(kx.^2 + kz.^2).*my_hat - iu.*ky.*(iu*kx.*mx_hat + iu.*kz.*mz_hat);  
            % RHS part of 00 modes
            mx00_hatwk = mx_hat(:,1,1);
            mz00_hatwk = mz_hat(:,1,1);
	
	    % RK present substep terms ome, phi, tem and mx00_hat and mz00_hat: 
	    phi_hat  =  phi_hat + gama(rkstep)*dt.*RHSphi_hat;
	    ome_hat  =  ome_hat + gama(rkstep)*dt.*RHSome_hat;
	    mx00_hat = mx00_hat + gama(rkstep)*dt.*mx00_hatwk;
	    mz00_hat = mz00_hat + gama(rkstep)*dt.*mz00_hatwk;
            % Evolve temperature
	    tem_hat  =  tem_hat + gama(rkstep)*dt.*RHStem_hat + xi(rkstep)*dt.*RHStem_hatp;
	
	
	    % Implicit Euler for psi: 
            %evolve now tem
            %qqq now has phys of T(i+1)
	    %tem has last temperature in physical domain
	    kkk = gama(rkstep)*dt.*RHStem_hat + xi(rkstep)*dt.*RHStem_hatp;
            kkk = real(ifftn(kkk,ndim));
            qqq  = real (ifftn(tem_hat,ndim)); %qqq=phys of T(i+1)
	    kkk = -kkk./(tem.*qqq); %kkk=phys of  
	    %kkk = 1.0./qqq - 1.0./tem; %kkk=phys of  
            %drho = 1/T(i+1) - 1/T(i)
	    RHSpsi_hat=fftn(kkk,ndim); RHSpsi_hat = RHSpsi_hat.*dealias; 
            psi_hat=-RHSpsi_hat./(-k2o)./(dt.*beta(rkstep))-alpha(rkstep)/beta(rkstep).*psi_hat;
	        
        end
     %------------END OF RKSTEP LOOP-------------------------------------%
	t = t + dt;
end

% NOTES ABOUT CFL AND dt. FROM LESOC3:
%
% The formula below is given (not proved) by u. schumann (1975),
% linear stability of finite difference equations for
% three-dimensional flow problems, journal of computational
% physics 18, 465-470
%                                   1
%    dt <=   ----------------------------------------------
%             u     v      w          1         1       1
%            -- +  --  +  --  + 2k ( ---    +  ---  +  ---  )
%            dx    dy     dz         dx**2    dy**2   dz**2
%

