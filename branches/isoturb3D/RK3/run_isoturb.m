function [] = run_isoturb(N,r,Re,CFL,T)
%N = 128;
nx = N*3/2;
ny = N*3/2;
nz = N*3/2;
Lx = 2*pi;
Ly = Lx;
Lz = Lx;
ii = 0;
path0 = '/share/drive/toni/isoturb3D/';
%path0 = 'Re500s80/';
%path0 = 'results_Re200/';
%CFL = 0.2
b2 = [-1.0]
%Re = 300
%Re = 350
%Re = 160
%T = 2.0 %Max time
%First run incompressible to check if it works
disp(['Running constant density case'])
%p0 = 0.0
p1 = -1.0
basename0 = [path0 'case_inc']
init_isoturb3D(basename0,N,Lx,Ly,Lz,ii)
isoturb3D(basename0,ii,Re,1.0,Lx,Ly,Lz,nx,ny,nz,CFL,T,p1)
%isoturb3D(basename0,ii,Re,1.0,Lx,Ly,Lz,nx,ny,nz,CFL,T,p0,p1)
%Run now the real cases
%r = 8.0
for i=1:length(b2)
	%p0 = a1(i)
	p1 = b2(i)
	disp(['Running case '  num2str(i)])
	disp(['beta2  =  '  num2str(p1)])
	basename = [path0 'case_' int2str(i)]
	%no need of call init_istoturb3D again!
	copyfile([basename0 '.0.mat'],[basename '.0.mat']);
	isoturb3D(basename,ii,Re,r,Lx,Ly,Lz,nx,ny,nz,CFL,T,p1)
	%isoturb3D(basename,ii,Re,r,Lx,Ly,Lz,nx,ny,nz,CFL,T,p0,p1)
end

