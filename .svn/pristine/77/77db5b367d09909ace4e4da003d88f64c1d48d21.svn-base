N = 96;
%N = 128;
nx = N*3/2;
ny = N*3/2;
nz = N*3/2;
Lx = 2*pi;
Ly = Lx;
Lz = Lx;
ii = 0;
path0 = '/data2/toni/work/isoturb3D/';
%path0 = 'results_Re200/';
CFL = 0.1
a1 =  [0.1 ,0.1,0.0,29.0/96.0,-0.2,     2.0,-1.0]
b2 =  [0.28,0.3,0.3, 5.0/24.0,0.4 ,-2.0/3.0,-2.0]
%a1 = [-1.0,-1.0,-1.0, 2.0,4.0,0.0,0.0,1.0,0.1]
%b2 = [-1.0, 0.3, 3.0,-2.0,4.0,-2.0,4.0,1.0,0.28]
%Re = 300
Re = 280
%Re = 160
T = 2.0 %Max time
%First run incompressible to check if it works
p0 = 0.0
p1 = -2.0
basename0 = [path0 'case_inc']
%init_isoturb3D(basename0,N,Lx,Ly,Lz,ii)
%disp(['Running constant density case'])

%isoturb3D(basename0,ii,Re,1.0,Lx,Ly,Lz,nx,ny,nz,CFL,T,p0,p1)
%Run now the real cases
r = 8.0
for i=5:length(a1)
	p0 = a1(i)
	p1 = b2(i)
	disp(['Running case '  num2str(i)])
	disp(['alpha1 =  '  num2str(p0)])
	disp(['beta2  =  '  num2str(p1)])
	basename = [path0 'case_' int2str(i)]
	%no need of call init_istoturb3D again!
	copyfile([basename0 '.0.mat'],[basename '.0.mat']);
	isoturb3D(basename,ii,Re,r,Lx,Ly,Lz,nx,ny,nz,CFL,T,p0,p1)
end

