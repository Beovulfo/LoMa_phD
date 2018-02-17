%aalmagro (Mar 2015)
%---------------------------------------------------------%
% script to create nice isosurfaces from image files
% from tofis
%
%PROCEDURE:
%  1)first run with lee=1, cut =0 => we read time
%  2)check dm on time read using stats,modify dm variable.
%  3)run with lee=0, cut=1 in order to cut data
%  4)run with lee=0, cut=0 in order to make plots
%  Make sure to change variables:
%   name, path, num, ymax,uref,var, dm
%lee=1;
%cut=1;
lee=1;
cut=0;
%lee=0;cut=0;
%name='s10';
%name='s10_011';
%path='/data/toni/domain3s10Re160_';
%path='/data/toni/Re160s40my1101e_';
%path='/data2/toni/turb/s40SSmov_';
path='/share/drive/toni/Re160s80/case1/SS/s80aSS_';
%path='/data/toni/domain3s10Re160_';
%where to save or read .mat (matlab space)
for num=1:16
        name=['s80aSS_' num2str(num,'%03d') ];
        matname =['/share/drive/toni/Re160s80/case1/SS/field_' name]
%path='/data/toni/domain3s10Re160_';
	%num=1;
	var='upxz';
	uref=0.0;
	ymax=30.0;
	savefields(path,num,var,uref,ymax,matname);
	%CUT
	z1=1;z2=288;x1=1;x2=768;
	dm0=1.0
	cutfields
        %if space not available...
        load([ matname '_cut.mat']);
	sigma = 4.0
	falpha = 0.3
        %momentum thickness read from stats
       	%dm = 3.8
       	%dm = 6.1
       	dm = 1.0
	x = x./dm;
	z = z./dm;
	posy = posy./dm;
	yy = yy./dm;
	ylim1= -30.0 %already scaled
	ylim2= 30.0  %already scaled
	pintafields
	savefig(h0,name)
end


