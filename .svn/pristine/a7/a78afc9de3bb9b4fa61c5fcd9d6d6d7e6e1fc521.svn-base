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
lee=0;
cut=0;
%lee=0;cut=0;
%name='s10';
name='s10_011';
%name='s40';
%path='/data/toni/domain3s10Re160_';
%path='/data/toni/Re160s40my1101e_';
path='/data/toni/domain3s10Re160_';
%where to save or read .mat (matlab space)
matname =['/data/toni/field_' name]
if lee==1 
%path='/data/toni/domain3s10Re160_';
	num=11;
	var='upxz';
	uref=0.0;
	ymax=20.0;
	savefields(path,num,var,uref,ymax,matname);
else
%to cut fields
	if cut ==1
		z1=1;z2=576;x1=1;x2=768;
		dm0=1.0
		cutfields
	else
                 %if space not available...
          	load([ matname '_cut.mat']);
		sigma = 4.0
		falpha = 0.3
        %momentum thickness read from stats
        	%dm = 3.8
        	%dm = 6.1
        	dm = 5.5
		x = x./dm;
		z = z./dm;
		posy = posy./dm;
		yy = yy./dm;
		ylim1= -3.0 %already scaled
		ylim2= 3.0  %already scaled

		pintafields
		savefig(h0,name)
	end
end


