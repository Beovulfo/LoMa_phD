%AALMAGRO
%Script for matlab,
%Using fields obtained by "tofis" this script plots some relevant variables from it.
%They will be comparable with Rogers&Moser.
close all
clear all

iplane=1;
%fileroot='IC3D';
%fileroot='ICOBIN';
%fileroot='ICml08turb';
%fileroot='ICml07turb'
%fileroot='ICml08turb';
%fileroot='ICml10turb';
fileroot='mlpantano17_001';
%fileroot='sample';
%fileroot='mlpantano17_035';
%fileroot='mlpantano01_040';

%fileroot='ml05turb_050'
%fileroot='ICml06turb';
%fileroot='ml08turb_028'
%filename='IC2D.ozyx';
%filename='IC2ndmode.ozyx';
%filename='IC3D.vpyx';
%fileroot='ml01loroll_014';
%scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/1.5])
tamfuente=14;

filename=strcat(fileroot,'.ozyx');
[time,x,y,wk1r]=readfieldyx(filename);
%subplot(2,3,1);
wz=figure;
pcolor(x,y,wk1r); shading flat
%contour(x,y,wk1r,10);
colormap(redblue);
xlabel('x(longitudinal)','FontSize',tamfuente)
ylabel('y(vertical)','FontSize',tamfuente)
axis equal
axis tight
%ylim([-3 3])
%ylim([5 10])
%caxis([-0.1 0.1])
titlefig=['Vor z'];colorbar;
filename=[fileroot 'vorz_plyx'];
title(titlefig,'FontSize',tamfuente)
set(gca,'fontsize',tamfuente);
print(wz,'-depsc',filename);
%!mv *.eps /data/toni/mixing/work/
%title(titlefig)

filename=strcat(fileroot,'.upyx');
[time,x,y,wk1r]=readfieldyx(filename);
%subplot(2,3,1);
wz=figure;
pcolor(x,y,wk1r); shading flat
%contour(x,y,wk1r,10);
colormap('redblue');
xlabel('x(longitudinal)','FontSize',tamfuente)
ylabel('y(vertical)','FontSize',tamfuente)
axis equal
axis tight
%ylim([-3 3])
%ylim([5 10])
%caxis([-0.5 0.5])
titlefig=['Vel U'];colorbar;
filename=[fileroot 'u_plyx'];
title(titlefig,'FontSize',tamfuente,'Interpreter','tex')
set(gca,'fontsize',tamfuente);
print(wz,'-depsc',filename);
%!mv *.eps /data/toni/mixing/work/



filename=strcat(fileroot,'.oxyx');
[time,x,y,wk1r]=readfieldyx(filename);
%subplot(2,3,2);
wz2=figure;
pcolor(x,y,wk1r); shading flat
%contour(x,y,wk1r,15);
%caxis([-1 1])
xlabel('x(longitudinal)','FontSize',tamfuente)
ylabel('y(vertical)','FontSize',tamfuente)
colormap(redblue);
axis equal
axis tight
%ylim([-3 3])
%colorbar
%caxis([-0.5 0.5])
%ylim([-2.5 2.5])
titlefig=['Vor x'];colorbar;
filename=[fileroot 'vorx_plyx'];
title(titlefig,'FontSize',tamfuente)
set(gca,'fontsize',tamfuente);
print(wz2,'-depsc',filename);
%!mv *.eps /data/toni/mixing/work/
%title(titlefig,'FontSize',15)


% filename=strcat(fileroot,'.oxyz');
% %filename='IC2ndmode.oxyz';
% %filename='IC3D.vpyz';
% [time,y3,z3,wk3r]=readfieldyz(filename);
% %subplot(2,3,3);
% wz3=figure;
% %pcolor(z3,y3,wk3r'); shading flat
% contour(z3,y3,wk3r',10);
% xlabel('z')
% ylabel('y')
% colorbar
% %ylim([-5 5]);
% %caxis([-1 1])
% titlefig=['Vorticidad x en plano yz, t=0'];
% filename=['vorx_plyzMP_t10.eps'];
% %print(wz3,'-depsc',filename);
% title(titlefig)
% 
% 
% %------------------------------------------------%
% figure(101)
% 
filename=strcat(fileroot,'.wpxz');
%filename='IC2ndmode.oxxz';
%filename='IC3D.vpxz';
[time,x2,posy,z2,wk2r]=readfieldxz(filename);
%subplot(2,3,4);
wxdata=wk2r(iplane,:,:);
wz3=figure;
pcolor(x2,z2,squeeze(wxdata)'); shading flat
posy(end)
xlabel('x(longitudinal)','FontSize',tamfuente)
ylabel('z(trans.)','FontSize',tamfuente)
colormap(redblue);
axis equal
axis tight
filename=[fileroot 'velw_plxz'];
colorbar
%ylim([-3 3])
%caxis([-0.2 0.2])
titlefig=['Vel W'];
title(titlefig,'FontSize',tamfuente)
set(gca,'fontsize',tamfuente);
print(wz3,'-depsc',filename);

%titlefig=['vor z en plano xz'];
% 
% title(titlefig)
% 
% 
% figure(102)
% filename=strcat(fileroot,'.wpxz');
% %filename='IC2ndmode.oxxz';
% %filename='IC3D.vpxz';
% [time,x2,posy,z2,wk2r]=readfieldxz(filename);
% %subplot(2,3,5);
% 
% wxdata=wk2r(iplane,:,:);
% pcolor(x2,z2,squeeze(wxdata)'); shading flat
% posy(end)
% xlabel('x')
% ylabel('z')
% colorbar
% %caxis([-1 1])
% titlefig=['Vel w en plano xz, t=0'];
% %filename=['vorx_plxz_IC3D.eps'];
% %print(wx,'-depsc',filename);
% title(titlefig)
% 
% 
% figure(103)
% filename=strcat(fileroot,'.upxz');
% %filename='IC2ndmode.oxxz';
% %filename='IC3D.vpxz';
% [time,x2,posy,z2,wk2r]=readfieldxz(filename);
% %subplot(2,3,6);
% wxdata=wk2r(iplane,:,:);
% pcolor(x2,z2,squeeze(wxdata)'); shading flat
% posy(end)
% xlabel('x')
% ylabel('z')
% %caxis([-1 1])
% colorbar
% titlefig=['Vel u en plano xz, t=0'];
% %filename=['vorx_plxz_IC3D.eps'];
% %print(wx,'-depsc',filename);
% title(titlefig)






