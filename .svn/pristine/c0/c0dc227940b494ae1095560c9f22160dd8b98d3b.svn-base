%opengl software
%name of the fig:
%name ='s10'
%sigma = 4.0
%falpha = 0.3

h0=figure(1)
set(gcf,'visible','off')
clf
[f,v] = isosurface(x,z,posy,u,uref);
h1 = patch('Faces',f,'Vertices',v)
isonormals(x,z,posy,u,h1);
isocolors(x,z,posy,yy,h1);
set(h1,'facecolor','interp','edgecolor','none','FaceAlpha',1);
%colormap jet
colormap(redblue)
%set(h1,'facecolor','interp','edgecolor','none');
axis equal 
%axis([x(1) x(end) posy(1) posy(end) z(1) z(end)])
%caxis([-10/dm0 10/dm0]);
view([0 20])
cl1 = camlight;
lighting gouraud
xlabel('x/\delta_m^0','fontname','times','fontsize',16)  
ylabel('y/\delta_m^0','fontname','times','fontsize',16)
zlabel('z/\delta_m^0','fontname','times','fontsize',16)

set(gca,'projection','perspective');
set(gca,'fontname','times','fontsize',16)
set(gcf,'position',[ 194   429      1649      543])   
set(gca,'Position',[.08 .08 .8 .9]) 
grid on
set(gcf,'visible','on')
%figname = [num2str(num,'%03d'), var]
%print('-dpng',figname,'-r300')



ens = o1.^2+o2.^2+o3.^2;
ensref = std(ens(:)); %third dimension is y
ensref = ensref*(1+ sigma) %ref = 2 means mean + 2 sigma
ens = single(ens);
[f2,v2] = isosurface(x,z,posy,ens,ensref);
h2 = patch('Faces',f2,'Vertices',v2);
isonormals(x,z,posy,u,h2);
isocolors(x,z,posy,u,h2);
set(h2,'facecolor','green','edgecolor','none','FaceAlpha',falpha);
view([0 20]);


set(gcf,'visible','on');
axis equal;
axis([x(1) x(end) z(1) z(end) ylim1 ylim2]);
caxis([ylim1 ylim2]);
set(gca,'Position',[.08 .08 .8 .9]);

figname = [num2str(num,'%03d'), 'iso_u']
print(h0,'-dpng',figname,'-r300')

