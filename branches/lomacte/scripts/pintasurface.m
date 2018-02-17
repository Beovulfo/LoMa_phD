function pintasurface(lee,path,num,var,uref)
%opengl software
if lee==1
  %filename = 'turb/3DincRe160_045.upxz'; 
   filename = [path,num2str(num,'%03d'),'.',var]
  [time,x,y,z,posy,u]=readfieldxz(filename);
  keyboard
  dm0 = 1  
  z1=1;z2=256;x1=1;x2=256;
  x = x(x1:x2)/dm0;
  y = y/dm0;
  z = z(z1:z2)/dm0;
  %z = z/dm0;
  posy = posy/dm0;

  u=single(u); 
  jj = find(abs(posy)<15/dm0); 
  posy=posy(jj); 
  u = u(jj,x1:x2,z1:z2); 
  u = permute(u,[3 2 1]);
  [zz,xx,yy] = ndgrid(z,x,posy); 
  clear xx zz
end

us = smooth3(u); 

figure(1)
set(gcf,'visible','off')
clf

[f,v] = isosurface(x,z,posy,us,uref);
h1 = patch('Faces',f,'Vertices',v)
isonormals(x,z,posy,us,h1);
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

set(gca,'projection','perspective');
set(gca,'fontname','times','fontsize',16)
set(gcf,'position',[ 194   429      1649      543])   
set(gca,'Position',[.08 .08 .8 .9]) 
grid on

set(gcf,'visible','on')
%figname = [num2str(num,'%03d'), var]
%print('-dpng',figname,'-r300')


%xlabel('x/\delta_m^0','fontname','times','fontsize',16)  
%ylabel('y/\delta_m^0','fontname','times','fontsize',16)
%zlabel('z/\delta_m^0','fontname','times','fontsize',16)
