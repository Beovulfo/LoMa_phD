function pintaens(path,num,ref,falpha,name)
h0 = figure(1)
set(gcf,'visible','off')

filename = [path,num2str(num,'%03d'),'.','o1xz']
[time,x,y,z,posy,o1]=readfieldxz(filename);
dm0 = 1  ;
x = x(x1:x2)/dm0;
y = y/dm0;
z = z(z1:z2)/dm0;
posy = posy/dm0;

o1=single(o1); 
jj = find(abs(posy)<15/dm0); 
posy=posy(jj); 
o1 = o1(jj,x1:x2,z1:z2); 
o1 = permute(o1,[3 2 1]);
[zz,xx,yy] = ndgrid(z,x,posy); 


filename = [path,num2str(num,'%03d'),'.','o2xz']
[time,x2,y2,z2,posy2,o2]=readfieldxz(filename);
time
o2=single(o2); 
%jj = find(abs(posy)<15); 
%posy=posy(jj); 
o2 = o2(jj,x1:x2,z1:z2); 
o2 = permute(o2,[3 2 1]);

filename = [path,num2str(num,'%03d'),'.','o3xz']
[time,x2,y2,z2,posy2,o3]=readfieldxz(filename);
o3=single(o3); 
%jj = find(abs(posy)<15); 
%posy=posy(jj); 
o3 = o3(jj,x1:x2,z1:z2); 
o3 = permute(o3,[3 2 1]);

%Calculate enstrophy
u = o1.^2+o2.^2+o3.^2;
uref = std(u(:)); %third dimension is y
uref = uref*(1+ ref) %ref = 2 means mean + 2 sigma
u = single(u);
[f2,v2] = isosurface(x,z,posy,u,uref);
h2 = patch('Faces',f2,'Vertices',v2);
isonormals(x,z,posy,u,h2);
isocolors(x,z,posy,u,h2);
set(h2,'facecolor','green','edgecolor','none','FaceAlpha',falpha);
view([0 20]);


set(gcf,'visible','on');
axis equal;
%axis([-172/dm0 172/dm0 -43/dm0 43/dm0 -25/dm0 25/dm0]);
%caxis([-10/dm0 20/dm0]);
set(gca,'Position',[.08 .08 .8 .9]);

figname = [num2str(num,'%03d'), 'ens',name,num2str(falpha*10,'%02d')]
print(h0,'-dpng',figname,'-r300')
colorbar;

