function savefields(path,num,var,uref,ymax,matname)
%opengl software
  %filename = 'turb/3DincRe160_045.upxz'; 
   filename = [path,num2str(num,'%03d'),'.',var]
  [time,x,y,z,posy,u]=readfieldxz(filename);
  dm0 = 1  
  y = y/dm0;
  %ymax = 20;
 %z = z/dm0;
  posy = posy/dm0;
  time

  jj = find(abs(posy)<ymax/dm0); 
  posy=posy(jj); 
  u=single(u(jj,:,:)); 
  [zz,xx,yy] = ndgrid(z,x,posy); 

  save([matname '.mat'],'x','z','y','yy','posy','u','uref','time');
  clear xx zz u
%Other fields
  filename = [path,num2str(num,'%03d'),'.','o1xz']
  [time,xv,yv,zv,posy2,o1]=readfieldxz(filename);
  o1=single(o1(jj,:,:));
  save([matname '.mat'],'o1','-append');
  clear xv yv zv o1

 
  filename = [path,num2str(num,'%03d'),'.','o2xz']
  [time,xv,yv,zv,posy2,o2]=readfieldxz(filename);
  o2=single(o2(jj,:,:)); 
  save([matname '.mat'],'o2','-append');
  clear xv yv zv o2

   filename = [path,num2str(num,'%03d'),'.','o3xz']
  [time,xv,yv,zv,posy2,o3]=readfieldxz(filename);
  o3=single(o3(jj,:,:)); 
  save([matname '.mat'],'o3','-append');
  clear xv yv zv o3



