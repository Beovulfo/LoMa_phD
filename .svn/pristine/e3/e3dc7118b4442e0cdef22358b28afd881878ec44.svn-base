load([matname '.mat']);
x = x(x1:x2)/dm0;
z = z(z1:z2)/dm0;
%posy is already broken
[zz,xx,yy] = ndgrid(z,x,posy); 
clear xx,zz; 

u = u(:,x1:x2,z1:z2); 
u = permute(u,[3 2 1]);

o1 = o1(:,x1:x2,z1:z2); 
o1 = permute(o1,[3 2 1]);
%jj = find(abs(posy)<15); 
%posy=posy(jj); 
o2 = o2(:,x1:x2,z1:z2); 
o2 = permute(o2,[3 2 1]);
%jj = find(abs(posy)<15); 
%posy=posy(jj); 
o3 = o3(:,x1:x2,z1:z2); 
o3 = permute(o3,[3 2 1]);

save([matname '_cut.mat'],'x','z','y','yy','posy','u','o1','o2','o3','uref','time');
display('Cut fields saved!')

