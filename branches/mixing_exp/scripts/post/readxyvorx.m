%clear all
close all
nfiles =30;
string1='ml3dloroll_';
files=1:nfiles;
ind=1;
fig=figure(1);
for ii=files
    if ii<10
        string2=strcat('00',num2str(ii));
    elseif ii<100
        string2=strcat('0',num2str(ii));       
    else 
        string2=num2str(ii);
    end
    filename=strcat(string1,string2,'.oxyx');
    [time,x,y,temp]=readfieldyx(filename);
    data(:,:,ind) = temp;
    timev(ind)=time;
    %contour(x,y,temp,50)
    
    pcolor(x,y,temp),shading flat;
    title (num2str(time))
    colorbar;
    axis([-3.5 3.5 -3 3])
    pause(0.1)
    M(ind)=getframe
    %image(F(ind).cdata)
    %colormap(F(ind).colormap)
    clear temp
    ind=ind+1;
end
save(string1,'y','timev','data')
figure(2)
movie(M,1)
movie2avi(M,'vorx_evo.avi', 'compression', 'none');


%choose whan time
nstep=14
timev(nstep)
figure(3)
pcolor(x,y,data(:,:,nstep));shading flat; colorbar;
%colormap('bone');
%contour(x,y,data(:,:,nstep),50);
axis([-3.5 3.5 -3 3])
title (num2str(timev(nstep)))
%choose whan time
nstep=20
timev(nstep)
figure(4)
pcolor(x,y,data(:,:,nstep));shading flat;colorbar;
%colormap('bone');
%contour(x,y,data(:,:,nstep),50);
axis([-3.5 3.5 -3 3])
title (num2str(timev(nstep)))




