clear all
close all
nfiles =60;
string1='mlpantano17_';
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
    filename=strcat(string1,string2,'.oz3d');
    [time,x,y,temp]=readfieldyx(filename);
    data(:,:,ind) = temp;
    timev(ind)=time;
    %contour(x,y,temp,50)
%    filename=strcat(string1,string2,'.oz3d');
%    [time,x,y,temp2]=readfieldyx(filename);
%    data2(:,:,ind) = temp2;
%    timev(ind)=time;
    

    title (num2str(time))
    pcolor(x,y,temp),shading flat;
    %fig=contour(x,y,temp,20);colorbar;
    %axis([-3.5 3.5 -3 3])
    caxis([-1 1]); 
    M(ind)=getframe
    %image(F(ind).cdata)
    %colormap(F(ind).colormap)
    clear temp
    ind=ind+1;
end
save(string1,'y','timev','data')
figure(2)
movie(M,1)
movie2avi(M,'mlpantano17_vorz.avi', 'compression', 'none');


%choose whan time
nstep=7
timev(nstep)
figure(3)
%pcolor(x,y,data(:,:,nstep));shading flat; colorbar;
contour(x,y,data(:,:,nstep),10); colorbar;
caxis([-1 1]); 
%colormap('bone');
%contour(x,y,data(:,:,nstep),50);
%axis([-3.5 3.5 -3 3])
title (num2str(timev(nstep)))

%choose one time
nstep=7
timev(nstep)
figure(4)
pcolor(x,y,data(:,:,nstep));shading flat;colorbar;
%contour(x,y,data(:,:,nstep),10);colorbar;
colormap('bone');
%contour(x,y,data(:,:,nstep),50);
%axis([-3.5 3.5 -3 3])
caxis([-10 10]); 
title (num2str(timev(nstep)))




