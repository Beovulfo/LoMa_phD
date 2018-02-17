clear all
close all
%nfiles =30;
nfiles =60

%string1='ml3dloroll_';
%string1='ml09turb_';
%string1='mlpantano13_';
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
    filename=strcat(string1,string2,'.ozyx');
    [time,x,y,temp]=readfieldyx(filename);
    data(:,:,ind) = temp;
    timev(ind)=time;
    pcolor(x,y,temp),shading flat;colormap(bone);
    title(num2str(floor(time)));
    xlabel('x','Fontsize',14,'Interpreter','Latex')
    ylabel('y','Fontsize',14,'Interpreter','Latex')
    axis equal
    axis tight
    %colorbar;
    %axis([-15 15 -5 5])
    caxis([-0.2 0.2])
    pause(0.5)
    M(ind)=getframe(fig)
    %image(F(ind).cdata)
    %colormap(F(ind).colormap)
    clear temp
    ind=ind+1;
end
save(string1,'y','timev','data')
%fgure(2)
fpsvalue=6;
%movie(M,2,fpsvalue)
movie2avi(M,'vorz_mlpantano17.avi', 'compression', 'none','quality',100,'fps',6);


%choose one time
nstep=30
timev(nstep)
wz=figure(3)
pcolor(x,y,data(:,:,nstep));shading flat; colorbar;
colormap('bone');
axis equal
axis tight
%contour(x,y,data(:,:,nstep),50);
%axis([-3.5 3.5 -7 7])
caxis([-0.1 0.1])
xlabel('x');ylabel('y');
filename=['vorz_plyx_t' num2str(timev(nstep)) '.eps'];
title (num2str(timev(nstep)))
print(wz,'-depsc',filename);


% % %choose one time
% % nstep=8
% % timev(nstep)
% % figure(3)
% % pcolor(x,y,data(:,:,nstep));shading flat; colorbar;
% % %colormap('bone');
% % %contour(x,y,data(:,:,nstep),50);
% % %axis([-3.5 3.5 -7 7])
% % title (num2str(timev(nstep)))
% % 
% % %choose whan time
% % nstep=20
% % timev(nstep)
% % figure(4)
% % pcolor(x,y,data(:,:,nstep));shading flat;colorbar;
% % %colormap('bone');
% % %contour(x,y,data(:,:,nstep),50);
% % %axis([-3.5 3.5 -7 7])
% % title (num2str(timev(nstep)))


% fig=figure(100);
% ind=1
% for ii=files
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end 
%     filename=strcat(string1,string2,'.oxyz');
%     [time,y,z,temp]=readfieldyz(filename);
%     %data2(:,:,ind) = temp;
%     %timev(ind)=time;
%     %contour(x,y,temp,20)
%     
%  
%     pcolor(z,y,temp'),shading flat;colormap('jet');
%     title(num2str(time));
%     colorbar;
%     %axis([-15 15 -5 5])
%     caxis([-0.5 0.5])
%     pause(0.5)
%     M(ind)=getframe
%     %image(F(ind).cdata)
%     %colormap(F(ind).colormap)
%     clear temp
%     ind=ind+1;
% end


