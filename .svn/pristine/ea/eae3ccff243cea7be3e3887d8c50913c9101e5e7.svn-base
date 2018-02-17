clear all
close all
%nfiles =30;
nfiles =30

% Omega Z 3D

% string1='mlpantano17_';
% files=1:nfiles;
% ind=1;
% fig=figure(1);
% for ii=files
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end 
%     filename=strcat(string1,string2,'.oz3d');
%     [time,x,y,temp]=readfieldyx(filename);
%     data(:,:,ind) = temp;
%     timev(ind)=time;
%     pcolor(x,y,temp),shading flat;colormap(redblue);
%     %title(num2str(floor(time)));
%     xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
%     %colorbar;
%     axis equal
%     axis tight
%     %colorbar;
%     %axis([-15 15 -5 5])
%     caxis([-0.2 0.2])
%     %pause(0.5)
%     M(ind)=getframe(fig)
%     %image(F(ind).cdata)
%     %colormap(F(ind).colormap)
%     clear temp
%     ind=ind+1;
% end
% save(string1,'y','timev','data')
% %fgure(2)
% fpsvalue=6;
% %movie(M,2,fpsvalue)
% movie2avi(M,'vorz3d_mlpantano17br.avi', 'compression', 'none','quality',100,'fps',6);
%%-----------------------------------------------------------------------%


% %VELOCIDAD Up YX+XZ
% %---------------------%
% string1='ml01loroll_';
% files=1:nfiles;
% ind=1;
% fig=figure(1);
% for ii=files
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end 
%     filename=strcat(string1,string2,'.upyx');
%     [time,x,y,temp]=readfieldyx(filename);
%     filename=strcat(string1,string2,'.upxz');
%     [time2,x2,posy,z2,temp2]=readfieldxz(filename);
%     timev(ind)=time;
%     
%     
%     subplot(2,1,1)
%     pcolor(x2,z2,squeeze(temp2(1,:,:))'),shading flat;colormap(redblue);
%     xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('z(transversal)','Fontsize',14,'Interpreter','Latex')
%     %colorbar;
%     %axis equal
%     xlim([min(x(:)) max(x(:))])
%     %axis tight
%     %colorbar;
%     %axis([-15 15 -5 5])
%     caxis([-0.5 0.5])
%     
%     
%     subplot(2,1,2)
% 
%     pcolor(x,y,temp),shading flat;colormap(redblue);
%     %title(num2str(floor(time)));
%     xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
%     %colorbar;
%     xlim([min(x(:)) max(x(:))])
%     %axis equal
%     %axis tight
%     %colorbar;
%     %axis([-15 15 -5 5])
%     ylim([-3 3])
%     xlim([min(x(:)) max(x(:))])
%     caxis([-0.1 0.1])
%     %pause(0.5)
%     M(ind)=getframe(fig)
%     %image(F(ind).cdata)
%     %colormap(F(ind).colormap)
%     clear temp
%     clear temp2
%     ind=ind+1;
% end
% %save(string1,'y','timev','data')
% %axis equal
% %fgure(2)
% fpsvalue=6;
% %movie(M,2,fpsvalue)
% movie2avi(M,'up_ml01lorollbr.avi', 'compression', 'none','quality',100,'fps',6);
% %------------------------------------------------------------------------%
% 
% 

% %Vorticidad Z (3D) XY + XZ
% string1='mlpantano17_';
% files=1:nfiles;
% ind=1;
% fig=figure(1);
% for ii=files
%     if ii<10
%         string2=strcat('00',num2str(ii));
%     elseif ii<100
%         string2=strcat('0',num2str(ii));       
%     else 
%         string2=num2str(ii);
%     end 
%     filename=strcat(string1,string2,'.upyx');
%     [time,x,y,temp]=readfieldyx(filename);
%     filename=strcat(string1,string2,'.upxz');
%     [time2,x2,posy,z2,temp2]=readfieldxz(filename);
% %     filename=strcat(string1,string2,'.ozyz');
% %     [time3,y3,z3,temp3]=readfieldyz(filename);
%      timev(ind)=time;
% 
%     
%     
%     g=subplot_tight(2,1,1,[0.1 0.1])
%     p = get(g,'position');
%     %p(4) = p(4)*1;  % Add 10 percent to height
%     %p(2) = p(2)+0.1*p(2);
%     %set(g, 'position', p);
%    
%     pcolor(x2,z2,squeeze(temp2(1,:,:))'),shading flat;colormap(redblue);
%     %xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('z(transversal)','Fontsize',14,'Interpreter','Latex')
%     %colorbar;
%     xlim([min(x(:)) max(x(:))])
%     ylim([min(z2(:)) max(z2(:))])
%     axis equal
%     %axis tight
%     %axis tight
%     %xlim([min(x(:)) max(x(:))])
%     %axis tight
%     %axis([-15 15 -5 5])
%     caxis([-0.5 0.5]); colorbar
%     
%     set(gca,'XTick',[]);
%   
%     %grid on
%     
%     g2=subplot_tight(2,1,2,[0.1 0.1])   
%     p2 = get(g2,'position');
%     %p2(4) = p2(4)*1/0.7;  % Add 10 percent to height
%     %p2(1)=p(1); p2(2)=p2(2)-0.1*p(2);
%     %set(g2, 'position', p2)
%     pcolor(x,y,temp),shading flat;colormap(redblue);
%     %title(num2str(floor(time)));
%     xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
%     xlim([min(x(:)) max(x(:))])
%     ylim([min(z2(:)) max(z2(:))])
%     %axis equal
%     %axis tight
%     %axis tight
%     %axis([-15 15 -5 5])
%     %ylim([-3 3])
%     %xlim([min(x(:)) max(x(:))])
%     caxis([-0.5 0.5])
%     axis equal
%     colorbar
%     %grid on
%     %axis tight
%     %pause(0.5)
%  
%     
% %     g3=subplot(2,2,2)   
% %     p3 = get(g3,'position');
% %     %p2(4) = p2(4)*1.5;  % Add 10 percent to height
% %     %set(g2, 'position', p2)
% %     pcolor(z3,y3,temp3'),shading flat;colormap(redblue);
% %     %title(num2str(floor(time)));
% %     xlabel('z(transversal)','Fontsize',14,'Interpreter','Latex')
% %     ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
% %     %colorbar;
% %     xlim([min(x(:)) max(x(:))])
% %     axis equal
% %     %axis tight
% %     %colorbar;
% %     %axis([-15 15 -5 5])
% %     ylim([-3 3])
% %     xlim([min(z3(:)) max(z3(:))])
% %     %caxis([-1 0])
% %     %pause(0.5)
%     
%     
%     
%     
%     M(ind)=getframe(fig)
%     %image(F(ind).cdata)
%     %colormap(F(ind).colormap)
%     clear temp
%     clear temp2
%     clear temp3
%     ind=ind+1;
% end
% %save(string1,'y','timev','data')
% %axis equal
% %fgure(2)
% fpsvalue=6;
% %movie(M,2,fpsvalue)
% %movie2avi(M,'oz_mlpantano17br.avi', 'compression', 'none','quality',100,'fps',6);
% movie2avi(M,'up_mlpantano17.avi', 'compression', 'none','quality',100,'fps',6);




%Vorticidad Z (3D) XY + XZ ML01LOROLL
string1='/data2/toni/work/mlBIG01_';
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
    filename=strcat(string1,string2,'.z3xz');
    [time2,x2,posy,z2,temp2]=readfieldxz(filename);
    filename=strcat(string1,string2,'.ozyz');
    [time3,y3,z3,temp3]=readfieldyz(filename);
    timev(ind)=time;

    
    
%    g=subplot_tight(2,1,1,[0.1 0.1])
%    p = get(g,'position');
    %p(4) = p(4)*1;  % Add 10 percent to height
    %p(2) = p(2)+0.1*p(2);
    %set(g, 'position', p);
   
    pcolor(x2,z2,squeeze(temp2(1,:,:))'),shading flat;colormap(redblue);
    %xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
    ylabel('z(transversal)','Fontsize',14,'Interpreter','Latex')
    %colorbar;
    axis equal
     xlim([min(x(:)) max(x(:))])
     ylim([-3 3])
    %axis equal
    %axis tight
    %axis tight
    %axis([-15 15 -5 5])
    caxis([-1 1])
    set(gca,'XTick',[]);
    colorbar
    
%    g2=subplot_tight(2,1,2,[0.1 0.1])   

 %   p2 = get(g2,'position');
    %p2(4) = p2(4)*1/0.7;  % Add 10 percent to height
    %p2(1)=p(1); p2(2)=p2(2)-0.1*p(2);
    %set(g2, 'position', p2)
    pcolor(x,y,temp),shading flat;colormap(redblue);
    %title(num2str(floor(time)));
    xlabel('x(longitudinal)','Fontsize',14,'Interpreter','Latex')
    ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
    xlim([min(x(:)) max(x(:))])
    axis equal
    %axis tight
    %axis tight
    %axis([-15 15 -5 5])
    ylim([-3 3])
    xlim([min(x(:)) max(x(:))])
    caxis([-1 1])
    colorbar
    %pause(0.5)
    
%     g3=subplot(2,2,2)   
%     p3 = get(g3,'position');
%     %p2(4) = p2(4)*1.5;  % Add 10 percent to height
%     %set(g2, 'position', p2)
%     pcolor(z3,y3,temp3'),shading flat;colormap(redblue);
%     %title(num2str(floor(time)));
%     xlabel('z(transversal)','Fontsize',14,'Interpreter','Latex')
%     ylabel('y(vertical)','Fontsize',14,'Interpreter','Latex')
%     %colorbar;
%     xlim([min(x(:)) max(x(:))])
%     axis equal
%     %axis tight
%     %colorbar;
%     %axis([-15 15 -5 5])
%     ylim([-3 3])
%     xlim([min(z3(:)) max(z3(:))])
%     %caxis([-1 0])
%     %pause(0.5)
    
    
    
    saveas(fig,sprintf('figure_%d.png',ind),'png'); 
%    M(ind)=getframe(fig)
    %image(F(ind).cdata)
    %colormap(F(ind).colormap)
    clear temp
    clear temp2
    clear temp3
    ind=ind+1;
end
%save(string1,'y','timev','data')
%axis equal
%fgure(2)
%fpsvalue=6;
%movie(M,2,fpsvalue)
%movie2avi(M,'oz_mlpantano17br.avi', 'compression', 'none','quality',100,'fps',6);
%movie2avi(M,'oz3d_ml01loroll.avi', 'compression', 'none','quality',100,'fps',7);


