clear all
close all
nfiles =24;
string1='ml';
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
    filename=strcat(string1,string2,'.upyx');
    [time,x,y,temp]=readfieldyx(filename);
    data(:,:,ind) = temp;
    timev(ind)=time;
    pcolor(x,y,temp),shading flat;
    pause(0.01)
    M(ind)=getframe
    %image(F(ind).cdata)
    %colormap(F(ind).colormap)
    clear temp
    ind=ind+1;
end
save(string1,'y','timev','data')
figure(2)
movie(M,1)
close all
%movie2avi(M,'vorz_evolution.avi', 'compression', 'none');
