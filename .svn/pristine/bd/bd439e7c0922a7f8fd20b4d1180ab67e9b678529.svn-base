%clear all
%close all

%File 01 limits and name
%string0='/data2/toni/mldiffSCAL01';
%string0='/data2/toni/mlturbSCAL02';
string0='/data2/toni/mlRMhiroll01';
 
%Re=160;
Re=500;
nu=1/Re;
ibeg=1;
iend =60;

string1=strcat(string0,'_01_');

for ii=1:iend
    ind=ii;
    if ii<10
        string2=strcat('00',num2str(ii));
    elseif ii<100
        string2=strcat('0',num2str(ii));       
    else 
        string2=num2str(ii);
    end
    filename=strcat(string1,string2,'.sta');
    [y,time,temp]=lee_sta(filename);
    if ind == 1;
        my =size(temp,1);
        nvar = size(temp,2);
        data=zeros(my,nvar,iend-ibeg);
    end
    data(:,:,ind) =temp;
    timev(ind)=time;
    clear temp
end
save('y','timev','data')

figure(2)
hold on
for i=ibeg:iend
%plot u00 evolution
   plot(y,data(:,1,i))
   pause(0.1)
end
