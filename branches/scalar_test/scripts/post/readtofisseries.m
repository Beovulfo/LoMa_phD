clear all
close all

%nfiles =30;
nfiles =60;


%string1='ml3dloroll_';
%string1='ml02hiroll_';
string1='mlpantano17_';
variable='ozyx';

files=1:nfiles;


%string1='mlEXPre500_01_';
ind=1;
for ii=files
    if ii<10
        string2=strcat('00',num2str(ii));
    elseif ii<100
        string2=strcat('0',num2str(ii));       
    else 
        string2=num2str(ii);
    end
    filename=strcat(string1,string2,'.',variable);
    [time,y,z,temp]=readfieldyz(filename);
%write in hdf5 for variable "variable"
    if (ii==1) 
      hdf5name=strcat(string1,'.h5');
      group=strcat('/',variable);
      dataset=strcat(group,'/',string2)
      if (exist(hdf5name)==0) %if hdf5 file do not exist...
         hdf5write(hdf5name,dataset,temp);
      else 
         hdf5write(hdf5name,dataset,temp,'WriteMode','append');
      end

    else
    dataset=strcat(group,'/',string2)
    hdf5write(hdf5name,dataset,temp,'WriteMode','append');
    end
    %if ii == 1;
    %    my =size(temp,1);
    %    nvar = size(temp,2);
    %    data=zeros(my,nvar,nfiles);
    %end
    data(:,ind) = temp(:);
    [my,mgalz]=size(temp);

    meandata(:,ind)=mean(temp,1);
    mindata(ind)=min(meandata(:,ind));
    %mindata(ind)=min(temp(:))
    timev(ind)=time;
    ind=ind+1;
    clear temp
end

 %[t,ener,dm]=read_cf('ml3dloroll_01.cf')
% [t,ener,dm]=read_cf('ml02hiroll_01.cf');;
 [t,ener,dm]=read_cf('mlpantano17_01.cf');;
 
 save(string1,'y','timev','data','mindata')
 %plot(timev,-mindata'./dm(1),'r*')
 plot(timev,-mindata','g.')
 title('Evoluci�n temporal del m�nimo de vorticidad z en el mid-braid')
 xlabel('t')
 ylabel('-\omega_b')
 hold on
% wbdig=load('wbdig.csv');
% plot(wbdig(:,1),wbdig(:,2),'r*')
