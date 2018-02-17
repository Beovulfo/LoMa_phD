nfiles =102;
string1='mlBCNeuman_';
for ii=1:nfiles
    if ii<10
        string2=strcat('00',num2str(ii));
    elseif ii<100
        string2=strcat('0',num2str(ii));       
    else 
        string2=num2str(ii);
    end
    filename=strcat(string1,string2,'.sta');
    [y,temp]=lee_sta(filename);
    if ii == 1;
        my =size(temp,1);
        nvar = size(temp,2);
        data=zeros(my,nvar,nfiles);
    end
    data(:,:,ii) =temp;
    clear temp
end
