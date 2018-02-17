%===================READING FROM FILE==========================%
%
my=513 %mlpantano13y2
ymax=10; %half heigth of the box
y=create_mesh1D(ymax,'manl',1 ,my); 
%mgalx=768;
mgalx=192;
%mgalz=288;
mgalz=192;
mx=mgalx*2/3;
mz=mgalz*2/3-1

%
filename0='/data2/toni/mlRMhiroll01sc001.';
ii=20

  if ii<10
     string2=strcat('00',num2str(ii));
  elseif ii<100
     string2=strcat('0',num2str(ii));
  else
     string2=num2str(ii);
  end
  filename=strcat(filename0,string2)

  finp=fopen(filename,'r','l');
nplanes=mx/2;
ntotr=2*my*mz;

 for (i=1:nplanes)
     dummy3=fread(finp,1,'int');
     wk1(i,:)=fread(finp,ntotr,'float');
     dummy4=fread(finp,1,'int');
 end

 fclose(finp)

 scalar=zeros(1,my);
 for (j=1:my)
    scalar(j)=wk1(1,2*j-1,1);
 end

%start transformation
scal_phys=zeros(mgalx,my,mgalz);
%1)Sum real part and imaginary part in order to create the complex number
scal_complex=wk1(:,2*(1:my)-1,:)+sqrt(-1)*wk1(:,2*(1:my),:);
%2)For each XZ plane we are going to Transfrom from FOU2FIS
for j=1:my
   tempFF=squeeze(scal_complex(:,j,:));
   %2.1) transform ifft for Z
   tempFP=ifft(tempFF,mgalz,2);
   %2.2) transform ifft for x
   tempPP=ifft(tempFP,mgalx,1,'symmetric');
   scal_phys(:,j,:)=tempPP;
end

figure(1)
plot(y,scalar,'b-');
