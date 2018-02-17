file11='cfttest_single.txt';        
file12='rffttest_single.txt';        
file21='cfttest_double.txt';        
file22='rffttest_double.txt';        

fileID = fopen(file11);        
C = textscan(fileID,'%f %f %f');
fclose(fileID);
phys_4=C{:,1};refou_4=C{:,2};imfou_4=C{:,3};

fileID = fopen(file12);        
R = textscan(fileID,'%f %f');
fclose(fileID);
phys_r4=R{:,1};refou_r4=R{:,2};

fileID = fopen(file21);        
C = textscan(fileID,'%f %f %f');
fclose(fileID);
phys_8=C{:,1};refou_8=C{:,2};imfou_8=C{:,3};

fileID = fopen(file22);        
R = textscan(fileID,'%f %f');
fclose(fileID);
phys_r8=R{:,1};refou_r8=R{:,2};

sizez=size(phys_4)
sizex=size(phys_r4)
mgalz=sizez(1)
mgalx=sizex(1)
vx=linspace(-pi,pi,mgalx)
vz=linspace(-pi,pi,mgalz)

figure(1)
plot(vz,abs(refou_4-refou_8))
plot(vz,abs(imfou_4-imfou_8))

figure(2)
plot(vx,abs(refou_r4-refou_r8))

