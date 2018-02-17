close all
clear all
itmax=100;
implicit=load('mlIMP_01_.mat');
explicit=load('mlEXP_01_.mat');
y=implicit.y;
time_imp=implicit.timev;
time_exp=explicit.timev;
for (i=1:itmax)
 dm_imp(i)=0.25*trapz(y,1-implicit.data(:,1,i).^2)
 dm_exp(i)=0.25*trapz(y,1-explicit.data(:,1,i).^2)
end
figure(1)
hold on
plot(time_imp,dm_imp,'b')
plot(time_exp,dm_exp,'r')



