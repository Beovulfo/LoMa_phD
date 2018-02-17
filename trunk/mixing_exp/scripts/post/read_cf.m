%Reads the data form the .cf file
function [time,ener,reynota]=read_cf(filename)
%clear all
%close all
%clc
%PONER EL NOMBRE DEL ARCHIVO A LEER
dataCF=dlmread(filename)
%Info of the CF file:
%time,-1.*Wz0,WzL,sqrt(reynota),ener,massu,massw,MPI_WTIME()+iter_time,commtimer-ctm2
%ener: 
%ener(1)=sqrt(u1^2)
%ener(2)=sqrt(u2^2)
%ener(3)=sqrt(u3^2)
%ener(4)=sqrt(u1*u2)
%ener(5)=sqsrt(u1*u3)
%ener(6)=sqsrt(u2*u3)
%ener(7,8,9)=vor1,vor2,vor3
%
ener=dataCF(:,5:13);
reynota=dataCF(:,4);
time=dataCF(:,1);

%figure(1)
%plot(time,ener(2))
%lambda=0.004151284675632;%Real de ee(1) primer autovalor                                       
%semilogy(time,ener(:,2),time,ener(1,2)*exp(3*lambda/2*time),'o')
%figure(1)
%plot(time,ener(:,2),time,ener(1,2)*exp(3*lambda/2*time),'o') 
%plot(time,ener(:,2),'o') 
%title('u2 energy')
%xlabel('time')
%ylabel('sqrt(energy v)')


