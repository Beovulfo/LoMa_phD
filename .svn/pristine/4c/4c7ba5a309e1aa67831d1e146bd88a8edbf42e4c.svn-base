
function []=plotframes(filename,ibeg,iend,nx,nz)

figure(1)
hold on
figure(2)
hold on
for i=ibeg:iend
	[y,scalar,mode]=read_sca_out(filename,i,nx,nz)
        figure(1)
        plot(y,mode)
        figure(2)
        plot(y,scalar,'r')
        pause(1)
end

