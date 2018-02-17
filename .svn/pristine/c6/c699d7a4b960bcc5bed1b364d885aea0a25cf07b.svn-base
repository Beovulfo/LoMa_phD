%a05 = load('caca/caca.stats'); 
%[ksp05 tsp05 sp05] = readsp('caca/caca.spe_u');
 
%a05 = load('results/s40CFL025pi96.stats'); 
a05 = load('results/s80pi96Re500.stats'); 
%[ksp05 tsp05 sp05] = readsp('results/s40CFL025pi96.spe_u'); 
[ksp05 tsp05 sp05] = readsp('results/s80pi96Re500.spe_u'); 

%a10 = load('r1_48c_Re1000/r1_48c_Re1000.stats'); 
%[ksp10 tsp10 sp10] = readsp('r1_48c_Re1000/r1_48c_Re1000.spe_u'); 
%a10 = load('caca/caca.stats'); 
%[ksp10 tsp10 sp10] = readsp('caca/caca.spe_u'); 

%a10 = load('results/s80pi96Re500.stats'); 
%[ksp10 tsp10 sp10] = readsp('results/s80pi96Re500.spe_u'); 
%a20 = load('results/s80pi96Re500.stats'); 
%a20 = load('results/s40CFL025pi96.stats'); 
%[ksp20 tsp20 sp20] = readsp('results/s80pi96Re500.spe_u'); 
%[ksp20 tsp20 sp20] = readsp('results/s40CFL025pi96.spe_u'); 


figure(1); % compare disipation to kinetic energy evolution: 

plot(0.5*(a05(1:end-1,1)+a05(2:end,1)),diff(0.5*sum(a05(:,6:8),2))./diff(a05(:,1)),'.b', a05(:,1),-a05(:,10),'-b'); hold on; 

%plot(0.5*(a10(1:end-1,1)+a10(2:end,1)),diff(0.5*sum(a10(:,6:8),2))./diff(a10(:,1)),'.r', a10(:,1),-a10(:,10),'-r'); hold on

%plot(0.5*(a20(1:end-1,1)+a20(2:end,1)),diff(0.5*sum(a20(:,6:8),2))./diff(a20(:,1)),'.k', a20(:,1),-a20(:,10),'-k'); hold off

figure(2); clf; 

ii05 = find(tsp05>tref(1) & tsp05 < tref(2)); 
for i=ii05'; 
    loglog(ksp05*a05(i,12),sp05(i,:)./a05(i,10).^(2/3),'k'); 
    hold on
    %pause
end     

%ii10 = find(tsp10>tref(1) & tsp10 < tref(2)); 
%for i=ii10'; 
%    loglog(ksp10*a10(i,12),sp10(i,:)./a10(i,10).^(2/3),'r'); 
%    hold on
    %pause
%end     

%ii20 = find(tsp20>tref(1) & tsp20 < tref(2)); 
%for i=ii20'; 
%    loglog(ksp20*a20(i,12),sp20(i,:)./a20(i,10).^(2/3),'k'); 
%    hold on
%    %pause
%end     

plot(ksp05*a05(i,12),(ksp05*a05(i,12)).^(-5/3),'--k'); 
%plot(ksp10*a10(i,12),(ksp10*a10(i,12)).^(-5/3),'--k'); 
hold off
