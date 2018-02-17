%This program loads some simulation results in order to compare them
%Use .cf files directly

fileroot='/data2/toni/';

fname='inc02cte_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,4),'r-')
listlabel=[{'inc02cte'}];

fname='inc02nocte_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,4),'r*')
listlabel=[listlabel,{'inc02nocte'}];

fname='inc02ctenewdy_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,4),'ro')
listlabel=[listlabel,  {'inc02ctenewdy'}];
%listlabel=strjoin({listlabel, 'inc02ctenewdy'}, ',');

fname='inc02noctenewdy_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,4),'rs')
listlabel=[listlabel, {'inc02noctenewdy'}];

%LOMA RESULTS IN BLUE

fname='inc03cte_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,3),'b-')
listlabel=[listlabel, {'inc03cte'}];

%This one is not working
%----------------------------
%fname='inc03nocte_01.cf';
%filename=[fileroot,fname];
%cf=load(filename);
%figure(1);hold on;
%plot(cf(:,1),cf(:,3),'b*')
%listlabel=[listlabel,{'inc03nocte'}];

fname='inc03ctenewdy_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,3),'bo')
listlabel=[listlabel,  {'inc03ctenewdy'}];
%listlabel=strjoin({listlabel, 'inc02ctenewdy'}, ',');

fname='inc03noctenewdy_01.cf';
filename=[fileroot,fname];
cf=load(filename);
figure(1);hold on;
plot(cf(:,1),cf(:,3),'bs')
listlabel=[listlabel, {'inc03noctenewdy'}];


legend(listlabel);
xlabel('time');
ylabel('dm');

