% PROGRAM aimed to create mesh for LISOuc3 %
% with my points creates stretched mesh    %
% A. ALMAGRO -  UC3M  (May 29 - 2013)      %
%------------------------------------------%
% function "create_mesh1D"                 %
%------------------------------------------%
function y = create_mesh1D(type,delta,my)

%delta=3; %strecthing factor
if (type=='tanh')
	eps=linspace(0,1,my);
	y=tanh(delta*(eps-1/2))/tanh(delta/2);
	%y=L*(1+tanh(delta*(eps-1/(2*L)))/tanh(delta/(2*L)));
end

if (type=='line')
       %1< delta < 2 for this mesh
        eps=linspace(0,1,my);
	L=2 %dimension of the channel
	A=-delta;	
	yc=L/2;
	y=L*eps+A*(yc-L*eps).*(1-eps).*eps;
        %this goes from 0 to L
        y=y-L/2; %we want from -1 to 1
end
%hyper mesh
if (type=='hype')
       %delta should be around 1.02
	eps=linspace(-1,1,my);
        C=delta;
        K=atanh(1/C);
        L=2;
        y=L/2*(1+C*tanh(K*eps));
        y=y-L/2; %we want from -1 to 1
end

	%y=L*(1+tanh(delta*(eps-1/(2*L)))/tanh(delta/(2*L)));
%if (type=='sin1')
%      delta=sin(0.33*eps);
%      y=2*sin(delta.*eps);
%end
