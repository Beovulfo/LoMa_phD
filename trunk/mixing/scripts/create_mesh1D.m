% PROGRAM aimed to create mesh for LISOuc3                      %
% with my points creates stretched mesh                         %
% A. ALMAGRO -  UC3M  (May 29 - 2013)                           %
%---------------------------------------------------------------%
% function "create_mesh1D"                                      %
%---------------------------------------------------------------%
% input: length -> half length of box                           %
% type: type of mesh you want to use                            %
% delta: parameter for mesh definition                          %
% my: number of points in mesh                                  %
%---------------------------------------------------------------%
function y = create_mesh1D(length,type,delta,my)

%tanh mesh
%delta=3; %strecthing factor
if (type=='tanh')
	eps=linspace(0,1,my);
	y=tanh(delta*(eps-1/2))/tanh(delta/2);
    % amplify mesh to lenght
       y=length/max(y)*y;
      
	%y=L*(1+tanh(delta*(eps-1/(2*L)))/tanh(delta/(2*L)));
end

%linear mesh
if (type=='line')
%not checked with lenght yet
       %1< delta < 2 for this mesh
        eps=linspace(0,1,my);
	L=2 %dimension from bot to top
	A=-delta;	
	yc=1;
	y=L*eps+A*(yc-L*eps).*(1-eps).*eps;
        %this goes from 0 to L
        y=y-L/2; %we want from -L/2 to L/2

    % amplify mesh to length
        for j=1:my
            y(j)=length*y(j)
        end
end

%hyper mesh
if (type=='hype')
%not checked with lenght yet
	eps=linspace(-1,1,my);
        y=tanh(eps/delta);
    % amplify mesh to length
        y=length/max(y).*y;
end

%mixing layer mesh
if (type=='mixl')
	eta=delta;%0<eta<0.5
	y=1:my; %init
	y(5:(my-4))=linspace(-eta*pi,eta*pi,my-8);
	gamma=1/(5*(tan(y(6))-tan(y(5)))+tan(y(my-5)));
	y=gamma*tan(y);
	tau=y(6)-y(5);
	y(1)=-1;
	y(my)=1;
	for j=1:5
	    y(j+1)=y(j)+tau;
	end

	for j=(my-4):my
	    y(j)=y(j-1)+tau;
	end
%amplify mesh for lenght
        for j=1:my
            y(j)=length*y(j);
        end

end	

end
