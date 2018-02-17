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
        y=linspace(-length,length,my);

    % amplify mesh to length
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
if (type=='beta')
%y=tanh(y/y0)
	eta=linspace(-length,length,my);
        y0=delta; %final layer thickness
        y=length*sinh(eta/y0);
end


if (type=='manl')
       %1< delta < 2 for this mesh
        y=linspace(-1,1,my);
        a=[1,0,28.5,-47,6,12];%Dy(0)=Dy0,Dy(1)=Dy0/2;Dy(0.6)=3Dy0;
%integrating Dy(Dy is a polynomial of order 5) gives y2
        y2=a(1)*y+a(2)*y.^2/2+a(3)*y.^3/3+a(4)*y.^4/4+a(5)*y.^5/5+a(6)*y.^6/6;
%now create the symmetric part
        if (rem(my,2)==0)
        y2(1:floor(my/2))=-y2(end:-1:floor(my/2)+1);
       else
        y2(1:floor(my/2))=-y2(end:-1:floor(my/2)+2);
       end
             %y2(j)=-(a(1)*abs(y(j))+a(2)*abs(y(j))^2/2+a(3)*abs(y(j))^3/3+a(4)*abs(y(j))^4/4+a(5)*abs(y(j))^5/5+a(6)*abs(y(j))^6/6);
%scale y
        y=y2/max(y2)*length;
       

    % amplify mesh to length
end
%
if (type=='file')
     y = load('mesh.txt')
     y = squeeze(y')
end





end
