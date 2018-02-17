%Create new mesh
%choose my
%my = 531;Ly= 60;
%alpha = 0.5;beta=0.6; gg=0.1;

y2=linspace(-1,1,my);
Dy=(1+alpha)+alpha*tanh((y2-beta)./gg);
for j=1:my/2
     Dy(j) = Dy(end-j+1);
end
y2 = cumtrapz(y2,Dy);
ynew = 2*y2*Ly/y2(end) - Ly;
                %np.savetxt(filename, ynew)
plot(ynew(2:end)/dw,diff(ynew),'.')