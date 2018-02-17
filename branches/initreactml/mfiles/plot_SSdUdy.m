clear all; close all; 

%vcasesLF10 = {'S01LF10','S05LF10','S10LF10','S15LF10'};
%vcasesLF20 = {'S01LF20','S05LF20','S10LF20','S15LF20'};
%vcasesLF03 = {'S01LF03','S05LF03','S10LF03','S15LF03'};
%vcases = vcasesS05;
%figname = 'SSdUdy_S10';
%figname = 'SS_S10';
%vcasesS01 ={'S01LF03','S01LF10','S01LF20'};
%vcasesS05 ={'S05LF03','S05LF10','S05LF20'};
%vcasesS10 ={'S10LF03','S10LF10','S10LF20'};
%hig  = {'q40S9sigma00','q40S9sigma05','q20S4sigma00','q0'}


%INPU
%----------------------------------------
figname = 'SSdUdy_S01';
vcases ={'S01LF03','S01LF10','S01LF20'};
fsave = false;



count = 1
%vcases=vcases1;

n=length(vcases);
%list=jet(n);
list=hsv(n);
style=['b-','g--','r..'];

for fname1 = vcases
    %Convert to physical for a giving y vector
    fpath = './'; fname=fname1{1}; 
    fSSphys = [fpath 'SSphys_' fname '.mat']
    load(fSSphys)
    Ts = 1./rhos;
%y1 = 200*exp(-0.05*x).*sin(x);
%y2 = 0.8*exp(-0.5*x).*sin(10*x);
    rhodUdy = -rhos.*dUdy;
    
    figure(69);
    if count ==1
        [AX,H1,H2] = plotyy(y,U,y,-rhos.*dUdy,'plot');
         set(get(AX(1),'Ylabel'),'String','U') 
         set(get(AX(2),'Ylabel'),'String','-\rho dU/dy') 
         xlabel('y') 
         %title('Velocity and deriv') 
         set(H1,'LineStyle','--')
         set(H2,'LineStyle',':')
         %legend(AX,fname)
     end

    axes(AX(1))
    xlim([-3,3])
    ylim([-1.1,1.1])
    hold on
    [K(count)]=plot(y,U,style(count))
    set(K,'LineStyle','--')
    %plot(y,U,'col',list(count,:))
    %plot(y,rhodUdy,'col',list(count,:),'LineStyle','--')

    axes(AX(2))
    xlim([-3,3])
    %ylim([1,6])
    hold on
    L(count) = plot(y,-rhos.*dUdy,style(count))

       % plot(y,Ts,'col',list(count,:))

     count = count + 1;

end
    legend(L,vcases)
    xlim([-3,3]); grid on
    print([figname],'-dpng','-r0')
