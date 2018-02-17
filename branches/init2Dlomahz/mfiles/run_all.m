clear all; close all;
startup; %set environment and paths
%Set dominium
%Chose if start from previous resultso or not.
vcases0 = {'S01LF10','S01LF20','S01LF03'};
vcases1 = {'S05LF10','S05LF20','S05LF03'};
vcases2 = {'S10LF10','S10LF20','S10LF03'};
vcases3 = {'S15LF10','S15LF20','S15LF03'};
%vcases3 = {'S15LF20','S15LF03'};
%vdom    = [20,30,30];
%vcases3 = {'S15LF10'};
%vcasesLF10 = {'S01LF10','S05LF10','S10LF10','S15LF10'};
%vcasesLF20 = {'S01LF20','S05LF20','S10LF20','S15LF20'};
%vcasesLF20 = {'S05LF20','S10LF20','S15LF20'};
%vcasesLF03 = {'S05LF03','S10LF03','S15LF03'};
%vcasesg6  = {'S01LF20g6','S15LF20g6'};
%vcasesLF03 = {'S01LF03','S05LF03','S10LF03','S15LF03'};
%vcasesLF03 = {'S01LF03','S05LF03','S10LF03','S15LF03'};
%vhiguera ={'q40S9sigma00','q40S9sigma05','q20S4sigma00','q0'}
%vcases = {'S01LF11'}
%vcasesLF03 = {'S05LF03','S10LF03','S15LF03'};
%vcasesLF03 = {'S01LF03','S05LF03','S10LF03','S15LF03'};
%vdom    = [10,10,10];

%vcases = {'S01LF10g2','S01LF10g3','S05LF10g1','S05LF10g2','S05LF10g3'};
%vcases  = {'S015LF10g1','S015LF10g2','S015LF10g3','S015LF10g0'}
%vcases = {'S015LF10g60','S015LF10g01','S015LF10g70','S015LF10g80','S010LF10g01'}
%vcases = {'q0','S10LF10','S10LF03'};
%vcases = {'S40LF03g65'}
vcases = {'S01LF10g6','S01LF20g6','S01LF03g6'};




%vcases = vcasesg6;

%initialize
icase = 1
for defcase = vcases
    %vectb= [50,100, 200,500,1E3];
    %vectb= [200,500,1E3];
    vectb= [20,50,100,200,300,500];
    fguess=0;
    for beta = vectb
        %Load parameters of case within loop for beta
        defcase{1}
        eval(defcase{1});
        domlen     = 15; 
        %domlen     = vdom(icase); 
        bet1 = beta
        params; %Calculate Combustion parameters from input Params
        tuneChebfun; %Settings for the solver
%Call solver
        mixlayer2
        %mixlayer
        pause(10)
        fguess=1;
    end
    fpath = './'; fname = defcase{1}
    cheb2physB
    icase = icase + 1;
end
exit;
