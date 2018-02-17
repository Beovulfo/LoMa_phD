clear all; close all;
startup; %set environment and paths
tuneChebfun

%vcasesLF03={'S01LF03','S05LF03','S10LF03','S15LF03'};
%vdom    = [20,30,30];
%vcases3 = {'S15LF20','S15LF03'};

vcasesLF10 = {'S01LF10','S05LF10','S10LF10','S15LF10'};
vcasesLF20 = {'S01LF20','S05LF20','S10LF20','S15LF20'};
vcasesLF03 = {'S01LF03','S05LF03','S10LF03','S15LF03'};
vcaseshig  = {'q20S4sigma00'}
vcaseshig  = {'q40S9sigma00','q40S9sigma05','q20S4sigma00','q0'}
%vcases =  [vcasesLF10,vcasesLF20,vcasesLF03];
%vcases=;
%vcases = [vcasesLF10, vcasesLF20, vcasesLF03];
%vcases = vcaseshig
%
vcases = vcaseshig;
%vcases = {'S01LF10g05'}
%vcases={'q0','S10LF10','S10LF03'};
vcases={'S015LF10g01'};
%vcases = {'q0'}


for fname1 = vcases
    %Convert to physical for a giving y vector
    fpath = './'; fname=fname1{1} 
    cheb2physB
    %cheb2phys
end
