%% Analysis of MP only - TRB paper - Run experiments:
clear
close all
global X1 SP1 M2
clc

%prepare input variables:
[indata,DBLplan,demandGroup2] = prepareInput();
load scenarios.mat
per = [0.1 0.25 0.50 0.75]; %percentage of all nodes for MP
p = per(4); 
%load MP_0 scenario settings (all nodes listed)
load(strcat('input_',scenarios{2},'.mat'))


parfor r = 1:10 %random choices
    
    sel = []; %selected
    while 1
        r_n = ceil(rand()*length(MPnodes));
        if ~ismember(r_n,sel)
            sel = [sel r_n];
        end
        if length(sel) >= round(p*length(MPnodes))
            break
        end
    end
    
    rMPnodes = MPnodes(sel);
    fname = strcat('output_MP_4_',num2str(p*100),'_',num2str(r),'.mat');
    SaF_3(demandGroup2,indata,rMPnodes,PC,fname);
    
end
