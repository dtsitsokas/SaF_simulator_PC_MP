%% Analysis of MP only - TRB paper - Run experiments:
clear
close all
clc

%prepare input variables:
load indata.mat
path = 'Y:\common\Dimitris\Max Pressure\hEART Results\Random Node Sel Results\Random node sets\';
per = [0.05 0.1 0.15 0.2 0.25 0.30]; %percentage of all nodes for MP
%p = per(1);
%load MP_0 scenario settings (all nodes listed)
load 'input_MP_0.mat'

%% Create and save random sets of variables 
for j=1:1%length(per)
    p = per(j)
    for r = 1:10 %random choices

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
        save(strcat(path,'RandomNodeSet_',num2str(p*100),'_',num2str(r)),'rMPnodes');

    end
end

%fname = strcat(path,'output_MP_R_',num2str(p*100),'_',num2str(r),'.mat');
%SaF_3(indata,MPnodes,PC,fname);