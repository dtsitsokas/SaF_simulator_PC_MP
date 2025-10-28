%% Main program to run tests for the combined MP+PC controller
% Run random sets 
clear
close all
%global X1 SP1 M2
clc

%prepare input variables (important if any parameter change or 
% if Initialization is changed ):
%[indata] = prepareInput();

load indata.mat

%% Call simulator for random sets 
%path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\fourth set results\Random Node Sel Results\';
path0 = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Random node sets\';
path1 = 'Y:\common\Dimitris\Max Pressure\hEART Results\Random Node Sel Results\';
p = 0.30; 
PC.mode=0; 
for r = 1:10
    fname = strcat(path1,'output_MP_R_',num2str(p*100),'_',num2str(r),'.mat');
    load(strcat(path0,'RandomNodeSet_',num2str(p*100),'_',num2str(r),'.mat'),'rMPnodes');
    SaF_3(indata,rMPnodes,PC,fname);
    disp(strcat('Replication no ',num2str(r),' finished.'))
end


