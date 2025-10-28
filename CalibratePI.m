% Calibrate PI controller ---- 
clear
clc
% Prepare PC input - set parameters 

%prepare input variables:
[indata,~,demandGroup2] = prepareInput();

% load scenario settings
load('input_PC.mat')

%based on weights from waiting vehicles over time 
testK_I{1} = [1.0000 -1.2007 0;...
    0   -1.5103    1.0000;...
   -0.8329    1.0000         0;...
   0    1.0000    0.6621];
scale = 0.1;
test_set{1} = [6724; 8252; 4453]*0.8; %critical
%%
all_results = {};
for i=1:length(testK_I)
    
    %decide how to set Kp Ki
    
    PC.K_I = testK_I{i}*scale;
    % PC.K_P = [];
    PC.n_set = test_set{i};   
    
    % Call simulator for all scenarios
    tic
    [all_results{i,j}] = SaF_3forPICalibration(demandGroup2,indata,PC);
    toc 
    
end
 
%% FInd the best combination (minimum TTT)
x = [];
for j=1:5
    x = [x all_results{1,j}.TTT];
end
plot(x)