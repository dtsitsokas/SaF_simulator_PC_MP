%From NC (or any) scenario - information about regions 
clear 
clc
%load results of the scenario 
load('output_NC.mat')
load('FinalInput.mat')
% - In which region cars queue for longer? (veh*hours of waiting per
% region)

weight = zeros(1,3);
for j=1:3 
    weight(j) = sum(sum(outdata.w(reInd{j},:)));
end
weight = weight*10^(-08);

%Construct K_I based on the weights: \

%1-2 / 3-2 / 2-1 / 2-3 
PC.K_P = [1 -weight(2)/weight(1) 0; 0 -weight(2)/weight(3) 1; -weight(1)/weight(2) 1 0; 0 1 weight(3)/weight(2)];
