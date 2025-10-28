function [accumulations,VHT] = runSim4calib(parameters,s_kp_ki)
%RUNSIM4CALIB Executes one simulation run for given PC parameters and
%   returns total travel time and regional accumulations at specific times
%   Function to run SaF traffic model with Perimeter Control for specific
%   values of parameters K_p, K_i and set-point n* and return traffic
%   performance measures: total travel time (VHT) and regional
%   accumulations n_i in specific time intervals 
%
%   Inputs: 
%   ------
%   parameters: PI controller parameteres & set-point [K_p K_I setpoint
%   activ Criterion] 
%
%   Outputs: 
%   --------
%   accumulations: array of regional accumulations at specific times 
%   VHT: total vehicle hours travelled of the simulation (objective) 
%   -----------------------------------------------------------------
%

%update PC struct based on parameters 
load('input_PC_0','PC'); 
PC = updatePCparameters(PC,parameters,s_kp_ki);
load('indata.mat','indata')
%call simulation 
%[VHT, accumulations] = SaF_3_4calib_UpdtTRsubset_2(indata,PC); 
[VHT, accumulations] = SaF_3_4calib_UpdtTR_2(indata,PC); 

end

