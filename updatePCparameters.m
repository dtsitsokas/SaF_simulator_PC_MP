function [PC] = updatePCparameters(PC,parameters, s_kp_ki)
%PC = struct containing all information for perimeter control 
%parameters = [Kp, Ki, setpoint, activationCriterion] 
%s_kp_ki = size of Kp, Ki matrices 

s_setP = 3; %number of regions  

PC.K_P = zeros(s_kp_ki);
PC.K_I = zeros(s_kp_ki);

%setpoint 
PC.n_set = transpose(parameters(end-s_setP-1:end-2)); 

%parameters 
k = 1;
for i=1:s_kp_ki(1)
    PC.K_I(i,:) = parameters(k:k+s_kp_ki(2)-1);
    PC.K_P(i,:) = parameters(k+s_kp_ki(1)*s_kp_ki(2):k+s_kp_ki(1)*s_kp_ki(2)+s_kp_ki(2)-1);
    k = k + s_kp_ki(2);
end

%activation criterion 
%from 0 to 0.20 (activ when at least one region is less than 20% of the setpoint) 
PC.actCrit = parameters(end-1);
PC.deactCrit = parameters(end); 
%
return