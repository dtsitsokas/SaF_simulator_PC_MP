
function [defac] = demandFactor(DT, kmax, incr_factor)
%demandFactor : Creates a mat file that will be loaded in the model with
%the demand factor in every step 
% Input: 
% -----------------------------------
%          DT: The time step of the simulation 
%        kmax: The max number of time-steps that will run 
% incr_factor: The desired scaling factor for the simulation (1 = Aimsun)
% -----------------------------------
% Output: 
% defac.mat: A mat file with variable defac (demand factor) which contains the values of 
% the demand factor for each time step 

% ==== Aimsun Pattern for San Francisco (DBL paper) ======================
% m: time points of change / incr: values of the demand factor 
% The demand pattern for cars is also applied for the buses demand 
% We consider zero demand after 3 hours of simulation 
% m = [0.25 0.5 0.75 1 2 2.25 2.5 2.75 3 100];  
% incr = [1 1.25 1.5 1.75 2 1.75 1.5 1.25 1 0]*4*incr_factor; %to x4 einai gia na diorthwsw to la8os me to matlab kai na eimai consistent me to aimsun
% 
% defac = zeros(1,kmax);
% 
% for k=1:kmax
%     for j=1:length(m)
%         if k*DT <= m(j)
%             defac(k) = incr(j);
%             break
%         end
%     end    
% end


%% ===== Aimsun pattern for Barcelona (MaxPressure)
% 15 mins warm-up, 105 mins full, rest 4.5 hours with zero 
% Considering the same OD matrix (needs modification for different ones)
m = [0.25 2 100];
incr = [1 2.02 0]*incr_factor; 

%m = [0.25 0.5 0.75 1 2 2.25 2.5 2.75 3 100];  
%incr = [1 1.25 1.5 1.75 2.02 1.75 1.5 1.25 1 0]*4*incr_factor; %to x4 einai gia na diorthwsw to la8os me to matlab kai na eimai consistent me to aimsun

defac = zeros(1,kmax);

for k=1:kmax
    for j=1:length(m)
        if k*DT <= m(j)
            defac(k) = incr(j);
            break
        end
    end    
end

save('defac.mat','defac');
% More patterns from previous cases in the DBL folder end

