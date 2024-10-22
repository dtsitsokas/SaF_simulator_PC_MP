
function [defac] = demandFactor(DT, kmax, incr_factor)
%demandFactor : Creates a mat file that will be loaded in the model with
%the dynamic multiplier for the OD demand matrix, which will eventually
%create the demand profile (considering a global increase/decrease factor
%and the desired relative dynamic profile) 

% Input: 
% -----------------------------------
%          DT: The time step of the simulation 
%        kmax: The max number of time-steps that will run 
% incr_factor: The desired scaling factor for the simulation (1 = Aimsun)
% -----------------------------------
% Output: 
% defac.mat: A mat file with variable defac (demand factor) which contains the values of 
% the demand factor for each time step 


% ===== Dynamic demand pattern for Barcelona case study (MP / PC) 

% Considering the same OD matrix (needs modification for different ones)
m = [0.25 2 100]; %duration of time-windows for different multiplier [hours],
% for the last one you can put a high number > max sim time, for flexibility (e.g. if you run 6 or 8 hours until network empty)  - 
% here: % 15 mins warm-up, 105 mins full, rest 4.5 hours with zero 
incr = [1 2.02 0]*incr_factor; 

%create the multiplier vector (demand factor) 
defac = zeros(1,kmax);
for k=1:kmax
    for j=1:length(m) %for every time window 
        if k*DT <= m(j)
            defac(k) = incr(j);
            break
        end
    end    
end

save('defac.mat','defac');

