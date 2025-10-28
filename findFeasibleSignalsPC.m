function [out_1,out_2] = findFeasibleSignalsPC(agg_u, greens_P, avg_x ,or_index, gn_R, PC, ind, twin_ind, minG)

% adjustSignalsPC: Translate greens calculated from the PC controller to
% feasible greens for all intersections before they are assigned  
%
%      agg_u: Control variable produced from the PI controller for an i-j
%             movement (scalar) 
%   greens_P: [Cell] Current durations of all stages of all nodes in ind (1
%               cell array per node)
% greens_Pco: [Cell] Current durations of all stages of the twin nodes 
%      avg_x: Average queues of the network at the end of the node's cycles  
%   or_index: index of origin links of all approaches 
%         PC: Perimeter Control struct (all info)
%       gn_R: Max green time fluctuation
%        ind: Indices of all PC intersections of the i-j movement in PC structure
%   twin_ind: Indices of all PC twin intersections
%       minG: Minimum green time imposed to all phases 

% Returns 
% feasibleGreens: the new feasible traffic signal settings (all phases) for the intersections 
% of the approach i-j

% Detailed explanation goes here: 
% Takes the calculated green time for the interregional approach i-j from the
% controller and produces feasible green times for all intersections of the 
% approach i-j that fit the constraints of: 
% 1. min/max green per phase
% 2. max fluctuation between consecutive control cycles 
% 3. Cycle remains constant (adjust all phases) 
% 4. Time be integer  
% 5. Queue balancing in the objective function 
%
% Make sure agg_u is between min/max values e.g. (0,1) or (PC.u_min,
% PC.u_max) - Impose min or max
% indu = agg_u > PC.u_max;
% agg_u(indu) = PC.u_max;
% indu = agg_u < PC.u_min;
% agg_u(indu) = PC.u_min;

%The current (previous) greens of the involved stages (first the main phases(all), then the secondary phases (all))
greensPC_prev = zeros(1,2*length(ind)); %initialize previous greens array - 2 phases per node
queues = zeros(size(greensPC_prev)); %prepare queues array
for i=1:length(ind) %for every PC node
    %previous greens (main and secondary phase) 
    greensPC_prev(i) = greens_P{i}(PC.stagesInvolved(ind(i),1)); %main phase green 
    greensPC_prev(i+length(ind)) = greens_P{i}(PC.stagesInvolved(ind(i),2)); %secondary phase green
    
    %approaches of main phase of node -- 
    approaches = PC.stageApproaches{ind(i),PC.stagesInvolved(ind(i),1)}; %the approaches involved at the main stage of node ind  
    
    %add approaches of main phase of twin node 
    if twin_ind(i)>0
        approaches = [approaches PC.stageApproaches{twin_ind(i), PC.stagesInvolved(twin_ind(i),1)}]; %add the approaches of the twin node 
    end
    
%     queue_set = zeros(1,length(approaches));
%     for j=1:length(approaches) 
%         queue_set(j) = avg_x(or_index(approaches(j)));%for every approach, add the last recorded average queue value
%     end
%     queues(i) = sum(queue_set);
    
    queues(i) = sum(avg_x(or_index(approaches)));
    
    
    %queues of secondary phase -- 
    approaches = PC.stageApproaches{ind(i),(PC.stagesInvolved(ind(i),2))};
    %Add queues of twin approaches 
    if twin_ind(i)>0 %if there is twin node 
        approaches = [approaches PC.stageApproaches{twin_ind(i), PC.stagesInvolved(twin_ind(i),2)}]; %add the approaches of the twin node 
    end
    
%     queue_set = zeros(1,length(approaches));
%     for j=1:length(approaches) 
%         queue_set(j) = avg_x(or_index(approaches(j)));  
%     end
%     queues(i+length(ind)) = sum(queue_set); 

    queues(i+length(ind)) = sum(avg_x(or_index(approaches)));
end

%Minimum greens per phase (same for all) 
minGreens = ones(size(greensPC_prev))*minG;

%Maximum allowed fluctuation per phase (same for all)
maxFluct = ones(size(greensPC_prev))*gn_R;

%Formulate the optimization problem  
yalmip('clear')
%Define variables: the feasible greens of the phases involved 
x = intvar(1,length(greensPC_prev)); %the feasible greens to be assigned after solution

%the total available assignable green per node
sumgreens = greensPC_prev(1:length(ind))+greensPC_prev(length(ind)+1:end); %the sum of greens - this can also be read from the data (check) 

%New greens derived by the controller (all equal to agg_u) ------ 
greensPC = zeros(size(greensPC_prev));
greensPC(1:length(ind)) = agg_u; %the new green of the primal stage 
greensPC(length(ind)+1:end) = sumgreens-greensPC(1:length(ind)); %the new green of the secondary stage  

%Define Constraints: greens > min, fluctuation, sum of greens constant
%(cycle maintained) 
Const = [x-minGreens >= 0, -maxFluct <= x-greensPC_prev <= maxFluct,...
    x(1:length(ind))+x(length(ind)+1:end)-sumgreens==0]; 

%Queues for queue balancing: Q = diagonal matrix with the sum of queues of
%every approach: Average queues of the phase? (+ queues of the twins???) 


Q = diag(queues); 

%Define and objective: 
Obj = (x-greensPC)*Q*(x-greensPC)';

% Set some options for YALMIP and solver
opt = sdpsettings('verbose',1,'solver','gurobi','verbose',0);

% Solve the problem
sol = optimize(Const,Obj,opt);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 newGreens = value(x);
else
 disp('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

%Assign all greens + twin nodes
out_1 = newGreens(1:length(ind)); %greens of main phase of PC nodes
out_2 = newGreens(length(ind)+1:end); %greens of secondary phase of PC nodes 

end

