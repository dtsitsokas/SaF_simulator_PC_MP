% Script to automatize node selection process for MP in incremants
%
% repeats the following loop:
% - load reference case
% - set selection target (5%) - the increment
% - excludes nodes previously selected (if any)
% - run the selection process until target is reached
% - evaluate the new layout of MP nodes -> new reference case

% it will alternate between nodeSelection and simulation running

clear
close all
clc

% There are two ways of increasing percentage of nodes:
% 1) decrease congestion threshold of step 1 (increase the subset size of
% the preselection)
% 2) decrease congestion + variance thresholds up to some minimum values
% before further decresing congestion threshold, in order to include

%In this analysis we will go with (1): decrease accepted congestion
%threshold / Need to try the second one too

% Start from zero: Reference case is No Control
% set file name (full path) of the reference case results

path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Current results\ResultsMediumOD_capacityDrop\';
refname = strcat(path,'output_set5_MedDemand_NC');
p_step = 5; %percentage of first case
load('FinalInput','MP')
load('indata.mat','indata')
num_MPall = sum(MP.splan==1)-length(indata.specialInt);  %total number of MP eligible nodes
%Define peak period as starting and ending time step (2 hours: 1h00-3h00)
k_s = 41*indata.c_int; %starting step of the peak interval
k_e = 120*indata.c_int; %ending step of the peak interval

%Exclude special intersections
MPspInt = zeros(1,length(indata.specialInt)); %ind of special inters in MP
for i=1:length(indata.specialInt)
    %indices of special intersections in MP
    MPspInt(i) = find(MP.nodeID == indata.specialInt(i));
end

%upper limit applies for this case (only critical, not medium)
lim_var{1} = [0.08 0.08]; %region 1
lim_occ{1} = [0.2 0.29];
lim_var{2} = [0.05 0.09]; %region 2
lim_occ{2} = [0.25 0.29];
lim_var{3} = [0.05 0.06]; %region 3
lim_occ{3} = [0.1 0.14];

% - Set initial congestion threshold for peak period of 2 hours/80 steps
cong_threshold = 70; %number of steps node remains congested (out of 80)
PC.mode=0;
MPnodes = []; 
for s = 1:6 %cases to evaluate (5-10-15-20-25-30%) 
    
    scenario = strcat('MP_IS_',num2str(s));

    %- load reference case
    load(refname,'indata','outdata')
       
    while 1
        tMPnodes = [];
        for reg=1:3
            % Calculate for every node for the peak-period of 2 hours (0.5,2.5) = 80 cycles:
            % 1. How many cycles each node is "congested" - keep those that are congested for more than 70/80
            %    Congested node: at least one incoming link more than 80% of its capacity
            % 2. Put the remaining nodes in two graphs:
            %   - std of queues vs. mean x/c of all inc links of the two hours
            %   - std of queues vs. number of cycles a link remains congested in
            %   the two hours
            % 3. Start selecting nodes by scanning the figures from uppp and right
            % until I reach the wished node percentage (10% for C1, 25% for C2)
            
            for i=1:length(indata.nodereg{reg})
                %scan all nodes of the region
                ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
                %calculate mean of the average queue of all incoming links over
                %time - normalize over storage capacity
                if MP.approaches{ind}>0
                    MP_nodes{reg}(i) = ind;
                    %             node_occ = outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}));
                    node_occ = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
                    node_queues{reg}(i) = mean(mean(node_occ));
                    %calculate mean of the variance of queues of all incoming links
                    %over time
                    %             node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
                    node_var{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
                    % count number of cycles
                    node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
                end
            end
                              
            %decide on nodes - critical: over limits 2
            ind_1 = (node_queues{reg}>lim_occ{reg}(2)); %ind of nodes with occupancy above theshold (congestion)
            ind_2 = (node_var{reg}>lim_var{reg}(2)); %ind of nodes with variance above threshold
            ind_3 =  (node_congested{reg}>cong_threshold); %extra threshold of congestion
            %ind_4 = ind_1+ind_2 == 2;
            ind_4 = ind_1 + ind_2 +ind_3 == 3; %nodes satisfying all criteria simultaneously (AND operator)
            tMPnodes = [tMPnodes MP_nodes{reg}(ind_4)];  %store here the temporary MaxPressure nodes (ind in MP)
            
            
            not_congested{reg} = (~ind_3);
        end
        
        p_1 = length(unique([MPnodes tMPnodes]))/num_MPall*100
        
        if p_1 < s*p_step 
            disp(strcat('solution rejected: node percentage below target of',num2str(s*p_step))) 
            if cong_threshold >= 32 
                cong_threshold = cong_threshold - 1;
                
            else
                %reduce limits 
                lim_var{1}(2) = lim_var{1}(2)-0.005; %region 1
                lim_occ{1}(2) = lim_occ{1}(2)-0.005; 
                lim_var{2}(2) = lim_var{2}(2)-0.005; %region 2
                lim_occ{2}(2) = lim_occ{2}(2)-0.005;
                lim_var{3}(2) = lim_var{3}(2)-0.005; %region 3
                lim_occ{3}(2) = lim_occ{3}(2)-0.005;
            end
            tMPnodes = []; 
        else
            %keep the layout 
            MPnodes = [MPnodes tMPnodes];
            tMPnodes = []; 
            break 
        end 
        
    end
           
    % Exclude special intersections (MP cannot be applied due to signal
    % settings) 
    MPnodes = setxor(MPnodes, intersect(MPspInt,MPnodes));

    % Save selected set of nodes:
    save(strcat('MPnodesCase_',scenario,'.mat'),'MPnodes','p_1','lim_occ','lim_var','cong_threshold')

    % Generate the input_scenario file
    prepareExperimentfunct(scenario, scenario, strcat('Incremental Adding, step_',num2str(s), ' percentage ',num2str(p_1)))
    
    % Run simulation for the new scenario 
           
    fname = strcat(path,'output_',scenario,'_',num2str(s*5),'.mat');
    
    SaF_3(indata,MPnodes,PC,fname);
    
    disp(strcat('Simulation of scenario ', scenario ,' is finished.'))
    
    % Update reference case 
    refname = fname; 
    
end

