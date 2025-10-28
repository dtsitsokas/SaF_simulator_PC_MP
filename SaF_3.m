%SaF_3: Simulator with virtual links for all centroids
%in the network. Single queues per link. Max-Pressure local control in
%selceted intersections and perimeter control options 

% INPUTS:
% indata: all network-related input data (struct)
% MPnodes: IDs of nodes where MP will be applied 
% PC: perimeter control related settings 
% fname: full path file-name for result file 
% savefull: binary indicator for saving (1) or not (0) result file 

%OUTPUTS: 
% r3: Total Travel Time in [veh x hours]
% *generation/saving of results file with all numerical results (queues
% etc)

function [outdata] = SaF_3(indata,MPnodes,PC,fname,savefull)

tic
% some additional settings: 
minSatFlowGating = 15; %lower bound of percentage of sat flow for external gating 

%% Initialize matrices ----------------------------
outdata.x = zeros(indata.NLinks2,indata.kmax);                               % Queues: (queue in time 1, ... , queue in time k)
outdata.m = zeros(indata.NLinks2,indata.kmax);                               % Moving part
outdata.w = zeros(indata.NLinks2,indata.kmax);                               % Waiting part
outdata.u = zeros(indata.NLinks2,indata.kmax);                               % Outflows
outdata.q = zeros(indata.NLinks2,indata.kmax);                               % Inflows
outdata.a_z = zeros(indata.NLinks2,indata.kmax);                             % Q arriving flows
outdata.rho_l = zeros(indata.NLinks2,indata.kmax);                           % min steps required to travel the distance until the end of the queue in constant v_ff
outdata.l_ff = zeros(indata.NLinks2,indata.kmax);                            % current length of the moving part of the road
tripsCompleted = zeros(length(indata.group3),indata.kmax);                   % Flows entering links with destination centroids / trips completed [veh/DT]
outdata.cumtripsCompleted=zeros(length(indata.group3),indata.kmax);          % Cummulative Flows entering links with destination centroids / trips completed
% tripsStarted=zeros(length(indata.group2),indata.kmax);                     % Flows entering origin links [veh/DT]
% cumtripsStarted=zeros(length(indata.group2),indata.kmax);                  % Cummulative Flows entering origin links
outdata.upair=zeros(length(indata.junct2.origin),indata.kmax);               % Outflow of the specific approach indata.junct.origin to indata.junct.destination at time step k
t_winSteps = ceil(indata.t_win/indata.DT);                                   % No of steps for the time window to update turn rates
speeds = zeros(size(indata.LinksP,1),ceil(indata.kmax/t_winSteps));          % speeds of all links per interval
outdata.node_heterog = zeros(size(indata.Nodes,1),ceil(indata.kmax/indata.c_int));  % Heterogeneity of nodes =
outdata.avg_q_cycle = zeros(size(indata.Links2,1),ceil(indata.kmax/indata.c_int));  % stores the average queues during the last cycle of the downstream node
outdata.greentimes2 = zeros([size(indata.greentimes2),ceil(indata.kmax/indata.c_int)]);

% --PC variables --------------------
if PC.mode==1
    % initialize matrices 
    % -- regional accumulations
    % three-region PC
    agg_n = zeros(indata.no_reg,ceil(indata.kmax/indata.c_int));                  % aggregated regional accumulations per control cycle
    
    %only external perimeter (all nodes)
    %agg_ng = zeros(1,ceil(indata.kmax/indata.c_int));
    
    % -- control variables
    % PC without external perimeter (only interregional)
    % agg_u = zeros(indata.no_adjReg,ceil(indata.kmax/indata.c_int));              % mean green time calculated by the controller per cycle
    
    % PC with external perimeter + selected nodes
    agg_u = zeros(indata.no_adjReg+3,ceil(indata.kmax/indata.c_int));              % mean green time calculated by the controller per cycle
    
    % only external perimeter PC (as in one region network)
    % agg_ug = zeros(1,ceil(indata.kmax/indata.c_int));
    
    activFlag = zeros(1,ceil(indata.kmax/indata.c_int)+1);                       % Flag about the activation status of the controller (0/1)
    actCheck = zeros(1,3);                                                       % Initialize activation check for the 3 regions
    greensToApply = zeros(length(PC.nodeID),2);                                  % Green settings ready to apply to the main approach of PC nodes when their cycle is completed
    applied_u = zeros(indata.no_adjReg+3,ceil(indata.kmax/indata.c_int));        % The applied values of u after projection to feasible space
    %applied_ug = zeros(1,ceil(indata.kmax/indata.c_int));
    init_greentimes = indata.greentimes2;                                        % The fixed-time signal settings initially given as input
    init_stageDur = PC.stageDur;                                                 % Initial fixed-time settings (in form of duration of phases) - same info
    init_stageDurMP = indata.MP.stageDur;

    %Initial values of applied_u (interregional + selected external)
    %Control variables = green time between adjacent regions 
    
    greens_P = [];
    for i=1:indata.no_adjReg %regions + external perimeter PC
        ind = PC.indices(i,PC.indices(i,:)>0); %all nodes that belong to this approach i-j (non-zero) - the k indices to refer to each node in PC structs
        for j = 1:length(ind)
            greens_P = [greens_P PC.stageDur{ind(j)}(PC.stagesInvolved(ind(j),1))]; %the greens of the main phase of all nodes controlling the i-j
        end
        applied_u(i,1) = mean(greens_P); %set the initial values for the controller
    end
    
    %initial values for VQ flows = 100 so that sat flow rate = 100/100 = 1
    for i=1+indata.no_adjReg:indata.no_adjReg+3 %regions + external perimeter PC
        applied_u(i,1) = 100; %set the initial values for the controller
    end
end

%indices of group3 pairs (exits) in junct2 per region 
indG3r = cell(1,indata.no_reg);
init_ExitLinksLanesR = cell(1,indata.no_reg);
for jk=1:indata.no_reg
    indG3r{jk} = []; 
    set1 = intersect(indata.reInd{jk},indata.group3);
    for j = 1:length(set1)
        indG3r{jk} = [indG3r{jk} find(set1(j) == indata.junct2.dest_index)];
    end
    init_ExitLinksLanesR{jk} = indata.junct2.lanesd(indG3r{jk});

end
    
k = 1; %time step index
clock = indata.DT*3600; %[sec]

% I use junct2 everywhere (also connections with virtual links)
sig = zeros(length(indata.junct2.origin),1);

% Check traffic signal lights
t = 1-indata.junct2.splan(:); %for every movement, 1 if there is no traffic signal, 0 if there is signal plan
sig(t==0) = signal(clock,indata.junct2.cycle(t==0),indata.junct2.offset(t==0),indata.greentimes2(t==0,:));
sig(t==1) = 1; %To make sure the movements without signal plan always get sig=1 (always open) 

% Outflows of all physical approaches
t_index = 1; %(when turn rate update is on) 
%t_index = min(ceil(k*indata.DT/indata.t_win),size(indata.junct2.turn,2));  %turnings change every t_win [hours] / Keeps the last column if the run is longer

%initialize VQs
for i=1:length(indata.group2)    
    outdata.w(indata.group2(i),1) = indata.DT*(indata.demandGroup2(i,k)*indata.defac(1));
end

outdata.upair(:,k) = sig.*OutflowS_2(indata.junct2.lanesu(:), indata.junct2.lanesd(:), indata.Links2(indata.junct2.or_index(:),2),...
    indata.capacity(indata.junct2.dest_index(:)), outdata.w(indata.junct2.or_index(:),k), outdata.x(indata.junct2.dest_index(:),k),indata.junct2.turn(:,t_index),indata.DT);

t_MP = ceil(indata.applyMPfreq*indata.MP.cycle/(3600*indata.DT)); %how many steps is a cycle for every node in MP

%Calculate the total outflow and inflow of the upstream links of every
%movement in indata.junct2 in time step t = 1 (includes vqueues-downstream
%and upstream-vlinks)
for i=1:indata.NPairs2
    %total outflows of all links going to dowstream links
    
    outdata.u(indata.junct2.or_index(i),k) = outdata.u(indata.junct2.or_index(i),k) + outdata.upair(i,k);
    %total inflows of all links coming from upstream links
    
    outdata.q(indata.junct2.dest_index(i),k) = outdata.q(indata.junct2.dest_index(i),k) + outdata.upair(i,k);
end


%Update VQ dynamics
for i=1:length(indata.group2)   
    outdata.w(indata.group2(i),k+1) = outdata.w(indata.group2(i),k) + indata.DT*(indata.demandGroup2(i,k)*indata.defac(1)- outdata.u(indata.group2(i),k));
end

%Calculate the arriving flow at the end of the queue: Based on the
%current queue of vehicles w - l_ff [NLinks*kmax] matrix -
%calculation for all links by using vectors
% Free flow length = ONLY for intermediate links (not for vlinks)

% l_ff(indata.group1,k)=max(0,(capacity(indata.group1)-w(indata.group1,k))*indata.vehlength./(indata.Links2(indata.group1,2)-buslanes(indata.group1,2)));
l_ff(indata.group1,k)=max(0,(indata.capacity(indata.group1)-outdata.w(indata.group1,k))*indata.vehlength./(indata.Links2(indata.group1,2)));

% Minimum time steps required to reach the end of the queue in a CONSTANT free flow speed
% rho_l: time step up to which we add the inflows at the
% cummulative
rho_l(indata.group1,k)=max(0,k-ceil(l_ff(indata.group1,k)/(indata.DT*indata.v_ff)));

% Calculation of the arrival flows at the end of the queues
% Cummulative arriving/entering flows:
test = rho_l(indata.group1,k)==1;
sum_q(indata.group1,k) = outdata.q(indata.group1,1).*test;
outdata.a_z(indata.group1,k) = sum_q(indata.group1,k);

% Calculations for every group of links
% Group 1 : Intermediate links (No origin centroid / No destination centroid)

outdata.m(indata.group1,k+1) = max(0,outdata.m(indata.group1,k)+indata.DT*((outdata.q(indata.group1,k)-outdata.a_z(indata.group1,k))));
outdata.w(indata.group1,k+1) = max(0,outdata.w(indata.group1,k)+indata.DT*((outdata.a_z(indata.group1,k)-outdata.u(indata.group1,k))));
outdata.x(indata.group1,k+1) = outdata.m(indata.group1,k+1) + outdata.w(indata.group1,k+1);

%Group 3 : Links with Destination centroid as their ending node
%I don't keep track of the link dynamics / I just store the inflows which
%are considered as copleted trips

%Trips completed in every time step / Cummulative
tripsCompleted(:,k)=indata.DT*(outdata.q(indata.group3,k));
outdata.cumtripsCompleted(:,k)=tripsCompleted(:,k); %Cummulative trips completed

tripsStarted(:,k)=indata.DT*(outdata.u(indata.group2,k)); %outflows of virtual queues into the network
outdata.cumtripsStarted(:,k)=tripsStarted(:,k); %Cummulative trips started


if PC.mode == 1
    
    PC.uniqueNodes = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
    
    % --- PC at all VQs in the sense of an inflow control (same for all VQs
    % of every region) - can happen externaly (move?)
    %initialize/store initial no of lanes for VQs
    
    PC.extGateLinksReg = cell(1,3); %index of VQ virtual links per region
    PC.extGateConnectionsReg = cell(1,3); %index of approaches between VQs and ntw links per region (in junct2)
    for i=1:indata.no_reg %for every region
        PC.extGateLinksReg{i} = intersect(indata.group2,indata.reInd{i}); %VQ link indices of the region
        %store indices of connections in junct to easily modify the VQ
        %lanes = the VQ saturation flow
        PC.extGateConnectionsReg{i} = ismember(indata.junct2.or_index,PC.extGateLinksReg{i});
    end
    PC.initialGateLanes = indata.junct2.lanesu;  %storage of all initial values (total saturation flows)
end

k_c = 0; %index for the control intervals independently of the activation of the controller

flag_modified_signals = zeros(length(indata.MP.nodeID),1); %if zero: initial settings (array including every intersection)

for k=2:indata.kmax %k = time counter of the simulation
    
    clock=clock+indata.DT*3600; %[sec]
    
    % Signals
    sig(t==0)=signal(clock,indata.junct2.cycle(t==0),indata.junct2.offset(t==0),indata.greentimes2(t==0,:));
    sig(t==1)=1; %To make sure the movements without signal plan always get sig=1
        
    % outflows of approaches
    outdata.upair(:,k) = sig.*OutflowS_2(indata.junct2.lanesu(:), indata.junct2.lanesd(:),...
        indata.Links2(indata.junct2.or_index(:),2),indata.capacity(indata.junct2.dest_index(:)), outdata.w(indata.junct2.or_index(:),k),...
        outdata.x(indata.junct2.dest_index(:),k),indata.junct2.turn(:,t_index),indata.DT);
   
    %Calculate the total Outflow and inflow of the upstream links of every movement in indata.junctions in time step t = 1
    for i=1:indata.NPairs2
        %total outflows of all links going to dowstream links
        outdata.u(indata.junct2.or_index(i),k) = outdata.u(indata.junct2.or_index(i),k) + outdata.upair(i,k);
        %total inflows of all links coming from upstream links
        outdata.q(indata.junct2.dest_index(i),k) = outdata.q(indata.junct2.dest_index(i),k) + outdata.upair(i,k);
    end
    
    %Update VQ dynamics
    for i=1:length(indata.group2)
        outdata.w(indata.group2(i),k+1) = outdata.w(indata.group2(i),k) + indata.DT*(indata.defac(k)*indata.demandGroup2(i,1+(k>indata.WU))- outdata.u(indata.group2(i),k));
    end
    
    %Calculate the arriving flow at the end of the queue: Based on the
    %current queue of vehicles w - l_ff [NLinks*kmax] matrix -
    %calculation for all links by using vectors
    
    % Free flow length
    % l_ff(indata.group1,k)=max(0,(capacity(indata.group1)-w(indata.group1,k))*indata.vehlength./(indata.Links2(indata.group1,2)-buslanes(indata.group1,2)));
    l_ff(indata.group1,k)=max(0,(indata.capacity(indata.group1)-outdata.w(indata.group1,k))*indata.vehlength./(indata.Links2(indata.group1,2)));
    
    % Minimum time steps required to reach the end of the queue in a CONSTANT free flow speed
    % rho_l: time step up to which we add the inflows at the
    % cummulative
    rho_l(indata.group1,k)=max(rho_l(indata.group1,k-1),k-ceil(l_ff(indata.group1,k)/(indata.DT*indata.v_ff)));
    
    % Calculation of the arrival flows at the end of the queues
    % Cummulative arriving/entering flows:
    for i=1:length(indata.group1)
        outdata.a_z(indata.group1(i),k) = sum(outdata.q(indata.group1(i),rho_l(indata.group1(i),k-1)+1:rho_l(indata.group1(i),k)));
    end
    
    %update dynamics of network links
    outdata.m(indata.group1,k+1) = max(0,outdata.m(indata.group1,k)+indata.DT*((outdata.q(indata.group1,k)-outdata.a_z(indata.group1,k))));
    outdata.w(indata.group1,k+1) = max(0,outdata.w(indata.group1,k)+indata.DT*((outdata.a_z(indata.group1,k)-outdata.u(indata.group1,k))));
    outdata.x(indata.group1,k+1) = outdata.m(indata.group1,k+1) + outdata.w(indata.group1,k+1);
    
    %Trips completed in every time step / Cummulative
    tripsCompleted(:,k)=indata.DT*(outdata.q(indata.group3,k));
    outdata.cumtripsCompleted(:,k)=outdata.cumtripsCompleted(:,k-1)+tripsCompleted(:,k); %Cummulative trips completed
    
    %Trips started:
    tripsStarted(:,k)=indata.DT*(outdata.u(indata.group2,k));
    outdata.cumtripsStarted(:,k)=outdata.cumtripsStarted(:,k-1)+tripsStarted(:,k); %Cummulative trips started
    
    
    % update service capacity of exit links as a function of accumulation
    % ---------------------------------------------------------------------
    % (exit MFD) - per region 

    if indata.redCapacity == 1
        for i=1:indata.no_reg 

            p = ExitServiceFunction(outdata.x(intersect(indata.group1,indata.reInd{i}),k),indata.capacity(intersect(indata.group1,indata.reInd{i})));
            indata.junct2.lanesd(indG3r{i}) = p*init_ExitLinksLanesR{i}; %per region
        end
    end
    
    % Update turn rates every t_winSteps steps
    % ---------------------------------------------------------------------
    if mod(k,t_winSteps)==0 
        if indata.updateTR == 1 && k<indata.kmax %update every t_winSteps
            t_index = t_index + 1;

            %update turn rates
            %defac_avg = mean(indata.defac(k-t_winSteps+1:k)); %average factor

            % estimate speeds of all links of group1 during the last time window
            % add free flow to all links
            % speeds(:,t_index) = speeds(:,t_index) + indata.v_ff/1000; %km/h - v_ff for all vlinks

            % for actual network links take the average speed of the time window

            is_empty = sum(outdata.x(:,(k-t_winSteps+1):k),2);
            is_empty_ind = find(is_empty==0); %virtual links are always empty
            speeds(is_empty_ind,t_index) = indata.v_ff/1000; %km/h

            non_empty_ind = find(is_empty>0);
            % define cost fundtion per link for non empty links 
            speeds(non_empty_ind,t_index) = sum(outdata.u(non_empty_ind,(k-t_winSteps+1):k)*indata.DT,2).*indata.Links2(non_empty_ind,3)...
                ./(1000*sum(outdata.x(non_empty_ind,(k-t_winSteps+1):k),2)*indata.DT); %km/hour

            speeds(speeds(:,t_index)>indata.v_ff/1000,t_index) = indata.v_ff/1000;

            ind1 = speeds(:,t_index)<1;
            %ind2 = speeds(:,t_index)>0;
            speeds(ind1,t_index) = 1; %minimum speed to estimate travel time in the link
            speeds(speeds(:,t_index) == 0,t_index) = indata.v_ff/1000; %speeds for upper vlinks
            speeds(2202:end,t_index) = indata.v_ff/1000; %speeds for upper vlinks and ending vlinks (not sure if necessary) 
            if sum(speeds(:,t_index)<=0)>0
                disp('error speeds')
            end

            %         for z=1:size(indata.Links2,1) %all links (also virtual?)
            %             is_empty = sum(x(z,(k-t_winSteps+1):k))==0; %binary for the case of link remaining empty
            %             speeds(z,t_index) = max(3000, (1-is_empty)*min([indata.v_ff, sum(u(z,(k-t_winSteps+1):k)*indata.DT)*indata.Links2(z,3)/(sum(x(z,(k-t_winSteps+1):k))/(t_winSteps+1)*(t_winSteps*indata.DT))])+is_empty*indata.v_ff);
            %         end

            %1. First step: find shortest paths(link to link)

            % make connectivity matrix A with travel times (sparse matrix)
            A1 = indata.A./speeds(:,t_index); %travel time in hours
            if sum(A1<0)>0 
                error('connectivity matrix A contains negative values of time')
            end
            %call function:
            %[connCounters, indata.stack_OD, indata.junct2] = updateTurnAndExitRates(...
            [~, indata.stack_OD, indata.junct2] = updateTurnAndExitRates_1(...
                indata.OD_links_2,indata.stack_OD,indata.junct2,indata.t_win,...
                indata.LinksP, A1, indata.downstrP, t_index,indata.empty_connCounters,mean(outdata.u(:,(k-t_winSteps+1):k),2));
                %indata.LinksP, A1, indata.downstrP, t_index,indata.empty_connCounters,mean(outdata.u(:,(k-t_winSteps+1):k),2));
            
            %update OD
            %indata.demandGroup2 = createDemandGroup2(indata.OD_links,...
            %indata.group2,indata.downstrP,indata.upstrP,connCounters);
            %disp('update turn rates')
            %t_index
        else  
            % no turn rate update -- go to next interval 
            % time interval for turn ratios (for constant turn rates)
            t_index = min(ceil((k)*indata.DT/indata.t_win),size(indata.junct2.turn,2));  %turnings change every t_win [hours] / Keeps the last column if the run is longer
        end
    end
    %----------------------------------------------------------------------------------
    
    
    if mod(k,indata.c_int)==0 %if we are at the end of a control cycle
        k_c = k_c + 1; %index of past control cycle (for all calculations of the step)
        
        %print status of simulation run
        %fprintf('Completed: %8.2f \n ' ,k_c/ceil(indata.kmax/indata.c_int)*100)
        
        
        % Calculate heterogeneity for every node = variance of all incoming queues (w or
        % x) = mean of the variance of incoming links (over time)
        for i=1:size(indata.MP.nodeID,1)
            %outdata.node_heterog(i,k_c) = mean(var(outdata.w(indata.junct.or_index(indata.MP.approaches{i}),k-90+1:k)));
            outdata.node_heterog(i,k_c) = mean(var(outdata.x(indata.junct.or_index(indata.MP.approaches{i}),k-90+1:k)));
        end
        
        % Calculate average queue at the end of every control cycle (should
        %be equal to the node signal cycle to be correct)
        for i=1:length(indata.junct.origin)
            outdata.avg_q_cycle(indata.junct.or_index(i),k_c) = mean(outdata.x(indata.junct.or_index(i),max(1,k-floor(indata.junct.cycle(i)/(3600*indata.DT))):k));
        end
        
        % Perimeter Control: Apply the controller after every cycle and calculate
        % new green times to be applied for the next cycle (when each
        % signal reaches the end of its own cycle)
        if PC.mode == 1
            
            %Measure current aggregated accumulations of regions at time step k (excluding VQs outside the network)
            for i=1:indata.no_reg
                %the sum of queues of all region links at the end of the control cycle (snapshot)
                agg_n(i,k_c) = sum(outdata.x(intersect(indata.group1,indata.reInd{i}),k));
            end
            
            % --- Network accumulation for 1-region PC ----
            % agg_ng(k_c) = sum(outdata.x(indata.group1,k));
            
            % Check/update activation/deactivation of controller for the
            %next cycle - three regions
            for i=1:indata.no_reg
                actCheck(i) = -(agg_n(i,k_c)-PC.n_set(i))/PC.n_set(i);
            end
            %actCheck
            %one region PC (network)
            %actCheck = -(agg_ng(k_c)-PC.n_set)/PC.n_set
            
            if activFlag(max(1,k_c-2):k_c) == 0 
                %controller currently inactive (for all regions)
                
                %activFlag(k_c+1:end) = sum(actCheck < PC.actCrit)>0; %activated when at least one is x% away from the target
                activFlag(k_c+1:end) = sum(actCheck < PC.actCrit)>=2; %at least 2 regions reaching the target
                %activFlag(k_c+1:end) = actCheck < PC.actCrit; %(one region)
            else
                %controller currently active
                
                activFlag(k_c+1:end) = ~(sum(actCheck > PC.deactCrit)==3); %deactivate when all regions are more than x% far from target
                %activFlag(k_c+1:end) = ~(actCheck > PC.deactCrit); %(one region)
            end
            
            
            % Formula for only external PC: agg_ug (from 0 to 1) is the % of
            % the saturation flow of the virtual links (gates) - one region
            % agg_ug(k_c+1) = applied_ug(k_c) + PC.K_P*((agg_ng(k_c) - agg_ng(max(1,k_c-1)))./max_ntw) + ...
            %    PC.K_I*((agg_ng(k_c)- PC.n_set)./PC.n_set);
            
            % Proportional-Integral Controler: agg_u are the greens for
            % every intersection of each interregional approach (same time)
            % Formula with external perimeter gating and interregional PC
            agg_u(:,k_c+1) = applied_u(:,k_c) + PC.K_P*((agg_n(:,k_c) - agg_n(:,max(1,k_c-1)))./indata.max_n') + ...
                PC.K_I*((agg_n(:,k_c)- PC.n_set)./PC.n_set);
            
            
            
            %% actions taken by PC 
            if activFlag(k_c+1) || sum(activFlag(max(1,k_c-1):k_c+1))>0 % if controller is activated for next cycle
                %disp('controller active')
                
                % --apply control results to all external gates (one region PC)
                %flagact = 1;
                %indata.junct2.lanesu(extGlobalGateConnections) = min(max(7/90, agg_ug(k_c+1)),1);
                
                % --external PC: apply control results to all VQs (external and internal)
                for i=1:indata.no_reg
                    %intended new value
                    val = agg_u(indata.no_adjReg+i,k_c+1)-applied_u(indata.no_adjReg+i,k_c);
                    check = abs(val)>indata.gn_R; %check how large is the change 
                    
                    if check
                        %too big change ->
                        if val>0
                            %increase
                            applied_u(indata.no_adjReg+i,k_c+1) = min(100,max(minSatFlowGating,applied_u(indata.no_adjReg+i,k_c) + indata.gn_R));
                        else
                            %decrease
                            applied_u(indata.no_adjReg+i,k_c+1) = min(100,max(minSatFlowGating,applied_u(indata.no_adjReg+i,k_c) - indata.gn_R));
                        end
                    else
                        %change within the desired limits
                        applied_u(indata.no_adjReg+i,k_c+1) = min(100,max(minSatFlowGating,agg_u(indata.no_adjReg+i,k_c+1)));
                    end
                    
                    %apply gating as adjustment of no of lanes upstream
                    indata.junct2.lanesu(PC.extGateConnectionsReg{i}) = applied_u(indata.no_adjReg+i,k_c+1)/100*PC.initialGateLanes(PC.extGateConnectionsReg{i});
                end
                
                
                % -- internal PC: apply control results to interregional gates + added external
                greens_P = {};
                % for i=1:size(agg_u,1)  %for every interregional approach i-j + selected external gates
                for i=1:(size(agg_u,1)-3)  %for every interregional approach i-j I solve the feasibility problem
                    
                    ind = PC.indices(i,PC.indices(i,:)>0); %all nodes that belong to this approach i-j (non-zero) - the k indices to refer to each node in PC structs
                    for j=1:length(ind)
                        greens_P{j} = PC.stageDur{ind(j)}; %the greens of all the phases of the set of nodes controlling the i-j
                    end
                    
                    %the twin nodes
                    twin_ind = PC.indicesco(i,PC.indices(i,:)>0); %all twin nodes corresp. to the nodes before
                    
                    %[newGreens_1, newGreens_2] = findFeasibleSignalsPC(agg_u(i,k_c+1), greens_P, outdata.avg_q_cycle(:,k_c), indata.junct.or_index, indata.gn_R, PC, ind, twin_ind, indata.minStageDur);
                    [newGreens_1, newGreens_2] = findFeasibleSignalsPC_3(agg_u(i,k_c+1), greens_P, outdata.avg_q_cycle(:,k_c), indata.junct.or_index, indata.gn_R, PC, ind, twin_ind, indata.minStageDur,indata.satFlow,indata.alpha, indata.beta);

                    greensToApply(ind,1) = newGreens_1; %main phase
                    greensToApply(ind,2) = newGreens_2; %secondary phase
                    
                    %the twin nodes take the exact same green times of the
                    %main nodes
                    greensToApply(twin_ind(twin_ind>0),1) = newGreens_1(twin_ind>0);
                    greensToApply(twin_ind(twin_ind>0),2) = newGreens_2(twin_ind>0);
                    
                    %update u (applied) based on final greens of the primal stage (for next cycle of controller)
                    applied_u(i,k_c+1) = mean(newGreens_1);
                end
                
            else
                if activFlag(max(k_c-5,1):k_c+1)==0 %inactive for the last 7 steps
                    %disp('controller inactive for long')
                    applied_u(:,k_c+1) = applied_u(:,k_c);
                    
                else
                    %disp('controller inactive')
                    %calculate the applied greens from the current settings
                    greens_P = [];
                    for i=1:indata.no_adjReg
                        ind = PC.indices(i,PC.indices(i,:)>0); %all nodes that belong to this approach i-j (non-zero) - the k indices to refer to each node in PC structs
                        for j = 1:length(ind)
                            greens_P = [greens_P PC.stageDur{ind(j)}(PC.stagesInvolved(ind(j),1))]; %the greens of the main phase of all nodes controlling the i-j
                        end
                        applied_u(i,k_c+1) = mean(greens_P); %set the initial values for the controller
                    end
                    applied_u(indata.no_adjReg+1:indata.no_adjReg+3,k_c+1) = 100; % sat flow to 100% 
                    indata.junct2.lanesu = min(PC.initialGateLanes,indata.junct2.lanesu + 0.2*PC.initialGateLanes); %gradual increase of sat flow 
                    %indata.junct2.lanesu = PC.initialGateLanes;
                end

            end
        end

    end
    
    timeInCycle = rem((indata.MP.offset + clock)',(indata.applyMPfreq*indata.MP.cycle)');
    
     
    if ~isempty(MPnodes) %if max pressure is applied 

        % Update traffic signal settings for MP + PC intersections every x cycles
        % (offset is taken into account)
        
        to_update = (timeInCycle(MPnodes)==0); %result: index of nodes in MPnodes that are at the end of the update period (x cycles)
        
        if sum(to_update)>0 %a least 1 node is updated (is at the end of the cycle)
            to_update = MPnodes(to_update); % result: index of nodes in MP
            
           [outdata.greentimes2(:,:,max(1,k_c)),indata.MP] = maxPressure(to_update,indata.MP,indata.junct2,indata.capacity,indata.satFlow,indata.greentimes2,outdata.x,k,t_MP,t_index,indata.gn_R);
            %[outdata.greentimes2(:,:,k_c),indata.MP] = maxPressure(to_update,indata.MP,indata.junct2,indata.capacity,indata.satFlow,indata.greentimes2,outdata.w,k,t_MP,t_index,indata.gn_R);
            indata.greentimes2 = outdata.greentimes2(:,:,max(1,k_c));
            
        end
    end
    
    % For Perimeter Control (part 2/2)-------------------------------------
    if PC.mode == 1
        % Apply updates to traffic signals of PC intersections that are at
        % the end of their cycle 
        
        to_update = (timeInCycle(PC.uniqueNodes)==0); 
        if sum(to_update)>0
            to_update = PC.uniqueNodes(to_update); %indices in MP
   
            if activFlag(k_c+1)
                %Controller currently active - Apply the calculated greens at the intersections
                for i = 1:length(to_update) %for every intersection that reached the end
                    
                    %find index in PC struct
                    ind = find(PC.nodeID == indata.MP.nodeID(to_update(i)));
                    
                    % update stage durations in PC and MP structures (both
                    % are necessary)
                    PC.stageDur{ind}(PC.stagesInvolved(ind,1:2)) = greensToApply(ind,:);
                    % required:
                    indata.MP.stageDur{to_update(i)}(PC.stagesInvolved(ind,1:2)) = greensToApply(ind,:);
                    
                    %update greentimes - for every approach invlovled:
                    %all approaches were affected
                    appr = indata.MP.approaches{to_update(i)}; %all approaches of the intersection
                    for j=1:length(appr)
                        s = indata.junct.stageNo(appr(j),indata.junct.stageNo(appr(j),:)>0); %stages where the approach takes green
                        for ss=1:length(s)
                            indata.greentimes2(appr(j),ss*2) = sum(indata.MP.stageDur{to_update(i)}(1:s(ss)));  %end time of respective stages
                            indata.greentimes2(appr(j),ss*2-1) = sum(indata.MP.stageDur{to_update(i)}(1:s(ss)-1)); %start time of repsective stages
                        end
                    end
                    flag_modified_signals(to_update(i)) = 1;
                end
            else
                % Controller is deactivated (last cycle) - return signals
                % to pre-timed settings if they are not already there (to make incremental?)
                
                for i = 1:length(to_update)
                    if flag_modified_signals(to_update(i))==1
                        ind = find(PC.nodeID == indata.MP.nodeID(to_update(i)));
                        appr = indata.MP.approaches{to_update(i)}; %all approaches of the intersection
                        
                        %immediately return to initial settings
                        indata.greentimes2(appr,:) = init_greentimes(appr,:);
                        PC.stageDur{ind} = init_stageDur{ind};
                        indata.MP.stageDur{to_update(i)} = init_stageDurMP{to_update(i)};
                        %label as reconciled to initial settings
                        flag_modified_signals(to_update(i)) = 0;
                    end
                end                 
            end
        end
    end
end

outdata.notserviced = sum(outdata.x(:,end))+sum(outdata.w(indata.group2,end));
outdata.virtualqueues = outdata.w(indata.group2,:);
outdata.MPnodes = MPnodes;
outdata.PC = PC;
outdata.turn = indata.junct2.turn; %the resulting turn rates
if PC.mode==1
    %gather results related to PC 
    outdata.agg_u = agg_u;
    outdata.agg_n = agg_n;
    outdata.applied_u = applied_u;
    outdata.activFlag = activFlag;
    outdata.PC_K_I = PC.K_I;
    outdata.PC_K_P = PC.K_P;
    outdata.SP = PC.n_set;
end
% Calculation of Total Travel Time 
r1 = sum(sum(outdata.x*indata.DT)); %time in network 
r2 = sum(sum(outdata.w(indata.group2,:)*indata.ksi*indata.DT)); %time in VQs
r3 = r1 + r2 %total TT (veh x hours)
outdata.r3 = r3;
outdata.r2 = r2; 
outdata.r1 = r1; 
toc
% generating/saving results file 
if savefull==1
    save(strcat(fname),'outdata','indata', 'MPnodes','PC','-v7.3');
    save(fname,'outdata','indata', 'MPnodes','PC','-v7.3');
elseif savefull==2
    %light save (only TTT): 
    notserviced = outdata.notserviced; 
    save(fname,'r1','r2','r3','notserviced');
end

return