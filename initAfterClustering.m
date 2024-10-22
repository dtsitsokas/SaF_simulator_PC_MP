function [indata,PC,indG3r,init_ExitLinksLanesR,agg_n,agg_u,actCheck,greensToApply,applied_u, init_greentimes,init_stageDurMP] = initAfterClustering(indata,PC)
% initialize for PC 

% -- regional accumulations
% three-region PC
agg_n = zeros(indata.no_reg,ceil(indata.kmax/indata.c_int));                  % aggregated regional accumulations per control cycle

%only external perimeter (all nodes)
%agg_ng = zeros(1,ceil(indata.kmax/indata.c_int));

% -- control variables
% PC without external perimeter (only interregional)
% agg_u = zeros(indata.no_adjReg,ceil(indata.kmax/indata.c_int));              % mean green time calculated by the controller per cycle

% PC with external perimeter + selected nodes
agg_u = zeros(indata.no_adjReg+indata.no_reg,ceil(indata.kmax/indata.c_int));              % mean green time calculated by the controller per cycle

% only external perimeter PC (as in one region network)
% agg_ug = zeros(1,ceil(indata.kmax/indata.c_int));


actCheck = zeros(1,indata.no_reg);                                                       % Initialize activation check for the 3 regions
greensToApply = zeros(length(PC.nodeID),2);                                  % Green settings ready to apply to the main approach of PC nodes when their cycle is completed
applied_u = zeros(indata.no_adjReg+indata.no_reg,ceil(indata.kmax/indata.c_int));        % The applied values of u after projection to feasible space
%applied_ug = zeros(1,ceil(indata.kmax/indata.c_int));
init_greentimes = indata.greentimes2;                                        % The fixed-time signal settings initially given as input                                              % Initial fixed-time settings (in form of duration of phases) - same info
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

PC.uniqueNodes = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];


% PC External 
% --- PC at all VQs in the sense of an inflow control (same for all VQs
% of every region) - can happen externaly (move?)
%initialize/store initial no of lanes for VQs

PC.extGateLinksReg = cell(1,indata.no_reg); %index of VQ virtual links per region
PC.extGateConnectionsReg = cell(1,indata.no_reg); %index of approaches between VQs and ntw links per region (in junct2)
for i=1:indata.no_reg %for every region
    PC.extGateLinksReg{i} = intersect(indata.group2,indata.reInd{i}); %VQ link indices of the region
    %store indices of connections in junct to easily modify the VQ
    %lanes = the VQ saturation flow
    PC.extGateConnectionsReg{i} = ismember(indata.junct2.or_index,PC.extGateLinksReg{i});
end
PC.initialGateLanes = indata.junct2.lanesu;  %storage of all initial values (total saturation flows)

end