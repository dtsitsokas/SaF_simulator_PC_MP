%Script/or function: Prepare input file for exp settings (per case)

function [] = prepareExperiment(scenario,PC_input_filename,PCmode,MPmode)
load parameters.mat
load FinalInput.mat 

%Set caseName for specific Max-Pressure Case (selection of nodes must be run first)
caseName = scenario; 

p_1 = 0; 
lim_var = {}; 
lim_occ = {}; 
cong_threshold = 0; 

%% Perimeter control settings 
PC.mode = PCmode; %on/off for the specified intersections 
if PC.mode >= 1 && init_clustering == 1 %ONLY used for static PC with initial clustering 
    
    %Info of the controlled intersections for PC  
    %[PC] = PCcontrolledIntersectionsInfo(PC);
    [PC] = PCcontrolledIntersectionsInfo_2(PC); %gating to all incoming inflow from VQs (external and internal)  
    
    %Find node indices within MP - keep same structure (connect the 2
    %structures) 
    PC.nodes = PC.junctionsID; %initialize
    for i=1:size(PC.junctionsID,1)
        for j=1:size(PC.junctionsID,2)
            if PC.junctionsID(i,j)>0
                PC.nodes(i,j) = find(MP.nodeID == PC.junctionsID(i,j)); %indices of controlled intersections in MP structure
            end
        end
    end
    
    %Find twin node (co) indices to MP - keep same structure 
    PC.nodesco = PC.junctionsIDco; 
    for i=1:size(PC.junctionsIDco,1)
        for j=1:size(PC.junctionsIDco,2)
            if PC.junctionsIDco(i,j)>0
                PC.nodesco(i,j) = find(MP.nodeID == PC.junctionsIDco(i,j)); %indices of twin controlled intersections in MP structure        
            end
        end
    end
    
    % exclude PC intersections from MaxPressure set (add them to specialIntersections)
    specialIntersections = [specialIntersections unique(PC.junctionsID(PC.junctionsID>0))' unique(PC.junctionsIDco(PC.junctionsIDco>0))'];
    
    %Copy info from the MP structure to create PC struct only for the controlled
    %intersections (all of them - main and co's)
    k=0; 
    PC.indices = zeros(size(PC.nodes)); 
    PC.indicesco = zeros(size(PC.nodesco));
    PC.stageDur = {};
    PC.stageBaseDur = {}; 
    PC.stageApproaches = cell(39,size(MP.stageApproaches,2));
    %PC.stagesInvolved = {};
    
    for i=1:size(PC.nodes,1)
        for j=1:size(PC.nodes,2)
            if PC.nodes(i,j)>0
                k = k + 1; %index in PC structure 
                PC.indices(i,j) = k; %index to associate node(i,j) to a PC record k 
                PC.stageDur{k} = MP.stageDur{PC.nodes(i,j)};
                PC.initStageDur{k} = MP.stageDur{PC.nodes(i,j)};
                PC.stageBaseDur{k} = MP.stageBaseDur{PC.nodes(i,j)};
                for ll=1:MP.noOfStages(PC.nodes(i,j))
                    PC.stageApproaches{k,ll} = MP.stageApproaches{PC.nodes(i,j),ll};
                end
                PC.nodeID(k) = MP.nodeID(PC.nodes(i,j)); 
                PC.noOfStages(k) = MP.noOfStages(PC.nodes(i,j)); 
                PC.cycle(k) = MP.cycle(PC.nodes(i,j)); %or PC.cyclePC(i,j) - check 
                PC.offset(k) = MP.offset(PC.nodes(i,j)); %or PC.offsetPC(i,j) - check 
                PC.stagesInvolved(k,1:2) = [PC.stagePC(i,j) PC.stage2PC(i,j)]; %pair of stages that change with PC 
                PC.timeToAllocate(k) = PC.sum_greensPC(i,j); %or sum(MP.stageDur{PC.nodes(i,j)}(PC.stagesInvolved(k,:))) %check this with my settings 
                %PC.upRegion(k) = regionsUpDown(i,1); 
                %PC.downRegion(k) = regionsUpDown(i,2);
            end
        end
    end

    %Add the twin nodes to PC struct (at the end) 
    for i=1:size(PC.nodesco,1)
        for j=1:size(PC.nodesco,2) 
            if PC.nodesco(i,j)>0
                k = k + 1; 
                PC.indicesco(i,j) = k; 
                PC.stageDur{k} = MP.stageDur{PC.nodesco(i,j)};
                PC.initstageDur{k} = MP.stageDur(PC.nodesco(i,j));
                PC.stageBaseDur{k} = MP.stageBaseDur(PC.nodesco(i,j));
                for ll=1:MP.noOfStages(PC.nodesco(i,j))
                    PC.stageApproaches{k,ll} = MP.stageApproaches{PC.nodesco(i,j),ll};
                end
                PC.nodeID(k) = MP.nodeID(PC.nodesco(i,j)); 
                PC.noOfStages(k) = MP.noOfStages(PC.nodesco(i,j)); 
                PC.cycle(k) = MP.cycle(PC.nodesco(i,j)); %or PC.cyclePC(i,j) - check 
                PC.offset(k) = MP.offset(PC.nodesco(i,j)); %or PC.offsetPC(i,j) - check 
                PC.stagesInvolved(k,1:2) = [PC.stagecoPC(i,j) PC.stage2coPC(i,j)]; %pair of stages that change with PC 
                PC.timeToAllocate(k) = PC.sum_greenscoPC(i,j); %or sum(MP.stageDur{PC.nodesco(i,j)}(PC.stagesInvolved(k,:))) %check this with my settings 
                %PC.upRegion(k) = regionsUpDown(i,1); 
                %PC.downRegion(k) = regionsUpDown(i,2);
            end
        end
    end

    %% PC PARAMETERS
    % [EITHER] set PC parameters manyally 
    %
  %   PC.K_I = [15   -10     0;
  %    0    -5    10;
  %  -15    10     0;
  %    0     5   -10;
  %  -20     0     0;
  %    0   -20     0;
  %    0     0   -10];
  % 
  %   % 
  %   % 
  %    PC.K_P = [150  -100     0;
  %    0   -50   100;
  % -150   100     0;
  %    0    50  -100;
  % -200     0     0;
  %    0  -200     0;
  %    0     0  -100];
  %   % 
  %   PC.n_set = [8000 8000 8000]';
  %   PC.actCrit = 0.0;
  %   PC.deactCrit = 0.15;
  %   PC.minRegAct = 1; 
  %   PC.maxRegAct = 3; 

    % [OR] load the already set parameters of an existing PC scheme (input file)
    load(PC_input_filename,'PC')
    PC.minRegAct = 2; 
    PC.maxRegAct = 3; %deactivate when all regions fulfill the criterion
    PC.mode = PCmode;
end

%% Max-Pressure settings - translate selected node IDs to MPnodes input

% Scenario: ALL/NO MP nodes  --- SET SET SET !!!!
if MPmode == 100
    % ""-----Case all nodes----------""
    MP_nodeIDs = MP.nodeID(MP.splan==1);  %all controlled nodes
elseif MPmode == 0 
    % ----Case no nodes (no MP)---------
    MP_nodeIDs = []; %no nodes 
elseif MPmode == 1 
    % ----Case specific nodes / Load specific nodes for MP
    load(strcat('MPnodesCase_',caseName,'.mat'))
    MP_nodeIDs = MP.nodeID(MPnodes);
elseif MPmode == 2
    %for random MP nodes + PC 
    load(caseName,'rMPnodes')
    MPnodes = rMPnodes; 
    clear rMPnodes 
    MP_nodeIDs = MP.nodeID(MPnodes);
else 
    error('wrong MPcode')
end

%% ""----Exclude non-eligible nodes and PC nodes --------""
if ~isempty(MP_nodeIDs) && init_clustering > 0 
    %exclude non-eligible nodes and PC intersections if PC is applied!
    if PC.mode >= 1 
        check = intersect(intersect(MP_nodeIDs,specialIntersections), PC.nodeID); 
    else
        check = intersect(MP_nodeIDs,specialIntersections);
    end
    if ~isempty(check)
        disp('removing non-eligible intersections') 
        MP_nodeIDs = setxor(MP_nodeIDs,check);
    end
    %reconstruct MPnodes 
    MPnodes = zeros(1,length(MP_nodeIDs)); 
    for i=1:length(MP_nodeIDs)
        MPnodes(i) = find(MP.nodeID == MP_nodeIDs(i)); %indices
    end
else
    MPnodes = [];
end

save(strcat('input_',scenario,'.mat'),'MPnodes','PC','p_1','lim_var','lim_occ','cong_threshold')