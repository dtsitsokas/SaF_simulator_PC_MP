%Script/or function: Prepare Experiment sets for different cases
function [] = prepareExperimentfunct(scenario, caseName, notes)
% Function to prepare input files for simulation about PC and MP settings 
% 
% scenario = the name for the input file to generate (name of the scenario)
% caseName = the set of MPnodes to load from file (already existing)  
% notes = string with comment about the scenario (optional, otherwise '-') 
% PC settings are set in the code (not used now) 


%% Set Scenario Name here: 
% scenario = 'MP_S1_PC1';
% notes = 'MP_S1 plus PC';
%load settings and initialize network
% parametersetting() %adjusted
load parameters.mat
% Initialization() %adjusted
load FinalInput.mat

%Set caseName for specific Max-Pressure Case (selection must be run first)
% caseName = 'MP_S1';  %MP case to load on the scenario 
% caseName = scenario; 
p_1 = 0; 
lim_var = {}; 
lim_occ = {}; 
cong_threshold = 0; 

%% Perimeter control settings 
PC.mode = 0; %on/off for the specified intersections 
if PC.mode == 1
    
    %Info of the controlled intersections for PC  
    %[PC] = PCcontrolledIntersectionsInfo(PC);
    [PC] = PCcontrolledIntersectionsInfo_2(PC); %gating to all incoming inflow from VQs - equally(?) 
    
    %Find node indices to MP - keep same structure 
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
    
%% --- Proportional-Integral Gain parameters K_I, K_P (global controller) 
% %     % --- from parameters
%     %parameters = ones(47,1);
%     %parameters(43:47,1) = [5000 6500 4000  0.15  0.25]';
%     %parameters = [9.2110080e-01   6.5326093e-01   4.5044572e-01   1.1005502e+00   5.8095665e-01   8.3534193e-01   1.1418881e+00   1.0244275e+00   7.2366215e-01   1.0566756e+00   1.0789376e+00   5.1380544e-01   9.3231970e-01   8.8262665e-01   4.8386712e-01   8.1053186e-01   1.3465794e+00   6.3751672e-01   5.7424690e-01   9.8463376e-01   1.0762841e+00   7.8532532e-01   7.4664950e-01   9.7925799e-01   1.3097016e+00   7.7303944e-01   1.0610478e+00   8.8999990e-01   1.1398435e+00   1.0101732e+00   1.1024765e+00   1.0348122e+00   7.9762243e-01   9.2579577e-01   1.0831839e+00   1.0480795e+00   8.1107085e-01   7.5673445e-01   1.0142730e+00   8.0657809e-01   1.3137070e+00   9.0596967e-01   2.1293285e+03   2.8127596e+03   1.4611331e+03   1.6883504e-01   3.1250000e-01]'; %subsets ? 
%     %parameters = [1.1251892e+00   9.1645699e-01   3.3674671e-01   1.3308482e+00   7.8408147e-01   1.1760306e+00   9.7945352e-01   5.0778594e-01   6.5812856e-01   4.5483242e-01   1.1086563e+00   8.1055917e-01   3.9442663e-01   1.2063737e+00   1.1262555e+00   1.4138117e+00   8.7772352e-01   8.1157415e-01   1.2810296e+00   5.6051597e-01   7.6916821e-01   2.6379003e-01   9.8388052e-01   7.4175040e-01   9.3353958e-01   6.0431712e-01   9.7289337e-01   1.0528996e+00   1.0088231e+00   9.5482561e-01   4.3581690e-01   1.1781879e+00   5.4274946e-01   6.0540530e-01   6.7500336e-01   9.5029478e-01   1.1926906e+00   1.3667023e+00   9.5724484e-01   7.4592112e-01   1.3469107e+00   6.7056988e-01   6.2500000e+03   7.5040291e+03   3.3297379e+03   1.8750000e-01   2.8701089e-01]'; %no subsets (best found so far)
%     
%     %parameters = [1.1251892e+00   9.1645699e-01   3.3674671e-01   1.3308482e+00   7.8408147e-01   1.1760306e+00   9.7945352e-01   5.0778594e-01   6.5812856e-01   4.5483242e-01   1.1086563e+00   8.1055917e-01   3.9442663e-01   1.2063737e+00   1.1262555e+00   1.4138117e+00   8.7772352e-01   8.1157415e-01   1.2810296e+00   5.6051597e-01   7.6916821e-01   2.6379003e-01   9.8388052e-01   7.4175040e-01   9.3353958e-01   6.0431712e-01   9.7289337e-01   1.0528996e+00   1.0088231e+00   9.5482561e-01   4.3581690e-01   1.1781879e+00   5.4274946e-01   6.0540530e-01   6.7500336e-01   9.5029478e-01   1.1926906e+00   1.3667023e+00   9.5724484e-01   7.4592112e-01   1.3469107e+00   6.7056988e-01   6.2500000e+03   7.5040291e+03   3.3297379e+03   1.8750000e-01   2.8701089e-01]'; 
% %     %parameters = [9.0567031e-01   1.5694316e+00   7.6927212e-01   1.1791133e+00   4.7367938e-01   2.3895216e-01   8.8876022e-01   6.0111213e-01   9.6708310e-01   9.0304682e-01   7.7309484e-01   2.0124903e-01   8.6780618e-01   5.1490107e-01   1.9730401e+00   6.8165001e-01   7.7401779e-01   1.0782205e+00   6.4182063e-01   1.4261640e+00   8.1453834e-01  -3.5555716e-02   7.7816407e-01   3.4511886e-01   1.0059375e+00   1.6472033e+00   6.2836676e-01   3.8267434e-01   2.9813072e-01   8.6374152e-02   3.0308717e+00   8.3670984e-01   2.8864197e-01   2.6609864e-02   3.3228824e-01   1.0819999e+00   6.1131247e-01   1.1652637e+00  -2.0040562e-01   9.9734905e-01   1.2816058e+00  -4.9201253e-02   2.0355293e+03   1.4422175e+03   1.9606215e+03   2.2771106e-01   8.5967289e-02]';
% %     %parameters = [1.9080798e+00  -1.3428721e-01   2.4925412e+00   1.6471137e+00   3.4211361e+00   1.7418989e+00   1.8736754e+00   1.1192957e+00   1.5161829e+00   2.4925964e-01   9.2319690e-01   1.0622088e+00   1.9150107e+00   2.1023381e+00   8.7883963e-02   2.4557394e+00   1.4678396e+00   9.8024217e-02   4.6986313e-01   9.9440411e-01   5.3585555e-01   4.6724066e-01   1.4342693e+00   7.7501210e-01   4.1913436e-01   9.4336947e-01   1.6849691e+00   5.9238177e-01   2.4251059e+00   1.3623731e+00   6.4936014e-01   1.1605556e-02   1.2174199e+00   2.6168012e-01   2.0615708e+00   1.8303445e+00   2.6402056e+00   1.6703403e+00   3.3465182e-01   6.2706118e-01   1.1731729e+00   1.1137247e+00   4.8016533e+03   4.0723276e+03   3.4882558e+03  -1.5686446e-02   2.3410948e-01]';
%     %parameters = [5.8348679e-01   1.0269476e+00   8.9066552e-01   5.0801158e-01   1.5217774e+00   1.7687969e+00   1.2535970e+00   2.8255784e+00   1.2339931e+00   9.1759815e-01   2.4414568e+00   2.4397228e-01   1.2816550e+00   2.3199623e-02   2.0790280e+00   1.1293552e+00   2.6295926e+00   1.9849184e+00   1.0267630e+00   1.1984403e+00   1.4420583e+00   4.5164438e-01   2.4692428e-01   1.5627594e+00   3.8951437e-01   1.2587932e+00   9.4360091e-01   1.5329790e+00   4.7063297e-01   9.4584470e-01   6.4491589e-01   1.0634245e-01   1.2726966e+00   1.0277599e+00   1.6794696e+00   2.4974332e+00   1.7189575e+00   2.0947008e+00   6.9246993e-01   2.1964714e-01   3.5005806e+00   3.1258210e+00   3.1250000e+03   6.7757794e+03   1.7790491e+03   8.6182661e-02   1.3101308e-01]';
%     %parameters = [6.9297962e-01   1.7571990e+00   7.2362690e-01   1.2859344e+00   4.2231725e-01   2.2094187e-01   7.7219378e-01   6.0861685e-01   6.5182964e-01   1.0520675e+00   4.9827343e-01   1.8108151e-01   8.5296675e-01   4.5489650e-01   2.1238457e+00   8.4491062e-01   7.8585442e-01   1.0065919e+00   6.3240968e-01   1.4522051e+00   6.6870451e-01  -2.8218581e-02   8.4507176e-01   3.6814072e-01   1.2299586e+00   1.4903749e+00   6.7239128e-01   3.5175194e-01   2.6529400e-01   8.3966571e-02   1.6908953e+00   4.7759785e-01   3.3472650e-01   1.5828166e-02   2.3044878e-01   9.9824150e-01   4.1806711e-01   1.1330998e+00  -2.0256653e-01   1.0331374e+00   1.3378911e+00  -4.2736969e-02   1.7863107e+03   1.0821360e+03   2.4156146e+03   1.8037529e-01   2.3970407e-02]';
%     %parameters = [5.8348679e-01   1.0269476e+00   8.9066552e-01   5.0801158e-01   1.5217774e+00   1.7687969e+00   1.2535970e+00   2.8255784e+00   1.2339931e+00   9.1759815e-01   2.4414568e+00   2.4397228e-01   1.2816550e+00   2.3199623e-02   2.0790280e+00   1.1293552e+00   2.6295926e+00   1.9849184e+00   1.0267630e+00   1.1984403e+00   1.4420583e+00   4.5164438e-01   2.4692428e-01   1.5627594e+00   3.8951437e-01   1.2587932e+00   9.4360091e-01   1.5329790e+00   4.7063297e-01   9.4584470e-01   6.4491589e-01   1.0634245e-01   1.2726966e+00   1.0277599e+00   1.6794696e+00   2.4974332e+00   1.7189575e+00   2.0947008e+00   6.9246993e-01   2.1964714e-01   3.5005806e+00   3.1258210e+00   3.1250000e+03   6.7757794e+03   1.7790491e+03   8.6182661e-02   1.3101308e-0]';
%     %parameters = [6.9297962e-01   1.7571990e+00   7.2362690e-01   1.2859344e+00   4.2231725e-01   2.2094187e-01   7.7219378e-01   6.0861685e-01   6.5182964e-01   1.0520675e+00   4.9827343e-01   1.8108151e-01   8.5296675e-01   4.5489650e-01   2.1238457e+00   8.4491062e-01   7.8585442e-01   1.0065919e+00   6.3240968e-01   1.4522051e+00   6.6870451e-01  -2.8218581e-02   8.4507176e-01   3.6814072e-01   1.2299586e+00   1.4903749e+00   6.7239128e-01   3.5175194e-01   2.6529400e-01   8.3966571e-02   1.6908953e+00   4.7759785e-01   3.3472650e-01   1.5828166e-02   2.3044878e-01   9.9824150e-01   4.1806711e-01   1.1330998e+00  -2.0256653e-01   1.0331374e+00   1.3378911e+00  -4.2736969e-02   1.7863107e+03   1.0821360e+03   2.4156146e+03   1.8037529e-01   2.3970407e-02]';
%     parameters = [6.9583963e-01   4.1324841e-01   2.6724263e-01   2.3186051e-01   7.7690061e-01   2.4773216e+00   9.6243602e-01   4.1933223e+00   5.8277174e-01   1.5772142e-01   3.4706498e-01   8.7125082e-02   4.1683486e-01  -1.0156448e-02   5.5478306e-01   7.9393188e-01   3.4303230e+00   1.0010728e+00   2.4460983e-01   8.2160230e-01   1.3467454e+00   1.2105935e-01   1.0382620e-01   2.5037531e-01   1.1537877e-01   1.1336600e+00   5.0655002e-01   1.3931855e+00   1.5335917e-01   2.2049652e-01   2.0582329e-01   1.1634151e-01   4.1843305e-01   5.2361831e-01   1.7824729e+00   1.4219845e+00   5.3119009e-01   1.7873625e+00   2.2359919e-02   7.2436609e-02   2.7117574e+00   7.8480316e-01   1.3353067e+03   5.6291586e+03   5.0206712e+02   9.8873541e-02  -3.9961034e-04]';
%     s_setP = 3; %number of regions
%     s_kp_ki = [7 3];
%     PC.K_P = zeros(s_kp_ki);
%     PC.K_I = zeros(s_kp_ki);
%     
%     %setpoint
%     %PC.n_set = transpose(parameters(end-s_setP-1:end-2));
%     PC.n_set = parameters(end-s_setP-1:end-2);
%     %parameters
%     k = 1;
%     for i=1:s_kp_ki(1)
%         PC.K_I(i,:) = parameters(k:k+s_kp_ki(2)-1);
%         PC.K_P(i,:) = parameters(k+s_kp_ki(1)*s_kp_ki(2):k+s_kp_ki(1)*s_kp_ki(2)+s_kp_ki(2)-1);
%         k = k + s_kp_ki(2);
%     end
%     
%     %activation criterion
%     %from 0 to 0.20 (activ when at least one region is less than 20% of the setpoint)
%     PC.actCrit = parameters(end-1);
%     PC.deactCrit = parameters(end);
% %     
   %% ----- 
    
% % 
    %PC.K_I = zeros(7,3);
%     PC.K_P = zeros(7,3);
    PC.K_I = [ 10 -10 0;
    0 -10 10;
    -10 10 0;
    0 10 -10;
    -3   0   0;
    0  -3   0;
    0   0  -1];       
    
    PC.K_P = [100 -100 0;
        0 -100 100;
        -100 100 0;
        0 100 -100;
        -50 0 0;
        0 -50 0;
        0 0 -20];
    
    PC.n_set = [7000 7000 5000]';
    PC.actCrit = 0.05;
    PC.deactCrit = 0.10; 
   
    %% -- for external PC only 
    %PC.K_P = -75;
    %PC.K_I = -7;
    %PC.n_set = 21360; 
    
    %PC.nodeInfo{i}() %(nodeID, phase_to_change, upstream_area, downstream_area)    
    
    % i = PC.node,    
    
    
    %% --- PC parameters - independent controllers 
%     
%     PC.iK_P = cell(3,1);  
%     PC.iK_P{1} = [-100 -100]'; 
%     PC.iK_P{2} = [-100 -100 -100]';
%     PC.iK_P{3} = [-100 -100]'; 
%     PC.iK_I = cell(3,1); 
%     PC.iK_I{1} = [-35 -15]'; 
%     PC.iK_I{2} = [-45 -7.5 -7.5]';
%     PC.iK_I{3} = [-25 -15]';
%     PC.actCrit = 0.25;
%     PC.deactCrit = 0.25;
%     PC.n_set = [8000  9000  5600]';
    %% ---
   
    load('FinalInput.mat','MP')
    
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
    
    
    %Min and Max values for the u - Don't need, I apply minimum green for
    %every phase/stage 
    %PC.u_min = 0.2; 
    %PC.u_max = 1.0; 
        
end


%% Define Nodes for Max-Pressure 
% scen = 3; %Base scenario for the classification of nodes 
% [MPnodes,p_1,lim_var,lim_occ,cong_threshold] = nodeSelectionMP(scen); %TO FIX/REPAIR when I implement the combined scenario (inputs missing)


%% Max-Pressure settings - translate selected node IDs to MPnodes input

% Scenario: ALL/NO MP nodes  --- SET SET SET !!!!

% ""-----Case all nodes----------""
%MP_nodeIDs = MP.nodeID(MP.splan==1);  %all controlled nodes

%----Case no nodes (no MP)---------
MP_nodeIDs = []; %no nodes 

%""----Exclude non-eligible nodes and PC nodes --------""
if ~isempty(MP_nodeIDs)
    %exclude non-eligible nodes and PC intersections if PC is applied!
    if PC.mode == 1
        check = intersect(intersect(MP_nodeIDs,specialIntersections), PC.nodeID); 
    else
        check = intersect(MP_nodeIDs,specialIntersections);
    end
    if ~isempty(check)
        disp('removing non-eligible intersections') 
        MP_nodeIDs = setxor(MP_nodeIDs,check);
    end
    
    MPnodes = zeros(1,length(MP_nodeIDs)); 
    for i=1:length(MP_nodeIDs)
        MPnodes(i) = find(MP.nodeID == MP_nodeIDs(i)); %indices
    end
else
    MPnodes = [];
end

%----Case specific nodes / Load specific nodes for MP 
if ~isempty(caseName)
    %----Case specific nodes / Load specific nodes for MP 
    load(strcat('MPnodesCase_',caseName,'.mat'))
end

save(strcat('input_',scenario,'.mat'),'MPnodes','PC','notes','p_1','lim_var','lim_occ','cong_threshold')