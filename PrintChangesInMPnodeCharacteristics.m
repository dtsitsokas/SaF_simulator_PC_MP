% print graphs of change in node characteristics before/after MP
clc
clear
close all

%results of the reference case (NC)
%(medium demand)
% fname = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\output_set5_MedDemand_NC';
% high demand 
fname = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\output_set5_highDemand_NC';

%input file of NC case (same for both demands) 
load('C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\input_NC.mat')

% path for the location of the input files of the examined cases 
% (high demand) 
pathIN = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs high demand\final runs input files\';
% pathIN = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';

%file name themes of Mpnodes files (inputs) 
% filethemeMP0 = 'MPnodesCase_MP_'; 
% filethemeMP_MS = 'MPnodesCase_MP_grid_';

filethemeMP0 = 'MPnodesCase_MP_0'; 
%filethemeMP_MS = 'MPnodesCase_MP_cd_MS';
filethemeMP_MS = 'MPnodesCase_PC_50b_MPf_MS';

%path for location of results files 
% pathTest = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
pathTest = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';

%themes of filenames of result files 
% resfilethemeMP0 = 'output_set5_highDemand_MP_';
% resfilethemeMP_MS = 'output_set5_highDemand_MP_grid_high_eq2t_MS'; 

%resfilethemeMP0 = 'output_set5_MedDemand_MP_';
resfilethemeMP0 = 'output_set5_highDemand_MP_';

%resfilethemeMP_MS = 'output_set5_MedDemand_MP_cd_MS'; 
resfilethemeMP_MS = 'output_set5_highDemand_PC_50b_MPf_MS';

CyclesInPeak = 80; %number of cycles of the peak period (2 hours/ 90 sec = 80 cycles)

%load results of reference case
load(fname,'indata','outdata')

MP = indata.MP; 
num_MPall = sum(MP.splan==1)-length(indata.specialInt); %total number of MP eligible nodes
mrkSize = 15; %size for scattters
colMP = 'b'; %color for MP nodes
colNonMP = 'g'; %color for nonMP nodes 

cases = [2 3 4]; %case id (direct assignment)
Nc_target = 33; %should be equivalent to the one used in reference case 

casesb = (2:6); %case id (incremental assignment)

colReg{1} = '#0072BD';
colReg{2} = '#D95319';
colReg{3} = '#77AC30';

%for filtered selection analysis 

%print changes in (over peak hour) + (overall):
% 1. average node congestion
% 2. duration of high node congestion in peak hour
% 3. variance of node queues

% horizontal axis: node ID
% vertical axis: difference wrt NC (for every scenario different point color
% selected nodes: triangles
% non-selected nodes: circles

% one set: direct assignment
% other set: incremental assignment

%Define peak period as starting and ending time step for figures
k_s = 41*indata.c_int; %starting point of the peak interval
k_e = 120*indata.c_int; %ending point of the peak interval
%k_s = 1;
%k_e = 6*3600;

%%
f3 = figure('Name','Difference wrt NC');
for j = 1:length(cases)
    
    f3s((j-1)*3+1) = subplot(length(cases),3,(j-1)*3+1);
    xlabel('node ID')
    title('node occupancy change')
    hold on
    grid on
    
    f3s((j-1)*3+2) = subplot(length(cases),3,(j-1)*3+2);
    xlabel('node ID')
    %ylabel('change node queue variance')
    if j==1
        title(strcat(num2str(cases(j)*5+5),' \% MP nodes \newline node variance change'))
    else
        title(strcat(num2str(cases(j)*5+5),' \% MP nodes'))
    end
    hold on
    grid on
    
    f3s((j-1)*3+3) = subplot(length(cases),3,(j-1)*3+3);
    xlabel('node ID')
    ylabel('(cycles)')
    title('congestion duration change')
    hold on
    grid on
end

%% Initialize figures for PDFs/CDFs 
f6 = figure('Name','Histograms of differences wrt NC');
% 3 x 2 = 6 plots per case
for j = 1:length(cases)
    
    f6s((j-1)*6+1) = subplot(length(cases)*2,3,(j-1)*6+1);
    xlabel('$m_1 - m_1^{\star}$')
    %ylabel('no of nodes')
    title(strcat(num2str(cases(j)*5+5),' \% MP nodes, Selected'))
    % Selected nodes: dif in avg node occupancy
    hold on
    
    
    f6s((j-1)*6+2) = subplot(length(cases)*2,3,(j-1)*6+2);
    xlabel('$m_2 - m_2^{\star}$')
    %ylabel('no of nodes')
    % Selected nodes: dif in avg node variance
    hold on
    
    
    f6s((j-1)*6+3) = subplot(length(cases)*2,3,(j-1)*6+3);
    xlabel('dif. congested cycles')
    %ylabel('no of nodes')
    hold on
    % Selected nodes: dif in duration of high congestion
    
    f6s((j-1)*6+4) = subplot(length(cases)*2,3,(j-1)*6+4);
    title(strcat(num2str(cases(j)*5+5),' \% MP nodes, Not selected'))
    xlabel('$m_1 - m_1^{\star}$')
    %ylabel('no of nodes')
    % Not Selected nodes: dif in avg node occupancy
    hold on
    
    f6s((j-1)*6+5) = subplot(length(cases)*2,3,(j-1)*6+5);
    % Not Selected nodes: dif in avg node variance
    xlabel('$m_2 - m_2^{\star}$')
    %ylabel('no of nodes')
    hold on
    
    f6s((j-1)*6+6) = subplot(length(cases)*2,3,(j-1)*6+6);
    % Not Selected nodes: dif in duration of high congestion
    xlabel('dif. congested cycles')
    %ylabel('no of nodes')
    hold on
end

%% Initialize figure for the CDF plots 
f7 = figure('Name','CDFs of differences wrt NC');
%CDFs in common plot for every case (dif color for MP and nonMP) 
for j = 1:length(cases)
   f7s((j-1)*3+1) = subplot(length(cases),3,(j-1)*3+1);

   hold on

   f7s((j-1)*3+2) = subplot(length(cases),3,(j-1)*3+2);
   hold on 
   
   f7s((j-1)*3+3) = subplot(length(cases),3,(j-1)*3+3);
   hold on 
    
end
    


%% Store values of the reference case (NC)
nodeMetricsNC = [];
for reg=1:3
    for i=1:length(indata.nodereg{reg})
        %scan all nodes of the region
        ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
        %calculate mean of the average queue of all incoming links over
        %time - normalize over storage
        
        
        if MP.approaches{ind}>0
            MP_nodesNC{reg}(i) = ind;
            node_occNC = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
            node_queuesNC{reg}(i) = mean(mean(node_occNC));
            %calculate mean of the variance of queues of all incoming links
            %over time
            %             node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
            node_varNC{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
            % count number of cycles
            node_congestedNC{reg}(i) = sum(sum(node_occNC>0.8)>0)/(indata.c_int*CyclesInPeak); %number of cycles (aggregated) when node is "congested"
            
        end
        
    end
    
    %node ID (index in MP) - node avg norm queue - node avg var of queues -
    %node no of high congestion cycles
    nodeMetricsNC{reg} = [MP_nodesNC{reg}(:) node_queuesNC{reg}(:) ...
        node_varNC{reg}(:) node_congestedNC{reg}(:)];
 
end


%% Load and calculate values of the examined case
for j=1:length(cases)
    
    % load examined case results
    %load node information of the case: nodes included in the plan
    if cases(j) == 0 && j==1 
        load(strcat(pathIN,filethemeMP0,num2str(cases(j))),'MPnodes') %100 percent 
    else
        load(strcat(pathIN,filethemeMP_MS,num2str(cases(j))),'MPnodes')
    end
    % create new 4-subplots figure per case
    % before/after figures in groups of 4 subplots (2 before - 2 after per case)
    
    f5(j) = figure('Name','Changes in node variance and congestion before/after');
    f5s(j,1) = subplot(2,2,1);
    %xlabel('avg node variance ($m_2$)')
    xlabel('$m_2$')
    %ylabel('avg node congestion ($m_1$)')
    ylabel('$m_1$')
    title('FTC')
    xlim([0 0.6])
    ylim([0 1])
    hold on
    grid on 
    
    f5s(j,2) = subplot(2,2,2);
    title('FTC')
    %xlabel('avg node variance ($m_2$)')
    xlabel('$m_2$')
    ylabel('N_c')
    xlim([0 0.6])
    ylim([0 80])
    hold on
    grid on
    
    f5s(j,3) = subplot(2,2,3);
    xlabel('$m_2$')
    ylabel('$m_1$')
    title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
    xlim([0 0.6])
    ylim([0 1])
    hold on
    grid on
    
    f5s(j,4) = subplot(2,2,4);
    xlabel('$m_2$')
    ylabel('N_c')
    title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
    xlim([0 0.6])
    ylim([0 80])
    hold on
    grid on
    
    
    % create subplot figures for the classification analysis (2 subplots)
    % color or size of the point corresponds to performance 
    f8(j) = figure('Name','Changes in MP criteria bef/aft + performace indication');
    f8s(j,1) = subplot(2,1,1); 
    xlabel('$m_2$')
    ylabel('$m_1$')
    title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
    hold on 
    grid on 
    
    f8s(j,2) = subplot(2,1,2); 
    xlabel('$m_2$')
    ylabel('$N_c$')
    title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
    hold on 
    grid on 
    
    
    %update before/after figures of NC case by coloring MP nodes     
    node_queues = {};
    node_var = {};
    node_congested = {};
    nodeMetricsTest = [];
    MPnT = cell(1,3);
    nonMPnT = cell(1,3);
    
    %load case results after MP application
    if cases(j)==0 && j==1
        load(strcat(pathTest,resfilethemeMP0,num2str(cases(j))),'outdata')
    else
        load(strcat(pathTest,resfilethemeMP_MS,num2str(cases(j))),'outdata')
    end
    %%
    for reg=1:3
        figure(f5(j))
        for i=1:length(indata.nodereg{reg})
            %scan all nodes of the region
            ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
            %calculate mean of the average queue of all incoming links over
            %time - normalize over storage capacity
            MP_nodes{reg}(i) = ind; %index in MP struct / separate nodes by region 
            if MP.approaches{ind}>0
                %MP_nodes{reg}(i) = ind; %index in MP
                node_occ = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
                node_queues{reg}(i) = mean(mean(node_occ));
                %calculate mean of the variance of queues of all incoming links
                %over time
                %             node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
                node_var{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
                % count number of cycles
                node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/(indata.c_int*CyclesInPeak); %number of cycles (aggregated) when node is "congested"
                if ismember(ind, MPnodes)
                    MPnT{reg} = [MPnT{reg} i]; %index of MP nodes in nodereg{reg}
                else
                    nonMPnT{reg} = [nonMPnT{reg} i]; %index of nonMP nodes in nodereg{reg}
                end
            end
        end
        %node ID (index in MP) - node avg norm queue - node avg var of queues -
        %node no of high congestion cycles (m_1 - m_2 - N_c)
        nodeMetricsTest{reg} = [MP_nodes{reg}(:) node_queues{reg}(:) ...
            node_var{reg}(:) node_congested{reg}(:)];
        
        
        %update the before/after figures
        
        subplot(f5s(j,1))
        scatter(node_varNC{reg}(nonMPnT{reg}),node_queuesNC{reg}(nonMPnT{reg}),mrkSize,'filled','MarkerFaceColor',colNonMP)
        scatter(node_varNC{reg}(MPnT{reg}),node_queuesNC{reg}(MPnT{reg}),mrkSize,'filled','MarkerFaceColor',colMP)
        
        subplot(f5s(j,2))
        scatter(node_varNC{reg}(nonMPnT{reg}),node_congestedNC{reg}(nonMPnT{reg}),mrkSize,'filled','MarkerFaceColor',colNonMP)
        scatter(node_varNC{reg}(MPnT{reg}),node_congestedNC{reg}(MPnT{reg}),mrkSize,'filled','MarkerFaceColor',colMP)
        
        subplot(f5s(j,3))
        scatter(node_var{reg}(MPnT{reg}),node_queues{reg}(MPnT{reg}),mrkSize,'filled','MarkerFaceColor',colMP)
        scatter(node_var{reg}(nonMPnT{reg}),node_queues{reg}(nonMPnT{reg}),mrkSize,'filled','MarkerFaceColor',colNonMP)
        
        subplot(f5s(j,4))
        scatter(node_var{reg}(MPnT{reg}),node_congested{reg}(MPnT{reg}),mrkSize,'filled','MarkerFaceColor',colMP)
        scatter(node_var{reg}(nonMPnT{reg}),node_congested{reg}(nonMPnT{reg}),mrkSize,'filled','MarkerFaceColor',colNonMP)
        
    end
    
    %% Calculate differences wrt NC
    for reg=1:3
        metricsDif{reg} = (nodeMetricsTest{reg}(:,2:end)-nodeMetricsNC{reg}(:,2:end));
        %percentmetricsDif(:,1) = metricsDif(:,1)./nodeMetricsNC(:,2)*100;
        %percentmetricsDif(:,2) = metricsDif(:,2)./nodeMetricsNC(:,3)*100;
        %percentmetricsDif(:,3) = metricsDif(:,3)./nodeMetricsNC(:,4)*100;
    end
    
    metricsDifG_MPn = [];
    metricsDifG_nonMPn = [];
    
    for reg=1:3
        
        
        %% Create two 3d figure m1-m2-(impr) & Nc-m2-(impr) - visualize characteristics of node selection
        % impr: Nc-Nc* or m1-m1* -> reflect to total delay
        % impr = m1-m1* or impr = Nc - Nc*  
        % plotting the before m1,m2 and Nc,m2 values (of the FTC)
        figure(f8(j))
        subplot(f8s(j,1)) %m2-m1 + color/size in impr (separation of the selected?) 
        %v_mrkSize_nonMP =  metricsDif{reg}(nonMPnT{reg},1); %m_1 - m_1* as a criterion for performance before/after 
        v_mrkSize_nonMP =  metricsDif{reg}(nonMPnT{reg},2);
        %v_mrkSize_nonMP =  metricsDif{reg}(nonMPnT{reg},3); %N_c - N_c* as a criterion for performance before/after 
        
        size_scale = 500;
        scatter(node_varNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP<0)),node_queuesNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP<0)),v_mrkSize_nonMP((v_mrkSize_nonMP<0))*(-size_scale),v_mrkSize_nonMP(v_mrkSize_nonMP<0),'^','filled') %criterion decreased (impr) -> triangles for nonMP 
        scatter(node_varNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP>=0)),node_queuesNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP>=0)),v_mrkSize_nonMP((v_mrkSize_nonMP>=0))*size_scale+0.1,v_mrkSize_nonMP((v_mrkSize_nonMP>=0)),'x') %criterion increased -> reverde triangle for nonMP
        
        %v_mrkSize_MP =  metricsDif{reg}(MPnT{reg},1);
        v_mrkSize_MP =  metricsDif{reg}(MPnT{reg},2);
        %v_mrkSize_MP =  metricsDif{reg}(MPnT{reg},3);
        scatter(node_varNC{reg}(MPnT{reg}(v_mrkSize_MP<0)),node_queuesNC{reg}(MPnT{reg}(v_mrkSize_MP<0)),v_mrkSize_MP(v_mrkSize_MP<0)*(-size_scale),v_mrkSize_MP(v_mrkSize_MP<0),'o','filled') %criterion decreased (impr) -> circles filled for MP nodes
        scatter(node_varNC{reg}(MPnT{reg}(v_mrkSize_MP>=0)),node_queuesNC{reg}(MPnT{reg}(v_mrkSize_MP>=0)),v_mrkSize_MP(v_mrkSize_MP>=0)*size_scale+0.1,v_mrkSize_MP(v_mrkSize_MP>=0),'s','filled') %criterion increased (not improv) -> squares filled for MP
        clb = colorbar;
        clb.Label.Interpreter = 'latex';
        %clb.Label.String = '$N_c - N_c^{\star}$';
       % clb.Label.String = '$m_1 - m_1^{\star}$';
        clb.Label.String = '$m_2 - m_2^{\star}$';
        clb
        legend('non MP -', 'non MP +','MP -','MP +')
        
        subplot(f8s(j,2))
        %scatter(node_varNC{reg}(nonMPnT{reg}),node_congestedNC{reg}(nonMPnT{reg}),v_mrkSize_nonMP*50,'MarkerEdgeColor',colNonMP)
        %scatter(node_varNC{reg}(MPnT{reg}),node_congestedNC{reg}(MPnT{reg}),v_mrkSize_MP*50,'s','filled','MarkerFaceColor',colMP)
        scatter(node_varNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP<0)),node_congestedNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP<0)),v_mrkSize_nonMP((v_mrkSize_nonMP<0))*(-size_scale),v_mrkSize_nonMP(v_mrkSize_nonMP<0),'^','filled')
        scatter(node_varNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP>=0)),node_congestedNC{reg}(nonMPnT{reg}(v_mrkSize_nonMP>=0)),v_mrkSize_nonMP((v_mrkSize_nonMP>=0))*(size_scale)+0.1,v_mrkSize_nonMP(v_mrkSize_nonMP>=0),'x')
        scatter(node_varNC{reg}(MPnT{reg}(v_mrkSize_MP<0)),node_congestedNC{reg}(MPnT{reg}(v_mrkSize_MP<0)),v_mrkSize_MP(v_mrkSize_MP<0)*(-size_scale),v_mrkSize_MP(v_mrkSize_MP<0),'o','filled')
        scatter(node_varNC{reg}(MPnT{reg}(v_mrkSize_MP>=0)),node_congestedNC{reg}(MPnT{reg}(v_mrkSize_MP>=0)),v_mrkSize_MP(v_mrkSize_MP>=0)*size_scale+0.1,v_mrkSize_MP(v_mrkSize_MP>=0),'s','filled')
        clb = colorbar;
        clb.Label.Interpreter =  'latex';
        %clb.Label.String = '$N_c - N_c^{\star}$';
        %clb.Label.String = '$m_1 - m_1^{\star}$';
        clb.Label.String = '$m_2 - m_2^{\star}$';
        clb
        legend('non MP -', 'non MP +','MP -','MP +')
        
        
      
        %% print results in figures
        figure(f3)
        subplot(f3s((j-1)*3+1))
        scatter(MP_nodes{reg}(nonMPnT{reg}),metricsDif{reg}(nonMPnT{reg},1),mrkSize,'filled'...
            ,'MarkerFaceColor',colNonMP)
        scatter(MP_nodes{reg}(MPnT{reg}),metricsDif{reg}(MPnT{reg},1),mrkSize,'filled',...
            'MarkerFaceColor',colMP)
        ylim([-0.3 0.3])
        
        subplot(f3s((j-1)*3+2))
        scatter(MP_nodes{reg}(nonMPnT{reg}),metricsDif{reg}(nonMPnT{reg},2),mrkSize,'filled',...
            'MarkerFaceColor',colNonMP)
        scatter(MP_nodes{reg}(MPnT{reg}),metricsDif{reg}(MPnT{reg},2),mrkSize,'filled',...
            'MarkerFaceColor',colMP)
        ylim([-0.2 0.1])
        
        
        
        subplot(f3s((j-1)*3+3))
        scatter(MP_nodes{reg}(nonMPnT{reg}),metricsDif{reg}(nonMPnT{reg},3),mrkSize,'filled',...
            'MarkerFaceColor',colNonMP)
        scatter(MP_nodes{reg}(MPnT{reg}),metricsDif{reg}(MPnT{reg},3),mrkSize,'filled',...
            'MarkerFaceColor',colMP)
        ylim([-70 70])
        
        
        % gather results to make distributions 
        metricsDifG_MPn = [metricsDifG_MPn; metricsDif{reg}(MPnT{reg},:)];
        metricsDifG_nonMPn = [metricsDifG_nonMPn; metricsDif{reg}(nonMPnT{reg},:)];
    end
    
%     %% Rank nodes according to the relationship A = a m_1 + b m_2 after filtering for Nc 
%     
%     % Find a, b so that after ranking nodes A, the top x% of N_C - N_c*
%     % (or other obj function expression) is maximum
%     % assuming a, b are \in (0,1) 
%     % load all eligible MP nodes 
%     elig_MPnodes = load('MPnodesCase_MP_0','MPnodes'); %all eligible nodes 
%     
%     best_found = []; %to store results of tests a,b  
%     for a=-2:0.05:2
%         for b=-2:0.05:2
%             nod = []; %node IDs in MP
%             A = [];  %selection function
%             B = []; %Nc for filtering
%             sel = []; %selected
%             med = []; %metrics difference (m_i-m_i* etc.) 
%             for reg = 1:3
%                 nod = [nod; nodeMetricsTest{reg}(:,1)]; %node ind in MP (they are added by region)
%                 A = [A; a*(0.85-nodeMetricsTest{reg}(:,2))+b*(nodeMetricsTest{reg}(:,3)-0.03)/0.6];
%                 B = [B; nodeMetricsTest{reg}(:,4)]; %store all Nc values 
%                 %I need to exclude non-candidate nodes ?
%                 
%                 %- the highest A is the highest avg queue, var and duration of high congestion
%                 %- the highest A, the highest the improvement potential
%                 
%                 med = [med; metricsDif{reg}(:,3)]; %performance N_c-N_C* (column 3) -
%                 % the higher the difference - the higher the degradation /
%                 % for improvement this should be negative
%                 % sel = [sel; (length(sel) + MPnT{reg})']; %indices of MPnodes in med, nod etc.
%                 sel = [sel; ismember(nodeMetricsTest{reg}(:,1),MPnodes)]; %logical array to say if selected or not
%             end
%             numTop = min(0.25*length(elig_MPnodes.MPnodes),length(MPnodes)); %count the top 20% of the nodes
%             
%             %Exclude non-eligible MPnodes
%             elig_ind = ismember(nod,elig_MPnodes.MPnodes); %ind of eligible
%             %only keep the eligible nodes
%             A = A(elig_ind);
%             B = B(elig_ind);
%             %med = med(elig_ind);
%             nod = nod(elig_ind);
%             %sel = sel(elig_ind);
%             
%             med_elig = med(elig_ind);
%             
%             % Filter by Nc
%             filtered = B>=Nc_target;
%             
%             % filter by Nc
%             A = A(filtered);
%             med_elig = med_elig(filtered);
%             
%             % sort by descending values of A
%             [A,Ia] = sort(A,'descend'); %max to min selection function
%             
%             Rmed = med_elig(Ia);  %performance of nodes sorted according to A
%             if numTop > length(A)
%                 error('Filtered subset not large enough. Decrease Nc_target or increase percentage target')
%                 
%             else
%                 sum_top = sum(A(1:numTop)); %sum of A of the top x% in the ranking
%                 sum_perf_top = sum(Rmed(1:numTop)); %sum of observed performance of the resulting top percent
%                 sum_perf_selected = sum(med(sel==1)); %sum of observed performance of selected nodes
%                 best_found = [best_found; a b sum_perf_top sum_perf_selected];
%             end
%         end
%     end
%     %plot performance of top x% of A + performance of selected nodes from
%     %previous 
%     figure
%     plot(best_found(:,3:end),'DisplayName','best_found(:,3:end)')
%     [min_s,ind] = min(best_found(:,3)) %somethig is wrong with this - the max is the solution (counter intuitive) 
%     best_found(ind,:)
    
    %% Rank nodes according to the relationship A = a m_1 + b m_2 + (1-a-b) m_3
    
    % Find a, b so that after ranking nodes A, the top x% of N_C - N_c*
    % (or other obj function expression) is maximum
    % assuming a, b are \in (0,1) 
    % load all eligible MP nodes 
    %     elig_MPnodes = load('MPnodesCase_MP_0','MPnodes');
    %     I = [];
    %     best_found = []; %to store results of tests a,b
    %     for a=0:0.1:1
    %         for b=0:0.1:1
    %              %for c=0:0.1:2
    %                 nod = []; %node IDs in MP
    %                 A = [];  %selection function
    %                 sel = []; %selected
    %                 med = []; %metrics difference
    %                 for reg = 1:3
    %                     nod = [nod; nodeMetricsTest{reg}(:,1)]; %node ind in MP (they are added by region)
    %                    % A = [A; 1.2*a*(0.85-nodeMetricsTest{reg}(:,2))+1.2*b*(nodeMetricsTest{reg}(:,3)-0.03)+1.2*c*nodeMetricsTest{reg}(:,4)/80]; %selection function (to calibrate)
    %                     A = [A; a*(0.85-nodeMetricsTest{reg}(:,2))+b*(nodeMetricsTest{reg}(:,3)-0.03)/0.6-(1-a-b)*nodeMetricsTest{reg}(:,4)]; %selection function (to calibrate)
    %
    %                     %I need to exclude non-candidate nodes ?
    %
    %                     %A = [A; a*nodeMetricsTest{reg}(:,2) + b*nodeMetricsTest{reg}(:,3) ...
    %                     %    + c*nodeMetricsTest{reg}(:,4)/80]; %selection function (to calibrate)
    %
    %
    %                     %- the highest A is the highest avg queue, var and duration of high congestion
    %                     %- the highest A, the highest the improvement potential
    %                     med = [med; metricsDif{reg}(:,3)]; %performance N_c-N_C* (column 3) -
    %                     % the higher the difference - the higher the degradation /
    %                     % for improvement this should be negative
    %                     % sel = [sel; (length(sel) + MPnT{reg})']; %indices of MPnodes in med, nod etc.
    %                     sel = [sel; ismember(nodeMetricsTest{reg}(:,1),MPnodes)]; %logical array to say if selected or not
    %                 end
    %                 numTop = min(0.25*length(elig_MPnodes.MPnodes),length(MPnodes)); %count the top 20% of the nodes
    %
    %                 %Exclude non-eligible MPnodes
    %                 elig_ind = ismember(nod,elig_MPnodes.MPnodes); %ind of eligible
    %                 %only keep the eligible nodes
    %                 A = A(elig_ind);
    %                 [A,Ia] = sort(A,'descend'); %min to max selection function
    %                 %I = [I Ia];
    %                 med = med(elig_ind);
    %                 nod = nod(elig_ind);
    %                 sel = sel(elig_ind);
    %
    %                 Rmed = med(Ia); %performance of nodes sorted according to A
    %                 sum_top = sum(A(1:numTop)); %sum of A of the top x% in the ranking
    %                 sum_perf_top = sum(Rmed(1:numTop)); %sum of observed performance of the resulting top percent
    %                 sum_perf_selected = sum(med(sel==1)); %sum of observed performance of selected nodes
    %                 best_found = [best_found; a b sum_perf_top sum_perf_selected];
    %              %end
    %         end
    %     end
    %     %plot performance of top x% of A + performance of selected nodes from
    %     %previous
    %     figure
    %     plot(best_found(:,3:end),'DisplayName','best_found(:,3:end)')
    %     [min_s,ind] = max(best_found(:,3)) %somethig is wrong with this - the max is the solution (counter intuitive)
    %     best_found(ind,:)
     
    
    %% Histograms/CDF for the selected nodes (difference in metrics)
    for k=1:3
        figure(f6)
        subplot(f6s((j-1)*6+k))
        histogram(metricsDifG_MPn(:,k),'Normalization','PDF','FaceColor','#0072BD')
        mu = mean(metricsDifG_MPn(:,k));
        xline(mu,'Color','r','LineWidth',2);
        if k==1 
            xlim([-0.5 0.5])
        elseif k==2
            xlim([-0.3 0.3])
        else
            xlim([-60 60])
        end
         
    end
    
    % Hist for the non-selected nodes
    for k=4:6
        figure(f6)
        subplot(f6s((j-1)*6+k))
        histogram(metricsDifG_nonMPn(:,k-3),'Normalization','PDF','FaceColor','#D95319')
        mu = mean(metricsDifG_nonMPn(:,k-3))
        xline(mu,'Color','r','LineWidth',2)
        if k==4
            xlim([-0.5 0.5])
        elseif k==5
            xlim([-0.3 0.3])
        else
            xlim([-60 60])
        end
    end
    %% CDF plots 
    for k=1:3
        figure(f7)
        subplot(f7s((j-1)*3+k))
        %CDF plot for the selected nodes
        h1 = cdfplot(metricsDifG_MPn(:,k));
        h1.LineWidth = 2;
         
        %CDF plot for the non selected nodes 
        if cases(j)>0
            h2 = cdfplot(metricsDifG_nonMPn(:,k));
            h2.LineWidth = 2;
        end
        if k==2
            title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),' \% MP nodes'))
        end
        
        
        
        if k==1
            xlim([-0.5 0.5])
            xlabel('$m_1 - m_1^{\star}$')
            ylabel('$f$')
            title('')
        elseif k==2
            xlim([-0.3 0.3])
            xlabel('$m_2 - m_2^{\star}$')
            ylabel('$f$')
        else
            xlim([-1 1])
            xlabel('$N_c - N_c^{\star}$')
            ylabel('$f$')
            title('')
            
        end
    end
    legend('MP','non-MP','Location','southeast')
    
    
end



% %% Remake the figures but for the cases of incremental assignment 
% % Triple figure with absolute values of metrics and different colors per
% % batch of assigned nodes - 
% colMP = cell(1,6);
% colMP{1} = '#0072BD';
% colMP{2} = '#D95319';
% colMP{3} = '#EDB120';
% colMP{4} = '#7E2F8E';
% colMP{5} = '#4DBEEE';
% colMP{6} = '#A2142F';
% 
% node_queues = {};
% node_var = {};
% node_congested = {};
% nodeMetricsTest = [];
% MPnT = cell(3,length(casesb)+1);
% nonMPnT = cell(3,length(casesb)+1);
    
% Initialize figure 
% f4 = figure('Name','Selection metrics - Incremental');
% for j = 1:length(casesb)+2
%     
%     f4s((j-1)*3+1) = subplot(length(cases)+1,3,(j-1)*3+1);
%     xlabel('node ID')
%     title('node occupancy change')
%     hold on
%     grid on
%     
%     f4s((j-1)*3+2) = subplot(length(cases)+1,3,(j-1)*3+2);
%     xlabel('node ID')
%     %ylabel('change node queue variance')
%     if j==1
%         title('Case \hspace{0.2 cm} FTC - node variance change')
%     else
%         title(strcat('Case \hspace{0.2 cm} ',num2str((j-1)*5),' \% MP nodes (IS) '))
%     end
%     hold on
%     grid on
%     
%     f4s((j-1)*3+3) = subplot(length(cases)+1,3,(j-1)*3+3);
%     xlabel('node ID')
%     ylabel('(cycles)')
%     title('congestion duration change')
%     hold on
%     grid on
% end
% 
% % Print results of No Control in first row (first reference) + nodes of first selection in color 
% load('MPnodesCase_MP_S0.mat','MPnodes') %same case of 5% in both ways (ref is NC) 
%  
% % find new indices 
%  for reg=1:3
%         for i=1:length(indata.nodereg{reg})
%             %scan all nodes of the region
%             ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
%             %calculate mean of the average queue of all incoming links over
%             %time - normalize over storage capacity
%             
%             if MP.approaches{ind}>0
%                 MP_nodes{reg}(i) = ind; %index in MP
%                 if ismember(ind, MPnodes)
%                     MPnT{reg,1} = [MPnT{reg,1} i];
%                 else
%                     nonMPnT{reg,1} = [nonMPnT{reg,1} i];
%                 end
%             end
%         end
%  end
% % print No Control Results 
% j=1;
% for reg=1:3
%     % print results in figures
%     figure(f4)
%     subplot(f4s((j-1)*3+1))
%     scatter(MP_nodes{reg}(nonMPnT{reg,1}),nodeMetricsNC{reg}(nonMPnT{reg,1},2),mrkSize,'filled'...
%         ,'MarkerFaceColor',colNonMP)
%     scatter(MP_nodes{reg}(MPnT{reg,1}),nodeMetricsNC{reg}(MPnT{reg,1},2),mrkSize,'filled',...
%         'MarkerFaceColor',colMP{1})
%     %ylim([-0.3 0.3])
%     
%     subplot(f4s((j-1)*3+2))
%     scatter(MP_nodes{reg}(nonMPnT{reg,1}),nodeMetricsNC{reg}(nonMPnT{reg,1},3),mrkSize,'filled',...
%         'MarkerFaceColor',colNonMP)
%     scatter(MP_nodes{reg}(MPnT{reg,1}),nodeMetricsNC{reg}(MPnT{reg,1},3),mrkSize,'filled',...
%         'MarkerFaceColor',colMP{1})
%     %ylim([-0.2 0.1])
%     
%     
%     subplot(f4s((j-1)*3+3))
%     scatter(MP_nodes{reg}(nonMPnT{reg,1}),nodeMetricsNC{reg}(nonMPnT{reg,1},4),mrkSize,'filled',...
%         'MarkerFaceColor',colNonMP)
%     scatter(MP_nodes{reg}(MPnT{reg,1}),nodeMetricsNC{reg}(MPnT{reg,1},4),mrkSize,'filled',...
%         'MarkerFaceColor',colMP{1})
%     %ylim([-70 70])
%     
% end
% 
% j=2; 
% 
% % Print results of first node set (S0)
% %load next MPnode set 
% load(strcat('MPnodesCase_MP_ISm_',num2str(2)),'MPnodes')
% 
% % load new reference case 
% load(strcat(pathTest,'output_set4_MP_S0'),'outdata')
% for reg=1:3
%         for i=1:length(indata.nodereg{reg})
%             %scan all nodes of the region
%             ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
%             %calculate mean of the average queue of all incoming links over
%             %time - normalize over storage capacity
%             
%             if MP.approaches{ind}>0
%                 MP_nodes{reg}(i) = ind; %index in MP
%                 node_occ = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
%                 node_queues{reg}(i) = mean(mean(node_occ));
%                 %calculate mean of the variance of queues of all incoming links
%                 %over time
%                 %             node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
%                 node_var{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
%                 % count number of cycles
%                 node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
%                 if ismember(ind, MPnodes)
%                     MPnT{reg,j} = [MPnT{reg,j} i];
%                 else
%                     nonMPnT{reg,j} = [nonMPnT{reg,j} i];
%                 end
%             end
%             
%             
%         end
%         MPnT{reg,j} = setxor(MPnT{reg,j-1},MPnT{reg,j}); %keep only new batch of nodes
%         %node ID (index in MP) - node avg norm queue - node avg var of queues -
%         %node no of high congestion cycles
%         nodeMetricsTest{reg} = [MP_nodes{reg}(:) node_queues{reg}(:) ...
%             node_var{reg}(:) node_congested{reg}(:)];
%         
% end
% 
% 
% 
% %print results of new case
% for reg=1:3
%     %% print results in figures
%     figure(f4)
%     subplot(f4s((j-1)*3+1))
%     scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},2),mrkSize,'filled'...
%         ,'MarkerFaceColor',colNonMP)
%     for k=1:j
%         scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},2),mrkSize,'filled',...
%             'MarkerFaceColor',colMP{k})
%     end
%     ylim([0 1])
%     
%     subplot(f4s((j-1)*3+2))
%     scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},3),mrkSize,'filled',...
%         'MarkerFaceColor',colNonMP)
%     for k=1:j
%         scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},3),mrkSize,'filled',...
%             'MarkerFaceColor',colMP{k})
%     end
%     ylim([0 0.5])
%     
%     
%     subplot(f4s((j-1)*3+3))
%     scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},4),mrkSize,'filled',...
%         'MarkerFaceColor',colNonMP)
%     for k=1:j
%         scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},4),mrkSize,'filled',...
%             'MarkerFaceColor',colMP{k})
%     end
%     ylim([0 80])
%     
% end
% 
% % Load next reference case and print nodes in multiple colors
% j=3; 
% for kl=1:length(casesb)
%    
%     node_queues = {};
%     node_var = {};
%     node_congested = {};
%     nodeMetricsTest = [];
%     %load new set of nodes 
%     if kl<length(casesb)
%         load(strcat('MPnodesCase_MP_ISm_',num2str(casesb(kl+1))),'MPnodes')
%     end 
%     %load new reference case 
%     load(strcat(pathTest,'output_set4_MP_ISm_',num2str(casesb(kl))),'outdata')
%     
%     for reg=1:3
%         for i=1:length(indata.nodereg{reg})
%             %scan all nodes of the region
%             ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
%             %calculate mean of the average queue of all incoming links over
%             %time - normalize over storage capacity
%             
%             if MP.approaches{ind}>0
%                 MP_nodes{reg}(i) = ind; %index in MP
%                 node_occ = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
%                 node_queues{reg}(i) = mean(mean(node_occ));
%                 %calculate mean of the variance of queues of all incoming links
%                 %over time
%                 %             node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
%                 node_var{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
%                 % count number of cycles
%                 node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
%                 if kl<length(casesb)
%                     if ismember(ind, MPnodes)
%                         MPnT{reg,j} = [MPnT{reg,j} i];
%                     else
%                         nonMPnT{reg,j} = [nonMPnT{reg,j} i];
%                     end
%                 end
%             end
%             
%         end
%         
%         if kl<length(casesb)
%             for km=1:j-1
%                 MPnT{reg,j} = setxor(MPnT{reg,j-km},MPnT{reg,j}); %keep only new batch of node
%             end
%         end
%         %node ID (index in MP) - node avg norm queue - node avg var of queues -
%         %node no of high congestion cycles
%         nodeMetricsTest{reg} = [MP_nodes{reg}(:) node_queues{reg}(:) ...
%             node_var{reg}(:) node_congested{reg}(:)];
%         
%         if j==7
%             MPnT{reg,7} = MPnT{reg,6};
%             nonMPnT{reg,7} = nonMPnT{reg,6} ;
%             
%             % print results in figures
%             figure(f4)
%             subplot(f4s((j-1)*3+1))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},2),mrkSize,'filled'...
%                 ,'MarkerFaceColor',colNonMP)
%             for k=1:j-1
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},2),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 1])
%             
%             subplot(f4s((j-1)*3+2))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},3),mrkSize,'filled',...
%                 'MarkerFaceColor',colNonMP)
%             for k=1:j-1
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},3),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 0.5])
%             
%             
%             subplot(f4s((j-1)*3+3))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},4),mrkSize,'filled',...
%                 'MarkerFaceColor',colNonMP)
%             for k=1:j-1
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},4),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 80])
%             
%             
%         else
%             
%             % print results in figures
%             figure(f4)
%             subplot(f4s((j-1)*3+1))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},2),mrkSize,'filled'...
%                 ,'MarkerFaceColor',colNonMP)
%             for k=1:j
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},2),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 1])
%             
%             subplot(f4s((j-1)*3+2))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},3),mrkSize,'filled',...
%                 'MarkerFaceColor',colNonMP)
%             for k=1:j
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},3),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 0.5])
%             
%             
%             subplot(f4s((j-1)*3+3))
%             scatter(MP_nodes{reg}(nonMPnT{reg,j}),nodeMetricsTest{reg}(nonMPnT{reg,j},4),mrkSize,'filled',...
%                 'MarkerFaceColor',colNonMP)
%             for k=1:j
%                 scatter(MP_nodes{reg}(MPnT{reg,k}),nodeMetricsTest{reg}(MPnT{reg,k},4),mrkSize,'filled',...
%                     'MarkerFaceColor',colMP{k})
%             end
%             ylim([0 80])
%             
%         end
%         
%     end
%     
%     j = j + 1;
% end

%%
% for j=1:7
%     figure(f4)
%     subplot(f4s((j-1)*3+3))
%     ylim([0 200])
% end
