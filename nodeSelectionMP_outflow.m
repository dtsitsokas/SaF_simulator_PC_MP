function [MPnodes_R] = nodeSelectionMP_outflow(demandCode,a,b,case_j,equation,caseName)
% nodeSelectionMP_outflow: Performs the selection of nodes for application of Max-Pressure
% by starting from nodes that serve the highest outflow (simple heuristic). Requires NC scenario results.


%% Load reference case file name (NC or other) 
if demandCode == 1
    %demand = 'med';
    fname = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\output_set5_MedDemand_NC';
    
elseif demandCode == 2
    % fname = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Current results\ResultsMediumOD_capacityDrop\output_set5_MedDemand_PC_61';
    %demand = 'high';
    
    %Case NC as reference for MP selection
    %fname = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Current results\ResultsHighOD_capacityDrop\output_set4_test_NC';
    fname = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\output_set5_highDemand_NC\';
    
    %Case PC as reference for MP selection 
    %fname = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\output_set5_highDemand_PC_50b';
    %fname = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\output_set5_highDemand_PC_50b\';
end

load(fname,'indata','outdata') %load results only here, nodes from the input
% load(fname,'indata','outdata','MPnodes')  %load results + nodes from fname

load('FinalInput','MP')
%%

num_MPall = sum(MP.splan==1)-length(indata.specialInt); %total number of MP eligible nodes
node_queues = {};
node_var = {};
node_outflows = {}; 

MPnodes = [];
%nodesPreSelected = MPnodes;

% - Set congestion threshold for peak period of 2 hours/80 steps - check
% perctage at the end and see if need to increase or decrease
cong_threshold = 0;

% (to ignore peak hour selection, set cong_threshold to 0 and peak-hour from
% 0 to end )

%settings for selection method based on selection function A
mode_A = 1; %if not 1, TRB method applied 
mode_C = 0; %filter by Nc (if 1: change the expression for A)
if mode_A==1
    A = cell(1,3); %stores the result of a m_1 + b m_2 + (1-a-b) N_c for every region
    O = cell(1,3);
    B = []; %all values of A gathered (all regions)
    C = []; %store Nc for separate filtering
    R = []; %store values of outflows (all regions) 
    
    MP_nodesAll = [];
    target_p = 0.05+case_j*0.05;
else
    target_p = 0.05; %set target p (trb method only)
end

%Exclude special intersections
MPspInt = zeros(1,length(indata.specialInt));
for i=1:length(indata.specialInt)
    %indices of special intersections in MP
    MPspInt(i) = find(MP.nodeID == indata.specialInt(i));
end

%Initialize common figrues (after first filtration)

f1 = figure('Name','node variance vs. node congestion');
%title('peak period (2h)')
xlabel('$m_2^n$')
ylabel('$m_1^n$')
hold on

f2 = figure('Name','node variance vs. no of cycles is congested');
%title('peak period (2h)')
xlabel('$m_2^n$')
ylabel('$N_c^n$')
hold on

f3 = figure('Name','node variance vs. node congestion');
%title('peak period (2h)')
xlabel('$m_2^n$')
ylabel('$m_1^n$')
hold on

f4 = figure('Name','node variance vs. no of cycles is congested');
%title('peak period (2h)')
xlabel('$m_2^n$')
ylabel('$N_c^n$')
hold on

colReg{1} = '#0072BD';
colReg{2} = '#0072BD';%'#D95319';
colReg{3} = '#0072BD';%'#77AC30';
not_picked = {};
not_signalized = {};

%Define peak period as starting and ending time step
%k_s = 41*indata.c_int; %starting point of the peak interval
%k_e = 120*indata.c_int; %ending point of the peak interval
k_s = 1;
k_e = 6*3600;

CyclesInPeak = (k_e - k_s)/indata.c_int+1;

node_queuesAll = [];
node_varAll = [];
node_congAll = [];
node_outflowsAll = []; 
for reg=1:3
    % calculate for every node of the entire simulation period:
    % 1. Total outflow of the node (sum of outflow of all incoming links) 
    % 3. Start selecting nodes from those with highest total outflow 
    
    for i=1:length(indata.nodereg{reg}) %scan all nodes of the region
        
        ind = find(MP.nodeID==indata.nodereg{reg}(i)); %index in MP
        
        MP_nodes{reg}(i) = ind;
        if MP.approaches{ind}>0
                     
            node_occ = outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}));
            
            node_queues{reg}(i) = mean(mean(node_occ));
            
            node_outflows{reg}(i) = sum(sum(outdata.u(indata.junct.or_index(MP.approaches{ind}),k_s:k_e))); 
            %calculate mean of the variance of queues of all incoming links
            %over time
            % node_var{reg}(i) = mean(var(outdata.x(junct.or_index(MP.approaches{ind}),:)./indata.capacity(junct.or_index(MP.approaches{ind}))));
            node_var{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
           
            %node_var{reg}(i) = mean(var(outdata.x(x_rows,k_s:k_e)./indata.capacity(x_rows)));
            
            % count number of cycles
            node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
        end
    end
    
    node_queuesAll = [node_queuesAll node_queues{reg}];
    node_varAll = [node_varAll  node_var{reg}];
    node_congAll = [node_congAll node_congested{reg}];
    node_outflowsAll = [node_outflowsAll node_outflows{reg}];
    
    % Selection method according to the function a m_1 + b m_2 + (1-a-b) N_c
    % put all nodes in one pool (region after region) - decide in the end
    if mode_A == 1
        %first case (no filtering of Nc)
        if ~mode_C==1
            if equation==2
                %eq 2:
                A{reg} =  1.2*a*(0.15-node_queues{reg}) - b*1.2*(node_var{reg}-0.03)/0.6 - abs(1-a-b)/CyclesInPeak*node_congested{reg}; %80 = number of cycles in the defined peak hour
            elseif equation==1
                %eq 1:
                A{reg} = -1.2*a*(0.9-node_queues{reg}) -1.2*b*(node_var{reg}-0.03)/0.6 - abs(1-a-b)/CyclesInPeak*node_congested{reg}; %80 = number of cycles in the defined peak hour
            elseif equation==3 %(1 prev version)
                A{reg} = -1.2*a*(0.9-node_queues{reg})-1.2*b*(node_var{reg}-0.03) - abs(1-a-b)/CyclesInPeak*node_congested{reg}; %the one used
     
               % A{reg} = -1.2*a*(node_queues{reg}-0.3) -1.2*b*(node_var{reg}-0.03) - abs(1-a-b)/CyclesInPeak*node_congested{reg}; %trying new (not used yet)

            end
        else
            %second case (with filtering)
            A{reg} = - a*(0.85-node_queues{reg}) - b*(node_var{reg}-0.03)/0.6;
        end
        O{reg} = node_outflows{reg}; 
        B = [B A{reg}]; %all A values in one array
        R = [R O{reg}];  %R: reference nodes (only outflow considered) 
        
        MP_nodesAll = [MP_nodesAll MP_nodes{reg}]; %all indices in MP struct of all nodes of all regions in one vector
        C = [C node_congested{reg}];
    else
        
        %decide on nodes - critical: over limits 2
        ind_1 = (node_queues{reg}>lim_occ{reg}(2)); %ind of nodes with occupancy above theshold (congestion)
        ind_2 = (node_var{reg}>lim_var{reg}(2)); %ind of nodes with variance above threshold
        ind_3 =  (node_congested{reg}>cong_threshold); %extra threshold of congestion
        ind_41 = MP.splan(MP_nodes{reg})==1; %caution: the set of not signalized nodes (~ind_41) is still ploted in the figures
        
        ind_4 = ind_1 + ind_2 + ind_3 + ind_41 == 4; %nodes satisfying all criteria simultaneously (AND operator)
        MPnodes = [MPnodes MP_nodes{reg}(ind_4)];  %store here the MaxPressure nodes (ind in MP)
        MPnodes = unique(MPnodes);
        
        %   medium + critical: over limits 1
        %      ind_1 = (node_queues{reg}>lim_occ{reg}(1));
        %      ind_2 = (node_var{reg}>lim_var{reg}(1));
        %      ind_3 =  (node_congested{reg}>cong_threshold); %extra threshold of congestion
        %      ind_4 = ind_1 + ind_2 +ind_3 == 3;
        %      MPnodes = [MPnodes MP_nodes{reg}(ind_4)];  %store here the MaxPressure nodes (ind in MP)
        %      MPnodes = setxor(MPnodes, intersect(MPspInt,MPnodes));
        
        %previous figures per region
        figure()
        title(strcat('Region ',num2str(reg)))
        scatter(node_var{reg}(ind_4),node_queues{reg}(ind_4),'filled','MarkerFaceColor',colReg{reg}) %selected
        xlabel('$m_2$ (variance)')
        ylabel('$m_1$ (congestion)')
        hold on
        scatter(node_var{reg}(~ind_4),node_queues{reg}(~ind_4),'filled','y') %not selected
        scatter(node_var{reg}(~ind_3),node_queues{reg}(~ind_3),'filled','m') %not congested
        scatter(node_var{reg}(~ind_41),node_queues{reg}(~ind_41),'filled','MarkerFaceColor',[17 17 17]/255) %not signalized (make invisible)
        
        x = (0:0.01:round(max(node_var{reg}+0.04),1));
        y = (0:0.01:round(max(node_queues{reg}+0.04),1));
        xlim([0 round(max(node_var{reg}+0.04),1)])
        ylim([0 round(max(node_queues{reg}+0.04),1)])
        %plot(x,ones(1,length(x)).*lim_occ{reg}(1),'--','LineWidth',1.5,'Color','r')
        plot(x,ones(1,length(x)).*lim_occ{reg}(2),'--','LineWidth',1.5,'Color','r')
        %plot(ones(1,length(y)).*lim_var{reg}(1),y,'--','LineWidth',1.5,'Color','r')
        plot(ones(1,length(y)).*lim_var{reg}(2),y,'--','LineWidth',1.5,'Color','r')
        legend(strcat('region ',num2str(reg)),'not selected','not congested')
        %saveas(gcf,strcat('clas_reg',num2str(reg),'_',scenario,'.fig'));
        %print(strcat('clas_reg',num2str(reg),'_',scenario),'-depsc','-painters')
        
        %common figures
        figure(f1)
        scatter(node_var{reg}(ind_4),node_queues{reg}(ind_4),'filled','MarkerFaceColor',colReg{reg})
        not_picked{reg} = (~ind_4);
        not_signalized{reg} = (~ind_41);
        %scatter(node_var{reg}(~ind_4),node_queues{reg}(~ind_4),'filled','y')
        
        figure(f2)
        scatter(node_var{reg}(ind_4), node_congested{reg}(ind_4) ,'filled','MarkerFaceColor',colReg{reg})
        scatter(node_var{reg}(~ind_4), node_congested{reg}(~ind_4) ,'filled','y')
        
        not_congested{reg} = (~ind_3);
    end
end


%Selection by the function A
if mode_A == 1
    % load the eligible MP nodes
    elig_MPnodes = load('MPnodesCase_MP_0','MPnodes');
    
    %define percentage target (first x nodes)
    num_Top = round(target_p*num_MPall);
    
    % exclude non eligible nodes from selection -
    elig_ind = ismember(MP_nodesAll,elig_MPnodes.MPnodes); %1-0 indicator of all eligible nodes
    B = B(elig_ind);
    C = C(elig_ind);
    MP_nodesAll = MP_nodesAll(elig_ind);
    node_queuesAll = node_queuesAll(elig_ind);
    node_varAll = node_varAll(elig_ind);
    node_congAll = node_congAll(elig_ind);
    R = node_outflowsAll(elig_ind); 
    
    %filter by Nc
    if mode_C==1
        filtered = C >= cong_threshold;
        %exclude nodes with Nc < congested_threshold
        B = B(filtered);
        MP_nodesAll = MP_nodesAll(filtered);
        node_queuesAll_f = node_queuesAll(filtered);
        node_varAll_f = node_varAll(filtered);
        node_congAll_f = node_congAll(filtered);
        
        
        if length(B) < num_Top
            error('Filtered subset not large enough. Decrease Nc_target or increase percentage target')
        end
        
        %sort B (reduced)
        [B, Ib] = sort(B); %min to max
        
        MPnodes = MP_nodesAll(Ib(1:num_Top));  %store here the MaxPressure nodes (ind in MP)
        MPnodes = unique(MPnodes);
        
        %print selection according to the proposed method A 
        figure(f1)
        scatter(node_varAll_f(Ib(1:num_Top)),node_queuesAll_f(Ib(1:num_Top)),'filled') %selected
        scatter(node_varAll_f(Ib(num_Top+1:end)),node_queuesAll_f(Ib(num_Top+1:end)),'filled','MarkerFaceColor','y') %not selected
        scatter(node_varAll(~filtered),node_queuesAll(~filtered),'filled','MarkerFaceColor','y') %not congested (filtered out)
        %saveas(gcf,strcat('sel_',num2str(1*(mode_A==0)+2*(mode_A==1)+(mode_C==1)),'_',num2str(round(100*target_p)),'a.emf'));
        
        figure(f2)
        scatter(node_varAll_f(Ib(1:num_Top)), node_congAll_f(Ib(1:num_Top)),'filled')
        scatter(node_varAll_f(Ib(num_Top+1:end)), node_congAll_f(Ib(num_Top+1:end)),'filled','MarkerFaceColor','y')
        scatter(node_varAll(~filtered), node_congAll(~filtered),'filled','MarkerFaceColor','y')
        %saveas(gcf,strcat('sel_',num2str(1*(mode_A==0)+2*(mode_A==1)+(mode_C==1)),'_',num2str(round(100*target_p)),'b.emf'));
        
    else
        
        [~, Ib] = sort(B);
        [~, Ir] = sort(R,'descend'); 
        
        MPnodes = MP_nodesAll(Ib(1:num_Top));  %store here the MaxPressure nodes (ind in MP)
        MPnodes = unique(MPnodes);
        
        %MP node selection with the reference way (max outflow) 
        MPnodes_R = MP_nodesAll(Ir(1:num_Top));  %store here the MaxPressure nodes (ind in MP)
        MPnodes_R = unique(MPnodes_R);        
        
        figure(f1)
        scatter(node_varAll(Ib(1:num_Top)),node_queuesAll(Ib(1:num_Top)),'filled') %selected
        scatter(node_varAll(Ib(num_Top+1:end)),node_queuesAll(Ib(num_Top+1:end)),'filled','MarkerFaceColor','y') %not selected
        
        %saveas(gcf,strcat('sel_',num2str(1*(mode_A==0)+2*(mode_A==1)+3*(mode_C==1)),'_',num2str(round(100*target_p)),'_',demand,'_','a.fig'));
%         saveas(gcf,strcat('selectionMP',demand,num2str(round(100*target_p)),'_figm1m2.fig'));
%         saveas(gcf,strcat('selectionMP',demand,num2str(round(100*target_p)),'_figm1m2.emf'));
%         print(strcat('selectionMP',demand,num2str(round(100*target_p)),'_figm1m2.eps'),'-depsc','-painters')
        
        figure(f2)
        scatter(node_varAll(Ib(1:num_Top)), node_congAll(Ib(1:num_Top))/CyclesInPeak,'filled')
        scatter(node_varAll(Ib(num_Top+1:end)), node_congAll(Ib(num_Top+1:end))/CyclesInPeak,'filled','MarkerFaceColor','y')
        
%         saveas(gcf,strcat('selectionMP',demand,num2str(round(100*target_p)),'_figNcm2.fig'));
%         saveas(gcf,strcat('selectionMP',demand,num2str(round(100*target_p)),'_figNcm2.emf'));
%         print(strcat('selectionMP',demand,num2str(round(100*target_p)),'_figNcm2.eps'),'-depsc','-painters')

        figure(f3)
        scatter(node_varAll(Ir(1:num_Top)),node_queuesAll(Ir(1:num_Top)),'filled') %selected
        scatter(node_varAll(Ir(num_Top+1:end)),node_queuesAll(Ir(num_Top+1:end)),'filled','MarkerFaceColor','y') %not selected
        
        
        figure(f4)
        scatter(node_varAll(Ir(1:num_Top)), node_congAll(Ir(1:num_Top))/CyclesInPeak,'filled')
        scatter(node_varAll(Ir(num_Top+1:end)), node_congAll(Ir(num_Top+1:end))/CyclesInPeak,'filled','MarkerFaceColor','y')
        
    end
end


%Exclude special intersections (MP cannot be applied due to signal
%settings)
MPnodes = setxor(MPnodes, intersect(MPspInt,MPnodes));
MPnodes_R = setxor(MPnodes_R, intersect(MPspInt,MPnodes_R));

% Check the percentage of selected nodes of the present scheme
p_1 = length(MPnodes)/num_MPall*100
p_1r = length(MPnodes_R)/num_MPall*100
%should be around our target/if not, adjust limits%

if mode_A == 0 %figures for TRB case
    
    for reg=1:3
        figure(f1)
        scatter(node_var{reg}(not_picked{reg}),node_queues{reg}(not_picked{reg}),'filled','y')
        scatter(node_var{reg}(not_congested{reg}),node_queues{reg}(not_congested{reg}),'filled','y')
        scatter(node_var{reg}(not_signalized{reg}),node_queues{reg}(not_signalized{reg}),'MarkerEdgeColor','y','MarkerFaceColor','w')
        
        figure(f2)
        scatter(node_var{reg}(not_picked{reg}),node_congested{reg}(not_picked{reg}),'filled','y')
        scatter(node_var{reg}(not_congested{reg}),node_congested{reg}(not_congested{reg}),'filled','y')
        scatter(node_var{reg}(not_signalized{reg}),node_congested{reg}(not_signalized{reg}),'MarkerEdgeColor','y','MarkerFaceColor','w')
    end
    
    %figure putting all nodes in two common graphs
    figure(f1)
    legend('region 1', 'region 2', 'region 3','not selected','not congested')
    %figure(f1)
    %saveas(gcf,strcat('sel_',num2str(1*(mode_A==0)+2*(mode_A==1)+(mode_C==1)),'_',num2str(round(p_1)),'a.emf'));
    %print(strcat('classification1_',scenario),'-depsc','-painters')
    figure(f2)
    legend('region 1', 'region 2', 'region 3','not selected','not congested')
    %saveas(gcf,strcat('sel_',num2str(1*(mode_A==0)+2*(mode_A==1)+(mode_C==1)),'_',num2str(round(p_1)),'b.emf'));
    
    % legend('region 1','not selected', 'region 2','not selected', 'region 3','not selected')
    %figure(f2)
    %saveas(gcf,strcat('classification2_',scenario,'.fig'));
    %print(strcat('classification2_',scenario),'-depsc','-painters')
    
    %how many of the preselected nodes are listed
%     test = sum((ismember(nodesPreSelected,MPnodes)))
    save(strcat('MPnodesCase_',caseName,'.mat'),'MPnodes','p_1','lim_occ','lim_var','cong_threshold')
else
    %Save selected set of nodes:
    %save(strcat('MPnodesCase_',caseName,'.mat'),'MPnodes','p_1','a','b')
    MPnodes = MPnodes_R; 
    save(strcat('MPnodesCase_',caseName,'.mat'),'MPnodes','p_1','a','b')
end


%% [optional] Plot the selected nodes on map  (use different colors) [optional - uncomment to print map]

load LINKS %links info [id, start_node, 3nd_node, length, cluster (in 4 reg)]
load nodes3 %nodes (ids, x coord, y coord)
load clus_final %clustering for links
nodes = nodes3;
clear nodes3
figure;
hold on;
for i=1:size(clus_final,1)
    ind=clus_final(i,1);
    gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off

% Print set of nodes with specific color
nodesC = MPnodes_R;
for i=1:length(nodesC)
    plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',25)
end
%title(caseName)
% saveas(gcf,strcat('mapNodes',caseName,'.fig'));
% saveas(gcf,strcat('mapNodes',caseName,'.emf'));
% print(strcat('mapNodes',caseName,'.eps'),'-depsc','-painters')

%close all 
%Plot colormap for the delay savings per link for one scenario -
end


