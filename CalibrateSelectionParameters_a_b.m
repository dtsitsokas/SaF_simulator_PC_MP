% Code to calibrate a, b parameters of MP node selection - optimize results
% Algorithm:
% ----------
% Load NC case and the TRB selection of a specific node percentage (30%) 

clc
clear
close all

%Which OD? med or high? 
%Load the before case/actual situation (NC) - set the path for medium or high demand case (reference)  
path = 'E:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
%fname = strcat(path,'output_set4_test_NC'); %NC of high
fname = strcat(path,'output_set5_MedDemand_NC'); %NC of med
%input file of NC case - why? 
load(strcat('input_NC.mat'))

%% load results of NC reference case (before) 
load(fname,'indata','outdata')
load('FinalInput','MP')

%file name of reference MP case (after) 
fnameInput_ref = 'MPnodesCase_MP_cd_S4'; %the file with the node selection
fname_ref = 'output_set5_MedDemand_MP_S4'; %the results of this selection
perc = 0.25; %the percentage of the loaded file 

%total number of MP eligible nodes
num_MPall = sum(MP.splan==1)-length(indata.specialInt); 

k_s = 41*indata.c_int; %starting point of the peak interval
k_e = 120*indata.c_int; %ending point of the peak interval

mrkSize = 15; %size for scattters
colMP = 'b'; %color for MP nodes
colNonMP = 'g'; %color for nonMP nodes 

colReg{1} = '#0072BD';
colReg{2} = '#D95319';
colReg{3} = '#77AC30';
%% Store values of the reference case - NC
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
            
            %m_1:
            node_queuesNC{reg}(i) = mean(mean(node_occNC));
            
            %m_2
            %calculate mean of the variance of queues of all incoming links
            %over time
            node_varNC{reg}(i) = mean(var(outdata.x(indata.junct.or_index(MP.approaches{ind}),k_s:k_e)./indata.capacity(indata.junct.or_index(MP.approaches{ind}))));
            
            % count number of cycles N_c: 
            node_congestedNC{reg}(i) = sum(sum(node_occNC>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
            
        end
        
    end
    
    %node ID (index in MP) - node avg norm queue - node avg var of queues -
    %node no of high congestion cycles
    nodeMetricsNC{reg} = [MP_nodesNC{reg}(:) node_queuesNC{reg}(:) ...
        node_varNC{reg}(:) node_congestedNC{reg}(:)];
 
end

%% load examined case results
%load node information of the case: nodes included in the plan
cases = 5; %id of the selection case (0 for 100%, 1 for 10%, 2 for 20% etc) 
j = 1; %one case only 

load(fnameInput_ref,'MPnodes') %selection of nodes
node_queues = {};
node_var = {};
node_congested = {};
nodeMetricsTest = [];
MPnT = cell(1,3);
nonMPnT = cell(1,3);

%load case results after MP application
load(strcat(path,fname_ref),'outdata')

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
% title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
xlim([0 0.6])
ylim([0 1])
hold on
grid on

f5s(j,4) = subplot(2,2,4);
xlabel('$m_2$')
ylabel('N_c')
% title(strcat(num2str(cases(j)*5+5*(cases(j)>0)+100*(cases(j)==0)),'\% nodes'))
xlim([0 0.6])
ylim([0 80])
hold on
grid on

for reg=1:3
    figure(f5)
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
            node_congested{reg}(i) = sum(sum(node_occ>0.8)>0)/indata.c_int; %number of cycles (aggregated) when node is "congested"
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
end

metricsDifG_MPn = [];
metricsDifG_nonMPn = [];

%% Rank nodes according to the relationship A = a m_1 + b m_2 + (1-a-b) m_3

%Find a, b so that after ranking nodes A, the top x% of N_C - N_c*
%(or other obj function expression) is maximum
%assuming a, b are \in (0,1)
%load all eligible MP nodes
elig_MPnodes = load('MPnodesCase_MP_0','MPnodes');
I = [];
best_found = []; %to store results of tests a,b
%grid search 
k=0;
Rnod = {};
for a=0:0.2:2
    for b=0:0.2:2
        %for c=0:0.1:2
        k = k + 1; 
        nod = []; %node IDs in MP
        A = [];  %selection function
        sel = []; %selected
        med = []; %metrics difference
        for reg = 1:3
            nod = [nod; nodeMetricsTest{reg}(:,1)]; %node ind in MP (they are added by region)
            %A = [A; 1.2*a*(0.85-nodeMetricsTest{reg}(:,2))+1.2*b*(nodeMetricsTest{reg}(:,3)-0.03)+1.2*c*nodeMetricsTest{reg}(:,4)/80]; %selection function (to calibrate)
            A = [A; 1.2*a*(0.15-nodeMetricsTest{reg}(:,2))-1.2*b*(nodeMetricsTest{reg}(:,3)-0.03)/0.6-abs(1-a-b)/80*nodeMetricsTest{reg}(:,4)]; %selection function (to calibrate)

            %I need to exclude non-candidate nodes ?
            
            %A = [A; a*nodeMetricsTest{reg}(:,2) + b*nodeMetricsTest{reg}(:,3) ...
            %    + c*nodeMetricsTest{reg}(:,4)/80]; %selection function (to calibrate)
            
            
            %- the highest A is the highest avg queue, var and duration of high congestion
            %- the highest A, the highest the improvement potential
            med = [med;  metricsDif{reg}(:,3)]; %performance N_c-N_C* (column 3) 
            %med = [med;  metricsDif{reg}(:,2) + metricsDif{reg}(:,3)]; -
            % the higher the difference - the higher the degradation /
            % for improvement this should be negative
            % sel = [sel; (length(sel) + MPnT{reg})']; %indices of MPnodes in med, nod etc.
            sel = [sel; ismember(nodeMetricsTest{reg}(:,1),MPnodes)]; %logical array to say if selected or not
        end
        numTop = round(min(perc*length(elig_MPnodes.MPnodes),length(MPnodes))); %count the top 20% of the nodes
        
        %Exclude non-eligible MPnodes
        elig_ind = ismember(nod,elig_MPnodes.MPnodes); %ind of eligible
        %only keep the eligible nodes
        A = A(elig_ind);
        [A,Ia] = sort(A); %min to max
        %I = [I Ia];
        med = med(elig_ind);
        nod = nod(elig_ind);
        sel = sel(elig_ind);
        
        Rmed = med(Ia); %performance of nodes sorted according to A
        sum_top = sum(A(1:numTop)); %sum of A of the top x% in the ranking
        sum_perf_top = sum(Rmed(1:numTop)); %sum of observed performance of the resulting top percent
        sum_perf_selected = sum(med(sel==1)); %sum of observed performance of selected nodes
        deviation = (sum_perf_top - sum_perf_selected); 
        Rnod{k} = nod(Ia);
        %end
        best_found = [best_found; a b deviation sum_perf_top sum_perf_selected sum_top];
        %real test of the selected nodes "nod"       
    end
end
%plot performance of top x% of A + performance of selected nodes from
%previous
figure
hold on 
%plot(best_found(:,3))
plot(best_found(:,4))
% plot(best_found(:,5))
%plot(best_found(:,6))
%plot(best_found(:,3:end),'DisplayName','best_found(:,3:end)')
legend('predicted perf.','predicted perf. of reference','real perf. of case')
[min_s,ind] = min(best_found(:,3)) %somethig is wrong with this - the max is the solution (counter intuitive)
best_found(ind,:)
%%
ind_k = find(best_found(:,3)<1);
%ind_l = find(best_found(:,3)<10);
path = 'E:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
aa = best_found(best_found(:,3)>2000,3);
[ab,iab] = sort(aa,'descend');

% for i=41:60%length(ind_k)
%     fname = strcat(path,'test2_classification_max',num2str(i),'.mat');
%     savefull = 2; %save big result files
%     PC.mode = 0;
%     real_VHT = SaF_3(indata,Rnod{ind_k(iab(i))}(1:numTop),PC,fname,savefull);
% end

% for i=1:length(ind_l)
%     fname = strcat(path,'test_classification_min',num2str(i),'.mat');
%     savefull = 2; %save big result files
%     PC.mode = 0;
%     real_VHT = SaF_3(indata,Rnod{ind_l(i)}(1:numTop),PC,fname,savefull);
% end

%% 

% for i=1:60%length(ind_k)
%     load(strcat(path,'test2_classification_max',num2str(i),'.mat'))
%     real_VHT_max(i) = r3;
% end
% 
% % for i=1:length(ind_l)
% %     load(strcat(path,'test_classification_min',num2str(i),'.mat'))
% %     real_VHT_min(i) = r3; 
% % end
% 
% figure
% hold on 
% plot(real_VHT_max)
% yyaxis right 
% plot(ab(1:60))
% % plot(real_VHT_min)
% legend('real VHT','prediction function')
% %max classification leads to better results -> Why? // ind 17 and 20 are
% %best 






