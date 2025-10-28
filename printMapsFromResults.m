% Print colormaps of the network based on results
clear
clc
close all
load scenarios.mat

scen_IDAll = 33; %set the case index in scenariosf (high) 
%scen_IDAll = 2; 
for ii = 1:length(scen_IDAll)
    scen_ID = scen_IDAll(ii);
    % high: 2, 18, 19, 33, 55
    % med: 2, 15, 45, 14,44
    demand = 2; %1 = med, 2 = high
    %% Set the path where the result file is ---
    if demand ==1
        path_g = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\'; %medium OD
        %name of results file
        fname = strcat('output_set5_MedDemand_',scenarios{scen_ID},'.mat'); %medium OD
        %name of input file (with MP and PC info)
        pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';
        %load reference case (for calculating differences)
        fname0 = strcat('output_set5_MedDemand_',scenarios{1},'.mat'); %medium OD
        
    elseif demand == 2
        path_g = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\'; %high OD
        fname = strcat('output_set5_highDemand_',scenarios{scen_ID},'.mat'); %high OD
        pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs high demand\final runs input files\';
        fname0 = strcat('output_set5_highDemand_',scenarios{1},'.mat'); %high OD
    end
    
    MPnodeSize = 40;
    PCnodeSize = 12;
    %setTimeAll = [0.75 2.01 2.75 6]; %time points of snapshots (in hours)
    setTimeAll = [0.75 2.01 2.75 4 8]; %end of simulation for high demand
    setTimes = setTimeAll;
    scenarios{scen_ID}
    
    fnameIn = strcat(pathIn, 'input_',scenarios{scen_ID});
    load(fnameIn, 'PC','MPnodes')
    
    mkdir(strcat(path_g,scenarios{scen_ID},'\'));
    path = strcat(path_g,scenarios{scen_ID},'\');
    
    
    load(strcat(path_g,fname0),'outdata','indata')
    ref =  cumsum(outdata.u*indata.DT,2);
    ref2 =  cumsum(outdata.x*indata.DT,2);
    ref3 =  cumsum(outdata.w*indata.DT,2);
    
    % ----------------------
    % calculate average link speeds for time window of x cycles before the time point
    % based on u and x:
    x_cycles = 4;
    speeds = ones(size(outdata.x,1),length(setTimeAll))*indata.v_ff/1000;
    for i=1:length(setTimeAll)
        % for x cycles before the time point (time window)
        % k_start = ceil(setTimeAll(i)/indata.DT)-x_cycles*90;
        % k_stop = min(indata.kmax,ceil(setTimeAll(i)/indata.DT)); %1 cycle is 90 sec and dt=1 sec
        
        % calculate mean speed from start every time
        k_start = 1;
        k_stop = floor(setTimeAll(i)/indata.DT);
        
        is_empty = sum(outdata.x(:,(k_start):k_stop),2);
        is_empty_ind = find(is_empty==0); %virtual links are always empty
        speeds(is_empty_ind,i) = indata.v_ff/1000; %km/h
        
        non_empty_ind = find(is_empty>0);
        % define cost fundtion per link for non empty links
        speeds(non_empty_ind,i) = sum(outdata.u(non_empty_ind,k_start:k_stop),2).*indata.Links2(non_empty_ind,3)...
            ./(1000*sum(outdata.x(non_empty_ind,k_start:k_stop),2)); %km/hour
        
        speeds(speeds(:,i)>indata.v_ff/1000,i) = indata.v_ff/1000;
        
        if sum(speeds(:,i)<0)>0
            disp('error speeds')
        end
        
    end
    ref4 = speeds;
    %--------------------------------------------------------
    save('ref_med','ref','ref2','ref3','ref4')
    
    %load results file
    load(strcat(path_g,fname),'outdata','indata')
    
    %% calculate speeds from results
    speeds = ones(size(outdata.x,1),length(setTimeAll))*indata.v_ff/1000;
    for i=1:length(setTimeAll)
        %k_start = ceil(setTimeAll(i)/indata.DT)-x_cycles*90;
        %k_stop = min(indata.kmax,ceil(setTimeAll(i)/indata.DT)); %1 cycle is 90 sec and dt=1 sec
        
        % calculate mean speed from start every time
        k_start = 1;
        k_stop = floor(setTimeAll(i)/indata.DT);
        
        
        is_empty = sum(outdata.x(:,(k_start):k_stop),2);
        is_empty_ind = find(is_empty==0); %virtual links are always empty
        speeds(is_empty_ind,i) = indata.v_ff/1000; %km/h
        
        non_empty_ind = find(is_empty>0);
        % define cost fundtion per link for non empty links
        speeds(non_empty_ind,i) = sum(outdata.u(non_empty_ind,k_start:k_stop),2).*indata.Links2(non_empty_ind,3)...
            ./(1000*sum(outdata.x(non_empty_ind,k_start:k_stop),2)); %km/hour
        
        speeds(speeds(:,i)>indata.v_ff/1000,i) = indata.v_ff/1000;
        
        if sum(speeds(:,i)<0)>0
            disp('error speeds')
        end
        
    end
    
    %% Print maps
    %clear all
    %close all
    load LINKS %links info [id, start_node, end_node, length, cluster (in 4 reg)]
    load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids)
    load clus_final %clustering for links
    load control %controlled intersections
    nodes = nodes3;
    
    %% print map with 3 clusters in different colors
    
    % figure('Name','Node selection C1')
    % hold on;
    % for i=1:size(clus_final,1)
    %     ind=clus_final(i,1);
    %     gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    %     hold on;
    % end
    % axis off
    
    
    
    %% Print specified set of nodes with specific color (all same length)
    
    % %Nodes of Selection C1 (Max-Pressure)
    % load('nodesC1n.mat')
    % hold on
    % for i=1:length(nodesC)
    %     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#7E2F8E','Markersize',35)
    % end
    % title('MP Nodes - C1')
    % saveas(gcf,strcat(path,'mapWithNodesC1n.fig'));
    % print(strcat(path,'mapWithNodesC1n'),'-depsc','-painters')
    
    %% Map with total link outflow (links serving the highest demand overall)
    
    % links_u = sum(outdata.u,2);
    % figure('Name','Total Outflow per link');
    % hold on;
    % for i=1:size(LINKS,1)
    %     ind = links_u(i)/max(links_u);
    %     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    %     hold on;
    % end
    % colorbar
    % axis off
    % title('Total link outflow')
    % saveas(gcf,strcat(path,'mapOutflowPerLink.fig'));
    % % saveas(gcf,'mapOutflowPerLink.emf');
    % print(strcat(path,'mapOutflowPerLink'),'-dmeta','-painters')
    
    %% Total VHT in virtual queues as circles on the start nodes
    % Plus total VHT spent in links
    
    % VHT_tot = sum(outdata.x,2)*indata.DT;
    % VHT_tot(indata.group2) = sum(outdata.w(indata.group2,:),2)*indata.DT;
    % figure('Name','Map total VQ-VHT per centroid & VHT per link');
    % hold on;
    % for i=1:size(clus_final,1)
    %     ind = VHT_tot(i)/max(VHT_tot(1:size(clus_final,1)));
    %     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    %     hold on;
    % end
    % axis off
    %
    % VHT_VQs = sum(outdata.virtualqueues,2);
    % %figure('Name','Map total VQ-VHT per centroid');
    % % hold on;
    % % for i=1:size(clus_final,1)
    % %     ind=clus_final(i,1);
    % %     gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    % %     hold on;
    % % end
    % % axis off
    % %Print specified set of nodes with specific color (all same length)
    % nodesVQs = indata.Links2(indata.group2,5);
    % myColorMap = jet(256);
    % colormap(myColorMap)
    % for i=1:length(nodesVQs)
    %     plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),...
    %         'Marker','.','Color',myColorMap(max(1,ceil(VHT_VQs(i)/max(VHT_VQs)*256)),:),...
    %         'Markersize',30)
    %     %     plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i)...
    %     %           ,3),'Marker','.','Color','#A2142F','Markersize',30)
    % end
    % title('VHT spent in VQs')
    % colorbar
    % saveas(gcf,strcat(path,'mapVHT_VQ_PerLink.fig'));
    % % saveas(gcf,'mapVHT_VQ_PerLink.emf');
    % print(strcat(path,'mapVHT_VQ_PerLink'),'-depsc','-painters')
    
    
    %% Total outflow of virtual queues as circles on the start nodes
    %
    % tot_VQs = sum(outdata.u(indata.group2,:),2);
    % figure('Name','Map total VQ outflow per centroid');
    % hold on;
    % for i=1:size(clus_final,1)
    %     ind=clus_final(i,1);
    %     gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    %     hold on;
    % end
    % axis off
    %
    % %Print specified set of nodes with specific color (all same length)
    % nodesVQs = indata.Links2(indata.group2,5);
    % myColorMap = jet(256);
    % colormap(myColorMap)
    % for i=1:length(nodesVQs)
    %     plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),'Marker','.','Color',myColorMap(max(1,ceil(tot_VQs(i)/max(tot_VQs)*256)),:),'Markersize',30)
    % end
    % title('Total VQ outflow')
    % colorbar
    % saveas(gcf,strcat(path,'mapTotal_VQ_Nodes.fig'));
    % % saveas(gcf,'mapVHT_VQ_PerLink.emf');
    % print(strcat(path,'mapTotal_VQ_Nodes'),'-depsc','-painters')
    
    
    %% Total time spent overall per link (already printed before)
    
%     VHT_tot = sum(outdata.x,2)*indata.DT;
%     VHT_tot(indata.group2) = sum(outdata.w(indata.group2,:),2)*indata.DT;
%     figure('Name','Total VHT per link');
%     hold on;
%     for i=1:size(clus_final,1)
%         ind = VHT_tot(i)/max(VHT_tot(1:size(clus_final,1)));
%         gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%         hold on;
%     end
%     axis off
%     colorbar
%     title('VHT per link')
%    
%     saveas(gcf,strcat(path,'mapVHT_PerLink_1.fig'));
%     % saveas(gcf,strcat(path,'mapVHT_PerLink_1.emf'));
%     print(strcat(path,'mapVHT_PerLink_2'),'-depsc','-painters')
    
    %% Highest average density/occupancy
    % Avg_occ = mean(outdata.x,2)./indata.capacity;
    % figure('Name','Map Mean Occupancy per link');
    % hold on;
    % for i=1:size(clus_final,1)
    %     ind = Avg_occ(i)/max(Avg_occ);
    %     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    %     hold on;
    % end
    % axis off
    % colorbar
    % title('Mean occupancy per link')
    % saveas(gcf,strcat(path,'mapOccupancyPerLink.fig'));
    % % saveas(gcf,'mapOccupancyPerLink.emf');
    % print(strcat(path,'mapOccupancyPerLink'),'-depsc','-painters')
    %
    
    %% print specific snapshot for link occupancy
    
    % snapshot occupancy
    
    for j=1:length(setTimeAll)
        setTime = setTimeAll(j);
        createSnapshot_2(setTime,indata.DT,outdata.x,indata.capacity,'link occupancy $x/c$',LINKS,nodes3)
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
        
        saveas(gcf,strcat(path,'mapOccupancyPeak',num2str(j),'.fig'));
        print(strcat(path,'mapOccupancyPeak',num2str(j)),'-depsc','-painters')
    end
    
    %% print map with average link occupancy (from the beginning of time)
    
    % average peak occupancy from start
    
    cummean = zeros(size(outdata.x));
    for j=1:size(setTimeAll,2)
        setTime = setTimeAll(j);
        ind = floor(setTime/indata.DT);
        cummean(:,ind) = mean(outdata.x(:,1:ind),2);
    end
    for j=1:length(setTimeAll)
        setTime = setTimeAll(j);
        createSnapshot_2(setTime,indata.DT,cummean,indata.capacity,'mean link occupancy $x/c$',LINKS,nodes3)
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);

        saveas(gcf,strcat(path,'mapMeanOccupancyFromStart_Peak',num2str(j),'.fig'));
        print(strcat(path,'mapMeanOccupancyFromStart_Peak',num2str(j)),'-depsc','-painters')
    end
    
    %% print map with cumulative waiting time (from the beginning of time)
    %
    % cumulative waiting time
    
    cummean = zeros(size(outdata.w));
    for j=1:size(setTimeAll,2)
        setTime = setTimeAll(j);
        ind = floor(setTime/indata.DT);
        cummean(:,ind) = sum(outdata.w(:,1:ind),2);
    end
    
    % normUp = ones(size(setTimeAll))*65000;
    normUp = [67500 470000 697000 870000]; %for NC high
    a1min = 0;
    for j=1:length(setTimeAll)
        setTime = setTimeAll(j);
    
        %figure
        %boxplot(cummean(:,floor(setTime/indata.DT)));
    
        createSnapshot_2(setTime,indata.DT,cummean,ones(size(cummean,1),1)*normUp(j),'cumulative waiting time $w$',LINKS,nodes3)
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
    
    
        c = colorbar('TickLabelInterpreter' ,'latex');
        c.Label.String = 'cumulative waiting time w (veh h)'; %'link occupancy $x/c$';
        c.Label.Interpreter = 'latex';
        c.TickLabels = a1min:(normUp(j)-a1min)/10:normUp(j);
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);

        saveas(gcf,strcat(path,'mapCumWaitingTime_Peak',num2str(j),'.fig'));
        saveas(gcf,strcat(path,'mapCumWaitingTime_Peak',num2str(j),'.emf'));
        print(strcat(path,'mapCumWaitingTime_Peak',num2str(j)),'-depsc','-painters')
    end
    close all
    %% Print snapshots for difference in cummulative throughput (cumulative outflow of links)wrt NC + add MP/PC nodes
    
    %cumulative throughput difference
    
    normalFactors = ones(1,4)*700; %max values for normalization
%     normalFactors = [0.15 1.5 5 15]*10^6;
%     a1min = [-0.15 -1.5 -5 -15]*10^6; %min values (lower limit of colorbar)
    
    a1min = ones(1,4)*(-700); 
    load ref_med
    a1 = cumsum(outdata.u*indata.DT,2); %cumulative link outflows
    
    for j=1:length(setTimes)
        
        setTime = setTimes(j);
        
        %a1max = max(a1);
        %a1max = ones(size(a1,1),1)*a1max(floor(setTime/indata.DT)); %normalize
        %a1max = ones(size(a1,1),1)*normalFactors(3);
        
        a2 = a1 - ref;
        %a2 = cumsum(a1,2)-cumsum(ref,2);
        %     a1min = min(a2);
        
        min(min(a2))
        max(max(a2))
        a2 = a2 + ones(size(a2)).*abs(a1min(j)) + 0.01;
        
%         figure
%         a3 = cumsum(a1,2)-cumsum(ref,2);
%         boxplot(a3(:,floor(setTime/indata.DT)))
        
        
        %a1max = ones(size(a1,1),1)*max(a2(:,ceil(setTime/indata.DT)));
        a1max = ones(size(a1,1),1)*normalFactors(j) + ones(size(a1,1),1).*abs(a1min(j)) + 0.01;
        
        createSnapshot_2(setTime,indata.DT,a2,a1max,'diference in cumulative u (veh/link)',LINKS,nodes3)
        
        
        
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        
        c = colorbar('TickLabelInterpreter' ,'latex');
        c.Label.String =' diference in cumulative output (veh/link)'; %'link occupancy $x/c$';
        c.Label.Interpreter = 'latex';
        %c.TickLabels = round((0:normalFactors(3)/10:normalFactors(3)));
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        c.TickLabels = a1min(j):(normalFactors(j)-a1min(j))/10:normalFactors(j);
        
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
        
        saveas(gcf,strcat(path,'mapThroughputDifPeak',num2str(j),'.fig'));
        saveas(gcf,strcat(path,'mapThroughputDifPeak',num2str(j),'.emf'));
        print(strcat(path,'mapThroughputDifPeak',num2str(j)),'-depsc','-painters')
        
        
    end
    
    %% Print snapshots for difference in cummulative time spent (vehxhours) wrt NC + add MP/PC nodes
    %
    
    %cumulative time spent difference
    
    normalFactors = ones(1,4)*20; %max values for normalization
    load ref_med
    a1 = cumsum(outdata.x*indata.DT,2); %cumulative link outflows
    
    for j=1:length(setTimes)
        
        setTime = setTimes(j);
        
        %a1max = max(a1);
        %a1max = ones(size(a1,1),1)*a1max(floor(setTime/indata.DT)); %normalize
        %a1max = ones(size(a1,1),1)*normalFactors(3);
        
        a2 = a1-ref2;
        %     a1min = min(a2);
        a1min = -20; %min values (lower limit of colorbar)
        min(min(a2))
        max(max(a2))
        a2 = a2 + ones(size(a2)).*abs(a1min) + 0.01;
        %a1max = ones(size(a1,1),1)*max(a2(:,ceil(setTime/indata.DT)));
        a1max = ones(size(a1,1),1)*normalFactors(j) + ones(size(a1,1),1).*abs(a1min) + 0.01;
        
        createSnapshot_2(setTime,indata.DT,a2,a1max,'diference in cumulative time spent (veh hours /link)',LINKS,nodes3)
        
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        
        c = colorbar('TickLabelInterpreter' ,'latex');
        c.Label.String ='diference in cumulative time spent (veh hours /link)'; %'link occupancy $x/c$';
        c.Label.Interpreter = 'latex';
        %c.TickLabels = round((0:normalFactors(3)/10:normalFactors(3)));
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        c.TickLabels = a1min:(normalFactors(j)-a1min)/10:normalFactors(j);
        
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
        
        
        saveas(gcf,strcat(path,'mapTimeSpentDifPeak',num2str(j),'.fig'));
        saveas(gcf,strcat(path,'mapTimeSpentDifPeak',num2str(j),'.emf'));
        print(strcat(path,'mapTimeSpentDifPeak',num2str(j)),'-depsc','-painters')
        
%             figure
%             a2 = a1 - ref2;
%             boxplot(a2(:,floor(setTime/indata.DT)))
    end
    close all
    %% Print shapshot for difference in cummulative waiting time
    
    %cumulative waiting time difference
    
    %normalFactors = ones(1,4)*1000; %max values for normalization [set!]
    normalFactors = [10 20 50 50];
    a1min = [-10 -20 -50 -50]; %min values (lower limit of colorbar)
    load ref_med
    a1 = cumsum(outdata.w*indata.DT,2); %cumulative link outflows
    
    for j=1:length(setTimes)
        
        setTime = setTimes(j);
        
        %a1max = max(a1);
        %a1max = ones(size(a1,1),1)*a1max(floor(setTime/indata.DT)); %normalize
        %a1max = ones(size(a1,1),1)*normalFactors(3);
        
        a2 = a1 - ref3;
        %     a1min = min(a2);
        
        
        %     figure
        %     a3 = a1 - ref3;
        %     boxplot(a3(:,floor(setTime/indata.DT)))
        %
        min(min(a2))
        max(max(a2))
        a2 = a2 + ones(size(a2)).*abs(a1min(j)) + 0.01;
        %a1max = ones(size(a1,1),1)*max(a2(:,ceil(setTime/indata.DT)));
        a1max = ones(size(a1,1),1)*normalFactors(j) + ones(size(a1,1),1).*abs(a1min(j)) + 0.01;
        
        createSnapshot_2(setTime,indata.DT,a2,a1max,'diference in cumulative w (veh/link)',LINKS,nodes3)
        
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        
        c = colorbar('TickLabelInterpreter' ,'latex');
        c.Label.String =' diference in cumulative waiting time (veh/link)'; %'link occupancy $x/c$';
        c.Label.Interpreter = 'latex';
        %c.TickLabels = round((0:normalFactors(3)/10:normalFactors(3)));
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        %c.TickLabels = a1min(ceil(setTime/indata.DT)):(a1max(1)-a1min(ceil(setTime/indata.DT)))/10:a1max(1);
        c.TickLabels = a1min(j):(normalFactors(j)-a1min(j))/10:normalFactors(j);
        
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
        
        saveas(gcf,strcat(path,'mapCumWaitingTimeDifPeak',num2str(j),'.fig'));
        saveas(gcf,strcat(path,'mapCumWaitingTimeDifPeak',num2str(j),'.emf'));
        print(strcat(path,'mapCumWaitingTimeDifPeak',num2str(j)),'-depsc','-painters')
        
        
    end
    
    
    %% average speed since start / difference wrt FTC
    
    
    normalFactors = ones(1,length(setTimes))*5;
    a1min = ones(1,length(setTimes))*-5; %min values (lower limit of colorbar)
    load ref_med
    a1 = speeds; %cumulative link outflows
    
    for j=1:length(setTimes)
        
        setTime = setTimes(j);
        
        a2 = a1 - ref4;
        
        min(min(a2))
        max(max(a2))
        a2 = a2 + ones(size(a2)).*abs(a1min(j)) + 0.01;
        %a1max = ones(size(a1,1),1)*max(a2(:,ceil(setTime/indata.DT)));
        a1max = ones(size(a1,1),1)*normalFactors(j) + ones(size(a1,1),1).*abs(a1min(j)) + 0.01;
        
        createSnapshot_2(j,1,a2,a1max,'diference in mean link speed wrt FTC (km/h)',LINKS,nodes3)
        
        %MP nodes
        if ~isempty(MPnodes)
            nodesC = MPnodes;
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#106e3a','Markersize',MPnodeSize)
            end
        end
        %PC nodes
        if PC.mode==1
            nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
            for i=1:length(nodesC)
                plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',PCnodeSize)
            end
        end
        
        c = colorbar('TickLabelInterpreter' ,'latex');
        c.Label.String = 'diference in mean link speed wrt FTC (km/h)';
        c.Label.Interpreter = 'latex';
        c.TickLabels = a1min(j):(normalFactors(j)-a1min(j))/10:normalFactors(j);
        
        nSecs = setTime*3600;
        text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
        
        saveas(gcf,strcat(path,'mapMeanLinkSpeedDifPeak',num2str(j),'.fig'));
        saveas(gcf,strcat(path,'mapMeanLinkSpeedDifPeak',num2str(j),'.emf'));
        print(strcat(path,'mapMeanLinkSpeedDifPeak',num2str(j)),'-depsc','-painters')
        
    end
    
    close all
end

%% Nodes that we control with PC + total outflow per link
% to check if nodes are well selected

% load(strcat('input_','PC_0','.mat'),'PC')
% nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
% links_u = sum(outdata.u,2);
% figure('Name','Map PC Nodes & total link outflow');
% hold on;
% for i=1:size(LINKS,1)
%     ind = links_u(i)/max(links_u);
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% colorbar
% axis off
%
% for i=1:length(nodesC)
%     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',35)
% end
% title('PC nodes + total link outflow')
% saveas(gcf,strcat(path,'mapNodesPC_plus_LinkOutflow.fig'));
% % saveas(gcf,'mapNodesPC.emf');
% print(strcat(path,'mapNodesPC_plus_LinkOutflow'),'-depsc','-painters')
% % print('mapWithNodesC1n.eps','-depsc','-painters')

%% Mean link occupancy + PC controlled nodes

% load(strcat('input_','PC_0','.mat'),'PC')
% nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
% Avg_occ = mean(outdata.x,2)./indata.capacity;
% figure('Name','Map PC Nodes & average link occupancy');
% hold on;
% for i=1:size(LINKS,1)
%     ind = Avg_occ(i)/max(Avg_occ);
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% colorbar
% axis off
% %node printing
% for i=1:length(nodesC)
%     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',35)
% end
% title('PC nodes + average link occupancy')
% saveas(gcf,strcat(path,'mapNodesPC_plus_LinkOccupancy.fig'));
% % saveas(gcf,'mapNodesPC.emf');
% print(strcat(path,'mapNodesPC_plus_LinkOccupancy'),'-depsc','-painters')
% % print('mapWithNodesC1n.eps','-depsc','-painters')


%% Plot: colormap-> aggregated time during which every link overpasses 80%
% of capacity (spill-backs) - SAME figure as mean occupancy per link

% spillback_count = outdata.x(1:indata.NLinks,:)>0.8*indata.capacity(1:indata.NLinks);
% spillback_count = sum(spillback_count,2);
% figure('Name','Spill-back time per link');
% hold on;
% for i=1:indata.NLinks
%     ind = spillback_count(i)/max(spillback_count);
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% axis off
% colorbar
% title('Spill-back frequency per link')
% saveas(gcf,strcat(path,'mapSpillBacksPerLink.fig'));
% % saveas(gcf,'mapOccupancyPerLink.emf');
% print(strcat(path,'mapSpillBacksPerLink'),'-dmeta','-painters')

%% Results of individual run per region -------
%% Accumulation per region (including VQs) (common graph - single run only)
% figure(f16)
% for r=1:indata.no_reg
%
%     %network links (group1) - x
%     sum_1 = sum(outdata.x(intersect(indata.group1, indata.reInd{r}),:));
%     %virtual queues (group2) - w
%     sum_2 = sum(outdata.w(intersect(indata.group2, indata.reInd{r}),:));
%     sum_r = sum_1 + sum_2;
%
%     %perform aggregation of values
%     accumulation=0;
%     k_MFD=0;
%     for i=1:t:indata.kmax
%         k_MFD = [k_MFD i*indata.DT]; %time for the MFD [hours]
%         accumulation = [accumulation mean(sum_r(i:i+(t-1)))]; %[veh] - mean accumulation over the time window (VQ included)
%     end
%
%     plot(k_MFD, accumulation,'LineWidth',1.5)
% end
% legend('region 1','region 2','region 3')
% saveas(gcf,'AccumulationPerRegion.fig');
% % saveas(gcf,'mapOccupancyPerLink.emf');
% print('AccumulationPerRegion','-dmeta','-painters')

%% Print figures with results of specific node to analyze

% Specify node

% print TS of green of duration of the main and secondary phase

% print TS of agg_u (calculated and applied) for the approach

% Print comparative graph of accumulation in both regions

% Print comparative graph of transfer flow


%% Print Demand pattern (maps)

% % print color table (origin-destination based on demand)
% load('FinalInput.mat')
% arrangedOD = zeros(size(OD));
% or_centroids = OD(2:end,1);
% dest_centroids = OD(1,2:end);
%
% %re-arrange the centroids of every region
% %classify origin centroids in regions
% groupsOrCentroids = zeros(1,length(or_centroids));
% for i=1:length(or_centroids)
%    i1 = find(indata.Links2(:,4)==or_centroids(i),1);
%    groupsOrCentroids(i) =  indata.Links2(i1,end);
% end
% %classify destination centroids in regions
% groupsDestCentroids = zeros(1,length(dest_centroids));
% for i=1:length(dest_centroids)
%    i1 = find(indata.Links2(:,5)==dest_centroids(i),1);
%    if ~isempty(i1)
%      groupsDestCentroids(i) =  indata.Links2(i1,end);
%    end
% end
% k = 1;
% for i=1:no_reg
%     %origins
%     a = or_centroids(groupsOrCentroids==i);
%     for j=1:length(a)
%         %centroid ID - origins
%         arrangedOD(k+1,1) = a(j);
%         k = k + 1;
%     end
%
% end
% k = 1;
% for i=1:no_reg
%     %destinations
%     b = dest_centroids(groupsDestCentroids==i);
%     for j=1:length(b)
%         %centroid ID
%         arrangedOD(1,k+1) = b(j);
%         k = k + 1;
%     end
% end
% % assign the demand value
% for i=1:length(or_centroids)
%     r = find(arrangedOD(:,1)==or_centroids(i));
%     for j=1:length(dest_centroids)
%        c = find(arrangedOD(1,:)==dest_centroids(j));
%        arrangedOD(r,c) = OD(i+1,j+1);
%     end
% end
%
% % plot heatmap of OD matrix with
% % - different color intensity depending on the value of the demand
%
% figure('Name','Heatmap OD')
% a = arrangedOD(2:end,2:end);
% heatmap(a)


