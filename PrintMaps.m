% Print map of the Barcelona network
clear
%close all 

load LINKS %links info [id, start_node, end_node, length, cluster (in 4 reg)]
load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids)
load clus_final %clustering for links
load control %controlled intersections
nodes = nodes3;

%% read new clustering data: 
l1 = readmatrix("links_clusters_n.txt"); 
n1 = readmatrix("nodes_clusters_n.txt"); 
not_f = []; 
for i=1:size(LINKS,1)
    in = LINKS(i,1)==l1(:,1); 
    if sum(in)==0
        %disp('link not found')
        not_f = [not_f LINKS(i,1)]; %links not found
    else
        LINKS(i,6) = l1(in,2);  
    end
end
clus_final = LINKS(:,6);
clus_final = clus_final + 1; %eliminate 0 because it causes problems 
control = readmatrix("PCjunctionsID_n.txt"); 

%% load control %controlled intersections 
%path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
%load(strcat(path,'indata'))
load indata
%% print network map and nodes for one network layout  

figure('units','normalized','outerposition',[0 0 1 1],'color','white');
axis off
hold on
for i=1:size(clus_final,1)
    ind=clus_final(i,1); 
    gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end

% Print set of nodes with specific color 

%Print PC nodes 
fname = 'C:\Users\tsitsoka\Desktop\SaF3 Simulator LUTS\inputPCscenarios\input_PC_50b';
load(fname,'PC')
% nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
% for i=1:length(nodesC)
%    plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',30)
% end

Col{1}='#A2142F';
Col{2}='#7E2F8E';
Col{3}='#77AC30';
Col{4}='#4DBEEE';
Mark{1} = '^';
Mark{2} = 'o';
Mark{3} = 'square';
Mark{4} = 'diamond';

PCnodeSize = 12;
for i=1:4
    for j=1:size(PC.nodes,2)
        nodesC = PC.nodes(i,j);
        if nodesC>0
            plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC),3),...
                'Marker',Mark{i},'MarkerFaceColor',Col{1},'MarkerEdgeColor',Col{1},'Markersize',PCnodeSize)
        end
    end
end
ylim([4.5805*10^6 4.5836*10^6])
xlim([4.277*10^5 4.33*10^5])

path = 'F:\Max Pressure files\fifth set results\paper results figures\';
%saveas(gcf,strcat(path,'mapRegionsPCnodes','.fig'));
print(strcat(path,'mapRegionsPCnodes'),'-depsc','-painters')

%Print MP nodes 

% load('nodesC1n.mat')
% hold on
% for i=1:length(nodesC)
%     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',25)
% end

% saveas(gcf,'mapWithNodesC1n.fig');
% print('mapWithNodesC1n.eps','-depsc','-painters')

%% Print multifigure of different node layouts 
% % Case: Direct node assignment 
% cases = [2 3 4];
% % Case: Incremental node assignment 
% % cases = [2 4 6]; 
% 
% for j=1:length(cases)
%     load(strcat('MPnodesCase_MP_cd_MS',num2str(cases(j))),'MPnodes'); %MED DEMAND
%     %load(strcat('MPnodesCase_MP_grid_high_eq2t_MS',num2str(cases(j))),'MPnodes'); %HIGH DEMAND 
%     %subplot(1,3,j)
%     figure()
%     hold on 
%     for i=1:size(clus_final,1)
%         ind=clus_final(i,1);
%         gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%         hold on;
%     end
%     title(strcat(num2str(cases(j)*5+5),' \% MP nodes'))    %direct
%     %title(strcat(num2str(cases(j)*5),' \% MP nodes (gradual)'))    %incremental
% 
%     %pbaspect([1 1 1]);
%     axis off 
%     box on
%     for i=1:length(MPnodes)
%         plot(nodes(nodes(:,1)==indata.MP.nodeID(MPnodes(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(MPnodes(i)),3)...
%             ,'Marker','.','Color','#A2142F','Markersize',25)
%     end
%     saveas(gcf,strcat('mapWithNodes_dir',num2str(j),'.fig'));
%     print(strcat('mapWithNodes_dir',num2str(j),'.eps'),'-depsc','-painters')
% end
% % magenda
% %'#7E2F8E'
% % dark red
% % '#A2142F'
% figname = 'mapRegionsPCnodes';
% % figname = 'mapOriginDensitieshigh';
% % figname = 'mapDestDensitiesmedium';
% % figname = 'mapDestDensitieshigh';
% 
% %save and print figure: 
% 
% %saveas(gcf,strcat(figname,'.fig'));
% %saveas(gcf,strcat(figname,'.emf'));
% %print(figname,'-depsc','-painters')

