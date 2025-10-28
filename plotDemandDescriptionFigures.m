% Script to generate figures for demand descrtiption

for demandCode = 1:2 %1 for med, 2 for high
    
    % load demand data (OD matrix)
    if demandCode == 1
        path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
        demandDesc = 'medium';
    elseif demandCode == 2
        path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
        demandDesc = 'high';
    end
    
    % path to save figures
    pathfig = 'F:\Max Pressure files\fifth set results\paper results figures\demadDescriptionFigures\';
    
    fname = strcat(path,'indata');
    load(fname)
    
    
    %***correction of the mistake in OD_links_2 (demand is wrong, but
    %simulations are correct because they run with demandGroup2)
    % OD_links_2 is corrected (need to rerun the cases to get the
    % corrections - for the input files of the paper we need this
    % correction: 
    
    if demandCode==1
       indata.OD_links_2(:,3) = indata.OD_links_2(:,3)*2; 
    else
       indata.OD_links_2(:,3) = indata.OD_links_2(:,3)*15;
    end
    
    
    demand = indata.OD_links_2; %this contains Origin_ID Dest_ID demand percntOfDemand
    origins = zeros(size(indata.group2));
    
    % Assign demand to entrance nodes
    for i=1:size(demand,1)
        ind = demand(i,1)==indata.group2;
        %nodes
        origins(ind) = origins(ind) + indata.defac(1000)*demand(i,3)*demand(i,4);
    end
    aa = [indata.group1; indata.group3];
    destinations = zeros(size(aa));
    
    %Assign demand to exit nodes (start node of exit Vlinks)
    for i=1:size(demand,1)
        while demand(i,2)>=2202
            demand(i,2) = find(indata.LinksP(:,1)==indata.upstrP(demand(i,2),3));
        end
        ind = demand(i,2)==aa;
        %nodes
        destinations(ind) = destinations(ind) + indata.defac(1000)*demand(i,3)*demand(i,4);
    end
    
    %% Create a map of origin densities (links with higher origin demand)
    
    load LINKS %links info [id, start_node, end_node, length, cluster (in 4 reg)]
    load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids)
    load clus_final %clustering for links
    load control %controlled intersections
    nodes = nodes3;
    
    %make simple map
    figure('Name','Origin density');
    hold on
    for i=1:size(clus_final,1)
        ind=clus_final(i,1);
        gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
        hold on;
    end
    
    %Print specified set of nodes with specific color (all same length)
    
    nodesVQs = indata.Links2(indata.group2,5); %all exit nodes of VQs (points of entry)
    myColorMap = jet(256);
    %myColorMap = autumn(256);
    colormap(myColorMap)
    for i=1:length(nodesVQs)
        plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),...
            'Marker','.','Color',myColorMap(max(1,ceil(origins(i)/max(origins)*256)),:),...
            'Markersize',30)%max(1,ceil(origins(i)/max(origins)*60)))
        %     plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i)...
        %           ,3),'Marker','.','Color','#A2142F','Markersize',30)
    end
    %title(strcat('origin densities ',demandDesc,' demand'))
    c = colorbar('Ticks',(0:0.1:1),...
        'TickLabels',round((0:max(origins)/10:max(origins))),...
        'TickLabelInterpreter' ,'latex')
    c.Label.String = 'demand (trips/h)';
    c.Label.Interpreter = 'latex';
    axis off
    %saveas(gcf,strcat(pathfig,'mapOriginDensities',demandDesc,'.fig'));
    %saveas(gcf,strcat(pathfig,'mapOriginDensities',demandDesc,'.emf'));
    %print(strcat(pathfig,'mapOriginDensities',demandDesc),'-depsc','-painters')
    
    %% Create a map for destinations density
    
    %make simple map
    figure('Name','Destination density');
    hold on
    for i=1:size(clus_final,1)
        ind=clus_final(i,1);
        gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
        hold on;
    end
    
    %Print specified set of nodes with specific color (all same length)
    nodesVQs = indata.Links2(aa(1:length(indata.group1)),5);
    nodesVQs = [nodesVQs; indata.Links2(aa(length(indata.group1)+1:end),4)];
    myColorMap = jet(256);
    %myColorMap = autumn(256);
    colormap(myColorMap)
    for i=1:length(nodesVQs)
        if destinations(i)>0
            plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),...
                'Marker','.','Color',myColorMap(max(1,ceil(destinations(i)/max(destinations)*256)),:),...
                'Markersize',30)%max(1,ceil(destinations(i)/max(destinations)*60)))
        end
    end
    %title(strcat('destination densities ',demandDesc,' demand'))
    c = colorbar('Ticks',(0:0.1:1),...
        'TickLabels',round((0:max(destinations)/10:max(destinations))),...
        'TickLabelInterpreter' ,'latex')
    c.Label.String = 'demand (trips/h)';
    c.Label.Interpreter = 'latex';
    axis off
    %saveas(gcf,strcat(pathfig,'mapDestDensities',demandDesc,'.fig'));
    %saveas(gcf,strcat(pathfig,'mapDestDensities',demandDesc,'.emf'));
    %print(strcat(pathfig,'mapDestDensities',demandDesc),'-depsc','-painters')
    
    
    %% figure for regional demand distribution
    
    ODRegDistr = zeros(3);
    %row = origin, column = destination
    
    for i=1:size(demand,1)
        r = indata.LinksP(demand(i,1),6); %origin
        c = indata.LinksP(demand(i,2),6); %destination
        ODRegDistr(r,c) = ODRegDistr(r,c) + indata.defac(1000)*demand(i,3)*demand(i,4);
    end
    
    figure('Name','Regional demand distribution')
    bar(ODRegDistr')
    legend('From 1', 'From 2', 'From 3')
    xlabel('To')
    ylabel('demand (veh/h)')
    saveas(gcf,strcat(pathfig,'RegDemandDistribution_',demandDesc,'.fig'));
    saveas(gcf,strcat(pathfig,'RegDemandDistribution_',demandDesc,'.emf'));
    print(strcat(pathfig,'RegDemandDistribution_',demandDesc),'-depsc','-painters')
    
    
    s = 0;
    g = 0;
    %Calculate trips
    for k=1:indata.kmax
        s = s + sum(indata.OD_links_2(:,3).*indata.OD_links_2(:,4))*indata.defac(k)*indata.DT;
        g = g + sum(indata.demandGroup2(:,1+(k>=901)))*indata.DT*indata.defac(k);
    end
    s
    g

end

