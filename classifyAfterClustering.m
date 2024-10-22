function [outputArg1,outputArg2] = classifyAfterClustering(inputArg1,inputArg2)

% Clustering with regions only

%Read clustering Reza 

% junctions3Is = dlmread('junctions_Is.txt');
% clusters3Is = dlmread('clusters4_Is.txt');

% no_reg = length(unique(junctions3Is(:,2))); %number of regions (1,2,3) 
% % no_adjReg = 4; is currently set from parametersetting 

% region2 = cell(1,no_reg); 
% for i=1:no_reg
%     region2{i} = clusters3Is(clusters3Is(:,2)==i,1);
% end
% 
% Links2 = [Links2 zeros(size(Links2,1),1)];
% 
% for j=1:length(region2)
%     for i=1:length(region2{j})
%         Links2(Links2(:,1)==region2{j}(i),6) = j;
%     end
% end
% 
% 
% reInd = cell(1,no_reg); 
% for i=1:no_reg
%     reInd{i} =  find(Links2(:,6)== i);
% end
% 
% nodereg = cell(1,no_reg);
% for i=1:no_reg
%     nodereg{i} = junctions3Is(junctions3Is(:,2)==i,1);
% 
% end
% 

%% classify virtual links to clusters: origins -> where the destination
% is / destinations -> where the origin is

% for i=size(Links,1)+1:size(Links2,1)
%     ind = junct2.origin==Links2(i,1);
%     if sum(ind)==0
%         ind = junct2.destination == Links2(i,1);
%         j = junct2.or_index(ind);
%     else
%         j = junct2.dest_index(ind);
%     end
% 
% 
%     for r=1:no_reg
%         if ismember(j,reInd{r})
%             reInd{r}(find(reInd{r},1,'last')+1) = i;
%             Links2(i,6)=r;
%         end            
%     end
% 
% 
% end

% LinksP(:,6) = zeros(size(LinksP,1),1);
% LinksP(1:length(Links2(:,end)),6) = Links2(:,end); 
% toReturn = []; 
% for i=size(Links2,1)+1:size(LinksP,1)
%     if LinksP(i,1) > 10^7 
%         toReturn = [toReturn i];  
%     else 
%         toSearch = [upstrP(i,3) downstrP(i,3)];
%         toSearch = toSearch(toSearch>0);
%         toSearch = toSearch(toSearch<10^7);   
%         LinksP(i,6) = LinksP(LinksP(:,1)==toSearch,6); 
%     end
% end
% for i = toReturn 
%     toSearch = [upstrP(i,3) downstrP(i,3)];
%     toSearch = toSearch(toSearch>0);
%     LinksP(i,6) = LinksP(LinksP(:,1)==toSearch,6);
% end


%% Calculate maximum no of vehicles per region (n_jam) 

% max_n = zeros(1,no_reg);
% for i=1:no_reg
%     max_n(i) = sum(capacity(reInd{i})); %add capacities of all links the region
% end
% % max_ntw = sum(max_n);  %jam accumulation of the entire network 
% 
% clear reInd_1 reInd_2 reInd_3
% clear nodereg_1 nodereg_2 nodereg_3 region_1 region_2 region_3


%% clustering for nodes --- (centroids not included in the clustering)

%    %% Test if clustering is done correctly - [ok]
%     load LINKS.mat %careful - different structure wrt Links2 or Links
%     load nodes.mat %careful - different than Nodes
%     %plot network with colors
%     figure;
%     hold on;
%     for i=1:size(LINKS,1)
%         ind=Links2(i,6);
%         gplotdc1([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%         hold on;
%     end



%% update service capacity of exit links as a function of accumulation
%     % ---------------------------------------------------------------------
%     % (exit MFD) - per region 
% 
%     if indata.redCapacity == 1
%         for i=1:indata.no_reg 
% 
%             p = ExitServiceFunction(outdata.x(intersect(indata.group1,indata.reInd{i}),k),indata.capacity(intersect(indata.group1,indata.reInd{i})));
%             indata.junct2.lanesd(indG3r{i}) = p*init_ExitLinksLanesR{i}; %per region
%         end
%     end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end