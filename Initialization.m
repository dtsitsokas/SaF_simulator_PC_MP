function []=Initialization()
%% Read all input files and create the FinalInput.mat to be used in a
%simulation

M2 = 4; %The maximum number of times an approach takes ROW in 1 cycle (I set 4 times) - NETWORK DEPENDENT
SP1 = 1; % [0=turnings.txt  1=turningsPic.txt  2:turning rates equal to all downstream links for the whole time]
centroids = dlmread('centroids.txt');
centroidsSeparation = dlmread('centroidsSeparation.txt');
downstr = dlmread('downstr.txt');
greens = dlmread('greens.txt'); %check if file is clean (Aimsun put garbage at the end)
nodes = dlmread('nodes.txt'); %check if file is clean (Aimsun put garbage at the end)- 692
OD = dlmread('OD.txt');
ODWU = dlmread('OD_WarmUp.txt'); %comment if doesn't exist
signals = dlmread('turningSignals.txt');
pythonNodes = dlmread('pythonNodes.txt');
pythonLinks = dlmread('pythonLinks.txt');
turnings = dlmread('turnings.txt');
turnings2 = dlmread('turningsPic.txt');
connectivity = dlmread('conectivity.txt');
rightLn = dlmread('rightLn.txt');

%parametersetting() %has been executed already
load('parameters.mat','minStageDur','vehlength','v_ff','init_clustering');

%% Exceptional adjustemnts for Barelona *****
% remove centroid 57008 (column 26 of OD, line 25 of centroids)

% Modify original OD to avoid certain local gridlocks 
[OD] = modifyOD(OD);

% ** To replace original OD matrix with new one (higher demand OD) uncomment: 
% load('newOD_high.mat','OD_1');
% OD = OD_1;
% clear OD_1

OD = [OD(:,1:25) OD(:,27:end)];
ODWU = [ODWU(:,1:25) ODWU(:,27:end)];
%OD(104,2:end) = OD(104,2:end)*0.5; %adjustment specific for Barcelona (to empty the network) 
centroids = [centroids(1:24,:); centroids(26:end,:)];
Links(402,5) = 69278; %correction for Bcn ntw
centroidsSeparation = [centroidsSeparation(1:24,:); centroidsSeparation(26:end,:)];
nodes = nodes(nodes(:,1)~=57008,:);

%also modifications at the end 

%% Create matrix: Links, Nodes
%load FirstInput
Links=pythonLinks;
Nodes=pythonNodes;
for i=1:length(centroids(:,1))
    index=find(centroids(i,2)==Links(:,1)); 
    Links(index, Links(index,:)==-1)=centroids(i,1);
    if ~isnan(centroids(i,3))
        index=find(centroids(i,3)==Links(:,1));
        Links(index,Links(index,:)==-1)=centroids(i,1);
    end
end
clear pythonLinks pythonNodes 

%% Plot the city (optional) ------------------------------------------------
% NLinks = size(Links,1);
% figure('Name','Map','units','normalized','outerposition',[0 0 1 1],'color','white')
% hold on;
% for i = 1:NLinks
%     S_node = Links(i,4);
%     E_node = Links(i,5);
%     if (S_node ~= -1 && E_node ~= -1)
%         Id_S = find (Nodes(:,1) == S_node);
%         Id_E = find (Nodes(:,1) == E_node);
%         gplotdc111([0 1;0 0],[Nodes([Id_S Id_E],2) Nodes([Id_S Id_E],3)],1)
%     end
% end
% pbaspect([1.5 1 1]); %NETWORK DEPENDENT
% axis off
% clear S_node E_node Id_S Id_E
% ------------------------------------------------------------------------

%% Create matrix : upstr (for every link the set of upstream links)

%Create File upstr with format: 
%Link ID - #ofUpstreamLinks - IDsofUpstreamLinks

NLinks = length(Links(:,1));
NNodes = length(Nodes(:,1));

upstr=downstr(:,1);
for i=1:NLinks
    ind=[];
    for j=3:size(downstr,2)
        a1 = find(downstr(:,j)==downstr(i,1));
        ind=[ind; a1];
    end
    upstr(i,2)=length(ind);
    l=1;
    while l<=upstr(i,2)
        upstr(i,l+2)=downstr(ind(l),1);
        l=l+1;
    end
end
clear a1 ind index


%% Construction of a structure to store junction data and turning ratios : junct

i=1; %row in nodes.txt
k=1; %row in the new table
step1 = length(nodes(:,1))/3; % -number of intersections included in "nodes"
step2 = length(turnings(:,1))/step1; % number of different turn ratio values (depends on total simulation time)
g_iter=1;
iter=g_iter;
if size(nodes,1)/3~=size(greens,1)
    disp('error in green/nodes dimensions')
    return
end
%extract pairs of links (approaches from nodes)
while i<= length(nodes(:,1))
    for j=1:nodes(i,3) %column in nodes.txt
        junct.node(k) = nodes(i,1); %intersection id of the approach 
        junct.node_index(k) = find(nodes(i,1)==Nodes(:,1));
        junct.origin(k) = nodes(i+1,j); %origin link id
        junct.destination(k) = nodes(i+2,j); %destination link id
        junct.splan(k) = nodes(i,2); %Control type of junction: <0: error, 0: fixed, 1: actuated, 2: external
        junct.lanesu(k) = connectivity(i+1,j); %no_of_lanes_up
        junct.lanesd(k) = connectivity(i+2,j); %no_of_lanes_down
        junct.rightu(k) = rightLn(i+1,j); %right_lane_included_in_lanes_up
        junct.rightd(k) = rightLn(i+2,j); %right_lane_included_in_lanes_down
        if junct.splan(k)<=0 % or -1 ? take the unsignalized intersections
            junct.goverc(k) = 1; %check for the T junctions/unsignalized
        else
            junct.goverc(k) = greens((i+2)/3,j); %g/C ratio for the approach (check approach 73054->3125)
        end
        for l=1:step2
            junct.turn(k,l) = turnings(iter,j); %no_of_vehicles_intending_to_turn/total_link_outflow (from snapshots every 20 sec) = averages every 900 sec (15') from AAPI of Aimsun
            iter = iter + step1;
        end
        k = k+1;
        iter = g_iter;
    end
    g_iter = g_iter + 1;
    iter = g_iter; %row in turnings.txt (go to next node)
    i = i+3;
end

%all signals are either uncotrolled or pre-timed
junct.splan(junct.splan>0) = 1; %controlled/pre-timed (any type, e.g. fixed, external etc)
junct.splan(junct.splan<0) = 0; %uncontrolled 

%% manually adjust signal plans (Barcelona) ***** 
junct.splan(2774) = 0;
junct.splan(129) = 0;

%for a small set of links that have downstream but not connected to a node
%(and not controlled) - Aimsun bugs (link-link without node) 
for i=1:NLinks
    check = sum(junct.origin == Links(i,1));
    if check == 0
        %link not listed as origin yet
        %check if it has downstream links
        if downstr(i,2)>0
            %disp('downstream links: ')
            %disp(downstr(i,2))
            if downstr(i,2)>1
                disp('error')
                return
            end
            %add pair to junct:
            junct.node(k) = -10; %no node - direct connection between links 
            junct.node_index(k) = -10; %no node - direct connection between links 
            junct.origin(k) = Links(i,1);
            junct.destination(k) = downstr(i,3);
            junct.splan(k) = 0; %uncontrolled
            junct.lanesu(k) = Links(i,2);
            junct.lanesd(k) = Links(Links(:,1)==downstr(i,3),2);
            junct.rightu(k) = 1;
            junct.rightd(k) = 1;
            junct.goverc(k) = 1;
            junct.turn(k,l:step2) = 1;
            k = k + 1;
        end
    end
end

% it should be ok now - (check)
% for i=1:NLinks
%     check = sum(junct.destination == Links(i,1));
%     if check==0
%        %link not listed as destination yet
%        %check if it has upstream links
%        if upstr(i,2)>0
%            disp('error upstream')
%        end
%     end
% end

stageDur = (-10)*ones(NNodes,11); %max no of phases per node = 11? [row = nodeID*same sequence as in nodes][minimum duration] 
stageVarDur = (-10)*ones(NNodes,11); %variable duration, if AR or <15 sec - variable = 0 
bigARnodes = []; 
bigARdur = [];

%% Read the file of traffic signals - pre-timed schedule
k = 1;
greentimes=zeros(length(junct.origin),M2*2); %one row per approach 
junct.offset=zeros(1,length(junct.origin)); %initialize
junct.cycle=zeros(1,length(junct.origin)); %initialize
junct.stageNo=zeros(length(junct.origin),5); %4 columns - max no of stages the same approach belongs 
junct.stageDur=zeros(length(junct.origin),5);
junct.stageVarDur=zeros(length(junct.origin),5);
while k <= size(signals,1) % k = first row of node block in signals
    if signals(k,2)>0 %node has declared signal settings in this file 
        %check if the file is correctly printed
        if sum(signals(k,1)==Nodes(:,1))==0
            disp('error in file turningSignals')
            disp('line k = ')
            disp(k)
            return
        end
        NofPhases=signals(k,2);
        offset=signals(k,3);
        cycle=sum(signals(k+1,:));
        cumtime=0; %cummulative time counter in the cycle to find time points to separate phases
        node_index = find(Nodes(:,1)==signals(k,1),1);
        for i=1:NofPhases %i=phase
            j=1; %j for the jth or-dest pair that has ROW in this phase i
            flag=0;
            stageDur(node_index,i) = signals(k+1,i);
            stageVarDur(node_index,i) = max(0,signals(k+1,i)-minStageDur); %if 0 - it does not change 
            while flag==0
                if signals(k+2*i,j)~=0 %There are approaches that DO have ROW at this phase (not ALL RED)
                    or_index=find(junct.origin==signals(k+2*i,j)); %index of link in junct.origin
                    dest_index=find(junct.destination==signals(k+2*i+1,j)); %index of link in junct.dest
                    index=intersect(or_index,dest_index); %index = index in the junct for this pair (approach)
                    %if isempty(index)
                    %    pause
                    %end
                    junct.cycle(index)=cycle;
                    junct.offset(index)=offset;
                    col = find(junct.stageNo(index,:),1,'last');
                    if isempty(col)
                        col = 1;
                    else
                        col = col + 1; 
                    end
                    junct.stageNo(index,col)=i; %no of stage (sequence) 
                    junct.stageDur(index,(find(junct.stageNo(index,:),1,'last'))+1)=signals(k+1,i); %initial duration of stage
                    junct.stageVarDur(index,(find(junct.stageNo(index,:),1,'last'))+1)=max(0,signals(k+1,i)-minStageDur); %variable duration=initial-minimum

                    % Find the green times in the cycle
                    % This pair of orig-dest in junct has green time in
                    % this phase (phase i) which lasts from
                    if i==1 %For the first phase only for each node
                        greentimes(index,1)=0;
                        greentimes(index,2)=signals(k+1,i);
                        
                    else
                        if sum(greentimes(index,:))==0 %if the row is empty
                            ind=0;
                        else
                            ind=find(greentimes(index,:)~=0,1,'last'); %find the last cell with a non-zero value to store in the next one
                        end
                        greentimes(index,ind+1)=cumtime;
                        greentimes(index,ind+2)=cumtime+signals(k+1,i); %lower-upper limits of green time
                        
                    end
                    j=j+1; %move to the next column (next pair of O-D)
                else
                    
                    if nodes(nodes(:,1)==signals(k,1),2)>0 %node is controlled
                        %mark nodes with large All-Red phases
                        if sum(signals(k+2*i,:)==0)==size(signals,2) && signals(k+1,i)>10
                            bigARnodes = [bigARnodes signals(k,1)];
                            bigARdur = [bigARdur signals(k+1,i)];
                            %disp('large AR')
                        end
                    end
                    %Case: No more approaches have ROW in this phase
                    %Update cummulative time
                    cumtime=cumtime+signals(k+1,i); %update time in cycle
                    flag=1; %continue to next phase
                end
            end
        end
        k=k+2*NofPhases+2; %Go to next node
    else
        k=k+1; %Go to next node
    end
end
%
clear cumtime dest_index flag NofPhases offsrt or_index
%fix mistakes of Aimsun --- (check all nodes of prob_approaches to see if 
%correction is required to all of them !!) 
prob_approaches = unique(junct.node(((junct.cycle==0)+(junct.splan==1))==2));
for i=1:length(prob_approaches)
    junct.splan(junct.node==prob_approaches(i))=0;
end
greentimes(greentimes==0)=1000;
greentimes(greentimes(:,1)==1000,1)=0; 
%% Fix the turning rates as they should be / normalize per link - add to junct

%Construct a matrix : turns
%Every row refers to a link as in Links(:,1) - Every column contains the no of turns
%for every destinatio ID from the same upstream link (of the row) (if there are
%downstream links) for every time interval of 15 mins (3 dimensions:
%1) Upstream link 2)Downstream links 3)time interval (1 to simulation_time/15mins(from APII))
%* problem with aimsun = some links are connected to each other without
%node between them - problem in reading junct origins-dest
turns=zeros(NLinks, length(downstr(1,:))-2, length(junct.turn(1,:)));
for i=1:NLinks
    a=find(Links(i,1)==junct.origin); %a is the index array in junct of all pairs link-downstream_link
    for j=1:downstr(i,2) %if it enters, the link has downstream links
        for k=1:length(junct.turn(1,:))
            turns(i,j,k)= junct.turn(a(j),k);
        end
    end
end

%Normalize turns
sum1 = sum(turns,2);
%turns2=turns; %Just in case I need the original values of turns
for i=1:length(turns(1,1,:))
    for j=1:length(turns(:,1,1))
        sum1=sum(turns(j,:,i));
        if sum1~=0
            for k=1:length(turns(1,:,1))
                turns(j,k,i)=turns(j,k,i)/sum1;
            end
        end
    end
end


for i=1:NLinks
    a=find(Links(i,1)==junct.origin); %a is the index array in junct of all pairs link-downstream_link
    for j=1:downstr(i,2) %if it enters, the link has downstream links
        for k=1:length(junct.turn(1,:))
            junct.turn(a(j),k)=turns(i,j,k);
        end
    end
end

clear step1 step2 k l j i  iter g_iter


%% Use the turningsPic file to define turning rates - based on the downstr file:
if SP1==1
    for i=1:NLinks
        if downstr(i,2)>0
            index1 = find(junct.origin==downstr(i,1));
            for j=1:downstr(i,2)
                index2 = find(junct.destination==downstr(i,2+j));
                a=intersect(index1,index2);
                for t=1:size(turnings2,1)/NLinks
                    junct.turn(a,t)=turnings2((t-1)*NLinks+i,j);
                end
            end 
        end
    end
    clear i index1 index2 a t j
elseif SP1==2
    % Make turning rates equal to all destinations from each upstream link
    % for the total simulation time
    for i=1:NLinks
        if downstr(i,2)>0
            index1 = find(junct.origin==downstr(i,1));
            for j=1:downstr(i,2)
                index2 = find(junct.destination==downstr(i,2+j));
                a=intersect(index1,index2);
                for t=1:size(turnings2,1)/NLinks
                    junct.turn(a,t)=1/downstr(i,2);
                end
            end
            
        end
    end
end

NPairs=length(junct.origin);

% Check data
%  a = Links(find(Links(:,2)>=2),1); %Links with 2 or more lanes
%  b = intersect(a,buses(:,4:end)); %Links with 2 or more lanes and at least 1 bus line traversing the link
%  c = Links(find(Links(:,2)>=3),1); %Links with 3 or more lanes
%  d = intersect(c,buses(:,4:end)); %Links with 3 or more lanes and at least 1 bus line traversing the link

clear turnings a1 bindex i ind index nodes %a b c d

% Sum the demand from origin centroids to any destination
OD(2:end,2:end) = OD(2:end,2:end)*60/105; %transform to [veh/hour]: it was veh/105'
ODWU(2:end,2:end) = ODWU(2:end,2:end)*60/15; %[veh/hour]: it was veh/15'
centroidDemand = OD(2:end,1); %=use the initial OD read from file
centroidDemand(:,2) = sum(OD(2:end,2:end),2);

centroidDemandWU = ODWU(2:end,1); %warm_up - comment if does not exist
centroidDemandWU(:,2)=sum(ODWU(2:end,2:end),2);
%load centroidDemand.mat %Here I load new centroid demand to match the mode shares for the mode choice model


%% Second grouping: links that are connected to centroids: origins + sinks

% Create matrix genAttract - specify connection of centroids to elements
% (nodes or links) 
genAttract = ones(size(centroids))*(-10);
genAttract(:,1) = centroids(:,1);
centroidsSeparation = [centroidsSeparation zeros(size(centroidsSeparation,1),1)];
%build generateTo and attractFrom, from centroidsSeparation
for i=1:size(centroidsSeparation,1)
    %find position of the separator and end_of_line [generatesTo 1
    %attractsFrom]
    sep = find(centroidsSeparation(i,:)==1,1);
    linend = find(centroidsSeparation(i,:)==0,1);
    %replace values before separator (generateTo)
    genAttract(i,2:sep-1)=1;
    %replace values after separator (attractFrom)
    genAttract(i,sep:linend-2)=-1;
end


% Create a matrix similar to genAttract but with the multipliers per
% centroid (divide flow equally to all connected elements) 
centroidMultiplier = zeros(size(genAttract));
centroidMultiplier(:,1) = centroids(:,1);
for i=1:size(genAttract,1)
    ind_gen = (genAttract(i,:)==1); %generate
    counter_o = sum(ind_gen);
    centroidMultiplier(i,ind_gen) = 1/counter_o;
    ind_att = (genAttract(i,:)==-1); %attract
    counter_o = sum(ind_att);
    centroidMultiplier(i,ind_att) = 1/counter_o;
end

%% Create a matrix similar to centroids with the indices of links in LINKS
%indexLinkToCentroid = replace all elements with the connected links 
%I have connections between centroids and links -> link indices

indexLinkToCentroid = [centroids(:,1) zeros(size(centroids))];
linksGenAttract = [centroids(:,1) zeros(size(centroids))];
linksCentroidMultiplier = [centroids(:,1) zeros(size(centroids))];
for i=1:size(centroids,1)
    for j=2:size(centroids,2)
        if centroids(i,j)>0
            %try to locate element in Links
            check = find(Links(:,1)==centroids(i,j)); %index of the links
            %find where the current last element is
            last_elem = find(indexLinkToCentroid(i,:)==0,1);
            %update dimensions if necessary
            if isempty(last_elem)  %case where row is full! Always need zero an the end of the row
                %increase number of columns in all matrices by 1
                indexLinkToCentroid = [indexLinkToCentroid zeros(size(indexLinkToCentroid,1),1)];
                linksCentroidMultiplier = [linksCentroidMultiplier zeros(size(linksCentroidMultiplier,1),1)];
                linksGenAttract = [linksGenAttract  zeros(size(linksGenAttract,1),1)];
                last_elem = find(indexLinkToCentroid(i,:)==0,1);
            end
            if isempty(check) %centroid connected to node
                if genAttract(i,j)==1 %generate
                    %check how many outgoing links start from this node
                    ind = Links(:,4)==centroids(i,j); %find links starting from this node/centroid
                    c1 = sum(ind); %no of links with the node as start
                elseif genAttract(i,j)==-1 %attract
                    %check how many incoming links end to this node
                    ind = Links(:,5)==centroids(i,j); %find links ending to this node/centroid
                    c1 = sum(ind); %no of links with the node as end
                else
                    disp('error - check')
                    break
                end
                
                %check dimensions again
                if last_elem+c1-1>size(indexLinkToCentroid,2)
                    %increase number of columns in all matrices
                    indexLinkToCentroid = [indexLinkToCentroid zeros(size(indexLinkToCentroid,1),c1-1)];
                    linksCentroidMultiplier = [linksCentroidMultiplier zeros(size(linksCentroidMultiplier,1),c1-1)];
                    linksGenAttract = [linksGenAttract  zeros(size(linksGenAttract,1),c1-1)];
                    last_elem = find(indexLinkToCentroid(i,:)==0,1);
                end
                if c1>=2
                    %                     if genAttract(i,j)==1 %generate
                    %                         %check if the node's outgoing links are destinations
                    %                         ign = intersect(group3,find(ind));
                    %                         ind = setxor(find(ind),ign);
                    %                     end
                    %                     if genAttract(i,j)==-1 %attract
                    %                         ign = intersect(group1,find(ind));
                    %                         ind = setxor(find(ind),ign);
                    %                     end
                    ind = find(ind);
                    indexLinkToCentroid(i,last_elem:last_elem+length(ind)-1) = Links(ind,1)';
                    linksCentroidMultiplier(i,last_elem:last_elem+length(ind)-1) = centroidMultiplier(i,j)/length(ind);
                    linksGenAttract(i,last_elem:last_elem+length(ind)-1) = genAttract(i,j); %copy the code for generation/attraction
                else
                    indexLinkToCentroid(i,last_elem) = Links(ind,1);
                    linksCentroidMultiplier(i,last_elem) = centroidMultiplier(i,j);
                    linksGenAttract(i,last_elem) = genAttract(i,j);
                end
            else
                %found in links - replace with the index of the link
                indexLinkToCentroid(i,last_elem) = Links(check,1);
                linksCentroidMultiplier(i,last_elem) = centroidMultiplier(i,j);
                linksGenAttract(i,last_elem) = genAttract(i,j);
            end
        end
    end
end

%transform to Link indices
for i=1:size(indexLinkToCentroid,1)
    for j=2:size(indexLinkToCentroid,2)
        if indexLinkToCentroid(i,j)>0
            indexLinkToCentroid(i,j) = find(Links(:,1)==indexLinkToCentroid(i,j));
        end
    end
end


%ANY link can be origin and destination at the same time - even while being an intermediate link (network-wise)

%% Check which intersections need separate queues for Max Pressure: when the same upstream link takes ROW in two phases (not identical) towards different destinations

% identify special intersections where Max Pressure will not be accurate
% with one queue per link 
specialIntersections = [];
specialLinks = [];
for i=1:NNodes
    %find node in signals
    indNode = find(signals(:,1)==Nodes(i,1),1);
    %check every incoming link
    noPhases = signals(indNode,2);
    inLinks = [];
    for j=1:noPhases
        inLinks = [inLinks signals(indNode+2*j,:)];
    end
    inLinks = unique(inLinks);
    inLinks = inLinks(inLinks~=0); %ids of inLinks
    for t=1:length(inLinks)
        %check if this link is connected to >1 downstream links
        noTimesAsOrigin = (junct.origin == inLinks(t));
        if sum(noTimesAsOrigin)>1
            %check if these two approaches take ROW simultaneously
            toCheck = greentimes(noTimesAsOrigin,:); %matrix - check if all rows are equal
            check = toCheck(1,:);
            k=1;
            while k+1<=size(toCheck,1)
                if (toCheck(k+1,:)==check)
                    k = k+1;
                else
                    %indicate as special case intersection
                    specialIntersections = [specialIntersections Nodes(i,1)];
                    specialLinks = [specialLinks inLinks(t)];
                    break
                end
            end
        end
    end    
end



%% origin centroids - indexes (rows) in centroids for origin centroids:
originCentroids = []; %specify which rows of "centroids" refer to origin centroids
for i=1:size(linksGenAttract,1)
    if sum(linksGenAttract(i,2:end)==1)>0
        originCentroids = [originCentroids; i];
    end
end


%% fixing problematic cases of intersections: 

for i=1:NPairs
    junct.or_index(i)=find(junct.origin(i)==Links(:,1));
    junct.dest_index(i)=find(junct.destination(i)==Links(:,1));
end

% !!Correct negative offsets 
ind = find(junct.offset<0);
junct.offset(ind) = junct.cycle(ind)+junct.offset(ind); 
% check ---------------
% ind = find(junct.offset<0,1);
% if ~isempty(ind)
%     disp('error offsets')
% end

%Identical phases listed as separate ones - fix offset accordingly 

%Fix no controlled intersections 

%% Create virtual links for ALL centroids (centroids always connected to virtual links
% and then v links connect to network)

connectLink2Centroid = zeros(size(Links,1),3); %matrix associating virtual links [1st col] to centroids [2nd col] (origin or dest)
greentimes2 = greentimes;
downstr2 = downstr;
Links2 = Links;
junct2 = junct;
kinit = size(Links,1);
jinit = length(junct.origin);

k = 0;
kd = 0;
ku = 0;
knodes = 80000; %numbers relating to virtual nodes
vlinksID = 90000; %numbers relating to virtual links
for i=1:size(indexLinkToCentroid,1) %row = centroid: one vlink as Origin, one vlink as Dest (highest level) 
    %scan every element that the centroid is connected to
    for j=2:size(indexLinkToCentroid,2) %all connections are to links
        if indexLinkToCentroid(i,j)>0 %nonZero element
            vlinksID = vlinksID + 1; %generate new virtual link connecting the centroid to the link
            ku = ku + 1;
            if linksGenAttract(i,j)==1 %centroid generates to link (they are all links) 
                k = k + 1;
                
                Links2(kinit + ku,1) = vlinksID; %virtual link ID > 90000
                Links2(kinit + ku,2) = Links(indexLinkToCentroid(i,j),2); %no of lanes of the downstream
                Links2(kinit + ku,3) = 0.01; %length
                Links2(kinit + ku,4) = centroids(i,1); %(start) the centroid
                Links2(kinit + ku,5) = Links(indexLinkToCentroid(i,j),4); % end = the start of the downstream link (the node)
                
                % update downstr
                downstr2(kinit + ku,1) = vlinksID;
                downstr2(kinit + ku,2) = 1;
                downstr2(kinit + ku,3) = Links(indexLinkToCentroid(i,j),1);
                
                %update virtual junct for connection: vlink -> link (origin
                %centroid)
                junct2.origin(jinit + ku) = vlinksID; %the new v.link
                junct2.destination(jinit + ku) = Links(indexLinkToCentroid(i,j),1); %the link
                junct2.splan(jinit + ku) = 0; %uncontrolled - everyone passes continuously
                junct2.lanesu(jinit + ku) = Links(indexLinkToCentroid(i,j),2); %no of lanes of downstream link of entry
                junct2.lanesd(jinit + ku) = Links(indexLinkToCentroid(i,j),2);
                junct2.rightu(jinit + ku) = 1; %all lanes - not important
                junct2.rightd(jinit + ku) = 1; %all lanes - not important
                junct2.goverc(jinit + ku) = 1;
                junct2.turn(jinit + ku,:) = ones(1,26); %everyone goes to the downstream link
                junct2.offset(jinit + ku) = 0; %not important since not controlled
                junct2.cycle(jinit + ku) = 0;
                junct2.or_index(jinit + ku) = kinit + ku;
                junct2.dest_index(jinit + ku) = indexLinkToCentroid(i,j);
                greentimes2(jinit + ku,:) = 0; 
                
                connectLink2Centroid(kinit + ku,1) = vlinksID;
                connectLink2Centroid(kinit + ku,2) = centroids(i,1);
                connectLink2Centroid(kinit + ku,3) = linksCentroidMultiplier(i,j)*linksGenAttract(i,j); %the multiplier of the node (origin)
                
                
                %this is important for demand generation -> this represents
                %the percentage of the flow that goes to this link from the
                %centroid - ENTRY rate from centroid - needs to be
                %estimated (now they are all equally divided) - from route
                %assignment this can change 
                
                %there are some centroids generating to multiple links
                %(similar to ones for nodes...) 
                
            elseif linksGenAttract(i,j)==-1 %centroid attracts from link
                %destination centroid in node - create virtual destination
                kd = kd + 1;
                %update downstr of the newly created link (include just the
                %link)
                downstr2(kinit + ku,1) = vlinksID;
                downstr2(kinit + ku,2:end) = 0; 
                
                
                %update downstr of the connected links
                row = find(downstr2(:,1)==Links(indexLinkToCentroid(i,j),1),1);
                column = find(downstr2(row,:)==0,1);
                if isempty(column)
                    downstr2 = [downstr2 zeros(size(downstr2,1),1)];
                    column = find(downstr2(row,:)==0,1);
                end
                if column==2
                    column = column + 1;
                end
                downstr2(row,column) = vlinksID;
                downstr2(row,2) = downstr2(row,2) + 1;
                
                Links2(kinit + ku,1) = vlinksID;
                Links2(kinit + ku,2) = Links(indexLinkToCentroid(i,j),2);
                Links2(kinit + ku,3) = 1000000; %length (to always receive outflow)
                Links2(kinit + ku,4) =  Links(indexLinkToCentroid(i,j),5); %the end of the upstream link
                Links2(kinit + ku,5) =  centroids(i,1); %the end (the centroid)
                
                
                %update virtual junct for connection: link -> vlink
                %(destination centroid)= sink
                junct2.origin(jinit + ku) =  Links(indexLinkToCentroid(i,j),1); %the link
                junct2.destination(jinit + ku) = vlinksID; %the new v.link
                copyTSindex = find(junct2.or_index==indexLinkToCentroid(i,j),1); 
                if isempty(copyTSindex)
                    copyTSindex = find(junct2.dest_index==indexLinkToCentroid(i,j),1);
                end
                junct2.splan(jinit + ku) = junct2.splan(copyTSindex); %same as its downstream or upstream
                junct2.lanesu(jinit + ku) = Links(indexLinkToCentroid(i,j),2); %no of lanes of downstream link of entry
                junct2.lanesd(jinit + ku) = Links(indexLinkToCentroid(i,j),2);
                junct2.rightu(jinit + ku) = 1; %all lanes - not important
                junct2.rightd(jinit + ku) = 1; %all lanes - not important
                junct2.goverc(jinit + ku) = 1;
                if downstr2(row,2) == 1
                    junct2.turn(jinit + ku,:) = 1; %only one downstream 
                else
                    %go to turningPics file and take the last recording
                    %every 15 mins 
                    row2 = indexLinkToCentroid(i,j); %index in downstr
                    for t=1:size(turnings2,1)/NLinks
                        junct2.turn(jinit + ku,t) = turnings2((t-1)*NLinks+row2,downstr(row2,2)+1);
                    end
%                     junct2.turn(jinit + ku,:) = 0; %node exits closed 
                end
                junct2.offset(jinit + ku) = junct2.offset(copyTSindex); %not important since not controlled
                junct2.cycle(jinit + ku) = junct2.cycle(copyTSindex);
                junct2.or_index(jinit + ku) = indexLinkToCentroid(i,j);
                junct2.dest_index(jinit + ku) = kinit + ku;
                greentimes2(jinit + ku,:) = greentimes2(indexLinkToCentroid(i,j),:);
                
                connectLink2Centroid(kinit + ku,1) = vlinksID;
                connectLink2Centroid(kinit + ku,2) = centroids(i,1);
                connectLink2Centroid(kinit + ku,3) = linksCentroidMultiplier(i,j)*linksGenAttract(i,j); %for attract this is not important / the turn rate is
                %this represents the percentage of the flow that the centroid may attract from this link (if all links were equivalent)
            end
        end
    end
end

NPairs2 = length(junct2.origin);
clear turnings2
%% Create matrix : upstr2 (for every link the set of upstream links) - CHECK THIS

%Link - Upstream Links : Create File upstr:
%Link ID - #ofUpstreamLinks - IDsofUpstreamLinks
NLinks2 = length(Links2(:,1));

upstr2=downstr2(:,1);
for i=1:NLinks2
    ind=[];
    for j=3:size(downstr2,2)
        a1 = find(downstr2(:,j)==downstr2(i,1));
        ind = [ind; a1];
    end
    upstr2(i,2)=length(ind);
    l=1;
    while l<=upstr2(i,2)
        upstr2(i,l+2) = downstr2(ind(l),1);
        l=l+1;
    end
end

%link 8770 has 2 destination centroids at its end: 57008 and 69278 -
%removed the one from the matrices

%% -----  Create new grouping of links:

%group1 = all intermediate links (all network real links - links with upstream and downstream links)

group1 = find(connectLink2Centroid(:,3)==0);

%group2 = all origin links (virtual) - all links starting from centroid
%(generating demand)

group2 = find(connectLink2Centroid(:,3)>0);

%group3 = all destination links (virtual) - all links ending at centroid
%(attracting demand)

group3 = find(connectLink2Centroid(:,3)<0);

%% Transform demand to be easily used for demand generation
% connectLink2CentroidIndex = connectLink2Centroid;
% % [vlinkID index_of_centroid_in_centroidDemand]
% for i=size(Links,1)+1:size(connectLink2CentroidIndex,1)
%     ind = (connectLink2CentroidIndex(i,2)==centroidDemand(:,1));
%     if sum(ind) > 0  %if the centroid is an origin
%         connectLink2CentroidIndex(i,2) = find(ind);
%     else
%         connectLink2CentroidIndex(i,2) = 0; %centroid is not listed in centroidDemand
%     end
% end

%% Prepare a matrix for shortest path checks
LinksP = Links2;
downstrP = downstr2;
LinksP(LinksP(:,3)==LinksP(1576,3),3)=0.01; %remove large distance of the destination links
centroidsP = centroids;
%add additional virtual links in LinksP and downstreamP for every centroid in nodes
%add one link connected to directly to centroid, then several links going
%to nodes (known turn rates), then to every link (unknown turn rates) 
kinit = 100000;

sinit = size(LinksP,1);
for i=1:size(centroids,1) %for every centroid = row
    if sum(centroids(i,2:end)>0)>1 || sum(centroids(i,2)==Links(:,1))==0  %centroid connected to multiple objects OR to at least one node  -> require highest vlink
        i_gen = genAttract(i,:)==1; %ind for origins
        i_att = genAttract(i,:)==-1; %ind for destinations 
        if sum(i_gen)>0 %centroid is origin for at least one object  
            %add  high level origin vlink 
            Ovlink = centroids(i,1)*1000; %the origin vlink id = centroid*1000  
            sinit = sinit + 1; %the index for LinksP matrix  
        
            downstrP(sinit,1) = Ovlink; %will add all its downstreams below
            downstrP(sinit,2:end) = 0; 
            
            %add to LinksP
            LinksP(sinit,1) = Ovlink;
            LinksP(sinit,2) = 1; %is not used in SaF, only in path detection
            LinksP(sinit,3) = 0.01; %is not used in SaF, only in path detection
            LinksP(sinit,4) = centroids(i,1); %start: the centroid
            LinksP(sinit,5) = centroids(i,1)*100; %end: some virtual node: centroid*100
                        
            % associate centroid with vlinkid in centroidsP
            centroidsP(i,1) = Ovlink; %first column of centroidsP = the highest origin link 
            
        end
        if sum(i_att)>0 %centroid is destination for at least one object 
            %add high level dest vlink 
            Dvlink = centroids(i,1)*1000 + 1; %the dest vlink id = centroid*1000+1  
            sinit = sinit + 1; 
            
            downstrP(sinit,1) = Dvlink;  
            downstrP(sinit,2:end) = 0; %no downsrteam
            
            %will be added as downstr to all upstream links below 
            
            %add to LinksP
            LinksP(sinit,1) = Dvlink;
            LinksP(sinit,2) = 1; %is not used in SaF, only in path detection
            LinksP(sinit,3) = 0.01; %is not used in SaF, only in path detection
            LinksP(sinit,5) = centroids(i,1); %end: the centroid
            LinksP(sinit,4) = centroids(i,1)*100+1; %start: some virtual node: centroid*10+1 
            
            % associate centroid with vlinkid in centroidsP
            centroidsP(i,11) = Dvlink; %11th column of centroidsP = the highest destination link 
            
        end
    else  %if centroid is connected only to one link = no high vlink : only one vlink (the one assigned before) 
        if genAttract(i,2)==1 %origin 
            centroidsP(i,1) = connectLink2Centroid(connectLink2Centroid(:,2)==centroids(i,1),1);
        elseif genAttract(i,2)==-1 %dest 
            centroidsP(i,11) = connectLink2Centroid(connectLink2Centroid(:,2)==centroids(i,1),1); %11th column of centroidsP = the highest destination link 
        end
    end
    for j=2:size(centroids,2) %for every link/node connected to centroid of row i 
        if centroids(i,j)>0 %if there is a connection to an element 
            if isempty(find(centroids(i,j)==Links(:,1),1)) %if not an original network link 
                %connected to network node --   
                if genAttract(i,j)==1
                    % origin centroid in node
                    
                    % create new node-vlink connecting to the higher level vlink of the centroid:
                    kinit = kinit + 1; 
                    sinit = sinit + 1;
                    
                   %add new vlink (of the node) as downstr of the high Ovlink of the centroid 
                    row = find(downstrP(:,1)==centroids(i,1)*1000,1); %row of the high Ovlink 
                    downstrP(row,2) = downstrP(row,2) + 1; %increase downstream links by 1 
                    column = find(downstrP(row,:)==0,1);
                    if isempty(column)
                        downstrP = [downstrP zeros(size(downstrP,1),1)];
                        column = find(downstrP(row,:)==0,1);
                    end
                    if column==2
                        column = column + 1;
                    end
                    downstrP(row,column) = kinit; %the node vlink is downstr of the high vlink  
                    
                    %add new vlink to downstrP
                    downstrP(sinit,1) = kinit; 
                    % find all prev vlinks that start from this centroid and end to
                    % this node
                    ind1 = (LinksP(:,4)==centroids(i,1)); %start from the centroid
                    ind2 = (LinksP(:,5)==centroids(i,j)); %end to the node
                    ind = find(ind1+ind2 == 2); %indices of rows in Links2
                    
                    %add the prev defined vlinks as downstr of the node
                    %vlink
                    downstrP(sinit,2) = length(ind);
                    for k = 1:length(ind)
                        downstrP(sinit,k+2) = LinksP(ind(k),1);
                    end
                    
                    % add new link to LinksP
                    LinksP(sinit,1) = kinit;
                    LinksP(sinit,2) = 1; %is not used in SaF, only in path detection
                    LinksP(sinit,3) = 0.01; %is not used in SaF, only in path detection
                    LinksP(sinit,4) = LinksP(row,5); %start = the end of the high vlink 
                    LinksP(sinit,5) = centroids(i,1)*10; %end = vnode 
                    
                    %Update start of vlinks connected to network links 
                    LinksP(ind,4) = LinksP(sinit,5); %start of links = end of node vlink 
                    
                    % replace node with vlinkid in centroidsP
                    centroidsP(i,j) = kinit; %row = refers to centroid, 1st col = high Ovlink, last = high Dvlink, interm. = node vlink
                    
                elseif genAttract(i,j)==-1
                    % destination centroid in node -- I don't need node
                    % vlinks here - just connect to the high Dvlink 
%                     kinit = kinit + 1;
%                     sinit = sinit + 1;
                    
                    % add Dvlink as downstr to all previously defined vlinks that start from the node : 
%                     downstrP(sinit,1) = kinit;
%                     downstrP(sinit,2) = 1;
%                     downstrP(sinit,3) = Dvlink;
                    
                    % add Dvlink as downstr of all links ending at the centroid
                    % find all links that end to this centroid and start from
                    % the node
                    ind1 = (LinksP(:,4)==centroids(i,j));
                    ind2 = (LinksP(:,5)==centroids(i,1));
                    ind = (ind1+ind2 == 2); %indices of rows in Links2
                    
                    
%                     for k = 1:length(ind)
%                         downstrP(ind(k),2) = downstrP(ind(k),2) + 1;
%                         column = find(downstrP(ind(k),:)==0,1);
%                         if isempty(column)
%                             downstrP = [downstrP zeros(size(downstrP,1),1)];
%                             column = find(downstrP(ind(k),:)==0,1);
%                         end
% %                         if column==2
% %                             column = column + 1;
% %                         end
%                         downstrP(ind(k),column) = Dvlink;
%                     end
                   
                    
                    downstrP(ind,2) = 1; 
                    downstrP(ind,3) = centroids(i,1)*1000 + 1; 
%                     
                    
                    %replace end node of these links = start of Dvlink 
                    LinksP(ind,5) = LinksP(LinksP(:,1)==centroids(i,1)*1000 + 1,4);
                    
                    % create new link:
                    % start = centroid, end=centroid
%                     LinksP(sinit,1) = kinit;
%                     LinksP(sinit,2) = 0; %is not used in SaF, only in path detection
%                     LinksP(sinit,3) = 0.01; %is not used in SaF, only in path detection
%                     LinksP(sinit,5) = centroids(i,1); %end = the centroid
%                     LinksP(sinit,4) = LinksP(LinksP(:,1)==Dvlink,5); %start = the end of the Dvlink node
                    
                    % replace node with vlinkid in centroidsP
                    centroidsP(i,j) = centroids(i,1)*1000 + 1; %not important as there is no node vlink for attraction centroids
                end
            else
                %Centroid connected to original network links -> vlinks are
                %already defined connecting links to the centroid (origin
                %node = centroid id) - I need to add those links as
                %downstream or upstream of the high vlinks 
                                
                if genAttract(i,j)==1
                    %for origins 
                    row = find(downstrP(:,1)==centroids(i,1)*1000,1); %row of the high Ovlink
                    downstrP(row,2) = downstrP(row,2) + 1; %increase downstream links by 1
                    column = find(downstrP(row,:)==0,1);
%                     if isempty(column)
%                         downstrP = [downstrP zeros(size(downstrP,1),1)];
%                         column = find(downstrP(row,:)==0,1);
%                     end
                    [a,~] = find(downstrP==centroids(i,j));
                    vlink = downstrP(a,1);
                    vlink = vlink(vlink>90000); 
                    downstrP(row,column:column-1+length(vlink)) = vlink; %the existing vlinks downstream of the high vlink
                    
                else
                    %for destinations = add the Dvlink as downstream of the
                    %existing vlinks (?) - not necessary 
                    row = find(downstrP(:,1)==centroids(i,j));
                    downstrP(row,2) = downstrP(row,2) + 1; %increase downstream links by 1
                    column = find(downstrP(row,:)==0,1); 
                    if isempty(column)
                        column = size(downstrP,2)+1;
                    end
                    downstrP(row, column) = centroids(i,1)*1000+1; 
                end
                    
                
                
            end
        end
    end
end
%% Create matrix : upstrP (for every link the set of upstream links) - CHECK THIS

%Link - Upstream Links : Create File upstr:
%Link ID - #ofUpstreamLinks - IDsofUpstreamLinks
NLinksP = length(LinksP(:,1));

upstrP=downstrP(:,1);
for i=1:NLinksP
    ind=[];
    for j=3:size(downstrP,2)
        a1 = find(downstrP(:,j)==downstrP(i,1));
        ind = [ind; a1];
    end
    upstrP(i,2)=length(ind);
    l=1;
    while l<=upstrP(i,2)
        upstrP(i,l+2) = downstrP(ind(l),1);
        l=l+1;
    end
end


%% 
OD_links = []; 
%make the OD pairs (centroids -> links)
for i=2:size(OD,1)
    orig = OD(i,1);
    if (sum(centroids(centroids(:,1)==orig,2:end)>0)>1) || (~ismember(centroids(centroids(:,1)==orig,2),Links(:,1)))  
        %multiple elements or first one node: Use the highest vlink 
%         O_links = centroidsP(centroids(:,1)==orig,1); %the highest level vlink
        
        %divide flow in distinct elements with flow multipliers
        ind_columns = [genAttract(centroids(:,1)==orig,:)==1 1==0];
        O_links = centroidsP(centroids(:,1)==orig,ind_columns); %take all links "generate"
        O_links_Mult = centroidMultiplier(centroids(:,1)==orig,ind_columns(1:end-1));
    else
        %Link: only one element - start directly from this link
        O_links = centroidsP(centroids(:,1)==orig,2);
        %start from the upstream of this link 
        O_links = upstrP(upstrP(:,1)==O_links,find(upstrP(upstrP(:,1)==O_links,:),1,'last'));
        O_links_Mult = 1; 
    end
    for j=2:size(OD,2)
        dest = OD(1,j);
        %create all possible pairs
        if sum(centroids(centroids(:,1)==dest,2:end)>0)>1 %multiple elements
            D_links = centroidsP(centroids(:,1)==dest,11); %highest dest vlink
        else
            D_links = centroidsP(centroids(:,1)==dest,2); %actual link/node? 
        end
        
        %O_mult = 1;
        %%or_centr %dest_centr %or_link %dest_link %or_link_index %dest_link_index %demand_multiplier %demand_WU %demand
        %OD_links = [OD_links; orig dest O_links D_links find(O_links==LinksP(:,1),1) find(D_links==LinksP(:,1),1) O_mult ODWU(i,j) OD(i,j)];
        
        % Alternative OD_links construction: flow from centroids is divided to the different elements by specific percentages 
        % When centroid is connected to multiple "generate" elements: OD start from
        % the elements (not the high level vlink - "attract" elements are all going
        % to the high level vlink of destination 
        for k=1:length(O_links) 
            %[origin_centroid, dest_centroid, Origin_link, destination_link,
            %origin_link_index, dest_link index, demand_split,
            %demand_warm_up, demand_regular]
            OD_links = [OD_links; orig dest O_links(k) D_links find(O_links(k)==LinksP(:,1),1)...
                find(D_links==LinksP(:,1),1) O_links_Mult(k) ODWU(i,j) OD(i,j)];   
        end
        
    end
end

%remove zeros
OD_links = OD_links(OD_links(:,9)~=0,:);

%remove impossible trips (Barcelona) 
% OD_links = OD_links([(1:1079) (1081:2637)],:);
% OD_links(1080,7) = 1; 

%% Prepare shortest path calculation - one time part (graph is made here) 

% Construct connectivity matrix A with distances (will be replaced with times) 
indA = LinksP(:,1);
A = zeros(size(LinksP,1)); %square matrix with all links (real+virtual)
for i=1:size(downstrP,1) %all links
    % ind1 = find(indA == downstr(i,1)); %this is i, the downstr + LINKS are
    % in the same order
    for j=1:downstrP(i,2) %all downstream links of the row
        ind2 = find(indA == downstrP(i,j+2)); %find the index of the downstream link
        A(i,ind2) = LinksP(i,3); %put the length [m] of the upstream link as dist
    end
end
A = A/1000; %[km]
A = sparse(A); 


%prepare initial bins for connection counters (empty) 
empty_connCounters = downstrP;
empty_connCounters(empty_connCounters==0) = -10; %cell not used
b = empty_connCounters(:,3:end);
b(b>0) = 0;
empty_connCounters(:,3:end) = b;
clear b
    
%% Calculate initial turn rates for all approaches based on distance shortest
%path for the first 15 mins of simulation: 

load('parameters','defac','t_win','DT')
defac_avg = mean(defac(1:ceil(t_win/DT))); %to fix = mean of 15 mins 
speeds = zeros(size(LinksP,1),1) + v_ff/1000; %(km/h)
A1 = A./speeds; %initial speeds 
%t_win = 15/60; %h %time window for the calculation 
t_ind = 1; %updates turn rates only on the first column 
%call function:                                     
[connCounters,stack_OD,newjunct2] = updateTurnAndExitRates(OD_links,defac_avg,junct2,t_win,LinksP, A1, downstrP, t_ind, empty_connCounters);
clear A1
%connCounters - correct if it is zero in all connections 

load('parameters','updateTR')
if updateTR~=0 
    junct2 = newjunct2; 
    %assign the calculated turn rates to all time intervals 
    % if dynamic update of turn rates is on: they will be replaced 
    for i=2:size(junct2.turn,2)
        junct2.turn(:,i) = junct2.turn(:,1);
    end
end
clear newjunct2;

%% Create connection between origins and virtual links of group2 for the
% %simulator: I need for every vlink in group2 -> total inflow (depends on
% %turning rates for node links) 
% 
[demandGroup2, OD_links_2] = createDemandGroup2(OD_links,group2,downstrP,upstrP,connCounters);

%correction for vlink 1572 going directly to exit link 2211 

ind = find(OD_links_2(:,1)==1572); 
OD_links_2(ind,1)=1571; 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

%correction for vlink 1606 

ind = find(OD_links_2(:,1)==1606); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1667); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1668); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1691); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1724); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1881); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1902); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1905); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

%Some destination of the follwing origins are non-reachable - remove entirely?  
ind = find(OD_links_2(:,1)==2060); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==2039); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==2180); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==2040); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

ind = find(OD_links_2(:,1)==1946); 
%remove demand from OD_links_2
OD_links_2 = OD_links_2(setxor(ind,1:size(OD_links_2,1)),:);

%% 
%indices of related approaches (origins/destinations) per link 
% for i=1:NLinks2
%     ind_or(i,:) = (junct2.or_index==i);
%     ind_or_pos(i) = sum(ind_or(i,:))>0;
%     ind_dest(i,:) = (junct2.dest_index==i);
%     ind_dest_pos(i) = sum(ind_or(i,:))>0;
% end

%% Modify demand in specific origin links to facilitate network emptying at
% the end of simulation time (6 hours preferrably) ***** 

%demandGroup2([147 214 285 291],:) = demandGroup2([147 214 285 291],:)*0.7;  

%% Correct turn rates that are read ---------------- ***

if updateTR==0 
    junct = correctTurnRates(junct,downstr);
    junct2 = correctTurnRates(junct2,downstr2);
end
%%
% check if all approaches have turn rates == 1
problematic_links = []; 
for i=1:NLinks2
    ind = (junct2.or_index == i);
    check = junct2.turn(ind,:); 
    check2 = sum(check,1);
    if sum(check2-1<abs(0.001))<length(check2)
        problematic_links = [problematic_links i];
    end
end
if ~isempty(problematic_links)
    % fix problematic links: links that end in nodes with multiple centroids 
    fix2 = [problematic_links(1:2) problematic_links(4:end)];
    fix3 = problematic_links(3);
    for i=1:length(fix2)
       ind = find(junct2.or_index==fix2(i));
       rates_ = junct2.turn(ind,:); 
       rates_ = rates_(end,:); 
       rates_ = rates_/2; 
       junct2.turn(ind(end-1),:) = rates_; 
       junct2.turn(ind(end),:) = rates_; 
    end
    for i=1:length(fix3)
       ind = find(junct2.or_index==fix3(i));
       rates_ = junct2.turn(ind,:); 
       rates_ = rates_(end,:); 
       rates_ = rates_/3; 
       junct2.turn(ind(end-1),:) = rates_; 
       junct2.turn(ind(end-2),:) = rates_;
       junct2.turn(ind(end),:) = rates_; 
    end
end
%% *** correction because of errors of Aimsun model for Barcelona: the above sections are
%listed uncontrolled but they have a traffic signal plan: 
fix_ = [1447 1450];
for i=1:length(fix_)
   ind = junct.or_index == fix_(i); 
   junct.splan(ind)=1; 
end

%% Create the structure that will be used for Max Pressure application: 
% For every intersection we need to know: 
% - # and duration of stages 
% % - all approaches and  take green in every stage 
MP.stageDur = {}; %initial (reference) stage duration 
MP.stageVarDur = {}; %variable part 
MP.stageBaseDur = {}; %constant part
MP.stagesInvolved = {}; %eligible stages (all small stages ignored)
MP.approaches = {}; %all intersection approaches
MP.stageApproaches = {};

%sn = size(junct.stageNo,2); 
for i=1:size(Nodes,1)
    MP.nodeID(i) = Nodes(i,1); %I can directly indicate indices in MP rather than nodeIDs 
    
    MP.noOfStages(i) = sum(stageDur(i,:)>0);
    MP.stageDur{i} = stageDur(i,stageDur(i,:)>0); 
    MP.stageVarDur{i} = stageVarDur(i,stageVarDur(i,:)>=0);
    MP.stageBaseDur{i} = MP.stageDur{i}-MP.stageVarDur{i}; %this always applies as minimum 
    MP.stagesInvolved{i} = find(MP.stageVarDur{i}>0); %# phases available for MP control 
    MP.timeToAllocate(i) = sum(MP.stageDur{i}(MP.stagesInvolved{i})); %the total green of the involved phases to assign (including the minimum) 
    MP.approaches{i} = find(junct.node_index==i); %indices of all approaches in junct that belong to this node 
    
    if isempty(MP.approaches{i})
        MP.approaches{i} = 0; 
        MP.splan(i) = -10;
        MP.cycle(i) = -10; 
        MP.offset(i) = -10;
    else
        MP.splan(i) = mode(junct.splan(MP.approaches{i}));
        MP.cycle(i) = mode(junct.cycle(MP.approaches{i}));
        MP.offset(i) = mode(junct.offset(MP.approaches{i}));
    end
    
    %stages involved: for every stage, associate the indices of the approaches involved 
    if MP.noOfStages(i)==0 || sum(MP.approaches{i})==0 
         MP.stageApproaches{i,1} = 0; 
    else
        for j=1:MP.noOfStages(i) 
            [ind,~] = find(junct.stageNo(MP.approaches{i},:) == j); %i=intersection, j=stage, MP.stageApproaches = [array of all approaches of the stage]
            if isempty(ind)
                MP.stageApproaches{i,j} = 0;
            else
                MP.stageApproaches{i,j} = MP.approaches{i}(ind);
            end
        end
    end
end


%fixing some signals 
ind = junct2.cycle==0; 
junct2.splan(ind)=0; 

% Capacities of links
capacity = Links2(:,2).*Links2(:,3)./vehlength;

%% Clustering with regions only

% MFDs - per region - clustering of ISTTT paper (Reza?) 
if init_clustering
    junctions3Is = dlmread('junctions_Is.txt');
    clusters3Is = dlmread('clusters4_Is.txt');
    no_reg = length(unique(junctions3Is(:,2))); %number of regions (1,2,3)
    % no_adjReg = 4; is currently set from parametersetting
    region2 = cell(1,no_reg);
    for i=1:no_reg
        region2{i} = clusters3Is(clusters3Is(:,2)==i,1);
    end

    Links2 = [Links2 zeros(size(Links2,1),1)];

    for j=1:length(region2)
        for i=1:length(region2{j})
            Links2(Links2(:,1)==region2{j}(i),6) = j;
        end
    end


    reInd = cell(1,no_reg);
    for i=1:no_reg
        reInd{i} =  find(Links2(:,6)== i);
    end

    nodereg = cell(1,no_reg);
    for i=1:no_reg
        nodereg{i} = junctions3Is(junctions3Is(:,2)==i,1);

    end


    %classify virtual links to clusters: origins -> where the destination
    %is / destinations -> where the origin is

    for i=size(Links,1)+1:size(Links2,1)
        ind = junct2.origin==Links2(i,1);
        if sum(ind)==0
            ind = junct2.destination == Links2(i,1);
            j = junct2.or_index(ind);
        else
            j = junct2.dest_index(ind);
        end


        for r=1:no_reg
            if ismember(j,reInd{r})
                reInd{r}(find(reInd{r},1,'last')+1) = i;
                Links2(i,6)=r;
            end
        end


    end

    LinksP(:,6) = zeros(size(LinksP,1),1);
    LinksP(1:length(Links2(:,end)),6) = Links2(:,end);
    toReturn = [];
    for i=size(Links2,1)+1:size(LinksP,1)
        if LinksP(i,1) > 10^7
            toReturn = [toReturn i];
        else
            toSearch = [upstrP(i,3) downstrP(i,3)];
            toSearch = toSearch(toSearch>0);
            toSearch = toSearch(toSearch<10^7);
            LinksP(i,6) = LinksP(LinksP(:,1)==toSearch,6);
        end
    end
    for i = toReturn
        toSearch = [upstrP(i,3) downstrP(i,3)];
        toSearch = toSearch(toSearch>0);
        LinksP(i,6) = LinksP(LinksP(:,1)==toSearch,6);
    end


    %Calculate maximum no of vehicles per region (n_jam)
    max_n = zeros(1,no_reg);
    for i=1:no_reg
        max_n(i) = sum(capacity(reInd{i})); %add capacities of all links the region
    end
    % max_ntw = sum(max_n);  %jam accumulation of the entire network

    clear reInd_1 reInd_2 reInd_3
    clear nodereg_1 nodereg_2 nodereg_3 region_1 region_2 region_3


    % clustering for nodes --- (centroids not included in the clutering)

    %% Test if clustering is done correctly - [ok]
    load LINKS.mat %careful - different structure wrt Links2 or Links
    load nodes.mat %careful - different than Nodes
    %plot network with colors
    figure;
    hold on;
    for i=1:size(LINKS,1)
        ind=Links2(i,6);
        gplotdc1([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
        hold on;
    end


    indG3r = cell(1,no_reg);
    init_ExitLinksLanesR = cell(1,no_reg);
    for jk=1:no_reg
        indG3r{jk} = [];
        set1 = intersect(reInd{jk},group3);
        for j = 1:length(set1)
            indG3r{jk} = [indG3r{jk} find(set1(j) == junct2.dest_index)];
        end
        init_ExitLinksLanesR{jk} = junct2.lanesd(indG3r{jk});

    end



end

%% Figure of generating demand based on OD matrix 

% ODRegDistr = zeros(3); 
% %row = origin, column = destination 
% 
% for i=1:size(OD_links,1)
%     r = LinksP(OD_links(i,5),6); %origin
%     c = LinksP(OD_links(i,6),6); %destination
%     ODRegDistr(r,c) = ODRegDistr(r,c) + OD_links(i,7)*OD_links(i,9);
% end
% 
% figure
% bar(ODRegDistr')
% legend('From 1', 'From 2', 'From 3')
% xlabel('To')
% ylabel('demand (veh/h)')
% saveas(gcf,'InterregionalDemandDistribution.fig')
% saveas(gcf,'InterregionalDemandDistribution.emf')
% close all 

clear t_win a c t_ind toSearch toReturn vlink a1 centroidsSeparation column cycle eli extralinks ind ind1 ind2 ind_r kd ku l len1 row sum1 vLinksID
clear c1 check counter i ind_att ind_gen ind indNode inLinks j k last_elem linend noPhases noTimesAsOrigin offset sep t toCheck
clear O_links O_links_Mult orig Ovlink sinit vlinksID counter_o D_links dest Dvlink jinit kinit knodes 
clear check2 col copyTSindex defac_avg fix2 fix3 fix_ i_att i_gen ind_columns indA node_index rates_ row2 
save('FinalInput');