% Create PC-suitable demand - Randomize demand for simulation tests 

%PC Suitable demand: 

%define ratios of intraregional demand (matrix 3x3) 

%separate centroids in groups per region 

%assign demand to pairs of regions x-y by using: 
% - mean value of demand*p(x,y) + random variance at every centroid 
% - verify later the ratios 
% - target: more demand in the city center from outside, less inside
% - all regional centroids (internal and external) treated in the same way 
% - change is made where? (be careful demand corrections happening in
% Initialization) 

%load current OD before all corrections
clear 
clc
OD = dlmread('OD.txt');
or_centroids = OD(2:end,1);
dest_centroids = OD(1,2:end); 

%separate centroids per region 
load('indata.mat')

load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids)
load clus_final %clustering for links

%classify origin centroids in regions
groupsOrCentroids = zeros(1,length(or_centroids)); 
for i=1:length(or_centroids)
   i1 = find(indata.Links2(:,4)==or_centroids(i),1); 
   groupsOrCentroids(i) =  indata.Links2(i1,end);
end

%classify destination centroids in regions 
groupsDestCentroids = zeros(1,length(dest_centroids)); 
for i=1:length(dest_centroids)
   i1 = find(indata.Links2(:,5)==dest_centroids(i),1); 
   if ~isempty(i1)
     groupsDestCentroids(i) =  indata.Links2(i1,end);
   end
end

% create new demand matrix 
OD_1 = OD;
a = OD(2:end,2:end);
tot_numb = sum(sum(a));
mean_orig = mean(sum(a,2)); %mean of total generating demand per centroid
mean_dest = mean(sum(a,1)); %mean of total arriving demand per centroid 

%set ratios i,j -> demand from i to j is x times a standard demand number 
rat = [1 10 1;2 10 2; 1 10 1];

baseDemand = (mean_orig/sum(sum(rat))); %766 veh/hour is the mean origin demand 
st_dev = baseDemand*0.25;  

%create deterministic demand
for i=1:indata.no_reg %origin region
    
    centr_or = find(groupsOrCentroids==i); %all origin centroids of the region
    
    for co=1:length(centr_or) %for every or centr of region i
        
        for cd = 1:length(dest_centroids) %for every dest centroid
            %define demand of every centroid of origin reg i to every centr of dest reg j
            switch groupsDestCentroids(cd) %select region of destination centroid
                case 1 %from region i to region (1)
                    OD_1(centr_or(co)+1,cd+1) = max(0,baseDemand/length(centr_or)*rat(i,1) + ...
                        (rand()-rand())*st_dev*rat(i,1));
                case 2
                    OD_1(centr_or(co)+1,cd+1) = max(0,baseDemand/length(centr_or)*rat(i,2) + ...
                        (rand()-rand())*st_dev*rat(i,2));
                case 3
                    OD_1(centr_or(co)+1,cd+1) = max(0,baseDemand/length(centr_or)*rat(i,3) + ...
                        (rand()-rand())*st_dev*rat(i,3));
            end
            
        end
    end
end

%Corrections for impossible trips 
cent = 69258; 
i = find(cent==OD_1(1,:));
OD_1(2:end,i)=0; %this centroid is problematic 

b = OD_1(2:end,2:end);
%global scaling
b = b*1.5; 
OD_1(2:end,2:end) = b; 

sum(sum(a))
sum(sum(b))

 


%plot the traffic pattern 
%old pattern
figure
surf(a)
%saveas(gcf,'old_demand_pattern.fig')
%new pattern 
figure 
surf(b) 
%saveas(gcf,'new_demand_pattern.fig')

%save the produced OD-Matrix (attention, it will replace the previous one)
%save('newOD.mat','OD_1')
%load('newOD.mat')
