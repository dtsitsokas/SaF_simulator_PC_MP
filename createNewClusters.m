% New clusters for improving the PC 

% for the existing cluster of 3 regions: create subsets of most congested
% links 

%[1]. First run result analysis for the NC scenario 
% - Identify links with average density of over 50% 
% - create subset of congested links belonging to this cluster 

subset_reg = cell(1,3); 
non_subset_reg = cell(1,3); 
limit_reg = [0.43 0.4 0.25]; 
for i=1:indata.no_reg 
    links_reg = Avg_occ(indata.reInd{i}); 
    subset_reg{i} = indata.reInd{i}(links_reg > limit_reg(i));
    non_subset_reg{i} = setxor(indata.reInd{i},subset_reg{i});
end
subset_reg

%plot the different subsets on the map 
figure('Name','Map subset of congested links per region');
hold on;
subset_final = clus_final; 
for i = 1:3 %set the region of which I need the subset 
    subset_final(subset_reg{i}(subset_reg{i}<=1570)) = 7; 
end
for i=1:size(clus_final,1)
    ind=subset_final(i,1); 
    gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off

ff1 = figure('Name','MFDs of subsets per region');
f1 = subplot(1,3,1);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 1')

f2 = subplot(1,3,2);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 2')

f3 = subplot(1,3,3);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 3')

%make MFDs based only on the subsets 
for i=1:3
        Rsum1(i,:) = sum(outdata.x(subset_reg{i},1:indata.kmax));  %region accumulation only
        Rsum5(i,:) = sum((outdata.m(subset_reg{i},1:indata.kmax)*indata.v_ff)./(indata.Links2(subset_reg{i},3)*ones(1,indata.kmax)),1); %total moving flow in all network links
        accumulation = 0;
        outflow = 0;
        k_MFD = 0;
        for j=1:t:indata.kmax
            k_MFD = [k_MFD j*indata.DT]; %time for the MFD [hours]
            accumulation = [accumulation mean(Rsum1(i,j:j+(t-1)))]; %[veh] - mean accumulation over the time window
            outflow = [outflow mean(Rsum5(i,j:j+(t-1)))]; %[veh/h]
        end
        if i==1
            subplot(f1)
        elseif i==2
            subplot(f2)
        else
            subplot(f3)
        end
        plot(accumulation,outflow,'Linewidth',2)
end

f1 = subplot(1,3,1);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 1')

f2 = subplot(1,3,2);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 2')

f3 = subplot(1,3,3);
xlabel('n (veh)')
ylabel('q (veh/h)')
title('region 3')

%save creation of subsets in new mat file 
save('subsets_regions.mat','subset_reg','non_subset_reg')
    