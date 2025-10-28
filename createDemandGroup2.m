function [demandGroup2,OD_links_2] = createDemandGroup2(OD_links,group2,downstrP,upstrP,connCounters)
%% Create connection between origins and virtual links of group2 for the
%simulator: I need for every vlink in group2 -> total inflow (depends on
%turning rates for node links) 
demandGroup2 = zeros(length(group2),2); %1st column = indx of link, 2nd = WU, 3rd = regular demand 
%this process must be performed for every new turn rates!! (to divide flow
%in node vlinks) 
OD_links_2 = []; %zeros(size(OD_links,1),4); %origin_link_ind(VQ of group2) - dest_link_ind - demand - multiplier 
k=1;
for i=1:size(OD_links,1) 
    row = OD_links(i,5); %index in LinksP
    if OD_links(i,3)>100000 %node vlink - need to generate the demand to its downstr (this doesnt exist)
        for j = 1:downstrP(row,2)
            ind = find(downstrP(:,1) == downstrP(row,j+2),1); 
            ind2 = find(group2==ind,1);
            demandGroup2(ind2,1) = demandGroup2(ind2,1) + OD_links(i,7)*OD_links(i,8)*connCounters(row,j+2); %warm-up
            demandGroup2(ind2,2) = demandGroup2(ind2,2) + OD_links(i,7)*OD_links(i,9)*connCounters(row,j+2); %regular    
            OD_links_2(k,1) = group2(ind2);
            OD_links_2(k,2) = OD_links(i,6);
            OD_links_2(k,3) = OD_links(i,9);
            OD_links_2(k,4) = OD_links(i,7)/downstrP(row,2);
            k = k + 1;
        end
    elseif OD_links(i,3)>90000  %link vlink
        ind2 = find(row==group2,1); 
        demandGroup2(ind2,1) = demandGroup2(ind2,1) + OD_links(i,7)*OD_links(i,8); %warm-up
        demandGroup2(ind2,2) = demandGroup2(ind2,2) + OD_links(i,7)*OD_links(i,9); %regular
        OD_links_2(k,1) = group2(ind2);
        OD_links_2(k,2) = OD_links(i,6);
        OD_links_2(k,3) = OD_links(i,7).*OD_links(i,9);
        OD_links_2(k,4) = 1;
        k = k + 1;
    else
        %original network vlink - need to generate demand in its upstr vlink
        ind = find(upstrP(row,find(upstrP(row,:)>90000,1,'last'))==downstrP(:,1),1);
        ind2 = find(group2==ind,1);
        demandGroup2(ind2,1) = demandGroup2(ind2,1) + OD_links(i,7)*OD_links(i,8); %warm-up
        demandGroup2(ind2,2) = demandGroup2(ind2,2) + OD_links(i,7)*OD_links(i,9); %regular
        OD_links_2(k,1) = group2(ind2);
        OD_links_2(k,2) = OD_links(i,6);
        OD_links_2(k,3) = OD_links(i,7).*OD_links(i,9);
        OD_links_2(k,4) = 1;
        k = k + 1;
    end
%     %origin link index
%     OD_links_2(i,1) = group2(ind2); 
end
%destination link index

%splitting of trips from every origin link of group2 to possible
%destinations based on the initial demand 
% OD_links_2(:,4) = OD_links_2(:,3); 
% for i=1:length(group2) 
%     ind = OD_links_2(:,1)==group2(i); 
%     aa = sum(OD_links_2(ind,4));
%     OD_links_2(ind,4) = OD_links_2(ind,4)/aa;
% end

%check if all multipliers are correct 
if sum(OD_links_2(:,3)==0)>0 
    error('zero multiplier for origin link! check!')
end

return


