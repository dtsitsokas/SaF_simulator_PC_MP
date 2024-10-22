function [OD] = redistCentroidDemand(OD,centroid,recipients,p_remove)

ind_c = OD(:,1) == centroid; 
ind = zeros(1,length(recipients));
%adding - distribute to recipients equally 
for i=1:length(recipients)
    ind(i) = find(recipients(i) == OD(:,1));
end 

%adding
OD(ind, 2:end) = OD(ind, 2:end) + ones(length(recipients),1)*...
    (1/length(recipients)*p_remove*OD(ind_c,2:end));

%removing 
OD(ind_c,2:end) = OD(ind_c,2:end)*(1-p_remove);
% fin_sum = sum(sum(OD(2:end, 2:end)))
end

