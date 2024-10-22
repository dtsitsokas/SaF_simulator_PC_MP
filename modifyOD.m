function [OD] = modifyOD(OD)
%Modify OD matrix to redistribute demand from certain centroids to
%neighboring centroids 

%Centroid 69087 (south-east) - long delays and high average
%density/spillbacks 
centroid = 69087; 
recipients = [69081 69083 69067];
p_remove = 0.20; %percentage of removal  
[OD] = redistCentroidDemand(OD,centroid,recipients,p_remove); 

%Centroid 69101 (west) 
centroid = 69101; 
recipients = [68951 68997 68994];
p_remove = 0.20; %percentage of removal  
[OD] = redistCentroidDemand(OD,centroid,recipients,p_remove); 

%Centroid 56508 (internal) 
centroid = 56508; 
recipients = [57796 68947 57104];
p_remove = 0.40; %percentage of removal  
[OD] = redistCentroidDemand(OD,centroid,recipients,p_remove); 

return 