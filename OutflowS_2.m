function [outflow] = OutflowS_2(or_lanes, dest_lanes, all_or_lanes, dest_capacity, w_orig, x_dest, turnRate,DT)
%OutflowS_2(or_lanes, dest_lanes, all_or_lanes, dest_capacity, w_orig, x_dest, turnRate,DT)
%Links: [LinkID	Number_of_Lanes	Length Starting_Node Ending_Node]
%green, cycle : minutes
%turns : veh/hour
%outflow : veh/hour
%satFlow : veh/hour
%DT : Model time step : hours
%x : vector (queues in this time step)

outflow = (-x_dest+dest_capacity > dest_lanes*1800*DT);
outflow(outflow==0) = outflow(outflow==0) + (x_dest(outflow==0)==0);
%Calculate the or-dest saturation flow s_{z,w} - Enhanced Version
outflow = outflow.*min([or_lanes*1800 dest_lanes*1800 all_or_lanes.*turnRate*1800 turnRate.*(w_orig)/DT],[],2); %(veh/h)
% satFlow=min([or_lanes dest_lanes],[],2)*1800; %(veh/h)

% satFlow=min([or_lanes dest_lanes],[],2)*1800;
%Case: Signals / DONT multiply by g/c - Previous version
%outflow = (x_dest<alpha.*dest_capacity); %Check downstream link
%outflow = outflow.*min([satFlow turnRate.*w_orig/DT],[],2); %[veh/h]
%-------------------------------------

%sumDesTurn = sum of turn rates from all incoming links of the junction to the dest_link if the approach k
%Case Signals / Enhanced SaF version / Dont multiply with g/c - Gia kapoio
%logo den doulevei - ta oxhmata menoun sta virtual queues
%outflow = min([satFlow turnRate.*w_orig/DT (dest_capacity-x_dest).*(turnRate./sumDesTurn)],[],2); %veh/h

% 2 - Apply the simple way
%satFlow = min([or_lanes dest_lanes],[],2)*1800; %(veh/h)
%outflow = (-x_dest+dest_capacity > dest_lanes*1800*DT);
%outflow(outflow==0) = outflow(outflow==0) + (x_dest(outflow==0)==0);
% if sum(outflow>1)>0
%     disp('error outflow')
%     a = find(outflow>1);
%     a
% end
    
%Check downstream link - full if there is space less than 1.2 vehicles
%outflow = outflow.*min([satFlow turnRate.*(w_orig)/DT],[],2); %[veh/h]
%problem with cases of shared lanes (2 directions) - need to make sure that
%the total outflow of the upstream link is <= Sat flow of the link 

end