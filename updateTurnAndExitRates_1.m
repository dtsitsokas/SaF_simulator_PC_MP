function [connCounters, stack_OD, junct2] = updateTurnAndExitRates_1(OD_links_2,...
    stack_OD,junct2,t_win,LinksP, A, downstrP, t_ind, connCounters,outflows)
%updateTurnAndExitRates: Calculates new turn and exit rates for a given
%time window for the current and new demand based on the average observed
%travel time of the last time window.

%  OD_links_2: the OD matrix of the current time interval (origin_link_ind,
%              dest_link_ind, demand (veh/hour), multiplier)
%    stack_OD: Ongoing trips [origin_link_index dest_link_index flow]
%       t_win: the time window duration for which we compute average turn rates [hours]
%           A: Connectivity matrix (link to link with travel time in hours) - all links real+virtual
%       t_ind: the index of 15-min interval of the simulation for storing the t
%               rate (1st, 2nd etc.)
%    outflows: simulated outflows of virtual queues (the trips starting) in veh/h
%       junct2
%       LinksP
%       downstrP
%    connCounters: bins as counters of connections between links


%1. New demand: for every OD pair:
%   - find shortest path (1 or more?)
%   - locate point to reach in tWindow (if trip is not finished in tW)
%   - update counters of links and connections
%   - save rest of path in stack for next interval

%2. Unfinished trips: For every path
%   - locate point to reach
%   - update counters
%   - save rest in stack if not finished

% create directed graph with values = travel times
Gr=digraph(A);

%1. Initialize counters - this can happen once outside - (input/load)

% link counters / can distibguish between defac_avg >/= 0! no new demand ->
% stack only ? (virtual queues?) -> I can take as new demand the generated
% inflow from virtual queues + ongoing trips?

%outflows: the actual outflows of all links
%transform outflows to O-Ds as in OD_Links

% initialize stack of ongoing trips (initially empty)
if isempty(stack_OD)
    origins = OD_links_2(:,1);
    destinations = OD_links_2(:,2);
    %flows = OD_links(:,7).*((t_ind>1)*OD_links(:,9)+(t_ind==1)*OD_links(:,8))*defac_avg*t_win;
    flows = outflows(OD_links_2(:,1)).*OD_links_2(:,4);
    spentInside = zeros(size(origins,1),1);
else
    %add stack O-D pairs of ongoing trips from previous cycles
    origins = [OD_links_2(:,1); stack_OD(:,1)];
    destinations = [OD_links_2(:,2); stack_OD(:,2)];
    %flows = [OD_links(:,7).*((t_ind>1)*OD_links(:,9)+(t_ind==1)*OD_links(:,8))*defac_avg*t_win; stack_OD(:,3)]; %veh in 15 mins (the time window)
    flows = [outflows(OD_links_2(:,1)).*OD_links_2(:,4); stack_OD(:,3)];
    spentInside = [zeros(size(OD_links_2,1),1); stack_OD(:,4)];
end
%probl_ODs = [];

%initialize stack_OD of the next interval
stack_OD = [];
P = cell(length(origins),1); %struct to store trips (not necessary to store)
for i=1:size(origins,1) %for every O-D
    if flows(i)>0 %if there is flow (new/waiting)
        %find all shortest paths of the OD matrix based on travel time
        s = origins(i);
        t = destinations(i);
        [P{i}, d] = shortestpath(Gr,s,t); %need to store?
        %P is cell with link indices in LinksP, d is time of trip
        if isempty(P{i})
            %probl_ODs = [probl_ODs; s t];
            %error('empty path - check')
            disp('empty path')
            LinksP(s,1)
            LinksP(t,1)
        else
            if t_win < d-spentInside(i) %trip exceeds the interval of calculation
                dur = 0;
                k = 0;
                %trace path for the interval
                while 1
                    k = k + 1; %index for path elements
                    dur = dur + A(P{i}(k),P{i}(k+1)) - (k==1)*spentInside(i); %current trip duration
                    if dur < t_win
                        % update connection counter of the next pair
                        if downstrP(P{i}(k),2)>1 %if the link k has several downstream links (if 1 we don't care - turn rate is 1)
                            % update connection counter (add the flow taking
                            % this direction)
                            col = downstrP(P{i}(k),:)==LinksP(P{i}(k+1),1);
                            if isempty(col)
                                error('error downstream link')
                            end
                            connCounters(P{i}(k),col) = connCounters(P{i}(k),col) + flows(i);
                            
                            % update upstream link's counter (add the exiting
                            % flow)
                            %linkCounters(P{i}(k)) = linkCounters(P{i}(k)) + flows(i);
                        end
                    else
                        %trip exceeded time window: save rest of path in stack
                        stack_OD = [stack_OD; P{i}(k) P{i}(end) flows(i) abs(t_win-dur)];
                        break %exit loop
                    end
                end
                
            else
                %the whole path is travelled in t_win
                for k=1:length(P{i})-1
                    if downstrP(P{i}(k),2)>1 %count turn rates only for links with several downstream links
                        % update connection counter
                        col = downstrP(P{i}(k),:)==LinksP(P{i}(k+1),1);
                        if isempty(col)
                            error('error downstream link')
                        end
                        connCounters(P{i}(k),col) = connCounters(P{i}(k),col) + flows(i);
                        
                    end
                end
            end
        end
    end
end
% if ~isempty(probl_ODs)
%     probl_ODs
%     %error('no path found for OD')    
%     pause
% end
% Calculate turn rates and assign them to the variables
alpha = 0.01; %all roads will always receive at least alpha

for i=1:size(connCounters,1) %number of real links of the network
    
    a_sum = sum(connCounters(i,3:3+connCounters(i,2)-1)); %sum of counted vehicles
    if a_sum>0 %case: vehicles are counted
        connCounters(i,3:3+connCounters(i,2)-1) = connCounters(i,3:3+connCounters(i,2)-1)./a_sum; %turn ratios
        if sum(connCounters(i,3:3+connCounters(i,2)-1)<=alpha)>0
            
            %adjust zeros - put at least "alpha" at every intersection
            
            sto = connCounters(i,3:3+connCounters(i,2)-1);
            connCounters(i,3:3+connCounters(i,2)-1) = (alpha + sto)/sum(sto + alpha);
            
        end
        %Assign to junct2:
        %(try average of last 3 ?)
        if downstrP(i,2)>1
            for j=1:downstrP(i,2)
                ind1 = junct2.origin == downstrP(i,1);
                ind2 = junct2.destination == downstrP(i,2+j);
                ind = find(ind1+ind2==2,1);
                %junct2.turn(ind,t_ind) = mean([junct2.turn(ind,max(1,t_ind-1):t_ind-1) connCounters(i,j+2)]);
                junct2.turn(ind,t_ind) = connCounters(i,j+2);
                if junct2.turn(ind,t_ind)<0
                    error('neg turn rate assigned 1')
                end
            end
        end
    else
        %case: no vehicles exited counted at this link
        %keep previous turn rate or divide equally if it is the first t_int:
        connCounters(i,3:3+(connCounters(i,2)-1)) = 1/connCounters(i,2); %divide equally
        for j=1:downstrP(i,2)
            %if downstrP(i,2)>1 && downstrP(i,1)<100000
                ind1 = junct2.origin == downstrP(i,1);
                ind2 = junct2.destination == downstrP(i,2+j);
                ind = find(ind1+ind2==2,1);
                if t_ind>1
                    junct2.turn(ind,t_ind) = junct2.turn(ind,t_ind-1); %keep the same as in the previous interval
                    %(problematic at the end of the simulation - vehicles remaining inside 
                    %junct2.turn(ind,t_ind) = ...
                     %   mean(junct2.turn(ind,max(1,t_ind-5):t_ind-1)); %take average of
                    if junct2.turn(ind,t_ind)<0
                        error('neg turn rate assigned 2')
                    end
                    %last 2
                else
                    junct2.turn(ind,t_ind) = connCounters(i,j+2); %assign the equally divided ratio
                    if junct2.turn(ind,t_ind)<0
                        error('neg turn rate assigned 3')
                    end
                end
            %end
        end
    end
end

return



