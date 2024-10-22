function [connCounters, stack_OD, junct2] = updateTurnAndExitRates(OD_links,defac_avg,junct2,t_win,LinksP, A, downstrP, t_ind, connCounters)
%updateTurnAndExitRates: Calculates new turn and exit rates for a given
%time window for the current and new demand based on the average observed
%travel time of the last time window.

%    OD_links: the OD matrix of the current time interval
%   defac_avg: the time-dependent multiplier of the OD matrix for the
%               interval
%    stack_OD: Ongoing trips [origin_link_index dest_link_index flow] 
%       t_win: the time window duration for which we compute average turn rates [hours]
%           A: Connectivity matrix (link to link with distances) - all links real+virtual
%      speeds: Vector of calculated speeds per link (for virtual links - free flow)
%       t_ind: the index of 15-min interval of the simulation for storing the t
%               rate (1st, 2nd etc.)
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

% initialize stack of ongoing trips (initially empty)
stack_OD = []; 
origins = OD_links(:,5);
destinations = OD_links(:,6);
flows = OD_links(:,8)*defac_avg*t_win;

P = cell(length(origins),1);
for i=1:size(origins,1) %for every O-D
    if flows(i)>0 %if there is flow (new/waiting)
        %find all shortest paths of the OD matrix based on travel time
        s = origins(i);
        t = destinations(i);
        [P{i}, d] = shortestpath(Gr,s,t); %need to store?
        %P is cell with link indices in LinksP, d is time of trip
        
        if t_win < d %trip exceeds the 15 mins interval 
            dur = 0;
            k = 0;
            %trace path for the 15 mins interval
            while 1
                k = k + 1; %index for path elements
                dur = dur + A(P{i}(k),P{i}(k+1)); %current trip duration
                if dur < t_win
                    % update connection counter of the next pair
                    if downstrP(P{i}(k),2)>1 %if the link k has several downstream links (if 1 we don't care - turn rate is 1) 
                        % update connection counter (add the flow taking
                        % this direction) 
                        col = downstrP(P{i}(k),:)==LinksP(P{i}(k+1),1);
                        connCounters(P{i}(k),col) = connCounters(P{i}(k),col) + flows(i);
                       
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
                    connCounters(P{i}(k),downstrP(P{i}(k),:)==LinksP(P{i}(k+1),1)) = connCounters(P{i}(k),downstrP(P{i}(k),:)==LinksP(P{i}(k+1),1)) + flows(i);
                    
                end
            end
        end
    end
end

% Calculate turn rates and assign them to the variables
alpha = 0.01; %all roads will always receive at least alpha 

for i=1:size(connCounters,1)
    %if linkCounters(i)>0
    a_sum = sum(connCounters(i,3:3+connCounters(i,2)-1));
     if a_sum>0 %case: vehicles are counted 
        connCounters(i,3:3+connCounters(i,2)-1) = connCounters(i,3:3+connCounters(i,2)-1)./a_sum; %turn ratios 
        if sum(connCounters(i,3:3+connCounters(i,2)-1)<=alpha)>0 
            %adjust zeros - put at least "alpha" at every intersection 
            sto = connCounters(i,3:3+connCounters(i,2)-1); 
            connCounters(i,3:3+connCounters(i,2)-1) = (alpha + sto)/sum(sto + alpha); 

        end
        %Assign to junct2:
        for j=1:downstrP(i,2)
            ind1 = junct2.origin == downstrP(i,1);
            ind2 = junct2.destination == downstrP(i,2+j);
            ind = find(ind1+ind2==2,1);
%             if t_ind < 2
                junct2.turn(ind,t_ind) = connCounters(i,j+2);
%             elseif t_ind > 10
%                 junct2.turn(ind,t_ind) = mean([connCounters(i,j+2) junct2.turn(ind,9:t_ind-1)]);
%             elseif t_ind > 16
%                 junct2.turn(ind,t_ind) = mean([connCounters(i,j+2) junct2.turn(ind,1:t_ind-1)]); 
%             else
%                 junct2.turn(ind,t_ind) = mean([connCounters(i,j+2) junct2.turn(ind,t_ind-4:t_ind-1)]);
%             end
            % Check maximum change 
            if junct2.turn(ind,t_ind)<0
                error('neg turn rate assigned 1')
            end
        end
    else 
        %case: no vehicles exited counted at this link
        %keep previous turn rate or divide equally if it is the first t_int:
        connCounters(i,3:3+(connCounters(i,2)-1)) = 1/connCounters(i,2); %divide equally 
        for j=1:downstrP(i,2)
            ind1 = junct2.origin == downstrP(i,1);
            ind2 = junct2.destination == downstrP(i,2+j);
            ind = find(ind1+ind2==2,1);
            if t_ind>1
                junct2.turn(ind,t_ind) = junct2.turn(ind,t_ind-1); %keep the same as in the previous interval
                if junct2.turn(ind,t_ind)<0
                    error('neg turn rate assigned 2')
                end
            else
                junct2.turn(ind,t_ind) = connCounters(i,j+2); %assign the equally divided ratio 
                if junct2.turn(ind,t_ind)<0
                    error('neg turn rate assigned 3')
                end
            end
        end
    end
end 

return



