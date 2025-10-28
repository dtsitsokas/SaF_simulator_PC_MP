function [greentimes2,MP] = maxPressure(to_update,MP,junct,capacity,satFlow,greentimes2,x,k,stepsInCycle_MP,t_index,gn_R)
% MAXPRESSURE %Calculates new greens and 
% Detailed explanation goes here

%update some traffic signal settigns
% to_update = MPnodes(to_update); %indices in MP struct of nodes to update -
%WHY? 
%update process
for jj = 1:length(to_update)%for every node that is to update
    i = to_update(jj); 
    st = max(1,k-stepsInCycle_MP(i)+1);  %start of the interval 
    %1.for every incoming link of every eligible
    %stage of the node Compute pressure
    node_pressures = []; %store the pressures of the node
    for s = MP.stagesInvolved{i} %s = eligible stage of node i in MP
        if MP.stageApproaches{i,s}~=0 %if the stage is not all red (for large AR stages that are eligible)
            stage_p = [];%store all the pressures of the current stage s
            incoming = unique(junct.or_index(MP.stageApproaches{i,s})); %unique Link_ind in Links() of all incoming links
            for inc = incoming %for every incoming link of the stage
                incInd = (junct.or_index(MP.stageApproaches{i,s})==inc);
                incInd = MP.stageApproaches{i,s}(incInd); %ind in junct of all approaches starting from the incoming link with index inc
                
                %calculate the pressure (either based on max or based on
                %average)
                
                %use max queues upstream and downstream (option 1)
                %stage_p = [stage_p (min(1,max(x(inc,st:k))/capacity(inc))-sum(indata.junct.turn(incInd,t_index).*min(1,max(x(indata.junct.dest_index(incInd),st:k))./capacity(indata.junct.dest_index(incInd)))))*satFlow(inc)];
                
                %use average queues upstream and downstream (option 2)
                stage_p = [stage_p max(0,(min(1,mean(x(inc,st:k))/capacity(inc))-sum(junct.turn(incInd,t_index).*min(1,mean(x(junct.dest_index(incInd),st:k),2)./capacity(junct.dest_index(incInd)))))*satFlow(inc))];
                % stage_p = [stage_p max(0,(min(1,max(x(inc,st:k))/capacity(inc))-sum(junct.turn(incInd,t_index).*min(1,max(x(junct.dest_index(incInd),st:k),2)./capacity(junct.dest_index(incInd)))))*satFlow(inc))];
                % stage_p = [stage_p max(0,(min(1,avg_q_cycle(inc)/capacity(inc))-sum(junct.turn(incInd,t_index).*min(1,mean(avg_q_cycle(junct.dest_index(incInd)),2)./capacity(junct.dest_index(incInd)))))*satFlow(inc))];

            end
            
            %define stage pressure (the max of the calculated pressures in the stage)
            node_pressures = [node_pressures max(0,sum(stage_p))];
        else
            node_pressures = [node_pressures 0];
        end
    end
    
    %distribute effective green to eligible phases
    greenP_p = MP.stageDur{i}(MP.stagesInvolved{i}); %the previously assigned greens of the involved stages
    %the real values of the greens of the stages involved
    greenP = (sum(node_pressures)>0)*round(node_pressures/sum(node_pressures)*MP.timeToAllocate(i)) + (sum(node_pressures)==0)*0;
    
    if sum(isnan(greenP))>0
        greenP(isnan(greenP)) = 0;
    end
    
    if sum(node_pressures)==0 %if all pressures zero - no change, keep prev settings
        greenP = greenP_p;
    else
        %the integer projection of the real values of greens (add it here and replace the code below)
        %(constraints: cycle maintained, minimums and max changes satisfied)
        
        greenP = feasibleGreensFromPressures(greenP,MP.stageBaseDur{i}(MP.stagesInvolved{i}),greenP_p,gn_R);
        
        %assign new stage durations
        MP.stageDur{i} = MP.stageBaseDur{i}; %the min for all (maybe unnecessary)
        MP.stageDur{i}(MP.stagesInvolved{i}) = greenP; %the newly calculated for the stages involved
        
        %update greentimes - for every approach invlovled:
        %all approaches were affected
        appr = MP.approaches{i}; %all approaches of the intersection
        for j=1:length(appr)
            s = junct.stageNo(appr(j),junct.stageNo(appr(j),:)>0); %stages where the approach takes green
            for ss=1:length(s)
                greentimes2(appr(j),ss*2) = sum(MP.stageDur{i}(1:s(ss)));  %end time of respective stages
                greentimes2(appr(j),ss*2-1) = sum(MP.stageDur{i}(1:s(ss)-1)); %start time of repsective stages
            end
        end
        
        
    end
 
end
end



