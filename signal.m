function c = signal(time,cycle,offset,greenTimes)
%SIGNAL: Binary function = 1 if the approach has ROW at time t, 0 otherwise
%--------------------------------------------------------------------------
%PARAMETERS: 
%time:        time (clock) that has passed since the beginning of the simulation [sec]
%cycle:       cycle time [sec]
%offset:      offest time [sec] - time in cycle in t=0
%greenTimes:  Set of pairs or time points within the cycle that the approach has
%             right of way | in the form [t1 t2 t3 t4 ...] which means 
%             (t1,t2) and(t3,t4) the approach has green light 
%--------------------------------------------------------------------------

timeInCycle = rem(offset+time,cycle)';
%timeInCycle(isnan(timeInCycle))=0; 

%Matrix <= 
matSE = timeInCycle <= greenTimes;
%Matrix >= 
matGE = timeInCycle > greenTimes;

c=zeros(size(greenTimes,1),1);
for i=1:2:size(greenTimes,2)
    c = c + matSE(:,i+1).*matGE(:,i);
end

end
