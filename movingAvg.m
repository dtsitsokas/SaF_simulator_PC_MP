function [x] = movingAvg(x,step)
%MOVINGAVG: Produce smoother vector by applying moving average 
%   Inputs
%       x: the non-smooth vector 
%    step: the number of steps before and after for the mov avg window 
%  Outputs
%       x: the smoothened vector

if 2*step+1>length(x)
    error('window for averaging too wide')
end
x_i = x; 
for i=1:length(x)
   x(i) = mean(x_i(max(1,i-step):min(length(x),i+step))); 
end

end

