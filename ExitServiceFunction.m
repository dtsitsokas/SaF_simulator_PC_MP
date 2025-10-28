function [p] = ExitServiceFunction(x, c)
%EXITSERVICEFUNCTION: Calculates the percentage of maximum exit
% service rate as a function of network current accumulation 
%
% x: the array with current link accumulations 
% c: array with link capacities 

% p: the percentage of available capacity 

occ = x./c;

if mean(occ) > 0.25
    p = max(0.10, 0.8-mean(occ));
else
    p = 1;
end

end

