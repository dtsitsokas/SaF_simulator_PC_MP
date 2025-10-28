%function to change demand values indata for running tests with values of
%OD as random values of normal distribution with mean the initial value and
%std a set percentage of it 

function [OD] = generateNormDistrOD(OD, pct) 

%Goal: replace OD mean values woth random variables drew from normal distr
%with variance as percentage of mean 

% OD = vector with mean demand per OD (veh/h) 
% pct = percentage of the mean to use as variance 

for i=1:length(OD)
    flag = 0; 
    while flag==0  
        newOD =  normrnd(OD(i),pct*OD(i));
        % accept only non-negative demand 
        if newOD >=0 
            OD(i) = newOD; 
            flag = 1; 
        end
    end
end

end