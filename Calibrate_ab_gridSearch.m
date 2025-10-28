clc
clear
close all
% Final a,b calibration / Grid Search 
% Simple method to find best a,b parameters for MP node selection by performing a GRID SEARCH - 

demandCode = 1; %1: med demand, 2: high demand
%ab_1 = zeros(121,3); 
k=0;
if demandCode == 1 %med
    pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
    %output file name to save after simulation is done
    path = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
    %fname = strcat(path,'output_set5_MedDemand_NC'); %NC of med
elseif demandCode == 2
    pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
    path = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';
    %fname = strcat(path,'output_set5_highDemand_NC');
end

% for a=0.1:0.1:1
%     for b=1.5:0.1:3.5
%         k = k + 1;
%         ab_1(k,1) = a; 
%         ab_1(k,2) = b;  
%     end 
% end


% for a=0.1:0.1:1
%     for b=0.1:0.1:1
%         k = k + 1;
%         ab_1(k,1) = a; 
%         ab_1(k,2) = b;  
%     end 
% end

for a=-2:0.3:2
    for b=-2.5:0.3:2.5
        k = k + 1;
        ab_1(k,1) = a; 
        ab_1(k,2) = b;  
    end 
end

caseName = 'ab_test38_oldCriteria_eq2_20prcnt_med';

% ab_1_med = ab_1;
%load('ab_test11_eq3_25prcnt_high.mat','ab_1')
for i=1:size(ab_1,1)
	load(strcat(pathInput,'indata.mat'),'indata')
    round(i/size(ab_1,1),2)
    fname = strcat(path,'ab_test38_',num2str(i),'.mat');
    savefull = 0; %save big result files
    %PC.mode = 0; %otherwise PC read from node selection 
    case_j = 3; %3 for 20%, 4 for 25% etc.  
    equation = 2;
    PCmode = 0; 
    MPmode = 1;
    ab_1(i,:)
    MPnodes = nodeSelectionMP(demandCode,ab_1(i,1),ab_1(i,2),case_j,equation,caseName);
    prepareExperiment(caseName,demandCode,PCmode,MPmode)
    load(strcat('input_',caseName,'.mat'))
    ab_1(i,3) = SaF_3(indata,MPnodes,PC,fname,savefull);
    close all 
end


%save('ab_1.mat','ab_1')
%save('ab_2.mat','ab_1')
save(strcat(caseName, '.mat'),'ab_1')

[~,sort_ind] = sort(ab_1(:,3));
sorted_ab_1 = ab_1(sort_ind,:);


%% Combine all values of a,b for a set of different runs (like in robust optimization) 

% Consider only tests for equation #3 (med/high, 10%/25%) -> 4 tests 

% %tests{1} = 'ab_test8_eq3_10prcnt_med.mat';
% %tests{2} = 'ab_test10_eq3_30prcnt_med.mat';
% tests{1} = 'ab_test9_eq3_10prcnt_high.mat';
% tests{2} = 'ab_test11_eq3_25prcnt_high.mat';
% path = 'F:\Max Pressure files\fifth set results\ab_calib_results\';
% 
% for i=1:2
%    %load the results 
%    load(strcat(path,tests{i}),'ab_1') 
%    ab_all(:,i) = ab_1(:,3); 
%    
% end

%normalize results 
% ab_all(:,1) = ab_all(:,1)/221400*0.25; 
% ab_all(:,3) = ab_all(:,3)/221400*0.25;
% ab_all(:,2) = ab_all(:,2)/484000*0.25;
% ab_all(:,4) = ab_all(:,4)/484000*0.25;


% ab_all(:,1) = ab_all(:,1)/221400; 
% ab_all(:,3) = ab_all(:,3)/221400;
% % ab_all(:,2) = ab_all(:,2)/484000;
% % ab_all(:,4) = ab_all(:,4)/484000;
% 
% %find the weighted sum of all cases 
% ab_sum = sum(ab_all,2);
%    
% %rank the sums and rearrange the ab matrix 
% [ranked_ab_sum, ranked_ind] = sort(ab_sum);
% ab = ab_1(ranked_ind,[1 2]); 
% 
% %best performing a,b 
% ab(1,:) 




