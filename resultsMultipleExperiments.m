%Script to produce Fig. 1: performance of MP only cases (random & selected
%cases) together with case 100%, wrt NC case (absolut values + percentage change)  
% Script for the case of FULL result files (MEDIUM OD)/ Boxplots (1)

clear 
clc
close all 

%Load scenario names 
load scenarios.mat

% path of NC and all MP selected result files 
pathNC = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
load(strcat(pathNC,'output_set5_MedDemand_NC.mat'),'outdata','indata') %NC file
pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';

% path of the result files of the random experiments 
path = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\RandomSetsResults\MP only\';
pathFig = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\Figures Results\';

% Horizontal axis - Node penetration rate 
perc = [5 10 15 20 25]; 

%% Initialize comparative figures 

%Create all figures first:
f1 = figure('Name','VHT - percent MP nodes');
%title('random MP node selection');
xlabel('MP node penetration rate (\%)');
ylabel('VHT');
box on
hold on

%%
f2 = figure('Name','VHT improvement - percent MP nodes');
%title('random MP node selection');
xlabel('MP node penetration rate (\%)');
ylabel('VHT change wrt FTC (\%)');
box on
hold on


%% No-Control Results: 1st column of the Table  
t = 90; %no of steps to aggregate for speed calculation
disp(strcat('Scenario: ',scenarios(1)))

%Total time spent in network and VQs -
%PHT_ntw
NCresults(1,1) = sum(sum(outdata.x*indata.DT));

%PHT_VQ
NCresults(2,1) = sum(sum(outdata.w(indata.group2,:)*indata.DT));

%Total travel time
%PHT_tot
NCresults(3,1) = NCresults(1,1) + NCresults(2,1);

%Vehicles remaining inside (not entirely serviced)
NCresults(4,1) = outdata.notserviced;

%TTT with penalty of remaining vehs 
NCresults(5,1) = NCresults(3,1) + NCresults(4,1)*0.25;

%Total waiting time in queues (network) 
NCresults(6,1) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(veh x hours)

%Mean queue over time and space 
NCresults(7,1) = mean(mean(outdata.x(indata.group1,:))); %mean queue over time and space 

%Mean of mean of Virtual Queues over time - mean virtual queue  
NCresults(8,1) = mean(mean(outdata.virtualqueues));

%Mean total link outflow (links and VQs) - mean link outflow  
NCresults(9,1) = mean(mean(outdata.u([indata.group1; indata.group2],:))); 

% Mean Space Mean Spead over time
SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
j=1;
for k=t:t:indata.kmax
    non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
    s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
    s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
        (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
    SMS(j) = mean(s);
    j = j + 1;
end
NCresults(10,1) = mean(SMS);


%% Load results of experiments - Create exper files: summary of results   

% TTT = zeros(10,6); %(no of replic, no of exper) 
% for j =1:length(perc) %No of experiment sets (different percentage of nodes for MP) 
%     
%     %load experiment input (indata is the same for all) 
%     load(strcat(path,'output_MP_R_',num2str(perc(j)),'_',num2str(1),'.mat'),'indata')
%     exper.percent = perc(j);
%     results = zeros(10,10); %adjust this - columns = no of replications
%     exper.id = j; 
%     for i=1:10 %id of replication 
% 
%         load(strcat(path,'output_MP_R_',num2str(exper.percent),'_',num2str(i),'.mat'),'outdata')
%         exper.nodes{i} = outdata.MPnodes;
%         
%         % Table of Results
% 
%         %Total time spent in network and VQs -
%         %PHT_ntw
%         results(1,i) = sum(sum(outdata.x*indata.DT));
% 
%         %PHT_VQ
%         results(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));
% 
%         %Total travel time
%         %PHT_tot
%         results(3,i) = results(1,i) + results(2,i);
% 
%         %Vehicles remaining inside 
%         results(4,i) = outdata.notserviced; 
% 
%         %Balanced PHT_tot (with penalty)
%         results(5,i) = results(3,i) + results(4,i)*0.25; 
%         
%         %Total waiting time in queues (network) 
%         results(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)
%         
%         %Mean queue over time and space 
%         results(7,i) = mean(mean(outdata.x(indata.group1,:)));
%         
%         %Mean of mean of Virtual Queues over time
%         results(8,i) = mean(mean(outdata.virtualqueues));
%         
%         %Mean total link outflow (links and VQs)
%         results(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));
%         
%         SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
%         l=1;
%         for k=t:t:indata.kmax
%             non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
%             s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
%             s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
%                 (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
%             SMS(l) = mean(s);
%             l = l + 1;
%         end        
%         results(10,i) = mean(SMS);
%     end
% 
%     %Statistics - overall results of experiment 
% 
%     %Mean (1st col) and Variance (2nd column) of all metrics  
%     for i=1:10 
%         exper.results(i,1) = mean(results(i,:)); %mean 
%         exper.results(i,2) = std(results(i,:)); %variance 
%     end
%     
%     [~,exper.best_repl_index(j)] = min(results(5,:));
%     
%     % Average change wrt NC scenario 
%     exper.improv = (exper.results(:,1)-NCresults(:,1))./NCresults(:,1)*100; %percentage of improv in TTT 
%     
%     exper.all_results = results; 
%     
%     TTT(:,j) = results(5,:)'; 
%     
%     
%     save(strcat(path,'exp_',num2str(exper.id),'_res.mat'), 'exper')
% 
% end

%% alternative way (time saving) of printing results after creating exper.mat files:

TTT2 = zeros(10,length(perc)); %(no of replic, no of exper) 
for j=1:length(perc)
   load(strcat(path,'exp_',num2str(j),'_res.mat'), 'exper') 
   TTT2(:,j) = exper.all_results(5,:)'; 
end


%% Print comparative figures 
% if exper.mat files are used comment the following 2 lines: 
%TTT2 = TTT;
%clear TTT

% TTT (boxplot) vs. percentage of nodes with MP 
figure(f1)
boxplot(TTT2) 
hold on 
ax = gca;
ax.TickLabelInterpreter = 'latex'; 
xticklabels({'$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$30\%$','$100\%$'})

%saveas(gcf,'TTT_boxplots.fig');
%%
figure(f2)
TTT_imp = (TTT2-NCresults(5,1))/NCresults(5,1)*100;
boxplot(TTT_imp) 
hold on 
%ax = gca;
%ax.TickLabelInterpreter = 'latex'; 
%xticklabels({'$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$30\%$'})
%saveas(gcf,'TTT_boxplots_percent.fig');
%print('TTT_boxplots_percent.eps','-depsc','-painters')

%% Load scenarios of selected nodes by selection algorithm
Sresults = zeros(10,length(perc));
for j=0:length(perc)-1
   %load results file - set file name convention 
   load(strcat(pathNC,'output_set5_MedDemand_MP_cd_MS',num2str(j)),'outdata') %the results
   load(strcat(pathIn,'input_MP_cd_MS',num2str(j)),'p_1') %the percentage
   %store VHT      
        % Table of Results
        i = j + 1; 
        %Total time spent in network and VQs -
        %PHT_ntw
        Sresults(1,i) = sum(sum(outdata.x*indata.DT));

        %PHT_VQ
        Sresults(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));

        %Total travel time
        %PHT_tot
        Sresults(3,i) = Sresults(1,i) + Sresults(2,i);

        %Vehicles remaining inside 
        Sresults(4,i) = outdata.notserviced; 

        %Balanced PHT_tot (with penalty)
        Sresults(5,i) = Sresults(3,i) + Sresults(4,i)*0.25; 
        
        %Total waiting time in queues (network) 
        Sresults(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)
        
        %Mean queue over time and space 
        Sresults(7,i) = mean(mean(outdata.x(indata.group1,:)));
        
        %Mean of mean of Virtual Queues over time
        Sresults(8,i) = mean(mean(outdata.virtualqueues));
        
        %Mean total link outflow (links and VQs)
        Sresults(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));
        
        SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
        l=1;
        for k=t:t:indata.kmax
            non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
            s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
            s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
                (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
            SMS(l) = mean(s);
            l = l + 1;
        end        
        Sresults(10,i) = mean(SMS);
        Sresults(11,i) = p_1; 

end



%% Load the 100% case

%load results file
load(strcat(pathNC,'output_set5_MedDemand_MP_0'),'outdata')
p_1 = 100;
%store VHT
% Table of Results
i = length(perc) + 1;
%Total time spent in network and VQs -
%PHT_ntw
Sresults(1,i) = sum(sum(outdata.x*indata.DT));

%PHT_VQ
Sresults(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));

%Total travel time
%PHT_tot
Sresults(3,i) = Sresults(1,i) + Sresults(2,i);

%Vehicles remaining inside
Sresults(4,i) = outdata.notserviced;

%Balanced PHT_tot (with penalty)
Sresults(5,i) = Sresults(3,i) + Sresults(4,i)*0.25;

%Total waiting time in queues (network)
Sresults(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)

%Mean queue over time and space
Sresults(7,i) = mean(mean(outdata.x(indata.group1,:)));

%Mean of mean of Virtual Queues over time
Sresults(8,i) = mean(mean(outdata.virtualqueues));

%Mean total link outflow (links and VQs)
Sresults(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));

SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
l=1;
for k=t:t:indata.kmax
    non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
    s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
    s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
        (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
    SMS(l) = mean(s);
    l = l + 1;
end
Sresults(10,i) = mean(SMS);
Sresults(11,i) = p_1;


%% Merge results

TTT_imp_s = (Sresults(5,:)-NCresults(5,1))/NCresults(5,1)*100;
% merge all improvements (medians of random and selections) 
med_imp = [median(TTT_imp) TTT_imp_s];
med_imp_perc = [perc Sresults(11,:)];
med_TTT = [median(TTT2) Sresults(5,:)];


%% Load scenarios of selected nodes by the benchmark case of max outflow

Sresults = zeros(10,length(perc));
for j=0:length(perc)-1
   %load results file - set file name convention 
   load(strcat(pathNC,'output_final_med_MP_cd_MSR',num2str(j)),'outdata') %the results
   %load(strcat(pathIn,'input_MP_cd_MSR',num2str(j)),'p_1') %the percentage
   p_1 = (0.05 + j*0.05); 
   %store VHT      
        % Table of Results
        i = j + 1; 
        %Total time spent in network and VQs -
        %PHT_ntw
        Sresults(1,i) = sum(sum(outdata.x*indata.DT));

        %PHT_VQ
        Sresults(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));

        %Total travel time
        %PHT_tot
        Sresults(3,i) = Sresults(1,i) + Sresults(2,i);

        %Vehicles remaining inside 
        Sresults(4,i) = outdata.notserviced; 

        %Balanced PHT_tot (with penalty)
        Sresults(5,i) = Sresults(3,i) + Sresults(4,i)*0.25; 
        
        %Total waiting time in queues (network) 
        Sresults(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)
        
        %Mean queue over time and space 
        Sresults(7,i) = mean(mean(outdata.x(indata.group1,:)));
        
        %Mean of mean of Virtual Queues over time
        Sresults(8,i) = mean(mean(outdata.virtualqueues));
        
        %Mean total link outflow (links and VQs)
        Sresults(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));
        
        SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
        l=1;
        for k=t:t:indata.kmax
            non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
            s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
            s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
                (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
            SMS(l) = mean(s);
            l = l + 1;
        end        
        Sresults(10,i) = mean(SMS);
        Sresults(11,i) = p_1; 

end
%% Merge results

TTT_imp_s = (Sresults(5,:)-NCresults(5,1))/NCresults(5,1)*100;

med_imp = [med_imp TTT_imp_s];
med_imp_perc = [med_imp_perc Sresults(11,:)];
med_TTT = [med_TTT Sresults(5,:)];


%%

figure(f2)
ylim([-25,15])
xlim([0,7])
scatter((1:length(perc)), med_imp(1:length(perc)),'*','SizeData',60) %median of random exp
scatter((1:length(perc)), med_imp(length(perc)+1:2*length(perc)),'v','filled','SizeData',60) %result of selections 
scatter((1:length(perc)), med_imp(2*length(perc)+2:end),'o','filled','SizeData',60) %result of selections
scatter(length(perc)+1, med_imp(2*length(perc)+1),'sq','filled','SizeData',60) % MP_0
ax = gca;
ax.TickLabelInterpreter = 'latex'; 
xticklabels({'$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$100\%$'})
yline(0,'Color','r','LineWidth',1.5)
title('medium demand')
legend('random MP node sets (median)','MP (proposed)','MP (max outflow)','MP to all nodes','Fixed-Time Control','Location','northeast')
saveas(gcf,strcat(pathFig,'Fig_1a_MP_MedDemand.fig'));
saveas(gcf,strcat(pathFig,'Fig_1a_MP_MedDemand.emf'));
print(gcf,strcat(pathFig,'Fig_1a_MP_MedDemand.eps'),'-depsc','-painters')

figure(f1)
xlim([0,7])
ylim([170000,250000])
scatter((1:length(perc)), med_TTT(1:length(perc)),'*','SizeData',60)
scatter((1:length(perc)), med_TTT(length(perc)+1:2*length(perc)),'v','filled','SizeData',60) %result of selections 
scatter((1:length(perc)), med_TTT(2*length(perc)+2:end),'o','filled','SizeData',60) %result of selections
scatter(length(perc)+1, med_TTT(2*length(perc)+1),'filled','SizeData',60)
yline(NCresults(5),'Color','r','LineWidth',1.5)
ax = gca;
ax.TickLabelInterpreter = 'latex'; 
xticklabels({'$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$30\%$','$100\%$'})
title('medium OD')
legend('random MP node sets (median)','MP (proposed)','MP (max outflow)','MP to all nodes','Fixed-Time Control','Location','northeast')

%% Add results of gradually selected nodes 
% ISm_results = zeros(size(Sresults,1),length(perc)-1);
% ISm_results(:,1) = Sresults(:,1);
% for i=2:length(perc)
%    %load results file 
%    load(strcat(pathNC,'output_set4_MP_ISm_',num2str(i)),'outdata')
%    load(strcat('input_MP_ISm_',num2str(i)),'p_1')
%    %store VHT      
%         % Table of Results
%  
%         %Total time spent in network and VQs -
%         %PHT_ntw
%         ISm_results(1,i) = sum(sum(outdata.x*indata.DT));
% 
%         %PHT_VQ
%         ISm_results(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));
% 
%         %Total travel time
%         %PHT_tot
%         ISm_results(3,i) = ISm_results(1,i) + ISm_results(2,i);
% 
%         %Vehicles remaining inside 
%         ISm_results(4,i) = outdata.notserviced; 
% 
%         %Balanced PHT_tot (with penalty)
%         ISm_results(5,i) = ISm_results(3,i) + ISm_results(4,i)*0.25; 
%         
%         %Total waiting time in queues (network) 
%         ISm_results(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)
%         
%         %Mean queue over time and space 
%         ISm_results(7,i) = mean(mean(outdata.x(indata.group1,:)));
%         
%         %Mean of mean of Virtual Queues over time
%         ISm_results(8,i) = mean(mean(outdata.virtualqueues));
%         
%         %Mean total link outflow (links and VQs)
%         ISm_results(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));
%         
%         SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
%         l=1;
%         for k=t:t:indata.kmax
%             non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
%             s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
%             s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3))./...
%                 (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
%             SMS(l) = mean(s);
%             l = l + 1;
%         end        
%         ISm_results(10,i) = mean(SMS);
%         ISm_results(11,i) = p_1; 
% 
% end
% TTT_imp_ISm = (ISm_results(5,:)-NCresults(5,1))/NCresults(5,1)*100;

% figure(f2)
% scatter((1:length(perc)), TTT_imp_ISm,'filled','SizeData',60)
% scatter(length(perc)+1, med_imp(end),'filled','SizeData',60)
% legend('random MP node sets (median)','selection MP nodes','incremental selection MP nodes','all MP nodes','Location','northeast')
% saveas(gcf,'TTT_impr_all.fig');
% print('TTT_impr_all.eps','-depsc','-painters')
