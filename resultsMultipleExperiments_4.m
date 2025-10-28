%Script to produce Fig. 1: performance of PC solo and PC+MP cases (random & selected
%cases) together with case 100%, wrt NC case (absolut values + percentage change)  
% Script for the case of REDUCED result files (HIGH OD) 

clear 
clc
% close all 

%Load scenario names 
load scenarios.mat

%% Variables to set: 
% path of NC and all MP selected result files (Medium or High Demand results)  

%pathNC = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Current results\ResultsHighOD_capacityDrop\';
pathNC = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';
load(strcat(pathNC,'output_set5_highDemand_NC.mat'),'outdata','indata') %NC file

% path of the result files of the random experiments (combo PC + MP random)
path = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\RandomMPPlusPC_SetResults\';
pathFig = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\Figures Results\';
pathIN = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs high demand\final runs input files\';

% Horizontal axis - Node penetration rate 
perc = [5 10 15 20 25]; 
demand = 'high'; 
figfName = 'Fig_2b_PCMP_HighDemand_a';
fname_PCMP100 = 'output_final_high_PC_501b'; 
fname_PC = 'output_final_high_PC_50b'; 
fname_prothema_PCMPR = 'output_PC_61_MP_R'; %format fname of random PC_MP results 
fnameinput_prothema_MPMS = 'input_PC_50b_MPf_MS';
% fnameinput_prothema_MPMS ='input_PC_50b_MPnc_MS'; %for the parameters of med
% fnameOutput_prothema_PCMPMS = 'output_set5_HighDemand_PC_50b_MPnc_MS'; %for the parameters of med to high
fnameOutput_prothema_PCMPMS = 'output_final_high_PC_50b_MP_MS';
fnameOutput_prothema_benchmark = 'output_final_high_PC_50b_MP_MSR';

%% Initialize comparative figures 

%Create all figures first:
f1 = figure('Name','VHT - percent MP nodes');
xlabel('MP node penetration rate (\%)');
ylabel('VHT');
box on
hold on
%%
f2 = figure('Name','VHT improvement - percent MP nodes');
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


%% Load reduced results of experiments 
TTT2 = zeros(10,length(perc)); %line = replication, column = percentage case 
for j =1:length(perc)
    
    
    for i=1:10 %replications 
        
        % load result file (indata is the same for all)
        load(strcat(path,fname_prothema_PCMPR,num2str(perc(j)),'_',num2str(i),'.mat'))
        
        % save TTT 
        TTT2(i,j) = r3 ; %+ notserviced*0.25; 
       
    end
end

%% alternative way (time saving) of printing results after creating exper.mat files:

% TTT2 = zeros(10,6); %(no of replic, no of exper) 
% for j=1:length(perc)
%    load(strcat(path,'exp_',num2str(j),'_res.mat'), 'exper') 
%    TTT2(:,j) = exper.all_results(5,:)'; 
% end
TTT2 = [zeros(size(TTT2,1),1) TTT2];
TTT2 = [TTT2 zeros(size(TTT2,1),1)];
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


%saveas(gcf,'TTT_boxplots.fig');
%%
figure(f2)
TTT_imp = (TTT2-NCresults(5,1))/NCresults(5,1)*100;
boxplot(TTT_imp) 
hold on 


%% Load scenarios of selected nodes by selection algorithm
Sresults = zeros(10,length(perc));
for j=0:length(perc)-1
   %load results file - set file name convention 
   load(strcat(pathNC,fnameOutput_prothema_PCMPMS,num2str(j)),'outdata')
   load(strcat(pathIN,fnameinput_prothema_MPMS ,num2str(j)),'p_1')
   
   %store VHT      
        % Table of Results
        i = j + 2; 
        
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

%% Load the PC solo case 

%load results file
load(strcat(pathNC,fname_PC),'outdata')
p_1 = 0;
%store VHT
% Table of Results
i = 1;
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


%% Load the PC + MP 100% case

%load results file
load(strcat(pathNC,fname_PCMP100),'outdata')
p_1 = 100;
%store VHT
% Table of Results
i = length(perc) + 2;
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


%%

TTT_imp_s = (Sresults(5,:)-NCresults(5,1))/NCresults(5,1)*100;
% merge all improvements (medians of random and selections) 
med_imp = [median(TTT_imp(:,2:end-1)) TTT_imp_s];
med_imp_perc = [perc Sresults(11,:)];
med_TTT = [median(TTT2(:,2:end-1)) Sresults(5,:)];


%% Load scenarios of selected nodes (benchmark, max outflow) 
Sresults = zeros(10,length(perc));
for j=0:length(perc)-1
   %load results file - set file name convention 
   load(strcat(pathNC,fnameOutput_prothema_benchmark,num2str(j)),'outdata')
   %load(strcat(pathIN,fnameinput_prothema_MPMS ,num2str(j)),'p_1')
   p_1 = 5 + j*5; 
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

%%

TTT_imp_s = (Sresults(5,:)-NCresults(5,1))/NCresults(5,1)*100;
% merge all improvements (medians of random and selections) 
med_imp = [med_imp TTT_imp_s];
med_imp_perc = [med_imp_perc Sresults(11,:)];
med_TTT = [med_TTT Sresults(5,:)];


%%
figure(f2)
scatter(1,med_imp(length(perc)+1),'^','filled','SizeData',60) %results PC solo
scatter((2:length(perc)+1), med_imp(1:length(perc)),'*','SizeData',60) %median of random exp
scatter((2:length(perc)+1), med_imp(length(perc)+2:2*length(perc)+1),'v','filled','SizeData',60) %result of PC + MP selections 
% scatter((2:length(perc)+1), med_imp(2*length(perc)+3:end),'o','filled','SizeData',60) %result of PC + MP benchmark
scatter(length(perc)+2, med_imp(2*length(perc)+1),'sq','filled','SizeData',60) % MP_0
ylim([-25,15])
% ylim auto
xlim([0,length(perc)+3])
ax = gca;
ax.TickLabelInterpreter = 'latex'; 
xticklabels({'$0\%$','$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$100\%$'})
yline(0,'Color','r','LineWidth',1.5)
title(strcat(demand,' demand'))
legend('single PC',' PC + MP to random node sets (median)','PC + MP (max outflow)','PC + MP to all nodes','Fixed-Time Control','Location','northeast')
%legend('single PC',' PC + MP to random node sets (median)','PC + MP (proposed)','PC + MP (max outflow)','PC + MP to all nodes','Fixed-Time Control','Location','northeast')

saveas(gcf,strcat(pathFig,figfName,'.fig'));
saveas(gcf,strcat(pathFig,figfName,'.emf'));
print(gcf,strcat(pathFig,figfName),'-depsc','-painters')

%%
% figure(f1)
% scatter(1,med_TTT(length(perc)+1),'filled','SizeData',60) %results PC solo
% scatter((2:length(perc)+1), med_TTT(1:length(perc)),'*','SizeData',60)
% scatter((2:length(perc)+1), med_TTT(length(perc)+2:end-1),'filled','SizeData',60)
% scatter(length(perc)+2, med_TTT(end),'filled','SizeData',60)
% ylim([140000,270000]);
% xlim([0,length(perc)+3])
% yline(NCresults(5),'Color','r','LineWidth',1.5)
% ax = gca;
% ax.TickLabelInterpreter = 'latex'; 
% xticklabels({'$0\%$','$5\%$','$10\%$','$15\%$','$20\%$','$25\%$','$30\%$','$100\%$'})
% title(strcat(demand,' demand'))
% legend('PC solo',' PC + random MP node sets (median)','PC + selection MP nodes','PC + all MP nodes','No Control','Location','northeast')
