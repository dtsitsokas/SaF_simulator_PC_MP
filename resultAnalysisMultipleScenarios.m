% Result analysis of any experiment (one replication) / Multiple scenarios  
% Print Figures  preparation 
clear 
clc
close all 
load scenarios.mat
cases = [1 2 4 6 8]; 
%% Set the path where the result file is --- 
path_g = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\fourth set results\';

% store the path of the scenario directory to save the reasulting figures (should exist already)
path = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\fourth set results\temporaryResults\';


%% Print table of resutls for different scenarios 

%load(strcat(path_g,'output_set4_',scenarios{1},'.mat'),'indata') 
load(strcat(path_g,'output_set4_noTRU_',scenarios{1},'.mat'),'indata') 

results = zeros(10+3*indata.no_reg,1);
t = 180; %time steps (sec) - time window for aggregation for speeds 
SMS = ones(ceil(indata.kmax/t),length(cases))*indata.v_ff/1000; %[km/h]

for i = 1:length(cases) 
    
        %load results of the specified scenario
%         load(strcat(path_g,'output_set4_',scenarios{cases(i)},'.mat'),'outdata','indata')
        load(strcat(path_g,'output_set4_noTRU_',scenarios{cases(i)},'.mat'),'outdata','indata')
        
        % Table of Results
        %Total time spent in network and VQs -
        %PHT_ntw
        results(1,i) = sum(sum(outdata.x*indata.DT));

        %PHT_VQ
        results(2,i) = sum(sum(outdata.w(indata.group2,:)*indata.DT));

        %Total travel time
        %PHT_tot
        results(3,i) = results(1,i) + results(2,i);

        %Vehicles remaining inside 
        results(4,i) = outdata.notserviced; 

        %Balanced PHT_tot (with penalty)
        results(5,i) = results(3,i) + results(4,i)*0.25; 
        
        %Total waiting time in queues (network) 
        results(6,i) = sum(sum(outdata.w(indata.group1,:)*indata.DT)); %(pas x hours)
        
        %Mean queue over time and space 
        results(7,i) = mean(mean(outdata.x(indata.group1,:)));
        
        %Mean of (sum of Virtual Queues) over time
        results(8,i) = mean(sum(outdata.virtualqueues));
        
        %Mean total link outflow (links and VQs) 
        results(9,i) = mean(mean(outdata.u([indata.group1; indata.group2],:)));
        
        %Mean SMS - only inside the physical network (excluding VQs) ------        
        j=1;
        for k=t:t:indata.kmax
            non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0;
            s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
            s(non_empty) = min(s(non_empty), (sum(outdata.u(non_empty,k-t+1:k)*...
                indata.DT,2).*indata.Links2(non_empty,3))./...
                (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
            SMS(j,i) = mean(s);
            j = j + 1;
        end
        
        %Mean SMS - considering VQ time 
        %         j=1;
        %         for k=t:t:indata.kmax
        %            
        %             VKT = sum(outdata.u([indata.group1; indata.group2],k-t+1:k)*indata.DT,2).*indata.Links2([indata.group1; indata.group2],3)/1000; 
        %             VHT = [sum(outdata.x(indata.group1,(k-t+1):k)*indata.DT,2); sum(outdata.w(indata.group2,(k-t+1):k)*indata.DT,2)];
        %             SMS(j,i)= sum(VKT) / sum(VHT);
        %             j = j + 1;
        %         end

        results(10,i) = mean(SMS(i,:));
        
        %VHT_per_region(3 regions)
        %network
        for j=1:indata.no_reg
            results(10+j,i) = sum(sum(outdata.x(intersect(indata.group1,indata.reInd{j}),:)*indata.ksi*indata.DT));
        end
        
        %VQs 
        for j=1:indata.no_reg
            results(10+indata.no_reg+j,i) = sum(sum(outdata.w(intersect(indata.group2,indata.reInd{j}),:)*indata.ksi*indata.DT));
        end
        
        %VHT total per region 
        for j=1:indata.no_reg
            results(10+2*indata.no_reg+j,i) = results(10+j,i) + results(10+indata.no_reg+j,i);
        end
end
disp('total VHT') 
results(3,:)
disp('Total VHT (adjusted):')
results(5,:)

%% Table tags 

% tableTag{1} = 'PHT network';
% tableTag{2} = 'PHT VQ';
% tableTag{3} = 'PHT total';
% tableTag{4} = 'Veh inside';
% tableTag{5} = 'PHT tot reconciled';
% tableTag{6} = 'Total waiting time in queues';
% tableTag{7} = 'Mean link occupancy';
% tableTag{8} = 'Mean virtual queue';
% tableTag{9} = 'Mean link outflow';
% tableTag{10} = 'Mean SMS (km/h)';

%% Print results (figures and table)

%Create figures
%Create all figures first:
% f1 = figure('Name','MFD n-q line');
% title('MFD: accumulation - flow');
% xlabel('accumulation (veh)');
% ylabel('flow (veh/h)');
% box on
% hold on

f2 = figure('Name','MFD n-p line');
title('MFD: accumulation - production');
xlabel('accumulation (veh)');
ylabel('production (veh km/h)');
box on
hold on

f2c = figure('Name','TS production');
xlabel('time (h)');
ylabel('production (veh km/h)');
box on
hold on

% 
% f2b = figure('Name','MFD n-o');
% title('MFD: accumulation - trip completion');
% xlabel('$n$ (veh)');
% ylabel('trip completion rate (veh/h)');
% grid on
% hold on
% box on

% f3 = figure('Name','MFDs region 1');
% title(strcat('MFD: region 1'));
% xlabel('accumulation (veh)');
% ylabel('flow (veh/h)');
% box on
% hold on
% %legend('region 1','region 2','region 3')
% 
% f3b = figure('Name','MFDs region 2');
% title(strcat('MFD: region 2'));
% xlabel('accumulation (veh)');
% ylabel('flow (veh/h)');
% box on
% hold on
% 
% f3c = figure('Name','MFDs region 3');
% title(strcat('MFD: region 3'));
% xlabel('n (veh)');
% ylabel('q (veh/h)');
% box on
% hold on

% f4 = figure('Name','Prod MFD per region');
% title('MFD: accumulation - total production');
% xlabel('n (veh)');
% ylabel('p (veh km/h)');
% %legend('region 1','region 2','region 3')
% hold on
% box on;

f5 = figure('Name','Accumulation n over time');
%title('Total accumulation over time')
box on
xlabel('time (h)')
ylabel('accumulation (veh)')
hold on

% f6 = figure('Name','Node heterogeneity over time');
% box on
% xlabel('time (h)')
% ylabel('mean std of node queues (veh)')
% hold on

f7 = figure('Name','Cumulative trip endings over time');
box on
xlabel('time (h)')
ylabel('cumulative trip endings (veh)')
hold on

% f8 = figure('Name','Network Occupancy - Time');
% title('Network occupancy');
% xlabel('time (h)');
% ylabel('occupancy');
% hold on
% box on
%

f9 = figure('Name','Congested Lane-kilometers - Time');
%title('Number of congested links');
xlabel('time (h)');
ylabel('congestged links (lane $\times$ km)');
hold on
box on

% f10 = figure('Name','Total flow - Time');
% %title('Total flow over time');
% xlabel('time (h)');
% ylabel('total flow (veh/h)');
% hold on
% box on

% f11 = figure('Name','Queue std - Time');
% %title('Network heterogeneity');
% xlabel('time (h)');
% ylabel('std of queues (veh)');
% hold on
% box on

% f12 = figure('Name','SMS entire ntw - time');
% %title('Space mean speed over time');
% xlabel('time (h)');
% ylabel('SMS (km/h)');
% hold on
% box on

f13 = figure('Name','VQ over time');
%title('Mean VQ over time');
xlabel('time (h)');
ylabel('total VQ (veh)');
hold on
box on

% f14 = figure('Name','Region 1 heterogeneity vs n_i');
% title('Heterogeneity region 1')
% xlabel('$n_i$ (veh)')
% ylabel('stdev of $n_i/c_i$')
% hold on
% box on 
% 
% f14b = figure('Name','Region 2 heterogeneity vs n_i');
% title('Heterogeneity region 2')
% xlabel('$n_i$ (veh)')
% ylabel('stdev of $n_i/c_i$')
% hold on
% box on 
% 
% f14c = figure('Name','Region 3 heterogeneity vs n_i');
% title('Heterogeneity region 3')
% xlabel('$n_i$ (veh)')
% ylabel('stdev of $n_i/c_i$')
% hold on
% box on 

f15 = figure('Name','Time-Series: lanes of Spill-backs');
xlabel('time (h)')
ylabel('spill-backs (lanes)')
hold on
box on 

colReg{1} = '#0072BD';
colReg{2} = '#D95319';
colReg{3} = '#77AC30';
colReg{4} = '#4DBEEE';
lineTypeScen{1} = ':';
lineTypeScen{2} = '-';
lineTypeScen{3} = '-.';
lineTypeScen{4} = '--';
lineTypeScen{5} = '--';

% f16 = figure('Name','Accumulation TS per region');
% hold on 
% xlabel('time (sec)')
% ylabel('accumulation (veh)')
% 
% f17 = figure('Name','Transfer Flows TS');
% hold on 
% xlabel('time')
% ylabel('flow (veh/h)')


% Time - Demand Multiplier
% figure('Name','Demand Multiplier - Time SaF_2')
% for k=1:indata.kmax
%     tot_demand(k) = indata.defac(k)*sum(indata.demandGroup2(:,1+(k>indata.WU)));
% end
% plot((1:indata.kmax)*indata.DT,tot_demand,'LineWidth',2)
% %title('Demand pattern');
% xlabel('time (h)');
% %ylim([0 1]);
% ylabel('total demand (veh/h)');
% % saveas(gcf,'demandPattern.fig');
% % saveas(gcf,'demandPattern.emf');
% % print('demandPattern.eps','-depsc','-painters')
%
for scen = 1:length(cases)
    
    %load results of the specified scenario
    load(strcat(path_g,'output_set4_',scenarios{cases(scen)},'.mat'),'outdata','indata')
  
    
    % Average Occupancy of every link over the total simulation time
    % oslink_avg_occupancy=mean(outdata.x(indata.group1,:),2)./indata.capacity(indata.group1);
    %save('link_avg_occupancy_1.mat','link_avg_occupancy');
    %links_axis=1:length(indata.group1);
    
    %         figure('Name','Histogram: Average link occupancy')
    %         histogram(link_avg_occupancy,10,'Normalization','probability')
    %         ylabel('probability');
    %         xlabel('average link occupancy')
    
    %         figure('Name','Bar chart: Average Link Occupancy')
    %         bar(links_axis,link_avg_occupancy)
    %         xlabel('link Index');
    %         ylabel('average occupancy');
    %         %saveas(gcf,'Res_AverageLinksOcc_S.fig');
    
    % VQs over time (one line per link in common graph)
%     figure('Name','Virtual Queues - Time SaF_2');
%     hold on
%     title(strcat('Virtual Queues - ',scenarios{scen}));
%     xlabel('Time [h]');
%     ylabel('Queue [veh]');
%     plot((1:size(outdata.virtualqueues(1,1:end-1),2))*indata.DT,outdata.virtualqueues(:,1:end-1));
%     %saveas(gcf,strcat(path,'TS_VirtualQsALL_',scenarios{scen},'.fig'));
%     %saveas(gcf,strcat(path,'TS_VirtualQsALL_',scenarios{scen},'.emf'));
    
    % MFDs: Entire network -
    t_MFD = 3/60; %(hours) time interval between measurements
    t = ceil(t_MFD/indata.DT); %The step / How many steps k*DT are included in the aggregation period
    
    sum1 = sum(outdata.x(indata.group1,1:indata.kmax),1); %sum1 = total network accumulation per step k [veh] (all original network links (group1))
    sum2 = sum(outdata.u,1); %sum2 = total outflow of all links and virtual queues in the network in every step k [veh/DT]
    sum4 = sum(outdata.a_z,1); %sum4 = total arriving flow at the end of the queue
    %sum5 = sum((outdata.m(indata.group1,1:indata.kmax)*indata.v_ff)./(indata.Links2(indata.group1,3)*ones(1,indata.kmax)),1); %sum5 = total moving flow in all network links
    %sum5b = sum(outdata.u(indata.group1,1:indata.kmax).*(indata.Links2(indata.group1,3)*ones(1,indata.kmax))./sum(indata.Links2(indata.group1,3)),1);
    sum3 = sum(outdata.u(indata.group2,:),1); %sum3 = total outflow of virtual queues (real inflow) in every step k [veh/DT]
    sum6 = sum(outdata.w(indata.group2,1:indata.kmax),1); %virtual queues 
    sum7 = sum(outdata.u(indata.group1,:).*(indata.Links2(indata.group1,3)/1000),1); % production = outflow of network links * link length
    %sum7 = sum(outdata.m(indata.group1,1:end-1)*indata.v_ff/1000,1); %production = moving*free_flow_speed
    sum8 = sum(outdata.q(indata.group3,:),1); %sum5 = Total Trips completed (inflow of group 3 links) every time step k [veh]
    %     tripCompletionRate=[0];
    production=0;
    accumulation=0;
    outflow=0;
    k_MFD=0;
    vqueues = 0; 
    tripCompletionRate = 0;
    for i=t:t:indata.kmax
        k_MFD = [k_MFD i*indata.DT]; %time for the MFD [hours]
        %accumulation = [accumulation mean(sum1(i-t+1:i)+sum6(i-t+1:i))]; %[veh] - mean accumulation over the time window (VQ included)
        accumulation = [accumulation mean(sum1(i-t+1:i))];
        vqueues = [vqueues mean(sum6(i-t+1:i))];
        %         outflow = [outflow mean(sum2(i-t+1:i))+mean(sum4(i:i+(t-1)))+mean(sum5(i-t+1:i))]; %[veh/h]
        outflow = [outflow mean(sum2(i-t+1:i))]; %[veh/h]
        % outflow = [outflow mean(sum2(i-t+1:i))]; %[veh/h]

        % outflow = [outflow mean(sum5(i-t+1:i))]; %[veh/h]
        %          outflow = [outflow mean(sum5b(i:i+(t-1)))]; %[veh/h] %only link outflow weighted with link lengths (ignores VQs and other flows)
         production = [production mean(sum7(i-t+1:i))]; %[veh*km/h]
         tripCompletionRate=[tripCompletionRate mean(sum8(i-t+1:i))]; %[veh/h]
    end
    
    % MFD Accumulation - Outflow
    %figure(f1)
    %plot(accumulation,outflow,'Linewidth',1.5)
    
%     legend(scenarios{1},scenarios{2},scenarios{3},scenarios{4})
%     legend('FTC','MP-0','MP-S1','MP-S2')
       
%       labels = strcat(labels, scenarios{scen}); 
%       legend(labels)
    %saveas(gcf,strcat(path,'MFD_n_q_ntw.fig'));
    %saveas(gcf,strcat(path,'MFD_n_q_ntw.emf'));
    %print(strcat(path,'MFD_n_q_ntw.eps'),'-depsc','-painters')
    
    
    % MFD: Accumulation [veh] - Total Travel Production [veh*km/h] (entire network)
        
    %production=production*2000/36;
    figure(f2)
    plot(accumulation,production,'Linewidth',1.5)
    
    figure(f2c)
    plot(k_MFD,production,'Linewidth',1.5)
    mean(production)
%    scatter(accumulation,production,'filled')
%     legend('FTC','MP-0','MP-S1','MP-S2')
%     legend(scenarios{1},scenarios{2},scenarios{3})
    %saveas(gcf,strcat(path,'MFD_n_p_ntw.fig'));
    %saveas(gcf,strcat(path,'MFD_n_p_ntw.emf'));
    %print(strcat(path,'MFD_n_p_ntw.eps'),'-depsc','-painters')


    % MFD: Accumulation [veh] - Trip Completion Rate [veh/h]
    %figure(f2b)
    %plot(accumulation,tripCompletionRate,'Linewidth',1.5)
    %saveas(gcf,'Res_MFD_Acc-TripCompletion_S.fig');
    
     % Flow MFDs per region - common graph for multiple runs
% 
%     for i=1:3
% %         Rsum1(i,:) = sum(outdata.x(intersect(indata.group1,indata.reInd{i}),1:indata.kmax),1) + sum(outdata.w(intersect(indata.group2,indata.reInd{i}),1:indata.kmax),1);  %region accumulation + VQs
%         Rsum1(i,:) = sum(outdata.x(intersect(indata.group1,indata.reInd{i}),1:indata.kmax));  %region accumulation only
%         Rsum2(i,:) = sum(outdata.u(indata.reInd{i},1:indata.kmax),1); %outflows of all links and virtual queues of the region
%         Rsum4(i,:) = sum(outdata.a_z(intersect(indata.reInd{i},indata.group1),1:indata.kmax),1); %total arriving flow at the end of the queue per region
%         Rsum5(i,:) = sum((outdata.m(intersect(indata.reInd{i},indata.group1),1:indata.kmax)*indata.v_ff)./(indata.Links2(intersect(indata.reInd{i},indata.group1),3)*ones(1,indata.kmax)),1); %total moving flow in all network links
%         Rsum3(i,:) = sum(outdata.u(intersect(indata.group2,indata.reInd{i}),1:indata.kmax),1); %total inflow of region
%         Rsum6(i,:) = sum(outdata.m(intersect(indata.group1,indata.reInd{i}),1:indata.kmax)*indata.v_ff/1000,1); %production = moving*free_flow_speed
%         %Rsum6(i,:) = sum(outdata.u(intersect(group1,reInd{i}),:).*(Links2(intersect(group1,reInd{i}),3)/1000),1); %production = outflow * link length
%         accumulation = 0;
%         outflow = 0;
%         k_MFD = 0;
%         for j=1:t:indata.kmax
%             k_MFD = [k_MFD j*indata.DT]; %time for the MFD [hours]
%             accumulation = [accumulation mean(Rsum1(i,j:j+(t-1)))]; %[veh] - mean accumulation over the time window
%             %outflow=[outflow sum(sum2(i:i+(t-1)))/t_MFD+sum(sum3(i:i+(t-1)))/t_MFD]; %[veh/h]
%             %                 outflow=[outflow mean(Rsum2(i,j:j+(t-1)))+mean(Rsum4(i,j:j+(t-1)))+mean(Rsum5(i,j:j+(t-1)))]; %[veh/h]
%             outflow=[outflow mean(Rsum5(i,j:j+(t-1)))]; %[veh/h]
%         end
% %         if i==1
% %             figure(f3)
% %         elseif i==2
% %             figure(f3b)
% %         else
% %             figure(f3c)
% %         end
% %         plot(accumulation,outflow,'Linewidth',1.5,'Color',colReg{i},'LineStyle',lineTypeScen{scen},'LineWidth', 2)
% %         scatter(accumulation,outflow,'filled','MarkerFaceColor',colReg{i})
%         
%         
%         % Heterogeneity Figures (per region): accumulation vs. stdev of queues x
%         %/check time - it takes too long / separate regions (?) or scenarios -
%         %too many lines
%                 
% %         ind = intersect(indata.group1,indata.reInd{i});
% %         occ = outdata.x(ind,1:indata.kmax);
% %         cap = zeros(length(ind),indata.kmax);
% %         for j = 1:indata.kmax
% %             cap(:,j) = indata.capacity(ind);
% %         end
% %         occ = occ./cap;
% %         Rstdev(i,:) = std(occ);
% %         stdev_x = 0;
% %         accumulation = 0;
% %         for j=1:t:indata.kmax
% %             accumulation = [accumulation mean(Rsum1(i,max(1,j):j+(t-1)))];
% %             stdev_x=[stdev_x mean(Rstdev(i,j:j+(t-1)))]; %[veh/h]
% %         end
% %         if i==1
% %             figure(f14)
% %         elseif i==2
% %             figure(f14b)
% %         else
% %             figure(f14c)
% %         end
% %         plot(accumulation, stdev_x,'LineWidth',1.5,'Color',colReg{i},'LineStyle',lineTypeScen{scen})
%         
%         
%     end

    %---- make the legend separately by running the command:
%    figure(f3)
%    legend('FTC','MP-0','MP-S1','MP-S2')
%     saveas(gcf,strcat(path,'MFD_n_q_reg_1.fig'));
%     saveas(gcf,strcat(path,'MFD_n_q_reg_1.emf'));
%     print(strcat(path,'MFD_n_q_reg_1.eps'),'-depsc','-painters')
 %   figure(f3b)
 %   legend('FTC','MP-0','MP-S1','MP-S2')
%     saveas(gcf,strcat(path,'MFD_n_q_reg_2.fig'));
%     saveas(gcf,strcat(path,'MFD_n_q_reg_2.emf'));
%     print(strcat(path,'MFD_n_q_reg_2.eps'),'-depsc','-painters')
   % figure(f3c)
  %  legend('FTC','MP-0','MP-S1','MP-S2')
%     saveas(gcf,strcat(path,'MFD_n_q_reg_3.fig'));
%     saveas(gcf,strcat(path,'MFD_n_q_reg_3.emf'));
%     print(strcat(path,'MFD_n_q_reg_3.eps'),'-depsc','-painters')
%     figure(f4)
%     % Production MFDs
%     for j=1:3
%         Rsum1(j,:) = sum(outdata.x(intersect(group1,reInd{j}),1:indata.kmax),1);
%         production=0;
%         accumulation = 0;
%         for i=1:t:indata.kmax
%             accumulation = [accumulation mean(Rsum1(j,i:i+(t-1)))]; %[veh] - mean accumulation over the time window
%             production=[production mean(Rsum6(j,i:i+(t-1)))]; %[veh*km/h]
%         end
%         plot(accumulation,production,'Linewidth',1.5,'Color',colReg{j},'LineStyle',lineTypeScen{scen})
%         %         scatter(accumulation, production,'filled')
%     end
    %---- make the legend separately by running the command:
    %figure(f4)
    %legend('NC-1','NC-2','NC-3','MP1-1','MP1-2','MP1-3')
    %legend off
    % -------
     
    % Time-Series: Total Network accumulation (+VQs included)
    figure(f5)
    plot(k_MFD, accumulation,'LineWidth',1.5)
    
    figure(f13)
    plot(k_MFD, vqueues,'LineWidth',1.5)
    
    %legend('network','virtual queues')
%      legend('FTC','MP-0','MP-S1','MP-S2')
%     legend(scenarios{1},scenarios{2},scenarios{3})
    %saveas(gcf,strcat(path,'TS_n_ntw.fig'));
    %saveas(gcf,strcat(path,'TS_n_ntw.emf'));
    %print(strcat(path,'TS_n_ntw.eps'),'-depsc','-painters')
    
    % Time-series: Number of lane*kms (% of the total lane-km of the network) that are over 80% capacity
    
    figure(f9)
    congested = outdata.x(:,1:end-1)>0.8*indata.capacity;
    congested = congested.*(indata.Links2(:,2).*indata.Links2(:,3)); 
    congested = sum(congested);
    m_congested = 0;
    %aggregate every 180 sec 
    for i=t:t:indata.kmax
        m_congested = [m_congested mean(congested(i-180+1:i))];
    end
%     plot(time, congested,'LineWidth',1.5)
    plot(k_MFD, m_congested,'LineWidth',1.5)
%     legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_congestedLinks.fig'));
    %saveas(gcf,strcat(path,'TS_congestedLinks.emf'));
    %print(strcat(path,'TS_congestedLinks.eps'),'-depsc','-painters')
    
    % Time-Series: Space mean speed (VKT/VHT)
%     figure(f12)
%     SMS = ones(1,ceil(indata.kmax/t))*indata.v_ff/1000;
%     j=1;
%     for k=t:t:indata.kmax
%        non_empty = sum(outdata.x(indata.group1,(k-t+1):k),2)>0; 
%        s = ones(length(indata.group1),1)*indata.v_ff/1000; %km/h
%        s(non_empty) = min(s(non_empty), sum(outdata.u(non_empty,k-t+1:k)*indata.DT,2).*indata.Links2(non_empty,3)./...
%            (1000*sum(outdata.x(non_empty,(k-t+1):k)*indata.DT,2))); %km/h
%        
%        SMS(j) = mean(s);
%        j = j + 1; 
%     end
%     results(10,scen) = mean(SMS); %mean SMS over time for the simulation
%     plot(k_MFD(2:end),SMS(:,scen),'LineWidth',1.5); 
%      legend('FTC','MP-0','MP-S1','MP-S2')
%     saveas(gcf,strcat(path,'TS_SMS.fig'));
%     saveas(gcf,strcat(path,'TS_SMS.emf'));
%     print(strcat(path,'TS_SMS.eps'),'-depsc','-painters')
    
    % Time-Series: Global network heterogeneity (std of all queues every time step)
    % entire queues or waiting queues 
%     figure(f11)
%     ntw_het = std(outdata.x(indata.group1,1:end-1));
%     het = 0; 
%     for i=t:t:indata.kmax
%         het = [het mean(ntw_het(i-t+1:t))];    
%     end
%     plot(k_MFD, het,'LineWidth',1.5)

%     legend('FTC','MP-0','MP-S1','MP-S2')
%     saveas(gcf,strcat(path,'TS_stdQueues_ntw.fig'));
%     saveas(gcf,strcat(path,'TS_stdQueues_ntw.emf'));
%     print(strcat(path,'TS_stdQueues_ntw.eps'),'-depsc','-painters')
%     
    % Time-Series: Total-outflow over time (outflows of network links and vqs)
    %figure(f10)
    %plot(k_MFD, outflow,'LineWidth',1.5)
%     legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_flow_ntw.fig'));
    %saveas(gcf,strcat(path,'TS_flow_ntw.emf'));
    %print(strcat(path,'TS_flow_ntw.eps'),'-depsc','-painters')    
%     
    
    % Time-Series: Mean/Sum VQ over time (real per time step / no aggregation)
    
    %figure(f13)
    %plot(time, mean(outdata.virtualqueues(:,1:end-1)),'LineWidth',1.5)
    %plot(time, sum(outdata.virtualqueues(:,1:end-1)),'LineWidth',1.5)
    %std VQ:
    % plot(time, std(outdata.virtualqueues(:,1:end-1)))
    %legend('mean VQ','sum VQ','st.dev VQ')
    %legend('mean VQ','st.dev VQ')  
%      legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_VQs_ntw.fig'));
    %saveas(gcf,strcat(path,'TS_VQs_ntw.emf'));
    %print(strcat(path,'TS_VQs_ntw.eps'),'-depsc','-painters')    
    
    % Time-series: Cumulative trips completed
    cumtripsAll = sum(outdata.cumtripsCompleted,1);
    p = cumtripsAll((t:t:indata.kmax));
    p = [0 p]; 
    figure(f7)
    plot(k_MFD, p,'LineWidth',1.5)  %all points
%     legend(scenarios{1},scenarios{2},scenarios{3})
%      legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_cumTripEndings_ntw.fig'));
    %saveas(gcf,strcat(path,'TS_cumTripEndings_ntw.emf'));
    %print(strcat(path,'TS_cumTripEndings_ntw.eps','-depsc','-painters')) 
%     
    % Time-series: node heterogeneity
    % Calculate hererogeneity of node first, then find mean of all nodes 
%     figure(f6)
%     tot_node_het = 0;
%     for k = t:t:indata.kmax
%         node_heterog = 0;
%         for i=1:size(indata.MP.nodeID,2)
%             if indata.MP.approaches{i}>0
%                 
%                 %Mean over time of the std of incoming queues around the
%                 %intersection every cycle
%                 node_heterog = [node_heterog mean(std(outdata.w(indata.junct.or_index(indata.MP.approaches{i}),k-t+1:k)))];
%             end
%         end
%         tot_node_het = [tot_node_het mean(node_heterog)];
%     end 
%     plot(k_MFD,tot_node_het,'LineWidth',1.5)
    
%      legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_std_nodes.fig'));
    %saveas(gcf,strcat(path,'TS_std_nodes.emf'));
    %print(strcat(path,'TS_std_nodes.eps','-depsc','-painters'))

    
    % Total Occupancy of the network 
%     figure(f8)
%     tot_occupancy = sum1/sum(indata.capacity(indata.group1));
%     time_int = (1:indata.kmax)*indata.DT;
%     plot(time_int,tot_occupancy(1:length(time_int)),'Linewidth',1.5)
%      legend('FTC','MP-0','MP-S1','MP-S2')
    %saveas(gcf,strcat(path,'TS_occupancy_ntw.fig'));
    %saveas(gcf,strcat(path,'TS_occupancy_ntw.emf'));
    %print(strcat(path,'TS_occupancy_ntw.eps'),'-depsc','-painters')     
    
  
   %     Movie - Showcase the congestion dynamic
    %     step = 5*36;% Time interval between snapshots for the video in # of ks - k = DT sec / time step = step * DT [sec]
    %     v=createmovie(outdata.x,indata.kmax,indata.capacity,Links,Nodes,DT,step,'S');
    %     Save some results
    %     accumulation=accumulation';
    %     outflow=outflow';
    %     save outflowS.mat outflow
    %     save accumulationS.mat accumulation
    
    % TIme-Series - Accumulation  
%     figure(f16)
%     for r=1:indata.no_reg
%        
%         %network links (group1) - x
%         sum_1 = sum(outdata.x(intersect(indata.group1, indata.reInd{r}),:));
%         %virtual queues (group2) - w
%         sum_2 = sum(outdata.w(intersect(indata.group2, indata.reInd{r}),:));
%         sum_r = sum_1 + sum_2;
%         
%         %perform aggregation of values
%         accumulation=0;
%         k_MFD=0;
%         for i=1:t:indata.kmax
%             k_MFD = [k_MFD i*indata.DT]; %time for the MFD [hours]
%             accumulation = [accumulation mean(sum_r(i:i+(t-1)))]; %[veh] - mean accumulation over the time window (VQ included)
%         end
%         if scen==1
%             plot(k_MFD, accumulation,'LineWidth',1.5,'Color',colReg{r}, 'LineStyle','--')
%         else
%             plot(k_MFD, accumulation,'LineWidth',1.5,'Color',colReg{r})
%         end
%     end
    
    

   % Time-Series of number of (lane)xkm of spillbacks (>85% of capacity)
    figure(f15)
    spillback_lanes = outdata.x(indata.group1,:)./indata.capacity(indata.group1);
    spillback_lanes = sum((spillback_lanes>=0.85).*indata.Links2(indata.group1,2));
    spillbacks = 0;
    for i=t:t:indata.kmax
        spillbacks = [spillbacks mean(spillback_lanes(i-t+1:i))]; %[veh] - mean accumulation over the time window (VQ included)
    end
    plot(k_MFD,spillbacks,'LineWidth',1.5)
    
%     %% Transfer flows between adjacent regions over time 
%     figure(f17) 
%     %indices of Links of every approach 
%     adjRg = [1 2;3 2;2 1;2 3]; 
%     for j=1:4
%         ind1 = ismember(indata.junct2.or_index, indata.reInd{adjRg(j,1)});
%         ind2 = ismember(indata.junct2.dest_index, indata.reInd{adjRg(j,2)});
%         ind = ind1+ind2 == 2; 
%         reg_outflow = sum(outdata.u(indata.junct2.or_index(ind),:));
%         regOut = 0;
%         for i=1:t:indata.kmax
%             regOut = [regOut mean(reg_outflow(i:i+(t-1)))]; %[veh] - mean accumulation over the time window (VQ included)
%         end
%         if scen==1
%             plot(k_MFD,regOut,'LineWidth',1.5,'Color',colReg{j},'LineStyle','--')
%         else
%             plot(k_MFD,regOut,'LineWidth',1.5,'Color',colReg{j})
%         end
%     end
%     legend('1-2 NC','3-2 NC','2-1 NC','2-3 NC','1-2','3-2','2-1','2-3')
end

%% Legends and save figures 



figure(f2)
legend('FTC','100\% MP nodes','10\% MP nodes','20\% MP nodes')
saveas(gcf,strcat(path,'MFD_n_p_ntw.fig'));
print(strcat(path,'MFD_n_p_ntw'),'-depsc','-painters')

figure(f2c)
legend('FTC','100\% MP nodes','10\% MP nodes','20\% MP nodes')
saveas(gcf,strcat(path,'TS_production.fig'));
print(strcat(path,'TS_production'),'-depsc','-painters')


figure(f5)
legend('FTC','100\% MP nodes','10\% MP nodes','20\% MP nodes')
saveas(gcf,strcat(path,'TS_accumulation.fig'));
print(strcat(path,'TS_accumulation'),'-depsc','-painters')

figure(f13)
legend('FTC','100\% MP nodes','10\% MP nodes','20\% MP nodes')
saveas(gcf,strcat(path,'TS_totVQs.fig'));
print(strcat(path,'TS_totVQs'),'-depsc','-painters')


figure(f7) 
legend('FTC','100\% MP nodes','10\% MP nodes','20\% MP nodes')
saveas(gcf,strcat(path,'cum_tripEndings.fig'));
print(strcat(path,'cum_tripEndings'),'-depsc','-painters')





%% Additional figures (maps) - Control related  



% Plot: colormap-> mean density per link over time (MP wrt NC) to see which links took
% more traffic 



%% Print map of the  Barcelona network
%clear all
%close all 
load LINKS %links info [id, start_node, end_node, length, cluster (in 4 reg)]
load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids) 
load clus_final %clustering for links
load control %controlled intersections 

%print map with 3 clusters in different colors 
nodes = nodes3; 
figure('Name','Node selection C1')
hold on;
for i=1:size(clus_final,1)
    ind=clus_final(i,1); 
    gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off


%% Print specified set of nodes with specific color (all same length) 

%Nodes of Selection C1 (Max-Pressure) 
load('nodesC1n.mat')
hold on
for i=1:length(nodesC)
    plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#7E2F8E','Markersize',35)
end
title('MP Nodes - C1')
saveas(gcf,strcat(path,'mapWithNodesC1n.fig'));
print(strcat(path,'mapWithNodesC1n'),'-depsc','-painters')

%% Map with total link outflow (links serving the highest demand overall)

links_u = sum(outdata.u,2);
figure('Name','Total Outflow per link');
hold on;
for i=1:size(LINKS,1)
    ind = links_u(i)/max(links_u); 
    gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
colorbar
axis off
title('Total link outflow')
saveas(gcf,strcat(path,'mapOutflowPerLink.fig'));
% saveas(gcf,'mapOutflowPerLink.emf');
print(strcat(path,'mapOutflowPerLink'),'-dmeta','-painters')

%% Total VHT in virtual queues as circles on the start nodes
% Plus total VHT spent in links 

VHT_tot = sum(outdata.x,2)*indata.DT;
VHT_tot(indata.group2) = sum(outdata.w(indata.group2,:),2)*indata.DT;
figure('Name','Map total VQ-VHT per centroid & VHT per link');
hold on;
for i=1:size(clus_final,1)
    ind = VHT_tot(i)/max(VHT_tot(1:size(clus_final,1)));
    gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off

VHT_VQs = sum(outdata.virtualqueues,2);
%figure('Name','Map total VQ-VHT per centroid');
% hold on;
% for i=1:size(clus_final,1)
%     ind=clus_final(i,1); 
%     gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% axis off
%Print specified set of nodes with specific color (all same length)  
nodesVQs = indata.Links2(indata.group2,5);
myColorMap = jet(256);
colormap(myColorMap)
for i=1:length(nodesVQs)
    plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),...
        'Marker','.','Color',myColorMap(max(1,ceil(VHT_VQs(i)/max(VHT_VQs)*256)),:),...
        'Markersize',30)
%     plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i)...
%           ,3),'Marker','.','Color','#A2142F','Markersize',30)
end
title('VHT spent in VQs')
colorbar
saveas(gcf,strcat(path,'mapVHT_VQ_PerLink.fig'));
% saveas(gcf,'mapVHT_VQ_PerLink.emf');
print(strcat(path,'mapVHT_VQ_PerLink'),'-depsc','-painters')


%% Total outflow of virtual queues as circles on the start nodes

tot_VQs = sum(outdata.u(indata.group2,:),2);
figure('Name','Map total VQ outflow per centroid');
hold on;
for i=1:size(clus_final,1)
    ind=clus_final(i,1); 
    gplotdc_([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off

%Print specified set of nodes with specific color (all same length)  
nodesVQs = indata.Links2(indata.group2,5);
myColorMap = jet(256);
colormap(myColorMap)
for i=1:length(nodesVQs)
    plot(nodes(nodes(:,1)==nodesVQs(i),2),nodes(nodes(:,1)==nodesVQs(i),3),'Marker','.','Color',myColorMap(max(1,ceil(tot_VQs(i)/max(tot_VQs)*256)),:),'Markersize',30)
end
title('Total VQ outflow')
colorbar
saveas(gcf,strcat(path,'mapTotal_VQ_Nodes.fig'));
% saveas(gcf,'mapVHT_VQ_PerLink.emf');
print(strcat(path,'mapTotal_VQ_Nodes'),'-depsc','-painters')


%% Total time spent overall per link 

VHT_tot = sum(outdata.x,2)*indata.DT;
VHT_tot(indata.group2) = sum(outdata.w(indata.group2,:),2)*indata.DT;
figure('Name','Total VHT per link');
hold on;
for i=1:size(clus_final,1)
    ind = VHT_tot(i)/max(VHT_tot(1:size(clus_final,1)));
    gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off
colorbar
title('VHT per link')

saveas(gcf,strcat(path,'mapVHT_PerLink_1.fig'));
% saveas(gcf,strcat(path,'mapVHT_PerLink_1.emf'));
print(strcat(path,'mapVHT_PerLink_2'),'-depsc','-painters')

%% Highest average density/occupancy 
Avg_occ = mean(outdata.x,2)./indata.capacity;
figure('Name','Map Mean Occupancy per link');
hold on;
for i=1:size(clus_final,1)
    ind = Avg_occ(i)/max(Avg_occ);
    gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
axis off
colorbar
title('Mean occupancy per link')
saveas(gcf,strcat(path,'mapOccupancyPerLink.fig'));
% saveas(gcf,'mapOccupancyPerLink.emf');
print(strcat(path,'mapOccupancyPerLink'),'-depsc','-painters')


%% print specific snapshot
setTime = 6; 
createSnapshot(setTime,indata.DT,outdata.x,indata.capacity,'occupancy',LINKS,nodes3)

%peak occupancy(2h)
setTime = 2; 
createSnapshot(setTime,indata.DT,outdata.x,indata.capacity,'occupancy',LINKS,nodes3)


%% Nodes that we control with PC + total outflow per link 
% to check if nodes are well selected 

% load(strcat('input_','PC_0','.mat'),'PC')
% nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
% links_u = sum(outdata.u,2);
% figure('Name','Map PC Nodes & total link outflow');
% hold on;
% for i=1:size(LINKS,1)
%     ind = links_u(i)/max(links_u); 
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% colorbar
% axis off
% 
% for i=1:length(nodesC)
%     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',35)
% end
% title('PC nodes + total link outflow')
% saveas(gcf,strcat(path,'mapNodesPC_plus_LinkOutflow.fig'));
% % saveas(gcf,'mapNodesPC.emf');
% print(strcat(path,'mapNodesPC_plus_LinkOutflow'),'-depsc','-painters')
% % print('mapWithNodesC1n.eps','-depsc','-painters')

%% Mean link occupancy + PC controlled nodes 

% load(strcat('input_','PC_0','.mat'),'PC')
% nodesC = [unique(PC.nodes(PC.nodes>0)); unique(PC.nodesco(PC.nodesco>0))];
% Avg_occ = mean(outdata.x,2)./indata.capacity;
% figure('Name','Map PC Nodes & average link occupancy');
% hold on;
% for i=1:size(LINKS,1)
%     ind = Avg_occ(i)/max(Avg_occ); 
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% colorbar
% axis off
% %node printing 
% for i=1:length(nodesC)
%     plot(nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),2),nodes(nodes(:,1)==indata.MP.nodeID(nodesC(i)),3),'Marker','.','Color','#A2142F','Markersize',35)
% end
% title('PC nodes + average link occupancy')
% saveas(gcf,strcat(path,'mapNodesPC_plus_LinkOccupancy.fig'));
% % saveas(gcf,'mapNodesPC.emf');
% print(strcat(path,'mapNodesPC_plus_LinkOccupancy'),'-depsc','-painters')
% % print('mapWithNodesC1n.eps','-depsc','-painters')


%% Plot: colormap-> aggregated time during which every link overpasses 80%
% of capacity (spill-backs) - SAME figure as mean occupancy per link 

% spillback_count = outdata.x(1:indata.NLinks,:)>0.8*indata.capacity(1:indata.NLinks);
% spillback_count = sum(spillback_count,2);
% figure('Name','Spill-back time per link');
% hold on;
% for i=1:indata.NLinks
%     ind = spillback_count(i)/max(spillback_count);
%     gplotdc_colormap([0 1;0 0],nodes([find(nodes(:,1)==LINKS(i,2)),find(nodes(:,1)==LINKS(i,3))],[2 3]),ind);
%     hold on;
% end
% axis off
% colorbar
% title('Spill-back frequency per link')
% saveas(gcf,strcat(path,'mapSpillBacksPerLink.fig'));
% % saveas(gcf,'mapOccupancyPerLink.emf');
% print(strcat(path,'mapSpillBacksPerLink'),'-dmeta','-painters')

%% Results of individual run per region ------- 
%% Accumulation per region (including VQs) (common graph - single run only) 
% figure(f16) 
% for r=1:indata.no_reg
%     
%     %network links (group1) - x
%     sum_1 = sum(outdata.x(intersect(indata.group1, indata.reInd{r}),:));
%     %virtual queues (group2) - w 
%     sum_2 = sum(outdata.w(intersect(indata.group2, indata.reInd{r}),:));
%     sum_r = sum_1 + sum_2; 
%     
%     %perform aggregation of values 
%     accumulation=0;
%     k_MFD=0;
%     for i=1:t:indata.kmax
%         k_MFD = [k_MFD i*indata.DT]; %time for the MFD [hours]
%         accumulation = [accumulation mean(sum_r(i:i+(t-1)))]; %[veh] - mean accumulation over the time window (VQ included)
%     end
%     
%     plot(k_MFD, accumulation,'LineWidth',1.5)
% end
% legend('region 1','region 2','region 3')
% saveas(gcf,'AccumulationPerRegion.fig');
% % saveas(gcf,'mapOccupancyPerLink.emf');
% print('AccumulationPerRegion','-dmeta','-painters')

%% Print figures with results of specific node to analyze 

% Specify node 

% print TS of green of duration of the main and secondary phase  

% print TS of agg_u (calculated and applied) for the approach 

% Print comparative graph of accumulation in both regions 

% Print comparative graph of transfer flow 