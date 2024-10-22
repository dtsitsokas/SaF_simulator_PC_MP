% Result analysis of any experiment (single replication) -
% Print Figures  preparation

function [] = resultAnalysis(S1,demand)
% clear
% clc
% close all

%type here the names of the cases you need to plot/get results or load
%scenarios:
%S1 = 5; %set the # of the case in scenarios - comparison with the NC case 
load('scenarios.mat','scenarios')

scen_ID = [1 S1] ; %choose from the scenarios the cases you need

%% Set the path where the result file is ---

%path for finding the results:
path_g = '/Users/tsitokas/Documents/GitHub/SaF-LUTS/results/';
% fnameProthem = 'output_set5_MedDemand_';
fnameProthem = strcat('output_final_',demand,'_'); %prothem of the file name for this demand case
fnameNC = strcat('output_final_',demand,'_NC'); %file for loading the for the indata
%path to save figures
pathFigJ = '/Users/tsitokas/Documents/GitHub/SaF-LUTS/results/';

xx = 5; %window for moving average
lineW = 2; %line width for all line plots 


%% Print table of resutls for different scenarios

%load the NC scenario 
load(strcat(path_g,fnameNC),'indata') %medium demand 
%% 
results = zeros(10+3*indata.no_reg,1);
t = 180; %time steps (sec) - time window for aggregation for speeds
SMS = ones(ceil(indata.kmax/t),length(scen_ID))*indata.v_ff/1000; %[km/h]
for jj=1:length(scen_ID)
    i = scen_ID(jj);
    % load results of the specified scenario
    load(strcat(path_g,fnameProthem,scenarios{i},'.mat'),'outdata','indata')
%     load(strcat(fnameProthem,scenarios{i},'_test.mat'),'outdata','indata')

    % incremental selection files 
    % load(strcat(path_g,'output_MP_IS_',num2str(i),'_',num2str(i*5),'.mat'),'outdata','indata')

    %load(strcat(path_g,'output_j_t_',scenarios{i},'_subsets.mat'),'outdata','indata')
    
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
    results(10,i) = mean(SMS(:, i));
    
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

% Print results (figures and table)

% Initialize figures
% -----------------------------

%% Figure 1: Multi figure MFD and n, VQ time-series --------------------
ff1 = figure('Name','Multi MFD plus n, VQ time-series');
%title('MFD: accumulation - flow');
%first row: MFDs n-e 
for i=1:3
    f1(1,i) = subplot(4,4,i); 
    xlabel('$n$ (veh)');
    ylabel('$e$ (veh/h)');
    title(strcat('region $',num2str(i),'$'))
    box on
    hold on
end

%entire network 
i=4;
f1(1,i) = subplot(4,4,i);
xlabel('$n$ (veh)');
ylabel('$e$ (veh/h)');
title('network')
box on
hold on


%second row: MFDs n-p
for i=5:7
    f1(2,i-4) = subplot(4,4,i); 
    xlabel('$n$ (veh)');
    ylabel('$p$ (veh/h)');
    title(strcat('region $',num2str(i-4),'$'))
    box on
    hold on
end

%entire network 
i=8;
f1(2,i-4) = subplot(4,4,i);
xlabel('$n$ (veh)');
ylabel('$p$ (veh/h)');
title('network')
box on
hold on



%third row: TS n

for i=9:11
    f1(3,i-8) = subplot(4,4,i); 
    xlabel('time (h)');
    ylabel('$n$ (veh)');
    title(strcat('region $',num2str(i-8),'$'))
    box on
    hold on
end

%entire network 
i=12;
f1(3,i-8) = subplot(4,4,i);
xlabel('time (h)');
ylabel('$n$ (veh)');
title('network')
box on
hold on



%fourth row: TS: VQ 

for i=13:15
    f1(4,i-12) = subplot(4,4,i); 
    xlabel('time (h)');
    ylabel('$VQ$ (veh)');
    title(strcat('region $',num2str(i-12),'$'))
    box on
    hold on
end

%entire network 
i=16;
f1(4,i-12) = subplot(4,4,i);
xlabel('time (h)');
ylabel('$VQ$ (veh)');
title('network')
box on
hold on

%Single figures
% 
% fig1 = figure('Name','MFD network');
% xlabel('$n$ (veh)');
% ylabel('$p$ (veh/h)');
% box on
% hold on
% 
% fig2 = figure('Name','TS n');
% xlabel('time (h)');
% ylabel('$n$ (veh)');
% box on
% hold on
% 
% fig3 = figure('Name','TS VQ');
% xlabel('time (h)');
% ylabel('$VQ$ (veh)');
% box on
% hold on
%% Figure 2: Multi-figure: traffic performance measures --------------
% all measures for the entire network (3x3 - 9 figures) 

ff2 = figure('Name','Multi-figure, performance measures');
f2(1) = subplot(3,3,1);
title('(a)');
ylabel('cong. links (lane $\times$ km)');
xlabel('time (h)');
box on
hold on

f2(2) = subplot(3,3,2);
title('(b)')
xlabel('time (h)')
ylabel('$n/c$')
box on
hold on

f2(3) = subplot(3,3,3);
title('(c)')
xlabel('time (h)')
ylabel('cum. trip endings (veh)')
box on
hold on

f2(4) = subplot(3,3,4);
title('(d)')
xlabel('time (h)')
ylabel('trip endings (veh/h)')
box on
hold on

f2(5) = subplot(3,3,5);
title('(e)')
xlabel('time (h)')
ylabel('$p$ (veh km/h)')
box on
hold on

f2(6) = subplot(3,3,6);
title('(f)')
xlabel('time (h)')
ylabel('SMS (km/h)')
box on
hold on

f2(7) = subplot(3,3,7);
title('(g)')
xlabel('time (h)')
ylabel('std nodes (veh)')
box on
hold on

f2(8) = subplot(3,3,8);
title('(h)')
xlabel('time (h)')
ylabel('std links (veh)')
box on
hold on

f2(9) = subplot(3,3,9);
title('(i)')
xlabel('$n$ (veh)')
ylabel('std of $n/c$')
box on
hold on

% fig4 = figure('Name','TS cum trip endings');
% xlabel('time (h)')
% ylabel('cum. trip endings (veh)')
% box on
% hold on
% 
% fig5 = figure('Name','TS std nodes');
% xlabel('time (h)')
% ylabel('std nodes (veh)')
% box on
% hold on


%% Figure 3: Transfer flows between regions + incoming flows from external VQs 
ff3 = figure('Name','Exit/Incoming Flows');
f3(1) = subplot(2,4,1); 
title('(a)');
ylabel('transfer flow 1 to 2');
xlabel('time (h)');
box on
hold on

f3(2) = subplot(2,4,2); 
title('(b)');
ylabel('transfer flow 3 to 2');
xlabel('time (h)');
box on
hold on

f3(3) = subplot(2,4,3); 
title('(c)');
ylabel('transfer flow 2 to 1');
xlabel('time (h)');
box on
hold on

f3(4) = subplot(2,4,4); 
title('(d)');
ylabel('transfer flow 2 to 3');
xlabel('time (h)');
box on
hold on


f3(5) = subplot(2,4,5); 
title('(e)');
ylabel('incoming flow 1');
xlabel('time (h)');
box on
hold on

f3(6) = subplot(2,4,6); 
title('(f)');
ylabel('incoming flow 2');
xlabel('time (h)');
box on
hold on

f3(7) = subplot(2,4,7); 
title('(g)');
ylabel('incoming flow 3');
xlabel('time (h)');
box on
hold on


color_b{1} = 'r';
color_b{2} = 'g';
color_b{3} = 'b';
color_b{4} = 'm';

regColor{1} = '#0072BD';
regColor{2} = '#D95319';
regColor{3} = '#77AC30';

%%


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

colReg{1} = '#0072BD';
colReg{2} = '#D95319';
colReg{3} = '#77AC30';
colReg{4} = '#4DBEEE';
lineTypeScen{1} = ':';
lineTypeScen{2} = '-';
lineTypeScen{3} = '-.';
lineTypeScen{4} = '--';
lineTypeScen{5} = '--';

%% TS - demand global multiplier (demand pattern) 
% ff0 = figure('Name','demand multiplier - time');
% for k=1:indata.kmax
%     tot_demand(k) = indata.defac(k)*sum(indata.demandGroup2(:,1+(k>indata.WU)));
% end
% plot((1:indata.kmax)*indata.DT,tot_demand,'LineWidth',2)
% xlabel('time (h)');
% ylabel('total demand (veh/h)');
Legend = cell(length(scen_ID),1);
for scen = 1:length(scen_ID) %print comparison figures for multiple scenarios 
    
    % store the path of the scenario directory to save the reasulting figures 
    %mkdir(strcat(path_g,scenarios{scen_ID(scen)},'\'));
    path = strcat(path_g,scenarios{scen_ID(scen)},'\');
    
    %load results of the specified scenario
    load(strcat(path_g,fnameProthem,scenarios{scen_ID(scen)},'.mat'),'outdata','indata')
    Legend{scen} = scenarios{scen_ID(scen)}; 
    %% VQs over time (one line per link in common graph)
%     figure('Name','Virtual Queues (all) - TS');
%     hold on
%     title(strcat('Virtual Queues - ',scenarios{scen_ID(scen)}));
%     xlabel('time [h]');
%     ylabel('V. Queue [veh]');
%     plot((1:size(outdata.virtualqueues(1,1:end-1),2))*indata.DT,outdata.virtualqueues(:,1:end-1));
    
    %% MFDs: Entire network -
    t_MFD = 3/60; %time interval between measurements
    t = ceil(t_MFD/indata.DT); %The step / How many steps k*DT are included in the aggregation period 
    sum1 = sum(outdata.x(indata.group1,1:indata.kmax),1); %total network accumulation
    sum2 = sum(outdata.u,1); %total outflow of all links and virtual queues in the network in every step k [veh/DT]
    sum5 = sum((outdata.m(indata.group1,1:indata.kmax)*indata.v_ff)./(indata.Links2(indata.group1,3)*ones(1,indata.kmax)),1); %sum5 = total moving flow in all network links
    sum6 = sum(outdata.w(indata.group2,1:indata.kmax),1); %virtual queues
    sum8 = sum(outdata.q(indata.group3,1:indata.kmax),1); %trip endings 
    sum7 = sum(outdata.u(indata.group1,:).*(indata.Links2(indata.group1,3)/1000),1); % production = outflow of network links * link length
    %sum7 = sum(outdata.m(indata.group1,1:indata.kmax)*indata.v_ff/1000,1); %production = moving*free_flow_speed
    
    congested = outdata.x(:,1:end-1)>0.85*indata.capacity; %links above 85% of capacity 
    congested = congested.*(indata.Links2(:,2).*indata.Links2(:,3)/1000); %(lane x km)
    congested = sum(congested);
    tot_occupancy = sum1/sum(indata.capacity(indata.group1));
    ntw_het = std(outdata.x(indata.group1,1:end-1));
    n_occ_std = std(outdata.x(indata.group1,1:indata.kmax)./indata.capacity(indata.group1)); 
    
    accumulation=0;
    outflow=0;
    k_MFD=0;
    vqueues = 0; 
    trip_endings = 0;
    production = 0;
    m_congested = 0; 
    t_occupancy = 0; 
    het = 0; 
    std_n_occ = 0; 
    for i=t:t:indata.kmax
        k_MFD = [k_MFD i*indata.DT]; %time for MFD [hours] - intervals of aggregation
        accumulation = [accumulation mean(sum1(i-t+1:i))]; %(no VQs included) 
        outflow = [outflow mean(sum2(i-t+1:i))]; %[veh/h]
        vqueues = [vqueues mean(sum6(i-t+1:i))];
        production = [production mean(sum7(i-t+1:i))]; %[veh*km/h]
        trip_endings = [trip_endings mean(sum8(i-t+1:i))]; %[veh/h]
        m_congested = [m_congested mean(congested(i-t+1:i))];
        t_occupancy = [t_occupancy mean(tot_occupancy(i-t+1:i))];
        het = [het mean(ntw_het(i-t+1:i))];
        std_n_occ = [std_n_occ mean(n_occ_std(i-t+1:i))];
    end
    
    % print the entire network figures 
    
    % trip endings 
    figure(ff1)
    subplot(f1(1,4))
    plot(accumulation,movmean(trip_endings,xx),'Linewidth',lineW)

    % production 
    subplot(f1(2,4))
    pp = movmean(production,xx);
    pp(1) = 0;
    plot(accumulation,pp,'Linewidth',lineW)
    
    % accumulation 
    subplot(f1(3,4))
    plot(k_MFD, accumulation,'Linewidth',lineW)
     
    %VQs
    subplot(f1(4,4))
    plot(k_MFD, vqueues,'Linewidth',lineW)
    
    %Single figures 
    
    % figure(fig1)
    % plot(accumulation,pp,'Linewidth',lineW)
    
    % figure(fig2)
    % plot(k_MFD, accumulation,'Linewidth',lineW)
    
    % figure(fig3)
    % plot(k_MFD, vqueues,'Linewidth',lineW)
    
    %performance measures - network 
    
    %congested links
    figure(ff2)
    subplot(f2(1))
    plot(k_MFD, m_congested,'LineWidth',lineW)
    
    %avg ntwrk occupancy 
    subplot(f2(2))
    plot(k_MFD,t_occupancy,'Linewidth',lineW)
    
    %cumulative trip endings
    cumtripsAll = sum(outdata.cumtripsCompleted,1);
    p = [0 cumtripsAll((t:t:indata.kmax))];
    subplot(f2(3))
    plot(k_MFD, p,'LineWidth',lineW)  
    
    %trip endings (rate)
    subplot(f2(4))
    tt = movmean(trip_endings,xx);
    tt(1) = 0; 
    plot(k_MFD,tt,'LineWidth',lineW)
    
    %production (rate)
    subplot(f2(5))
    pp = movmean(production,xx);
    pp(1) = 0; 
    plot(k_MFD,pp ,'LineWidth',lineW)
    
    %SMS 
    subplot(f2(6))
    plot(k_MFD(2:end), SMS(1:(length(k_MFD)-1),scen_ID(scen)),'LineWidth',lineW)
    
    %node std 
    subplot(f2(7))
    tot_node_het = 0;
    for k = t:t:indata.kmax
        node_heterog = [];
        for i=1:size(indata.MP.nodeID,2)
            if indata.MP.approaches{i}>0
                
                %Mean over time of the std of incoming queues around the
                %intersection every cycle
                node_heterog = [node_heterog mean(std(outdata.w(indata.junct.or_index(indata.MP.approaches{i}),k-t+1:k)))];
            end
        end
        tot_node_het = [tot_node_het mean(node_heterog)];
    end 
    plot(k_MFD,tot_node_het(1:length(k_MFD)),'LineWidth',lineW)
    
    %link queue std 
    subplot(f2(8))
    plot(k_MFD, (1:length(k_MFD)),'LineWidth',lineW)
    
    % n - std of n/c (network heterogeneity) 
    subplot(f2(9))
    plot(accumulation, std_n_occ,'LineWidth',lineW)     
    
    
    
    % figure(fig4)
    % %cumulative trip endings
    % plot(k_MFD, p,'LineWidth',lineW)  
    
    % figure(fig5)
    % plot(k_MFD,tot_node_het,'LineWidth',lineW)
    
    % print the regional figures 
    % Flow MFDs per region - common graph for multiple runs

    Rsum1 = zeros(3,indata.kmax); 
    Rsum2 = zeros(3,indata.kmax); 
    Rsum6 = zeros(3,indata.kmax); 
    Rsum7 = zeros(3,indata.kmax); 
    Rsum8 = zeros(3,indata.kmax); 
    figure(ff1) 
    for i=1:3
        
        % Rsum1(i,:) = sum(outdata.x(intersect(indata.group1,indata.reInd{i}),1:indata.kmax),1) + sum(outdata.w(intersect(indata.group2,indata.reInd{i}),1:indata.kmax),1);  %region accumulation + VQs
        Rsum1(i,:) = sum(outdata.x(intersect(indata.group1,indata.reInd{i}),1:indata.kmax));  %region accumulation only
        Rsum2(i,:) = sum(outdata.u(indata.reInd{i},1:indata.kmax),1); %outflows of all links and virtual queues of the region
        %Rsum4(i,:) = sum(outdata.a_z(intersect(indata.reInd{i},indata.group1),1:indata.kmax),1); %total arriving flow at the end of the queue per region
        %Rsum5(i,:) = sum((outdata.m(intersect(indata.reInd{i},indata.group1),1:indata.kmax)*indata.v_ff)./(indata.Links2(intersect(indata.reInd{i},indata.group1),3)*ones(1,indata.kmax)),1); %total moving flow in all network links
        %Rsum3(i,:) = sum(outdata.u(intersect(indata.group2,indata.reInd{i}),1:indata.kmax),1); %total inflow of region
        %Rsum7(i,:) = sum(outdata.m(intersect(indata.group1,indata.reInd{i}),1:indata.kmax)*indata.v_ff/1000,1); %production = moving*free_flow_speed
        Rsum7(i,:) = sum(outdata.u(intersect(indata.group1,indata.reInd{i}),1:indata.kmax).*(indata.Links2(intersect(indata.group1,indata.reInd{i}),3)/1000),1); %production = outflow * link length
        Rsum6(i,:) = sum(outdata.w(intersect(indata.group2,indata.reInd{i}),1:indata.kmax),1); %virtual queues
        Rsum8(i,:) = sum(outdata.q(intersect(indata.group3,indata.reInd{i}),1:indata.kmax),1); %trip endings 
        
        accumulation=0;
        outflow=0;
        vqueues = 0; 
        trip_endings = 0;
        production = 0;
        for j=t:t:indata.kmax
            accumulation = [accumulation mean(Rsum1(i,j-t+1:j))]; %[veh] - mean accumulation over the time window
            outflow = [outflow mean(Rsum2(i,j-t+1:j))]; %[veh/h]
            vqueues = [vqueues mean(Rsum6(i,j-t+1:j))];
            production = [production mean(Rsum7(i,j-t+1:j))]; %[veh*km/h]
            trip_endings = [trip_endings mean(Rsum8(i,j-t+1:j))]; %[veh/h]
        end
        %MFD n-e
        figure(ff1)
        subplot(f1(1,i))
        tt = movmean(trip_endings,xx);
        tt(1) = 0; 
        plot(accumulation,tt,'Linewidth',lineW)
        
        %MFD n-p
        subplot(f1(2,i))
        pp = movmean(production,xx);
        pp(1) = 0; 
        plot(accumulation, pp,'Linewidth',lineW)
        %plot(accumulation, production,'Linewidth',1.5,'LineStyle',lineTypeScen{scen})
        %scatter(accumulation, production,'filled','MarkerFaceColor',colReg{i})
        
        %TS: n 
        subplot(f1(3,i))
        plot(k_MFD, accumulation,'Linewidth',lineW)
        
        %TS: VQs
        subplot(f1(4,i))
        plot(k_MFD, vqueues,'Linewidth',lineW)
    
        %regional heterogeneity
        % Heterogeneity Figures (per region): accumulation vs. stdev of queues x
        %/check time - it takes too long / separate regions (?) or scenarios -
        %too many lines
        
        %         ind = intersect(indata.group1,indata.reInd{i});
        %         occ = outdata.x(ind,1:indata.kmax);
        %         cap = zeros(length(ind),indata.kmax);
        %         for j = 1:indata.kmax
        %             cap(:,j) = indata.capacity(ind);
        %         end
        %         occ = occ./cap;
        %         Rstdev(i,:) = std(occ);
        %         stdev_x = 0;
        %         accumulation = 0;
        %         for j=1:t:indata.kmax
        %             accumulation = [accumulation mean(Rsum1(i,max(1,j):j+(t-1)))];
        %             stdev_x=[stdev_x mean(Rstdev(i,j:j+(t-1)))]; %[veh/h]
        %         end
        %         if i==1
        %             figure(f14)
        %         elseif i==2
        %             figure(f14b)
        %         else
        %             figure(f14c)
        %         end
        %         plot(accumulation, stdev_x,'LineWidth',1.5,'Color',colReg{i},'LineStyle',lineTypeScen{scen})
        
        
    end
     
    
    %%     Movie - Showcase the congestion dynamic
    %     step = 5*36;% Time interval between snapshots for the video in # of ks - k = DT sec / time step = step * DT [sec]
    %     v=createmovie(outdata.x,indata.kmax,indata.capacity,Links,Nodes,DT,step,'S');
    %     Save some results
    %     accumulation=accumulation';
    %     outflow=outflow';
    %     save outflowS.mat outflow
    %     save accumulationS.mat accumulation
       
 
    %% Transfer flows between adjacent regions over time
    figure(ff3)
    %indices of Links of every approach
    adjRg = [1 2;3 2;2 1;2 3];
    for j=1:4
        ind1 = ismember(indata.junct2.or_index, indata.reInd{adjRg(j,1)});
        ind2 = ismember(indata.junct2.dest_index, indata.reInd{adjRg(j,2)});
        ind = ind1+ind2 == 2;
        reg_outflow = sum(outdata.u(indata.junct2.or_index(ind),:));
        regOut = 0;
        for i=t:t:indata.kmax
            regOut = [regOut mean(reg_outflow(i-t+1:i))]; %[veh] - mean accumulation over the time window (VQ included)
        end
        subplot(f3(j))
        plot(k_MFD,movmean(regOut,xx),'LineWidth',lineW)
    end
    
    %figures of starting trips per region (incoming flow - regulated by PC)
    for j=1:3
       %count outflows of virtual queues (origin links) 
       ind1 = intersect(indata.group2,indata.reInd{j});
       reg_inflow = sum(outdata.u(ind1,:));
       reg_in = 0; 
       for i=t:t:indata.kmax
          reg_in = [reg_in mean(reg_inflow(i-t+1:i))]; 
       end
       
       subplot(f3(j+4))
       plot(k_MFD,movmean(reg_in,xx),'LineWidth',lineW)  
       
    end
    
    if outdata.PC.mode == 1
        %reads from the last recorded clustering (needs to adapt to show
        %previous clusters if X and Y change) 
        load('fromToRegions.mat','X','Y') %-> change to .csv if needed
        

        %It prints applied_u for all rows of the junctionsID matrix (even
        %if no nodes are controlled - to be adjusted if needed) 
        figure('Name',strcat('Control Variables ',scenarios{scen_ID(scen)}))
        hold on
        for i=1:4
            %yyaxis left
            plot((1:length(outdata.applied_u))*indata.c_int*indata.DT,outdata.applied_u(i,:),'Color',color_b{i},'LineWidth',lineW)
        end
        for i=5:7
            %yyaxis right
            plot((1:length(outdata.applied_u))*indata.c_int*indata.DT,outdata.applied_u(i,:)/100*90,'Color',regColor{i-4},'LineWidth',lineW)
        end
        legend(strcat(num2str(X(1)),'-',num2str(Y(1))),strcat(num2str(X(2)),'-',num2str(Y(2))),...
            strcat(num2str(X(3)),'-',num2str(Y(3))),strcat(num2str(X(4)),'-',num2str(Y(4))),...
            'external 1','external 2','external 3','Location','westoutside')
        %legend('1-2','3-2','2-1','2-3','external 1','external 2','external 3','Location','westoutside')
    end
    
end

% %Legends and saving 
% figure(ff1)
% legend(Legend) 
% 
% figure(ff2) 
% legend(Legend)
% 
% figure(ff3) 
% legend(Legend)
% 
% figure(fig1)
% % legend('FTC','MP $100\%$','MP $25\%$','PC + MP $25\%$','location','southeast')
% legend('FTC','MP $100\%$','PC','PC + MP $25\%$','location','southeast')
% fname = 'MFDs_2b';
% saveas(gcf,strcat(pathFigJ,fname,'.fig'));
% saveas(gcf,strcat(pathFigJ,fname,'.emf'));
% print(strcat(pathFigJ,fname,'.eps'),'-depsc','-painters')
%     
% figure(fig2)
% % legend('FTC','MP $100\%$','MP $25\%$','PC + MP $25\%$')
% legend('FTC','MP $100\%$','PC','PC + MP $25\%$')
% fname = 'TSn_2b'; 
% saveas(gcf,strcat(pathFigJ,fname,'.fig'));
% saveas(gcf,strcat(pathFigJ,fname,'.emf'));
% print(strcat(pathFigJ,fname,'.eps'),'-depsc','-painters')
% 
% figure(fig3)
% % legend('FTC','MP $100\%$','MP $25\%$','PC + MP $25\%$')
% legend('FTC','MP $100\%$','PC','PC + MP $25\%$')
% fname = 'TSVQ_2b';
% saveas(gcf,strcat(pathFigJ,fname,'.fig'));
% saveas(gcf,strcat(pathFigJ,fname,'.emf'));
% print(strcat(pathFigJ,fname,'.eps'),'-depsc','-painters')
% 
% figure(fig4)
% % legend('FTC','MP $100\%$','MP $25\%$','PC + MP $25\%$','location','southeast')
% legend('FTC','MP $100\%$','PC','PC + MP $25\%$','location','southeast')
% fname = 'TScumTrips_2b';
% saveas(gcf,strcat(pathFigJ,fname,'.fig'));
% saveas(gcf,strcat(pathFigJ,fname,'.emf'));
% print(strcat(pathFigJ,fname,'.eps'),'-depsc','-painters')
% 
% figure(fig5)
% % legend('FTC','MP $100\%$','MP $25\%$','PC + MP $25\%$')
% legend('FTC','MP $100\%$','PC','PC + MP $25\%$')
% fname = 'TSstdNodes_2b';
% saveas(gcf,strcat(pathFigJ,fname,'.fig'));
% saveas(gcf,strcat(pathFigJ,fname,'.emf'));
% print(strcat(pathFigJ,fname,'.eps'),'-depsc','-painters')

%%  Plot: Maps 

% legend('FTC','PC','PC + MP $25\%$')
% saveas(gcf,strcat(pathFigJ,'multiMFD_2.fig'));
% print(strcat(pathFigJ,'multiMFD_2.eps'),'-depsc','-painters')

