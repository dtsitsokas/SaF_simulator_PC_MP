% Produce detailed figures for PC evaluation after simulation ends 
% (old - use the rsultAnalysis for all cases) 
clear
%clc
%close all
%path to retrieve results file 
% path = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';
path = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
load scenarios.mat 

%load results of No Control Case 
% load(strcat(path,'output_set5_highDemand_NC.mat'),'outdata','indata')
% load(strcat(path,'output_set5_MedDemand_NC_test.mat'),'outdata','indata')
load(strcat(path,'output_set5_MedDemand_MP_cd_MS4_test.mat'),'outdata','indata')

scen2comp = 45; %scenario to compare with NC (only one every time) 

%settings 
step_avg = 2*indata.c_int; %time window for the averaging of MFDs etc. [no of time steps]  
% fnameProthem = 'output_set5_highDemand_'; %(high OD) prothema of the name of the results file 
fnameProthem = 'output_set5_MedDemand_'; %(medium OD) 

% initialize figures 
%% MFDs - trip endings / travel production + time-series accumulation/VQs 
ff1 = figure('Name','Multifigure MFDs');
hold on 
%Accumulation - travel production/trip endings 
%Load and plot the data of the NC case in dashed lines
f1 = subplot(3,3,1);
xlabel('n (veh)')
ylabel('trip endings (veh/h)')
%ylabel('p (veh km/h)')
title('region 1')
hold on

f1b = subplot(3,3,2);
xlabel('n (veh)')
ylabel('trip endings (veh/h)')
%ylabel('p (veh km/h)')
title('region 2')
hold on

f1c = subplot(3,3,3);
xlabel('n (veh)')
%ylabel('p (veh km/h)')
ylabel('trip endings (veh/h)')
title('region 3')
hold on

f1d = subplot(3,3,4);
xlabel('n (veh)')
%ylabel('trip endings (veh/h)')
ylabel('p (veh km/h)')
title('region 1')
hold on

f1e = subplot(3,3,5);
xlabel('n (veh)')
%ylabel('trip endings (veh/h)')
ylabel('p (veh km/h)')
title('region 2')
hold on

f1f = subplot(3,3,6);
xlabel('n (veh)')
ylabel('p (veh km/h)')
%ylabel('trip endings (veh/h)')
title('region 3')
hold on

% Time-Series of n and n_Vq per region
f2 = subplot(3,3,7);
title('region 1')
xlabel('time (h)')
ylabel('n (veh)')
hold on

f2b = subplot(3,3,8);
title('region 2')
xlabel('time (h)')
ylabel('n (veh)')
hold on

f2c = subplot(3,3,9);
title('region 3')
xlabel('time (h)')
ylabel('n (veh)')
hold on

regColor{1} = '#0072BD';
regColor{2} = '#D95319';
regColor{3} = '#77AC30';

% ----------------------

% Figure for cumulative VHT over time (vs NC case)
ff5 = figure('Name','VHT over time');
f5a = subplot(1,4,1);
title('region 1')
xlabel('time (h)')
ylabel('VHT')
hold on

f5b = subplot(1,4,2);
title('region 2')
xlabel('time (h)')
ylabel('VHT')
hold on

f5c = subplot(1,4,3);
title('region 3')
xlabel('time (h)')
ylabel('VHT')
hold on

f5d = subplot(1,4,4);
title('total')
xlabel('time (h)')
ylabel('VHT')
hold on


% figures for total inflows/outflows per region
%Figure with TS of total inflows(VQ-transfer) and outflows(exit/transfer) per region
ff2 = figure('Name','TS Flows');
f4 = subplot(3,1,1); %region 1
title('flows region 1')
xlabel('time (h)')
ylabel('flow (veh/h)')
hold on

f4b = subplot(3,1,2); %region 2
title('flows region 2')
xlabel('time (h)')
ylabel('flow(veh/h)')
hold on

f4c = subplot(3,1,3); %region 3
title('flows region 3')
xlabel('time (h)')
ylabel('flow (veh/h)')
hold on

adjRg = [1 2;3 2;2 1;2 3]; %neighboring_regions
for j=1:4
    ind1 = ismember(indata.junct2.or_index, indata.reInd{adjRg(j,1)});
    ind2 = ismember(indata.junct2.dest_index, indata.reInd{adjRg(j,2)});
    nr_ind{j} = ind1+ind2 == 2;
    nr_outflow(j,:) = sum(outdata.u(indata.junct2.or_index(nr_ind{j}),:));
end

% print results of No Control Case 

for reg=1:indata.no_reg
    figure(ff1)
    ac = 0;
    pr = 0;
    vq = 0;
    nr_out = 0;
    nr_in = 0;
    r1 = 0;
    r3 = 0;
    tr_e = 0;
    
    % -- for MFDs with subsets of links (mean -> mean flow per link of the
    %subset / sum -> mean total flow of the subset)
    %load subsets of regions
    %load('subsets_regions.mat','subset_reg','non_subset_reg')
    %Rsum1(reg,:) = sum(outdata.x(subset_reg{reg},1:indata.kmax));  %subset accumulation only
    %Rsum5(reg,:) = mean((outdata.m(subset_reg{reg},1:indata.kmax)*indata.v_ff)./(indata.Links2(subset_reg{reg},3)*ones(1,indata.kmax)),1); %total moving flow in links
    %NRsum1(reg,:) = sum(outdata.x(non_subset_reg{reg},1:indata.kmax));  %non-subset accumulation only
    %NRsum5(reg,:) = mean((outdata.m(non_subset_reg{reg},1:indata.kmax)*indata.v_ff)./(indata.Links2(non_subset_reg{reg},3)*ones(1,indata.kmax)),1);
    
    subset_ac = 0;
    subset_outfl = 0;
    non_subset_ac = 0;
    non_subset_outfl = 0;
    
    for k=step_avg:step_avg:indata.kmax
        ac = [ac mean(sum(outdata.x(intersect(indata.group1,indata.reInd{reg}),k-step_avg+1:k),1))]; % accum in reg 1 (excluding VQ)
        
        %pr = [pr mean(sum(outdata.u(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)...
        %    .*(indata.Links2(intersect(indata.group1,indata.reInd{reg}),3)/1000),1))]; %production of real links in reg 1
        
        pr = [pr mean(sum(outdata.m(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)*indata.v_ff/1000 + ...
           outdata.u(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)...
          .*(outdata.w(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)*(indata.vehlength/1000)...
          ./indata.Links2(intersect(indata.group1,indata.reInd{reg}),3)/1000)))]; %production of real links in reg 2 (veh/km)
        
        vq = [vq mean(sum(outdata.w(intersect(indata.group2,indata.reInd{reg}),k-step_avg+1:k),1))];
        tr_e =[tr_e mean(sum(outdata.q(intersect(indata.group3,indata.reInd{reg}),k-step_avg+1:k),1))];
        nr_out = [nr_out mean(sum(nr_outflow(adjRg(:,1)==reg,k-step_avg+1:k),1))];
        nr_in = [nr_in mean(sum(nr_outflow(adjRg(:,2)==reg,k-step_avg+1:k),1))];
        r1 = [r1 mean(sum(outdata.u(intersect(indata.group2,indata.reInd{reg}),k-step_avg+1:k)))];
        r3 = [r3 mean(sum(outdata.q(intersect(indata.group3,indata.reInd{reg}),k-step_avg+1:k)))];
        %subset_ac = [subset_ac mean(Rsum1(reg,k-step_avg+1:k))]; %[veh] - mean accumulation over the time window
        %subset_outfl = [subset_outfl mean(Rsum5(reg,k-step_avg+1:k))]; %[veh/h]
        %non_subset_ac = [non_subset_ac mean(NRsum1(reg,k-step_avg+1:k))]; %[veh] - mean accumulation over the time window
        %non_subset_outfl = [non_subset_outfl mean(NRsum5(reg,k-step_avg+1:k))]; %[veh/h]
        
    end
    if reg==1
        subplot(f1)
    elseif reg==2
        subplot(f1b)
    else
        subplot(f1c)
    end
    hold on
    plot(ac,tr_e,'Color','k','LineStyle','--','LineWidth',1) %accumulation-trip endings
    
    if reg==1
        subplot(f1d)
    elseif reg==2
        subplot(f1e)
    else
        subplot(f1f)
    end
    hold on
    plot(ac,pr,'Color','k','LineStyle','--','LineWidth',1) %accumulation-production
    
    if reg==1
        subplot(f2)
    elseif reg==2
        subplot(f2b)
    else
        subplot(f2c)
    end
    hold on
    plot((0:step_avg:indata.kmax)*indata.DT,ac,'Color','#0072BD','LineStyle','--','LineWidth',1.5) %time-series accumulation
    plot((0:step_avg:indata.kmax)*indata.DT,vq,'Color','#D95319','LineStyle','--','LineWidth',1.5) %time-series virtual queues
    
    figure(ff2)
    if reg==1
        subplot(f4)
    elseif reg==2
        subplot(f4b)
    else
        subplot(f4c)
    end
    hold on
    
    %in/out-flows per region
    
    %total external inflow (outflows of all VQs)
    plot((0:step_avg:indata.kmax)*indata.DT, r1,'Color','r','LineStyle','--','LineWidth',1.5)
    %total exiting flow (trip endings = inflows of exit links)
    plot((0:step_avg:indata.kmax)*indata.DT, r3,'Color','g','LineStyle','--','LineWidth',1.5)
    %total tranfer outflow (to neighboring regions)
    plot((0:step_avg:indata.kmax)*indata.DT, nr_out,'Color','b','LineStyle','--','LineWidth',1.5)
    %total transfer inflow (from neighboring regions)
    plot((0:step_avg:indata.kmax)*indata.DT, nr_in,'Color','m','LineStyle','--','LineWidth',1.5)
    
    %VHT plot
    figure(ff5)
    
    if reg==1
        subplot(f5a)
    elseif reg==2
        subplot(f5b)
    else
        subplot(f5c)
    end
    hold on
    r1 = sum(outdata.x(intersect(indata.group1,indata.reInd{reg}),:)*indata.DT); %network
    r2 = sum(outdata.w(intersect(indata.group2,indata.reInd{reg}),:)*indata.DT); %VQ
    r3 = r1 + r2;
    plot((0:indata.kmax)*indata.DT, cumsum(r1),'LineWidth',1.5,'Color','#0072BD','LineStyle','-.') %ntw reg
    plot((0:indata.kmax)*indata.DT, cumsum(r2),'LineWidth',1.5,'Color','#D95319','LineStyle','-.') %VQ reg
    plot((0:indata.kmax)*indata.DT, cumsum(r3),'LineWidth',2,'Color','k','LineStyle','-.') %total region
    
    %     %MFDs for subset of links - first print NC results and then live
    %     %figures for the PC scenario
    %     figure(ff6)
    %     if reg==1
    %         subplot(f6a)
    %     elseif reg==2
    %         subplot(f6b)
    %     else
    %         subplot(f6c)
    %     end
    %     plot(subset_ac,subset_outfl,'Color',regColor{reg},'Linewidth',2,'LineStyle','--')
    %
    %     %MFDs for non-subsets of links
    %     figure(ff7)
    %     if reg==1
    %         subplot(f7a)
    %     elseif reg==2
    %         subplot(f7b)
    %     else
    %         subplot(f7c)
    %     end
    %     plot(non_subset_ac,non_subset_outfl,'Color',regColor{reg},'Linewidth',2,'LineStyle','--')
    
end

%total VHT (network)
r1 = sum(outdata.x((indata.group1),:)*indata.DT); %network
r2 = sum(outdata.w((indata.group2),:)*indata.DT); %VQ
r3 = r1 + r2;
figure(ff5)
subplot(f5d)
hold on
plot((0:indata.kmax)*indata.DT, cumsum(r3),'LineWidth',1.5,'Color','k','LineStyle','-.') %total ntw
plot((0:indata.kmax)*indata.DT, cumsum(r1),'LineWidth',1.5,'Color','#0072BD','LineStyle','-.') %ntw reg
plot((0:indata.kmax)*indata.DT, cumsum(r2),'LineWidth',1.5,'Color','#D95319','LineStyle','-.') %VQ reg

color_b{1} = 'r';
color_b{2} = 'g';
color_b{3} = 'b';
color_b{4} = 'm';
% Print results of controlled cases (1 or more) 

for j=1:length(scen2comp)
    %load results of No Control Case 
    load(strcat(path,fnameProthem,scenarios{scen2comp(j)},'_test.mat'),'outdata','indata')
    
    
    if outdata.PC.mode == 1
        % figure of control decision time-series
        f3 = figure('Name','TS of app_u for in/out time');
        xlabel('time (h)')
        ylabel('sec')
        %yyaxis right
        %ylabel('\%')
        hold on
    end

    for reg=1:indata.no_reg
        figure(ff1)
        ac = 0;
        pr = 0;
        vq = 0;
        nr_out = 0;
        nr_in = 0;
        r1 = 0;
        r3 = 0;
        tr_e = 0;

        % -- for MFDs with subsets of links (mean -> mean flow per link of the
        %subset / sum -> mean total flow of the subset)
        %load subsets of regions
        %load('subsets_regions.mat','subset_reg','non_subset_reg')
        %Rsum1(reg,:) = sum(outdata.x(subset_reg{reg},1:indata.kmax));  %subset accumulation only
        %Rsum5(reg,:) = mean((outdata.m(subset_reg{reg},1:indata.kmax)*indata.v_ff)./(indata.Links2(subset_reg{reg},3)*ones(1,indata.kmax)),1); %total moving flow in links
        %NRsum1(reg,:) = sum(outdata.x(non_subset_reg{reg},1:indata.kmax));  %non-subset accumulation only
        %NRsum5(reg,:) = mean((outdata.m(non_subset_reg{reg},1:indata.kmax)*indata.v_ff)./(indata.Links2(non_subset_reg{reg},3)*ones(1,indata.kmax)),1);

        subset_ac = 0;
        subset_outfl = 0;
        non_subset_ac = 0;
        non_subset_outfl = 0;
        
        for k=step_avg:step_avg:indata.kmax
            ac = [ac mean(sum(outdata.x(intersect(indata.group1,indata.reInd{reg}),k-step_avg+1:k),1))]; % accum in reg 1 (excluding VQ)

            %pr = [pr mean(sum(outdata.u(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)...
            %   .*(indata.Links2(intersect(indata.group1,indata.reInd{reg}),3)/1000),1))]; %production of real links in reg 1

            pr = [pr mean(sum(outdata.m(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)*indata.v_ff/1000 + ...
                outdata.u(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)...
                .*(outdata.w(intersect(indata.group1,indata.reInd{reg}),k-indata.c_int+1:k)*(indata.vehlength/1000)...
                 ./indata.Links2(intersect(indata.group1,indata.reInd{reg}),3)/1000)))]; %production of real links in reg 2 (veh/km)

            vq = [vq mean(sum(outdata.w(intersect(indata.group2,indata.reInd{reg}),k-step_avg+1:k),1))];
            tr_e =[tr_e mean(sum(outdata.q(intersect(indata.group3,indata.reInd{reg}),k-step_avg+1:k),1))];
            nr_out = [nr_out mean(sum(nr_outflow(adjRg(:,1)==reg,k-step_avg+1:k),1))];
            nr_in = [nr_in mean(sum(nr_outflow(adjRg(:,2)==reg,k-step_avg+1:k),1))];
            r1 = [r1 mean(sum(outdata.u(intersect(indata.group2,indata.reInd{reg}),k-step_avg+1:k)))];
            r3 = [r3 mean(sum(outdata.q(intersect(indata.group3,indata.reInd{reg}),k-step_avg+1:k)))];
            %subset_ac = [subset_ac mean(Rsum1(reg,k-step_avg+1:k))]; %[veh] - mean accumulation over the time window
            %subset_outfl = [subset_outfl mean(Rsum5(reg,k-step_avg+1:k))]; %[veh/h]
            %non_subset_ac = [non_subset_ac mean(NRsum1(reg,k-step_avg+1:k))]; %[veh] - mean accumulation over the time window
            %non_subset_outfl = [non_subset_outfl mean(NRsum5(reg,k-step_avg+1:k))]; %[veh/h]

        end
        hold on
        if reg==1
            subplot(f1)
        elseif reg==2
            subplot(f1b)
        else
            subplot(f1c)
        end
        hold on 
        %plot(ac,pr,'Color','k','LineStyle','--') %accumulation-production
        plot(ac,tr_e,'Color','r','LineWidth',1.5) %accumulation-trip endings
       
        
        if reg==1
            subplot(f1d)
        elseif reg==2
            subplot(f1e)
        else
            subplot(f1f)
        end
        hold on
        plot(ac,pr,'Color','r','LineWidth',1.5) %accumulation-production
        
           
        if reg==1
            subplot(f2)
        elseif reg==2
            subplot(f2b)
        else
            subplot(f2c)
        end
        hold on 
        plot((0:step_avg:indata.kmax)/3600,ac,'Color','#0072BD','LineWidth',1.5)
        plot((0:step_avg:indata.kmax)/3600,vq,'Color','#D95319','LineWidth',1.5)
        
        figure(ff2)
        if reg==1
            subplot(f4)
        elseif reg==2
            subplot(f4b)
        else
            subplot(f4c)
        end

        %in/out-flows per region

        %total external inflow (outflows of all VQs)
        plot((0:step_avg:indata.kmax)/3600, r1,'Color','r','LineWidth',1.5)
        %total exiting flow (trip endings = inflows of exit links)
        plot((0:step_avg:indata.kmax)/3600, r3,'Color','g','LineWidth',1.5)
        %total tranfer outflow (to neighboring regions)
        plot((0:step_avg:indata.kmax)/3600, nr_out,'Color','b','LineWidth',1.5)
        %total transfer inflow (from neighboring regions)
        plot((0:step_avg:indata.kmax)/3600, nr_in,'Color','m','LineWidth',1.5)

        %VHT plot
        figure(ff5)

        if reg==1
            subplot(f5a)
        elseif reg==2
            subplot(f5b)
        else
            subplot(f5c)
        end

        r1 = sum(outdata.x(intersect(indata.group1,indata.reInd{reg}),:)*indata.DT); %network
        r2 = sum(outdata.w(intersect(indata.group2,indata.reInd{reg}),:)*indata.DT); %VQ
        r3 = r1 + r2;
        plot((0:indata.kmax)*indata.DT, cumsum(r1),'LineWidth',1.5,'Color','#0072BD') %ntw reg
        plot((0:indata.kmax)*indata.DT, cumsum(r2),'LineWidth',1.5,'Color','#D95319') %VQ reg
        plot((0:indata.kmax)*indata.DT, cumsum(r3),'LineWidth',2,'Color','k') %total region

    %     %MFDs for subset of links - first print NC results and then live
    %     %figures for the PC scenario
    %     figure(ff6)
    %     if reg==1
    %         subplot(f6a)
    %     elseif reg==2
    %         subplot(f6b)
    %     else
    %         subplot(f6c)
    %     end
    %     plot(subset_ac,subset_outfl,'Color',regColor{reg},'Linewidth',2,'LineStyle','--')
    %     
    %     %MFDs for non-subsets of links
    %     figure(ff7)
    %     if reg==1
    %         subplot(f7a)
    %     elseif reg==2
    %         subplot(f7b)
    %     else
    %         subplot(f7c)
    %     end
    %     plot(non_subset_ac,non_subset_outfl,'Color',regColor{reg},'Linewidth',2,'LineStyle','--')
      
    
    end
    
    %total VHT (network)
    r1 = sum(outdata.x((indata.group1),:)*indata.DT); %network
    r2 = sum(outdata.w((indata.group2),:)*indata.DT); %VQ
    r3 = r1 + r2;
    figure(ff5)
    subplot(f5d)
    plot((0:indata.kmax)*indata.DT, cumsum(r3),'LineWidth',1.5,'Color','k','LineStyle','-') %total ntw
    plot((0:indata.kmax)*indata.DT, cumsum(r1),'LineWidth',1.5,'Color','#0072BD','LineStyle','-') %ntw reg
    plot((0:indata.kmax)*indata.DT, cumsum(r2),'LineWidth',1.5,'Color','#D95319','LineStyle','-') %VQ reg
    
    %Control decisions 
    if outdata.PC.mode == 1 
        figure(f3)
        hold on 
        for i=1:4
            %yyaxis left
            plot((1:length(outdata.applied_u))*indata.c_int*indata.DT,outdata.applied_u(i,:),'Color',color_b{i},'LineWidth',1.5,'LineStyle','-')
        end
        for i=5:7
            %yyaxis right
            plot((1:length(outdata.applied_u))*indata.c_int*indata.DT,outdata.applied_u(i,:)*90/100,'Color',regColor{i-4},'LineWidth',2,'LineStyle','-')
        end
        legend('1-2','3-2','2-1','2-3','Ext 1','Ext 2','Ext 3','Location','westoutside')

    end
end

% %% print snapshots of the network at specific times / select the metric 
% load LINKS %links info [id, start_node, 3nd_node, length, cluster (in 4 reg)]
% load nodes3 %nodes (ids, x coord, y coord) - updated file Nodes to Nodes3 (some missing centroids) 
% 
% snap_t = 6.00; %hours since simul start 
% metric = outdata.x(:,ceil(snap_t/indata.DT));
% metric_n = indata.capacity;
% metric_name = 'link occupancy';
% createSnapshot(snap_t,metric,metric_n,metric_name,LINKS,nodes3)
% 
% 
% %% save figures 
% 
% figure(ff1)
% %saveas(gcf,strcat(path,'\',scenarios{scen2comp},'\','multifigureMFD.fig'));
% 
% figure(ff2)
% legend()
% %saveas(gcf,strcat(path,'\',scenarios{scen2comp},'\','flowsTS.fig'));
% 
% figure(ff5)
% %saveas(gcf,strcat(path,'\',scenarios{scen2comp},'\','cumulativeVHT.fig'));
% 
% if outdata.PC.mode == 1
%     figure(f3)
%     legend('1-2','3-2','2-1','2-3','ext 1','ext 2','ext 3','Location','westoutside')
%     %saveas(gcf,strcat(path,'\',scenarios{scen2comp},'\','TS_controlVariables.fig'));
% end
% %%
% % Create link occupancy video from output files - normalize over a specific
% % value for all: to see the differences between cases 
% 
% start_time = 1/60; %[hours] 
% end_time = 6;
% 
% start_step = floor(start_time/indata.DT); % #first step 
% end_step = floor(end_time/indata.DT); % #last step
% 
% step = floor(indata.c_int/(3600*indata.DT));% no of steps between snapshots for the video in # 
% fname = strcat(scenarios{scen2comp},'_from',num2str(start_time*3600),...
%     'to',num2str(end_time*3600),'sec'); 
% 
% v = createmovie(outdata.x(:,(start_step:end_step)),step,start_step,indata.capacity...
%     ,LINKS,nodes3,indata.DT,fname);
%     

% Create difference maps (PC - NC): to observe differences between
% scenarios on the map 


%% fix the legends of the figures (to do in the end for the final figures) 
