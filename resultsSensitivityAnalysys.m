%figures for sensitivity analysis 

% No control case (x replications)
%clear 
%clc 

x = 20; 
demandCode = 2; 
j = 33; %scenario of comparison to NC  
pcntALL = [0.05 0.10 0.15 0.20]; 
load scenarios.mat

scenarios{33} = 'PC_50b_MP_MS4';

if demandCode == 1 %med
    pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
    %output file name to save after simulation is done
    path = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
    pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';
elseif demandCode == 2
    pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
    path = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';
end

VHTrepl = zeros(x,2); %2 columns: 1 is NC, 2 is the tested case 
% Read results files 
%indata: only network info used (constant) 
load(strcat(path,'ODrepliNo_',num2str(1),'_',num2str(100*pcntALL(1)),'_output_',scenarios{j},'.mat'),'indata')
%% VHTreplALL = cell(1,3);
for i=1:length(pcntALL)
    pcnt = pcntALL(i);
    for repl=1:x
        
        %load results file of NC
        fname = strcat(path,'ODrepliNo_',num2str(repl),'_',num2str(100*pcnt),'_output_',scenarios{1},'.mat');
        load(fname,'outdata')
        r1 = sum(sum(outdata.x*indata.DT));
        r2 = sum(sum(outdata.w(indata.group2,:)*indata.ksi*indata.DT));
        VHTrepl(repl,1) = r1 + r2; %1st column for the NC case
                
        %load results file of NC
        fname = strcat(path,'ODrepliNo_',num2str(repl),'_',num2str(100*pcnt),'_output_',scenarios{j},'.mat');
        load(fname,'outdata')
        r1 = sum(sum(outdata.x*indata.DT));
        r2 = sum(sum(outdata.w(indata.group2,:)*indata.ksi*indata.DT));
        VHTrepl(repl,2) = r1 + r2; %1st column for the NC case
    end
    VHTreplALL{i} = VHTrepl; 
    save(strcat('sensAnalResults_',num2str(demandCode),'_',num2str(100*pcnt),'.mat'),'VHTrepl')
    
end
save(strcat('sensAnalResultsALL_',num2str(demandCode), 'VHTreplALL'))

%% load Result files 

%load(strcat('sensAnalResultsALL_',num2str(demandCode)));

%% Create figures / boxplots 
%med
for i=1:4
   f = figure();
   boxplot(VHTreplALL{1,i})
   title(strcat(num2str(i*5),'\% std, (20 repl.)'))
   ax = gca;
   ax.TickLabelInterpreter = 'latex';
   %xticklabels({'FTC','MP 25\%'})
   xticklabels({'FTC','PC + MP 25\%'})
   ylabel('VHT')
end

%%
saveas(gcf,'sensAnalhigh20.fig');
print('sensAnalhigh20.eps','-depsc','-painters')

%% Create MFDs of the network for the FTC case 

ff1 = figure('Name', 'MFD n-p NC network');
hold on
xlabel('$n$ (veh)')
ylabel('$p$ (veh km/h)')
box on


ff2 = figure('Name', 'MFD n-e NC network');
hold on
xlabel('$n$ (veh)')
ylabel('$e$ (veh/h)')
box on
%
ff3 = figure('Name','Reg MFDs n-p');
hold on 
xlabel('$n$ (veh)')
ylabel('$p$ (veh km/h)')
box on 

xx = 5; %window for moving average
% read all files of 5% replications
colorReg{1} = 'r';
colorReg{2} = 'g';
colorReg{3} = 'b';
pcntALL = [0.05 0.10 0.15 0.20];
for pcntj = 1:length(pcntALL)
    pcnt = pcntALL(pcntj);
    for repl = 1:20
        %load results file of NC
        fname = strcat(path,'ODrepliNo_',num2str(repl),'_',num2str(100*pcnt),'_output_',scenarios{1},'.mat');
        load(fname,'outdata')
        % MFDs: Entire network -
        t_MFD = 3/60; %time interval between measurements
        t = ceil(t_MFD/indata.DT); %The step / How many steps k*DT are included in the aggregation period
        sum1 = sum(outdata.x(indata.group1,1:indata.kmax),1); %total network accumulation
        %sum2 = sum(outdata.u,1); %total outflow of all links and virtual queues in the network in every step k [veh/DT]
        %sum5 = sum((outdata.m(indata.group1,1:indata.kmax)*indata.v_ff)./(indata.Links2(indata.group1,3)*ones(1,indata.kmax)),1); %sum5 = total moving flow in all network links
        %sum6 = sum(outdata.w(indata.group2,1:indata.kmax),1); %virtual queues
        sum8 = sum(outdata.q(indata.group3,1:indata.kmax),1); %trip endings
        sum7 = sum(outdata.u(indata.group1,:).*(indata.Links2(indata.group1,3)/1000),1); % production = outflow of network links * link length
        %sum7 = sum(outdata.m(indata.group1,1:indata.kmax)*indata.v_ff/1000,1); %production = moving*free_flow_speed
        
        %congested = outdata.x(:,1:end-1)>0.85*indata.capacity; %links above 85% of capacity
        %congested = congested.*(indata.Links2(:,2).*indata.Links2(:,3)/1000); %(lane x km)
        %congested = sum(congested);
        %tot_occupancy = sum1/sum(indata.capacity(indata.group1));
        %ntw_het = std(outdata.x(indata.group1,1:end-1));
        %n_occ_std = std(outdata.x(indata.group1,1:indata.kmax)./indata.capacity(indata.group1));
        
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
            %         k_MFD = [k_MFD i*indata.DT]; %time for MFD [hours] - intervals of aggregation
            accumulation = [accumulation mean(sum1(i-t+1:i))]; %(no VQs included)
            %         outflow = [outflow mean(sum2(i-t+1:i))]; %[veh/h]
            %         vqueues = [vqueues mean(sum6(i-t+1:i))];
            production = [production mean(sum7(i-t+1:i))]; %[veh*km/h]
            trip_endings = [trip_endings mean(sum8(i-t+1:i))]; %[veh/h]
            %         m_congested = [m_congested mean(congested(i-t+1:i))];
            %         t_occupancy = [t_occupancy mean(tot_occupancy(i-t+1:i))];
            %         het = [het mean(ntw_het(i-t+1:i))];
            %         std_n_occ = [std_n_occ mean(n_occ_std(i-t+1:i))];
        end
        
        % print the entire network figures
        
        % trip endings
        figure(ff2)
        scatter([0 accumulation],[0 movmean(trip_endings,xx)],20,'o','filled','MarkerFaceColor','#0072BD')
        %plot(accumulation,movmean(trip_endings,xx))
        
        % production
        figure(ff1)
        pp = movmean(production,xx);
        pp(1) = 0;
        %scatter(accumulation(1:42),pp(1:42),20,'o','filled','b')
        
        %add extra poins by linear interpolation for 10 first points 
        a1 = pp(1:10);
        a2 = conv(a1,[0.5 0.5],'valid');
        a3 = accumulation(1:10);
        a4 = conv(a3,[0.5 0.5],'valid');
        accumFin = [a3 a4 accumulation(11:42)];
        ppFin = [a1 a2 pp(11:42)]; 
        scatter(accumFin, ppFin,20,'o','filled','b')
        
        %MFDs per region (in common figure)
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
%             %MFD n-e
%             figure(ff1)
%             subplot(f1(1,i))
%             tt = movmean(trip_endings,xx);
%             tt(1) = 0;
%             plot(accumulation,tt,'Linewidth',1.5)
            
            %MFD n-p
            figure(ff3)
            pp = movmean(production,xx);
            pp(1) = 0;
            %scatter(accumulation, pp,'o','filled',colorReg{i})
            %scatter(accumulation(1:42), pp(1:42),20,'o','filled',colorReg{i})
            
            a1 = pp(1:10);
            a2 = conv(a1,[0.5 0.5],'valid');
            a3 = accumulation(1:10);
            a4 = conv(a3,[0.5 0.5],'valid');
            accumFin = [a3 a4 accumulation(11:42)];
            ppFin = [a1 a2 pp(11:42)];
            
            scatter(accumFin, ppFin,20,'o','filled',colorReg{i})
            
            %plot(accumulation, production,'Linewidth',1.5,'LineStyle',lineTypeScen{scen})
            %scatter(accumulation, production,'filled','MarkerFaceColor',colReg{i})
            
        end
               
    end
    
end

figure(ff3)
legend('region 1','region 2','region 3','Location','southeast')



