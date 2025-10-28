%Print PC related figures 
%Perform an analysis: What happens to the system when PC is applied 
close all 
clear 
clc

%load input and output file of the specific case 
filename_0 = 'output_PC_0.mat'; 
% filename_0b = 'input_PC_0.mat';
% filename_0 = 'output_NC.mat'; 
filename_0b = 'input_PC_0.mat';

load(filename_0,'outdata','indata') %No control case
load(filename_0b,'PC'); %No control case

%the pretimed greens - 
u_0 = ones(size(outdata.applied_u)).*outdata.applied_u(:,1);  

%x-axis: time 
time = 1:size(outdata.x(1,1:end-1),2);%(1:indata.kmax);
time = time*indata.DT;


%%  Time series of the controller decisions during the last run 
% plot min/max, calculated/applied and pre-timed green over time

figure('Name','Controller decision');
hold on

x = time(1:90:end);
%y-axis: applied green times 
% the pretimed greens 
plot(x,u_0(:,2:end))
xlabel('time (h)')
ylabel('greens (sec)')
title('Min, Max and pre-timed greens')
a = PC.sum_greensPC;
a(a==0) = NaN;

% maximum greens
plot(x,ones(1,length(x)).*mean(a,2,'omitnan'),'LineStyle','--','LineWidth',1)

% minimum greens
plot(x,ones(1,length(x))*indata.minStageDur,'LineStyle','--','LineWidth',1)
legend('1-2','3-2','2-1','2-3','max 1-2','max 3-2','max 2-1','max 2-3','min')

% plot the applied greens from the simulation 
plot(x,outdata.applied_u(:,1:end-1),'LineStyle','-.','LineWidth',2.5)
% legend('1-2','3-2','2-1','2-3','max 1-2','max 3-2','max 2-1','max 2-3','min','u 1-2','u 3-2','u 2-1','u 2-3')

% The calculated greens from the controller 
plot(x,outdata.agg_u(:,1:end-1),'LineStyle',':','LineWidth',1.5)
legend('1-2','3-2','2-1','2-3','max 1-2','max 3-2','max 2-1','max 2-3','min','u 1-2','u 3-2','u 2-1','u 2-3','c 1-2','c 3-2','c 2-1','c 2-3')
saveas(gcf,'control_outputs.emf');
%% Check the inflows and outflows of each region specifically 

% Plot time series of: 
% - Inflows from Centroids (perimeter and inside)
% - Inflows from neighboring regions 
% - Outflows to centroids 
% - Outflows to neighboring regions
% - Accumulation of region 
% - Error (n-n_set) - distance from target
% - product 

i = 1; %region 1
figure('Name',strcat('TS Outflows region ',num2str(i)));
hold on

%inflow from centroids (new demand - not controllable)
q_cent = sum(outdata.u(intersect(indata.reInd{i},indata.group2),:)*indata.DT);
plot(time(1:90:end),q_cent(1:90:end),'LineWidth',2)
title(strcat('region ',num2str(i)))
xlabel('time (h)')
ylabel('flow (veh/sec)')

%inflow from neighboring region - sum upairs of connections i-j
ind_1 = indata.Links2(indata.junct2.or_index,6)==2; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==1; %destination = region 1
ind = ind_1 + ind_2 == 2; 
q_2_1 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), q_2_1(1:90:end),'LineWidth',2); 
legend('inflow centr','inflow 2-1')
saveas(gcf,strcat('reg_',num2str(i),'inflows.emf'));

figure('Name',strcat('TS Outflows region ',num2str(i)));
hold on
title(strcat('region ',num2str(i)))
%outflows to centroids (trip endings inside region) 
u_cent = sum(outdata.q(intersect(indata.reInd{i},indata.group3),:))*indata.DT; %inflows to destination links of the region 
plot(time(1:90:end),u_cent(1:90:end),'LineWidth',2,'LineStyle','-.');

%outflows to neighboring regions 
ind_1 = indata.Links2(indata.junct2.or_index,6)==1; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==2; %destination = region 1
ind = ind_1 + ind_2 == 2; 
u_1_2 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), u_1_2(1:90:end),'LineWidth',2,'LineStyle','-.'); 
legend('ourflow centr','outflow 1-2')
xlabel('time (h)')
ylabel('flow (veh/sec)')
saveas(gcf,strcat('reg_',num2str(i),'outflows.emf'));

%% Check variables for the controller 
figure('Name','TS accumulation');
%accumulation of region 

n_reg = sum(outdata.x(indata.reInd{i},:));
plot(time(1:90:end), n_reg(1:90:end-1)); 
xlabel('time (h)') 
ylabel('control error (\%)')

%control error 
err_i = n_reg - PC.n_set(i);
err_i = err_i/PC.n_set(i); 
plot(time(1:90:end), err_i(1:90:end-1))

%product 
yyaxis right 
prod_i = diff(n_reg); 
plot(time(1:90:end-1),prod_i(1:90:end))
ylabel('dn/dt')

title(strcat('region\hspace{0.2cm} ',num2str(i)))

saveas(gcf,strcat('reg_',num2str(i),'control_variables.emf'));

%% region 2 

i = 2; %region 1
figure('Name',strcat('TS Inflows region ',num2str(i)));
hold on
title(strcat('region ',num2str(i)))
%inflow from centroids (new demand - not controllable)
q_cent = sum(outdata.u(intersect(indata.reInd{i},indata.group2),:)*indata.DT);
plot(time(1:90:end),q_cent(1:90:end),'LineWidth',2)
title(strcat('region ',num2str(i)))
xlabel('time (h)')
ylabel('flow (veh/sec)')

%inflow from neighboring region - sum upairs of connections i-j
ind_1 = indata.Links2(indata.junct2.or_index,6)==1; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==2; %destination = region 1
ind = ind_1 + ind_2 == 2; 
q_1_2 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), q_1_2(1:90:end),'LineWidth',2); 

%inflow from neighboring region - sum upairs of connections i-j
ind_1 = indata.Links2(indata.junct2.or_index,6)==3; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==2; %destination = region 1
ind = ind_1 + ind_2 == 2; 
q_3_2 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), q_3_2(1:90:end),'LineWidth',2); 
legend('inflow centr','inflow 1-2','inflow 3-2')
saveas(gcf,strcat('reg_',num2str(i),'inflows.emf'));

f5 = figure('Name',strcat('TS Outflows region ',num2str(i)));
hold on
%outflows to centroids (trip endings inside region) 
u_cent = sum(outdata.q(intersect(indata.reInd{i},indata.group3),:))*indata.DT; %inflows to destination links of the region 
plot(time(1:90:end),u_cent(1:90:end),'LineWidth',2,'LineStyle','-.');
title(strcat('region ',num2str(i)))
xlabel('time (h)')
ylabel('flow (veh/sec)')

%outflows to neighboring regions 
ind_1 = indata.Links2(indata.junct2.or_index,6)==2; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==1; %destination = region 1
ind = ind_1 + ind_2 == 2; 
u_2_1 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), u_2_1(1:90:end),'LineWidth',2,'LineStyle','-.'); 

ind_1 = indata.Links2(indata.junct2.or_index,6)==2; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==3; %destination = region 1
ind = ind_1 + ind_2 == 2; 
u_2_3 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), u_2_3(1:90:end),'LineWidth',2,'LineStyle','-.'); 
legend('ourflow centr','outflow 2-1','outflow 2-3')
saveas(gcf,strcat('reg_',num2str(i),'outflows.emf'));
%% Check variables for the controller 
f6 = figure('Name','TS accumulation');
%accumulation of region 

n_reg = sum(outdata.x(indata.reInd{i},:));
plot(time(1:90:end), n_reg(1:90:end-1)); 
xlabel('time (h)') 
ylabel('control error (\%)')
title(strcat('region ',num2str(i)))
%control error 
err_i = n_reg - PC.n_set(i);
err_i = err_i/PC.n_set(i); 
plot(time(1:90:end), err_i(1:90:end-1))

%product 
yyaxis right 
prod_i = diff(n_reg); 
plot(time(1:90:end-1),prod_i(1:90:end))
ylabel('dn/dt')

title(strcat('region\hspace{0.2cm} ',num2str(i)))
saveas(gcf,strcat('reg_',num2str(i),'control_variables.emf'));
%% region 3 

i = 3; %region 1
figure('Name',strcat('TS Outflows region ',num2str(i)));
hold on

%inflow from centroids (new demand - not controllable)
q_cent = sum(outdata.u(intersect(indata.reInd{i},indata.group2),:)*indata.DT);
plot(time(1:90:end),q_cent(1:90:end),'LineWidth',2)
title(strcat('region ',num2str(i)))
xlabel('time (h)')
ylabel('flow (veh/sec)')

%inflow from neighboring region - sum upairs of connections i-j
ind_1 = indata.Links2(indata.junct2.or_index,6)==2; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==3; %destination = region 1
ind = ind_1 + ind_2 == 2; 
q_2_3 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), q_2_3(1:90:end),'LineWidth',2); 
legend('inflow centr','inflow 2-3')
saveas(gcf,strcat('reg_',num2str(i),'inflows.emf'));

figure('Name',strcat('TS Outflows region ',num2str(i)));
hold on
title(strcat('region ',num2str(i)))
%outflows to centroids (trip endings inside region) 
u_cent = sum(outdata.q(intersect(indata.reInd{i},indata.group3),:))*indata.DT; %inflows to destination links of the region 
plot(time(1:90:end),u_cent(1:90:end),'LineWidth',2,'LineStyle','-.');

%outflows to neighboring regions 
ind_1 = indata.Links2(indata.junct2.or_index,6)==3; %origin = region 2 
ind_2 = indata.Links2(indata.junct2.dest_index,6)==2; %destination = region 1
ind = ind_1 + ind_2 == 2; 
u_3_2 = sum(outdata.upair(ind,:))*indata.DT;
plot(time(1:90:end), u_1_2(1:90:end),'LineWidth',2,'LineStyle','-.'); 
legend('ourflow centr','outflow 3-2')
xlabel('time (h)')
ylabel('flow (veh/sec)')
saveas(gcf,strcat('reg_',num2str(i),'outflows.emf'));

%% Check variables for the controller 
f3 = figure('Name','TS accumulation');
%accumulation of region 

n_reg = sum(outdata.x(indata.reInd{i},:));
plot(time(1:90:end), n_reg(1:90:end-1)); 
xlabel('time (h)') 
ylabel('control error (\%)')

%control error 
err_i = n_reg - PC.n_set(i);
err_i = err_i/PC.n_set(i); 
plot(time(1:90:end), err_i(1:90:end-1))

%product 
yyaxis right 
prod_i = diff(n_reg); 
plot(time(1:90:end-1),prod_i(1:90:end))
ylabel('dn/dt')

title(strcat('region\hspace{0.2cm} ',num2str(i)))

saveas(gcf,strcat('reg_',num2str(i),'control_variables.emf'));



