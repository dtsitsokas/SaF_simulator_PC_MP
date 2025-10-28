%Create complementary figures based on a set of runs 
clear 
close all
clc
% Load settings 
% parametersetting() %adjusted
% load parameters.mat
% Initialization() %adjusted
 
load('parameters.mat','kmax','DT')
%prepare input variables:
[indata,DBLplan,demandGroup2] = prepareInput();

%% Load results of experiment
scenario = 'MP_ALL';
load(strcat('output_',scenario,'.mat'),'outdata')
load('input_PC.mat')
% x axis - Time 
x = (1:indata.kmax)*indata.DT;

 
%% calculate regional accumulations 
for i=1:3 
    y_i{i} = diff(sum(outdata.x(indata.reInd{i},:),1));
    figure
    plot(x,y_i{i},'LineWidth',1)
    title('Differences $n_i(k)-n_i(k-1)$ vs time')
    xlabel('time (h)');
    ylabel('difference (veh)');
    box on
    legend(strcat('reg ',num2str(i)))
end
%legend('reg 1','reg 2','reg 3') 
%% 
f11 = figure('Name','Differences n_i-n_set');
title('Differences $n_i(k)-n_i^*$ vs time')
xlabel('time (h)');
ylabel('divergence (veh)');
box on
hold on
x = (1:indata.kmax)*indata.DT;
% Calculate differences in n_i from set points 

for i=1:3 
    y_i = sum(outdata.x(indata.reInd{i},1:end-1));
    y_i = (y_i - PC.n_set(i));%/PC.n_set(i); 
    plot(x,y_i,'LineWidth',1)
end
legend('reg 1','reg 2','reg 3') 

%% Plot the values of green that the PI controller would give for these simulation results 
% Test different K_I, K_P parameter values
agg_u = zeros(4,402);
applied_u = zeros(4,402);
%calculate aggregated accumulations per region over time  
for j=1:3 
    agg_n(j,:) = sum(outdata.x(indata.reInd{j},:),1); 
end

%max accumulation per region 
max_n = zeros(1,indata.no_reg);
for i=1:indata.no_reg
    max_n(i) = sum(outdata.capacity(indata.reInd{i}));  
    weight_n(i) = max(agg_n(i,:))/max_n(i);
end

%calculate differences between consecutive cycles 
diff_n = zeros(3,401);
diff_sp = zeros(3,401);

for j=1:3
    i=0;
    for k_c = (indata.c_int+1):indata.c_int:indata.kmax
        i = i + 1;
        diff_n(j,i) = (agg_n(j,k_c)-agg_n(j,k_c-indata.c_int))/max_n(j);
        diff_sp(j,i) = (agg_n(j,k_c)-PC.n_set(j))/PC.n_set(j);
    end
end
aa = 30000;
PC.K_P = [10 -10 0; 0 -10 10; -10 10 0;0 10 -10]*100; 
PC.K_I = [aa -aa*0.8 -aa*0.33; aa*0.33 -aa*0.8 aa*0.4; -aa*0.8 aa 0; 0 aa*0.8 -aa*0.4].*weight_n; 

%% set the initial values for the controller
greens_P = []; 
for i=1:indata.no_adjReg
    ind = PC.indices(i,PC.indices(i,:)>0); %all nodes that belong to this approach i-j (non-zero) - the k indices to refer to each node in PC structs
    for j = 1:length(ind)
        greens_P = [greens_P PC.stageDur{ind(j)}(PC.stagesInvolved(ind(j),1))]; %the greens of all the phases of the set of nodes controlling the i-j
    end
    applied_u(i,1) = mean(greens_P); 
end

applied_u(:,2) = applied_u(:,1); 
agg_u(:,1:2) = applied_u(:,1:2); 

limit_ = 400;
a = PC.sum_greensPC;
a(a==0) = NaN;
max_greens = ones(4,400).*mean(a,2,'omitnan');
min_greens = ones(4,400)*indata.minStageDur;
time = []; 
i=1;
for k_c = indata.c_int:indata.c_int:limit_*indata.c_int %indata.kmax
    i = i + 1; 
    time = [time k_c*indata.DT];
    %here assume that the previous u is applied without modification 
    agg_u(:,i+1) = applied_u(:,i) + PC.K_P*(diff_n(:,i)) + PC.K_I*(diff_sp(:,i));
    applied_u(:,i+1) = max([min([agg_u(:,i+1) max_greens],[],2) min_greens],[],2) ;
end

%plot figures: 
figure('Name','Time-Series of calculated greens by PI')
hold on 
plot(time(3:limit_),applied_u(:,3:limit_))

%maximum 
plot(time(3:limit_),max_greens(:,3:limit_),'LineStyle','--','LineWidth',1.5)

%minimum
plot(time(3:limit_),min_greens(1,3:limit_),'LineStyle','-.','LineWidth',1.5)
legend('1-2','3-2','2-1','2-3','max 1-2','max 3-2','max 2-1','max 2-3','min')


%% plot the differences 
figure
plot(time(2:limit_),diff_n(:,2:limit_))
xlabel('time')
ylabel('Dn_i/Dt (veh)')
title('Derivative of n_i')
legend('1','2','3')

%plot the diff from SP
figure
plot(time(2:limit_),diff_sp(:,2:limit_))
xlabel('time')
ylabel('n_i-n_set (veh)')
title('Divergence from set point')
legend('1','2','3')

