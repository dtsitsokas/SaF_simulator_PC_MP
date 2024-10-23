%% Main program to run tests for the combined MP+PC controller

clear
close all
clc

%% --- Set the names of the scenarios (used for saving related file names/
% inputs and results)

% for MP scenarios = the number at the end should be indicating the
% percentage of nodes - 0 is 5%, 1 is 10%, etc. 
% NC = no adaptive control (fixed-time) 
scenarios{1} = 'NC';

% --- Max-Pressure control scenarios 
scenarios{2} = 'MP_0';    %max_pressure in all nodes 

%selection of MP with max outflow (benchmark) med&high
scenarios{3} = 'MP_cd_MSR00';
scenarios{4} = 'MP_cd_MSR01';
scenarios{5} = 'MP_cd_MSR02';
scenarios{6} = 'MP_cd_MSR03';
scenarios{7} = 'MP_cd_MSR04';
scenarios{8} = 'MP_cd_MSR05';

% selection of MP with the proposed method (m1,m2,Nc) / med&high
scenarios{9} = 'MP_cd_MS00';
scenarios{10} = 'MP_cd_MS01';
scenarios{11} = 'MP_cd_MS02';
scenarios{12} = 'MP_cd_MS03';
scenarios{13} = 'MP_cd_MS04';
scenarios{14} = 'MP_cd_MS05';

% ---perimeter control scenarios 

scenarios{15} = 'PC_50b'; %high PC
scenarios{16} = 'PC_501b'; %high (PC + MP100)
scenarios{17} = 'PC_62'; %med PC
scenarios{18} = 'PC_620'; %med (PC + MP100)

% --- Combined control scenarios 

%PC + MP selection - method- med OD (based on NC)
scenarios{19} = 'PC_62_MP_MS0';
scenarios{20} = 'PC_62_MP_MS1';
scenarios{21} = 'PC_62_MP_MS2';
scenarios{22} = 'PC_62_MP_MS3';
scenarios{23} = 'PC_62_MP_MS4';
scenarios{24} = 'PC_62_MP_MS5';

%PC + MP selection - benchmark - based on NC results - med OD
scenarios{25} = 'PC_62_MP_MSR00';
scenarios{26} = 'PC_62_MP_MSR01';
scenarios{27} = 'PC_62_MP_MSR02';
scenarios{28} = 'PC_62_MP_MSR03';
scenarios{29} = 'PC_62_MP_MSR04';
scenarios{30} = 'PC_62_MP_MSR05';

%PC + MP selection - based on NC results (high demand)
scenarios{31} = 'PC_50b_MP_MS00';
scenarios{32} = 'PC_50b_MP_MS01';
scenarios{33} = 'PC_50b_MP_MS02';
scenarios{34} = 'PC_50b_MP_MS03';
scenarios{35} = 'PC_50b_MP_MS04';
scenarios{36} = 'PC_50b_MP_MS05';

%PC + MP selection - based on PC results  - High OD
scenarios{37} = 'PC_50b_MP_MSR00';
scenarios{38} = 'PC_50b_MP_MSR01';
scenarios{39} = 'PC_50b_MP_MSR02';
scenarios{40} = 'PC_50b_MP_MSR03';
scenarios{41} = 'PC_50b_MP_MSR04';
scenarios{42} = 'PC_50b_MP_MSR05';

% MP 50% and 75% - benchmark selection
scenarios{43} = 'MP_cd_MSR09';
scenarios{44} = 'MP_cd_MSR14';

% MP 50% and 75% - method selection
scenarios{45} = 'MP_cd_MS09';
scenarios{46} = 'MP_cd_MS14';

% PC + MP 50% and 75% - method selection / medium
scenarios{47} = 'PC_62_MP_MS09';
scenarios{48} = 'PC_62_MP_MS14';

% PC + MP 50% and 75% - benchmark selection / medium
scenarios{49} = 'PC_62_MP_MSR09';
scenarios{50} = 'PC_62_MP_MSR14';

% same but no capacity drop 
scenarios{51} = 'MP_ncd_MS00';
scenarios{52} = 'MP_ncd_MS01';
scenarios{53} = 'MP_ncd_MS02';
scenarios{54} = 'MP_ncd_MS03';
scenarios{55} = 'MP_ncd_MS04';
scenarios{56} = 'MP_ncd_MS05';
scenarios{57} = 'MP_ncd_MS09';
scenarios{58} = 'MP_ncd_MS14';

scenarios{59} = 'PC_50b_MP_MS09';
scenarios{60} = 'PC_50b_MP_MS14';

scenarios{61} = 'PC_50b_MP_MSR09';
scenarios{62} = 'PC_50b_MP_MSR14';

scenarios{63} = 'ncd_NC';

scenarios{64} = 'PC_62_MP_ncd_MS00';
scenarios{65} = 'PC_62_MP_ncd_MS01';
scenarios{66} = 'PC_62_MP_ncd_MS02';
scenarios{67} = 'PC_62_MP_ncd_MS03';
scenarios{68} = 'PC_62_MP_ncd_MS04';
scenarios{69} = 'PC_62_MP_ncd_MS05';

% scenarios contains all scenario names 
save('scenarios.mat','scenarios')



%% Run Settings 

jj = 1; % IDs of scenarios to run 
demandCode = 2; %1 med, 2 high (defines which input files (demand, PC/MP settings etc.) 
PCcode = 0; %1 if PC is included, 0 else
MPcode = 0; %1 if specific MP is included, 0 if no MP, 100 if MP all, 2 if random MP  
MPnodeSelection = 1; %1 for the proposed methodology, 2 for the benchmark (ignore if no MP) 
dataChanged = 1; %set as 1 if changes require for indata to be modified/re-generated 
redCapacity = 1; % set 0 for no reduced exit rate, 1 for reduced exit rate (edit function!) 
savefull = 1; %1: save big result files, 0: don't save files 

% paths where to find input files (per demand scenario): 
%path_in_med = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
%path_in_high = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';

% for PC input (set filename of PC settings to load directly)
inputPCfile_med = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\runs high demand\final runs input files\input_PC_50b.mat';
inputPCfile_high = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\runs medium demand\input_PC_62.mat'; 


%*** load new indata (debugging) 
%*** load new indata (debugging) 
path_in_med = strcat(cd,'\'); 
path_in_high = strcat(cd,'\'); 
path_out_high = strcat(cd,'\'); 
path_out_med = strcat(cd,'\'); 

%set the file name of the results of the reference case for MP node
%selection (use NC case)
reffname_med = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\output_final_med_NC.mat';
reffname_high = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\output_final_high_NC.mat';

% paths where to save results: 
%path_out_med = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
%path_out_high = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';

%%

for k = 1:length(jj)
    j = jj(k); 
    scenarios{j} %display on screen the current scenario 
    
    if demandCode == 1
        if MPcode == 1
            % %settings for MP node selection 
            a = 0.5; %med
            b = 1.5; %med 
            equation = 3; %selection equation (3 for med, 2 for high)
            reffname = reffname_med; 
        end
        fprintf('medium demand')
        pathInput = path_in_med; 
        % PCfile (to load directly) 
        inputPCfile = inputPCfile_med; 
        % output file name to save after simulation is done
        fname = strcat(path_out_med,'output_final_med_',scenarios{j},'.mat');

    elseif demandCode == 2
        if MPcode == 1
            % %settings for MP node selection 
            a = 0.6; 
            b = 0.2; 
            equation = 2;
            reffname = reffname_high; 
        end
        fprintf('high demand')       
        pathInput = path_in_high; 
        inputPCfile = inputPCfile_high;
        fname = strcat(path_out_high,'output_final_high_',scenarios{j},'.mat');

    end
    
    % load network information (indata block) 
    if dataChanged == 1
        prepareInput();  %run when there changes in the network/demand etc. 
    end
    % load network input file 
    load(strcat(pathInput,'indata.mat'),'indata')
    indata.redCapacity = redCapacity;
    indata.demandCode = demandCode; 
    
    if MPcode==1
        case_j = str2double(scenarios{j}(end-1:end)); %filenames indicate node percentage at the end (00 is 5%, 01 is 10% etc)
        %MP node selection
        if MPnodeSelection == 1
            % proposed method 
            [~] = nodeSelectionMP(a,b,case_j,equation,scenarios{j},reffname); 
        elseif MPnodeSelection == 2
            % benchmark (max outflow based) 
            [~] = nodeSelectionMP_outflow(a,b,case_j,equation,scenarios{j},reffname); 
        end
    end
    
    % create input file with all experiment settings (MP, PC, specific scenario) 
    prepareExperiment(scenarios{j},inputPCfile,PCcode,MPcode)
        
    %load experiment settings
    %pathInp = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';
    %load(strcat(pathInp,'input_',scenarios{j},'.mat'))
    load(strcat('input_',scenarios{j},'.mat'))

    % execute simulation experiment / print TTT on screen 
    outdata = SaF_3(indata,MPnodes,PC,fname,savefull);
    %close all 

end

%resultAnalysis()

%% Call simulator for random sets 
% path0 = 'Y:\common\Dimitris\Max Pressure\fourth set results\ResultsMediumOD_noCapacityDrop\Random Node Sel Results\Random node sets\';  
% %path1 = 'Y:\common\Dimitris\Max Pressure\hEART Results\Random Node Sel Results\capacityDrop';
% %path1 = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\fifth set results\ResultsMediumOD_capacityDrop\RandomSetsResults\';
% path1 = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\Current results\ResultsHighOD_capacityDrop\RandomMPonly_Results\';
% 
% for p = 0.05:0.25:0.30
%     for r = 1:10
%         PC.mode=0; %case No PC  
%         fname = strcat(path1,'output_MP_R_',num2str(p*100),'_',num2str(r),'.mat');
%         load(strcat(path0,'RandomNodeSet_',num2str(p*100),'_',num2str(r),'.mat'),'rMPnodes')
%         SaF_3(indata,rMPnodes,PC,fname);
%         disp(strcat('Replication no ',num2str(r),' finished.'))
%     end
% end

%% Call sim for PC + random sets 
% % location of random sets
% % path0 = 'Y:\common\Dimitris\Max Pressure\fourth set results\ResultsMediumOD_noCapacityDrop\Random Node Sel Results\Random node sets\';  
% path0 = 'F:\Max Pressure files\fifth set results\Random node sets\';
% 
% %location to save results files 
% path1 = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\RandomSetsResults\MP plus PC - reduced results\';
% 
% demandCode = 1; %1 med, 2 high
% savefull = 2;
% livePlots = 0; 
% for p = 0.10:0.05:0.25
%     for r = 1:10
%         rMPcasename = strcat(path0,'RandomNodeSet_',num2str(p*100),'_',num2str(r),'.mat');
%         prepareExperiment_2(num2str(p*100),r,rMPcasename,demandCode)
%         
%         if demandCode == 1 %med
%             pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
%         elseif demandCode == 2
%             pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
%         end
%         
%         % indata = prepareInput();
%         load(strcat(pathInput,'indata.mat'),'indata')
% 
%         %load scenario settings
%         load(strcat('input_PC_62_MP_R',num2str(p*100),'_',num2str(r),'.mat'))
%     
%         fname = strcat(path1,'output_PC_62_MP_R',num2str(p*100),'_',num2str(r),'.mat');
%         SaF_3(indata,MPnodes,PC,fname,savefull,livePlots);
%         disp(strcat('Replication no ',num2str(r),'of case ',num2str(100*p),' finished.'))
%     end
% end

%% Run a specific case for different ODs (demand as random variable with normal distribution) 

%% Save the multiple demands to use the same for all experiments 

% jj = [1 33]; %which scenario cases to run
% for pcnt = 0.05:0.05:0.20
%     for k=1:length(jj)%length(scenarios)
%         j = jj(k);
%         scenarios{j}
%         
%         %prepare Experiment
%         demandCode = 2; %1 med, 2 high
%         %     a = 0.5;
%         %     b = 1.5;
%         %     case_j = j-jj(1);
%         %     equation = 3; %selection equation
%         %pcnt = 0.20; %percentage of mean as standard dev for OD demand
%         
%         if demandCode == 1 %med
%             pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\medium demand input\';
%             %output file name to save after simulation is done
%             path = 'F:\Max Pressure files\fifth set results\ResultsMediumOD_capacityDrop\';
%             pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs medium demand\';
%         elseif demandCode == 2
%             pathInput = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\input files\high demand input\';
%             path = 'F:\Max Pressure files\fifth set results\ResultsHighOD_capacityDrop\';
%             pathIn = 'C:\Users\tsitsoka\Documents\Dimitris Local Documents pc10\GitHub-My repositories\MaxPressure_Codes-\runs high demand\final runs input files\';
%         end
%         
%         %load scenario settings
%         load(strcat(pathIn,'input_',scenarios{j},'.mat'))
%         % indata = prepareInput();
%         for repl = 1:20
%             repl
%             load(strcat(pathInput,'indata.mat'),'indata')
%             
%             if j==1
%                 %generate OD (draw from normal)
%                 indata.demandGroup2(:,2) = generateNormDistrOD(indata.demandGroup2(:,2), pcnt);
%             else
%                 %OR
%                 %load OD already generated
%                 if demandCode == 1
%                     load(strcat('RepDemand_',num2str(100*pcnt),'.mat'))
%                 else
%                     load(strcat('RepDemand_high_',num2str(100*pcnt),'.mat'))
%                 end
%                 
%                 indata.demandGroup2(:,2) = RepDemand(:,repl);
%             end
%             
%             %MPnodes = nodeSelectionMP(demandCode,a,b,case_j,equation,scenarios{j});
%             %prepareExperiment(scenarios{j})
%             
%             fname = strcat(path,'ODrepliNo_',num2str(repl),'_',num2str(pcnt*100),'_output_',scenarios{j},'.mat');
%             savefull = 1; %save big result files
%             SaF_3(indata,MPnodes,PC,fname,savefull);
%             
%             % SaF_3subsets(indata,MPnodes,PC,fname); %with subsets of links
%         end
%         
%         
%         if j==1
%             % Save the generated demands
%             
%             for repl = 1:20
%                 load(strcat(path,'ODrepliNo_',num2str(repl),'_',num2str(pcnt*100),'_output_',scenarios{1},'.mat'),'indata')
%                 RepDemand(:,repl) = indata.demandGroup2(:,2);
%             end
%             
%             % save('RepDemand_5','RepDemand')
%             % save('RepDemand_10','RepDemand')
%             % save('RepDemand_15','RepDemand')
%             if demandCode == 1
%                 save(strcat('RepDemand_',num2str(pcnt*100)),'RepDemand')  
%             else
%                 save(strcat('RepDemand_high_',num2str(pcnt*100)),'RepDemand')
%             end
%         end
%         
%     end
% end
% 

