%% Main program to run tests for the combined MP+PC controller

clear
close all
clc

%% [1] Set here the names of the scenarios 
% (used for creating the related files/ inputs and results)

scenarios{1} = 'NC'; % NC = no adaptive control (fixed-time) 
% MP scenarios
scenarios{2} = 'MP_0';    %max_pressure in all nodes 
% PC scenarios
scenarios{3} = 'PC'; 
scenarios{4} = 'DynPC';

scenarios{5} = 'MP00'; 
scenarios{6} = 'DynPC_MP_0';

% scenarios contains all scenario names 
save('scenarios.mat','scenarios')


%% [2] Run Settings 

jj = [1 2 3]; % SET the scenarios ID that you want to run here, e.g. [1 3 5] for scenarios{1}, scenarios{3} and scenarios{5} 
demandCode = 1; %1 med, 2 high (defines which input files (demand, PC/MP settings etc.) 
PCcode = 0; %1 if PC static is included, 0 if no PC, 2 if PC dynamic 
MPcode = 100; %1 if specific MP is included, 0 if no MP, 100 if MP all, 2 if random MP  
MPnodeSelection = 1; %1 for the proposed methodology, 2 for the benchmark (ignore if no MP) 
dataChanged = 1; %set as 1 if changes require for indata to be modified/re-generated (to be on the safe side set it 1, 0 is to save on resources)
redCapacity = 0; % set 0 for no reduced exit rate, 1 for reduced exit rate (requires initial clustering) 
savefull = 1; %1: save big result files, 0: don't save files 
dim_demand = 0; %1 if the demand from dimitris' paper is used, 0 otherwise
%% [3] Set Paths for file loading/saving [Put your paths or for current directory: strcat(cd,'/'); ] 

% paths for existing input files (per demand scenario: med/high) ready to load: 
path_in_med = '/Users/tsitokas/Documents/GitHub/SaF-LUTS/scenarios dimitris/medium demand/'; 
path_in_high ='/Users/tsitokas/Documents/GitHub/SaF-LUTS/scenarios dimitris/high demand/';

% for PC parameters input files (set filename + path to load directly the PC settings)
inputPCfile_med = strcat(cd,'/inputPCscenarios/input_PC_50b.mat'); 
inputPCfile_high = strcat(cd,'/inputPCscenarios/input_PC_62.mat');

% Set the file name of the results of the reference case (base) for performing the MP node
% selection (we use the NC case)
reffname_med = strcat(cd,'/results/output_final_med_NC.mat');
reffname_high = strcat(cd,'/results/output_final_high_NC.mat');

% Set paths where to save results of the current simulation cases (already created the results folder): 
path_out_med = strcat(cd,'/results/');
path_out_high = strcat(cd,'/results/');

%% [4] Executing the requested runs

for k = 1:length(jj)
    j = jj(k); %the current scenario index
    scenarios{j} %display on screen the current scenario 
    inputPCfile = []; 
    if demandCode == 1
        if MPcode == 1
            % %settings for MP node selection algortihm
            a = 0.5; %med
            b = 1.5; %med 
            equation = 3; %selection equation (3 for med, 2 for high)
            reffname = reffname_med; 
        end
        fprintf('medium demand')
        pathInput = path_in_med; 
        
        if PCcode >=1 
            % PCfile (to load directly) 
            inputPCfile = inputPCfile_med; 
        end
        % output file name to save after simulation is done
        fname = strcat(path_out_med,'output_final_med_',scenarios{j},'.mat');

    elseif demandCode == 2
        if MPcode == 1
            %settings for MP node selection algortihm
            a = 0.6; 
            b = 0.2; 
            equation = 2;
            reffname = reffname_high; 
        end
        fprintf('high demand')       
        pathInput = path_in_high; 
        if PCcode >=1 
            inputPCfile = inputPCfile_high; 
            
        end
        %name of the results file 
        fname = strcat(path_out_high,'output_final_high_',scenarios{j},'.mat');

    end
    
    % load network information (indata block) 
    if dataChanged == 1
        prepareInput();  %run when there are changes in the network/demand etc. 
    end
    %else %load network input file 
    % load(strcat(pathInput,'indata.mat'),'indata')
    
 
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
      
    % load experiment settings

    load(strcat('indata.mat'),'indata'); %here it loads the indata of my previous scenarios 
    
    if dim_demand 
        % add missing elements - 
        % for NC I need the initial clustering equal to 1 (for the capacity
        % drop) 
        indata_dim = load(strcat(pathInput,'indata.mat'),'indata'); %here it loads the indata of my previous scenarios 
         
        indata.clustering = 1; 
        indata.demandGroup2 = indata_dim.indata.demandGroup2; 
        indata.defac = indata_dim.indata.defac; 
        indata.kmax = indata_dim.indata.kmax;  %8 hours simulation 
        indata.OD_links = indata_dim.indata.OD_links;
        indata.OD_links_2 =  indata_dim.indata.OD_links_2;

        clear indata_dim
    end

    indata.alpha = 0.4; 
    indata.beta = 0.9;
    indata.redCapacity = redCapacity;
    indata.demandCode = demandCode; 

    %for running the old indata files from prev paper
    % indata.clustering = 0; 
    load(strcat('input_',scenarios{j},'.mat'))
     
    % execute simulation experiment / print TTT on screen
    outdata = SaF_3_dynPC(indata,MPnodes,MPcode,PC,fname,scenarios{j},savefull);
    %outdata = SaF_3(indata,MPnodes,PC,fname,savefull);
    if demandCode == 1 
        demand = 'med';
    else
        demand = 'high';
    end
    %% [5] To produce some results compared to NC, run the following for each case: 
    resultAnalysis(j,demand) 

end



