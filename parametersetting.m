function [] = parametersetting()
%parametersetting: Creates a file 'parameters.mat' with all parameters of the run 

%% General Simulation parameters 
DT = 1/3600;                 % [hours] : time step for the simulation
Tmax = 6;                    % [hours] : Total simulation time
vehlength = 5;               % [metres]: average vehicle length for capacity computation
v_ff = 25*1000;              % [metres/h]: Speed that vehicles are supposed to have in the moving part of the road 
incr_factor = 1.5;           % Global Multiplier to scale demand OD matrices
ksi = 1;                     % Average number of passengers per private car [passengers] - (used to calculate VHT of cars)
t_win_cycles = 1;            % Time window for the calculation of speed in a link (in times of the maximum cycle in the network)- for the objective function of DBL problem 
warm_up = 15/60;             % Warm-up time (hours) - (must have an OD for warm-up / zero if not) 
updateTR = 1;                % (0 = use constant turn rates (read from file), 1 = update turn rates regularly based on shortest paths, 2 = use constant turn rates from shortest path of empty network) 
t_win = 15/60;               % [hours] Time interval between updating of turning rates (when we move from one column to the next in the junct.turn matrix) [hours]

% Create dynamic demand pattern (original OD matrix multiplier as function of time step)
kmax = ceil(Tmax/DT);        % Total number of time steps (index of last time step of simulation) (NO CHANGE)
defac = demandFactor(DT,kmax,incr_factor);  % Edit this function as preffered 


% Adaptive traffic signals
applyMPfreq = 1;                 % frequency of cycles to update traffic signal settings for MP (e.g. every 1 cycle) 
minStageDur = 3;                 % Minimum duration of stages to be considered for MP [sec] (3 if we consider AR stages for dynamic PC / 7 if not) 
gn_R = 5;                        % Maximum absolute change of green time of one phase between consecutive cycles 
c_cycle = 90/3600;               % Average control cycle [hours]
newPC_cycle = 30*60/3600;        % Interval for updating clustering and PC boundaries (for dynamic PC) [hours] 
c_int = ceil(c_cycle/DT);        % No of simulation steps per control cycle (both for MP and PC)
newPC_int = ceil(newPC_cycle/DT);% No of simulation steps per cycle of PC boundary and clustering udpate 
PCalpha = 0.4;                   % weight parameters for PC - queue balancing (1st term -> closeness to controller output)     
PCbeta = 0.9;                    % weight parameters for PC - queue balancing (2nd term -> assignment based on observed queues) 

%For PC clustering 
no_adjReg = 4;                % number of region adjacencies for first PC (1-2, 2-1, 2-3, 3-2) - (for dynamic PC, it will be read from the files, ignore here)
init_clustering = 1;            % 0 if not initial clustering, 1 if yes 

save('parameters')

end

