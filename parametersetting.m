function [] = parametersetting()
%parametersetting: Creates a file 'parameters.mat' with all parameters of the run 

%% General Simulation parameters 
DT = 1/3600;                 % [hours] : time step for the simulation
Tmax = 8;                    % [hours] : Total simulation time
vehlength = 5;               % [metres]: average vehicle length for capacity computation
v_ff = 25*1000;              % [metres/h]: Speed that vehicles are supposed to have in the moving part of the road 
incr_factor = 2.01;          % Global Multiplier to scale demand OD matrices (In accordance with Aimsun - to scale the OD matrix) - NETWORK DEPENDENT 
kmax = ceil(Tmax/DT);        % Total number of time steps (index of last time step of simulation)
ksi = 1;                     % Average number of passengers per private car [passengers]
t_win_cycles = 1;            % Time window for the calculation of speed in a link (in times of the maximum cycle in the network)- for the objective function of DBL problem 
warm_up = 15/60;             % warm-up time (hours) - must have an OD for warm-up / zero if not
updateTR = 1;                % (0 = use constant Aimsun turn rates, 1 = update turn rates regularly based on shortest paths, 2 = use constant turn rates from shortest path of empty network) 
t_win = 15/60;                % [hours] Time interval duration for each value of turning rates (when we move from one column to the next in the junct.turn matrix) [hours]
applyMPfreq = 1;             % frequency of cycles to update traffic signal settings for MP (e.g. every 1 cycle) 
minStageDur = 7;             % Minimum duration of stages to be considered for MP [sec]
gn_R = 5;                    % Maximum absolute change of green time of one phase between cosecutive cycles 
c_cycle = 90/3600;           % Control cycle (hours) 
c_int = ceil(c_cycle/DT);    % No of steps per control cycle (both for MP and PC)
PCalpha = 0.4;               % parameters for PC - queue balancing     
PCbeta = 0.9;                % parameters for PC - queue balancing 

%create demand pattern (OD multiplier per time step)
defac = demandFactor(DT,kmax,incr_factor); %Edit this function

% Bus parameters
% busStatus = 0;                % 0: no bus system considered, 1: bus system considered 
% vbus_ff = 25*1000;            % [metres/h]: Free flow speed of buses (in DBLs and in empty links) 
% busCapacity = 70;             % Passengers per single bus (Used for defining maximum shift to buses) 
% maxAllBoar_single_rate = 0.3; % Percentage of passengeres assumed to board and alight at every bus stop (average)  - ATTENTION: I need to further think about that 
% passengers_moving = 0;        % passengers changing mode at the beginning (from car to bus): (0, unless we check intermediate simulations)  
% tripL_c = 1.7194;             % Average car trip length [km] = (total_dist_travelled(aimsun)/no of cars
% avg_bus_trip_length = 1.4387; % Average bus trip length [km] = (total_bus_dist-aimsun)/no_of_buses : avg distance per single bus 
% tripL_b= 0.5*tripL_c;         % Average bus passenger trip length [km] - from calibration 

%Mode Shift 
% max_count_iter = 50;         % Max number of iterations for Mode choice adjustment 
% NoOfVehicles = 64545;        % No of vehicles - Used for the Mode choice adjustment 
% NoOfPass = 22589;            % No of passengers - Used for the Mode choice adjustment                            
% accepted_move_ratio = 0.001; % Accepted 0.1% of total passengers oscillating between car and bus in the iterative mode choice process - Used for the Mode choice adjustment

save('parameters')

end

