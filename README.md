# Store-and-Forward urban traffic simulator with enhanced capabilities and modules for multi-region perimeter control and decentralized Max Pressure control 

Simulator based on Store-and-Forward model with enhanced features, including possibility for Perimeter Control (with static boundaries) and local Max Pressure control. The code is written in MATLAB. It is used in our paper that is cited below. In case you find the code useful for yor research, please cite the following article in your prospective publication: 

Tsitsokas, D., Kouvelas, A., & Geroliminis, N. (2023). Two-layer adaptive signal control framework for large-scale dynamically-congested networks: Combining efficient Max Pressure with Perimeter Control. Transportation Research Part C: Emerging Technologies, 152, 104128. Full article [here](https://www.sciencedirect.com/science/article/pii/S0968090X23001171). 

```
@article{tsitsokas2023two,
  title={Two-layer adaptive signal control framework for large-scale dynamically-congested networks: Combining efficient Max Pressure with Perimeter Control},
  author={Tsitsokas, Dimitrios and Kouvelas, Anastasios and Geroliminis, Nikolas},
  journal={Transportation Research Part C: Emerging Technologies},
  volume={152},
  pages={104128},
  year={2023},
  publisher={Elsevier}
}
```

The code related to this work is partially funded by [Dit4Tram](https://dit4tram.eu/) “Distributed Intelligence & Technology for Traffic & Mobility Management” project from the European Union’s Horizon 2020 research and innovation program under Grant Agreement 953783.

Some brief instructions for using enhanced Store-and-Forward simulator 
----------------------------------------------------------------------
A] How to run a simulation experiment: 

1) Open 'parametersetting.mat' and set the desired parameters for the run 
	* make sure to check the demand dynamic profile in the function 'demandFactor'. The default one is a uniform 2-hour demand with a 15 mins warm-up period. 

2) Make sure all required input files are in the folder and up to date 

3) Open main SaF_3 file 

	* set filenames of the scenarios that you wish to run 
	* set all parameter values in section "Run Settings":

		- in jj set a vector of all "scenario" indices that you need to run, e.g. jj = [1 2 3] to run cases scenarios{1}, scenarios{2}, scenarios{3} 
		- demandCode: 1 for the medium demand case, 2 for the high demand case (to adjust certain parameters of PC/MP and set paths to read/write files)
		- PCcode: 0/1 for PC or not 
		- MPcode: 0/1 for MP or not 
		- MPnodeSelection: 1 or 2 (MP node selection method, ignore if no MP) 
		- dataChanged: if any model parameter or network data has changed (otherwise the "FinalInput" and "parameters" files can be loaded directly) 
		- reducedCapacity: 0/1 to set if the exit rate capacity of the network will be affected by the regional accumulation or not (to change the pattern of drop, edit the ExitServiceFunction)
		- saveFull: 0/1 to set if the results file needs to be generated and saved in the hard drive

	* set all paths/file names for required input files 
	* Run main_SaF3_MP_PC - result file will be saved at the predefined location
	* After simulation is finished, 'resultAnalysis' function will provide a table of results and some figures of the run, comparing with those of a benchmark case (by default the NC case) 
	
-------------------------------------------------------------------------
B] How to change the demand: 

	* Replace the O-D matrix file (OD.txt) for a completely different new demand or load a new OD matrix from mat file (uncomment lines 32-34 in Initialization) 
	* To increase/decrease while keeping the same OD, change the variable incr_factor in "parametersetting.mat" (multiplier of the demand) 
	* To change the dynamic demand profile, update the content of the function in "demandFactor" 
	* To use the demand cases from my paper, set the input file paths accordingly (uncomment lines 39-40 in main) and set dataChanged = 0 
	
C] How to maintain constant turn rates 

	* For constant turn ratios given from Aimsun set updateTR equal to 0 
	* For constant turn ratios based on distance-based shortest path, set updateTR equal to 2 
	* For cdynamic update of turn rates based on time-based shortest path, set updateTR equal to 1 

D] Basic simulation results  

	* Results are all gathered in struct "outdata", which is saved in the results file 
	* Run "ResultAnalysis" to obtain some standard figures describing the run: 

-----------------------------------------------------------------------------

Results file content: 

All results of a simulation experiment are included in the results file: 

* outdata struct: all model-related results of the experiment  
	- outdata.x: [rows: links - columns: time steps] number of cars in link (queues) in [veh] 
	- outdata.m: [same as x] number of moving vehicles in link [veh] 
	- outdata.w: [same as x] number of vehicles waiting in queue at link [veh] 
	- outdata.u: [same as x] outflow rate of link [veh/h] 
	- outdata.q: [same as x] inflow rate of link [veh/h] 
	- outdata.a_z: [same as x] vehicles arriving at the end of the queue of link [veh/h] - the transitioning flow from m to x 
	- outdata.cumtripsCompleted: [rows: exit links - columns: time steps] cumulative vehicles exiting network at exit links [veh] 
	- outdata.cumtripsStarted: [rows: entry links - columns: time steps] cumulative trips started at entry links [veh] 
	- outdata.upair [row: approaches - columns: time steps] outflow of approaches (pairs of upstream-downstream link) [veh/h] 
	- outdata.avg_q_cycle: [rows: links - columns: control cycles] average queue of link within a cycle [veh] 
	- outdata.greentimes2: [approaches x phases x control cycles] green times of all phases per cycle 
	- outdata.notserviced: vehicles remaining inside the network and in virtual queues after the end of the simulation [veh] 
	- outdata.virtualqueues: [rows: entry links - columns: time steps] virtual queues before the enry links [veh] 
	- outdata.MPnodes: indices of selected Max Pressure nodes 
	- outdata.PC: struct with all PC settings of the specific run 
	- outdata.turn: [rows: approaches - columns: time intervals of turn ratios update] the turn ratios used for the simulation (static or dynamic) 
	- outdata.r1: travel time spent inside the network (trip time) [veh x hours] 
	- outdata.r2: travel time spent in  virtual queues (waiting time) [veh x hours] 
	- outdata.r3: total travel time (sum of r1+r2) [veh x hours] 
	
* indata:  struct with all settings/parameters/demand and network information used for the specific experiment 
	- .ksi: number of pax/car (not used in this version, can be deleted) 
	- .vehlength: assumed lenght of vehicle (for storage capacity calculations and queue length estimations) 
	- .v_ff: assumed free flow speed in the moving part of links [km/h] 
	- .kax: number of time steps of one simulation experiment 
	- .DT: time step duration [hours] 
	- .defac: global demand factor multiplier per time step of simulation [vector: length equal to the number of time steps] 
	- .WU: number of time steps of the beginning of simulation where OD_WarmUp applies (duration of warm-up period) 
	- .t_win: duration of time window for turn ratios change [hours] 
	- .gn_R: maximum absolute change of green time of any phase between cosecutive cycles [sec] 
	- .specialInt: indices of intersections where MP cannot be applied in this form (due to phase/lane configuration) - they are exempt from MP 
	- .c_int: number of time steps that compose one control cycle 
	- .demandGroup2: [rows: entry links, col: demand cases = here 2: warm_up, regular] total trip demand per entry link [veh/h] - this one used in the simulator 
	- .MP: struct with several intersection details (signal plans, configureations etc.) used in MP and PC when updating traffic signals 
	- .applyMPfreq: interval between application of MP controller (no. of cycles)
	- .stageDur: FTC duration of phases per intersection [rows: approaches, col: phases (if -10: non applicable)]  
	- .stageVarDur: same as stageDur but only listing phases that can change in the adaptive control (the rest take zero) 
	- .updateTR: code to indicate update process of turn ratios 
	- .connCounters: bins to count vehicles passing between approaches (auxiliary) 
	- .OD_links: connection of ODs to demand (only for real links) - details in Initialization (auxiliary) 
	- .stack_OD: stack of ongoing trips (used in the turn ratio updates) 
	- .A: adjacency matrix for turn ratio update (shortest path calculation)  
	- .upstr2: [linkID - no of upstream links - IDs of upstream links] - includes virtual links (sources, sinks) 
	- .LinksP: [linkID - No of lanes - Length(m) - ID of start node - ID of end node - region where it belongs(clustering of 3)] - includes 2 layers of virtual links reaching to centroids 
	- .downstrP:  [linkID - no of downstream links - IDs of downstream links] - includes 2 layers of virtual links reaching to centroids 
	- .upstrP: same as upstr2 but including 2 layers of virtual links 
	- .empty_connCounters: bins for approaches for counting vehicles in turn ratio updates (empty) - for turn ratio updates calculation  
	- .OD_links_2: OD description matching demand to origin-destination links (including V links) - used in turn ratio update  
	- .no_reg: number of regions 
	- .no_adjReg: number of adjacent regions (for PC) 
	- .reInd: indices of links per region 
	- .nodereg: IDs of nodes per region 
	- .greentimes2: intervals of the cycle between which approach takes green [start_time_1 end_time_1 start_time_2 end_time_2 ...] (1000 means no more intervals) [rows: approaches]  
	- .group1: indices of all physical network links (intermediate links: not entries, not exits) - indices refer to the rows of LinksP or Links2 matrices 
	- .group2: indices of all entry links (virtual queues) 
	- .group3: indices of all exit links (virtual sinks) 
	- .junct: all information of approaches of junctions (intersections) - each variable refers to one approach (connection between origin-destination links) 
	- .junct2: same but including approaches between virtual and physical links as well (the one use in the simulation) 
	- .Links2:  same as Links but including one layer of virtual links 
	- .NLinks: number of physical links 
	- .Nlinks2: number of physical + one layer virtual links 
	- .Nodes: [nodeID - x coordinate - y coordinate] 
	- .NPairs2: number of connections(approaches) between links (including one layer vlinks) 
	- .satFlow: saturation flow per link [veh/h] including virtual links 
	- .capacity: link storage capacity including virtual links [veh] 
	- .max_n: maximum possible accumulation per region (regional capacity) [veh] 
	- .aplha: parameter theta_1 of 2nd layer of PC (for queue balancing) - term relating to being close to the controllers decision 
	- .beta: parameter theta_2 of 2nd layer of PC (for queue balancing) -  term relating to assigning green relative to the current queue 
---------------------------------------------------------------------------------

File description: 

*               main_SaF3_MP_PC(): main file to execute simulation experiments. Simulator is called and settings are made. 
*                  prepareInput(): function that prepares and formats all network-related input -> indata file 
*             prepareExperiment[]: function to create an input file with all experiment settings together -> input_... file 
*              parametersetting(): (called from prepareInput) - input modeling/simulation parameters -> parameters.mat
*                Initialization(): (called from prepareInput) - process network information -> FinalInput.mat 
*     			              SaF_3[]: function to execute one simulation experiment -> returns TTT, generates & saves results file 
*               nodeSelectionMP[]: selection of nodes for MP (proposed method) -> MPnodesCase_... (file wih selected MP nodes)
*       nodeSelectionMP_outflow[]: selection of nodes for MP based on max output -> MPnodesCase_... (file wih selected MP nodes) 
*                  demandFactor[]: (called from parametersetting): introduce dynamic demand profile 
*            createDemandGroup2[]: (called from Initialization): adjusts OD to usable form  
*              correctTurnRates[]: (called from Initialization): applies modifications to turn rates that are given as inputs for running with SaF
*               createPCinputFile: Creates and saves an input file with PC parameters of the PI controller 
*   			             modifyOD(): changes to OD to avoid some errors (obsolete) [called in Initialization]
*      				      createmovie(): creates a video with snapshots of the network maps (old) 
*     			     createSnapshot[]: creates a snapshot of the map of netowork with different colors, for any metric that we need to visualize 
*     				      addYalmipPath: for installing yalmip with matlab  
*     	   gplotdc (all versions): functions used in the creation of visual maps (slight variations between files, some may not be used any more)   
*                     movingAvg[]: function to calculate moving averages for smoother figures (probably not used, the built-in matlab function must be used) 
*   PCcontrolledIntersectionsInfo: information for PC intersections 
* PCcontrolledIntersectionsInfo_2: information for PC intersections (including external intersections) 
*    plotDemandDescriptionFigrues: Plot figures to visualize demand patterns of a specific OD (inter/intra regional, origin/dest centroid densities) 
* 						PrintMaps: produce network map with regional separation and PC nodes visualization 
* 				   resultAnalysis: produce some descriptive figures from the results files (MFDs, accumulation time-series etc, numerical results etc.)


Functions called within SaF_3: 
-                   OutflowS_2[]: Simulation outflows 
-      findFeasibleSignalsPC_3[]: calculate detailed signals based on results of PI controller & applying queue balancing (new) 
-   (old)findFeasibleSignalsPC[]: calculate detailed signals based on results of PI controller & applying queue balancing (old method)
-          ExitServiceFunction[]: calculates the drop in maximum possible exit rate at exit links in association to regional accumulations  
-       updateTurnAndExitRates[]: estimates turn rates based on current traffic conditions (used in empty network =  distance-based shortest paths) 
-     updateTurnAndExitRates_1[]: estimates new turn rates based on current traffic conditions (time-wise shortest paths) 
-                  maxPressure[]: 1st step of MP  
-  feasibleGreensFromPressures[]: 2nd step of MP (opt problem) 
- 						signal[]: checks if an approach has green or not according to the current control plan 
- 

--------------------------------------------------------------------------------

Notes: 

* For running scenarios with PC or MP, yalmip and Gurobi need to be installed in the computer 
* Several input data files are required for the runs (exist in the folder) 
* In folder "scenarios dimitris" there are the input files I used for my 2 demand scenarios (medium, high) - they can be loaded directly in main 
* PC settings for the PI controller calibrated for the demand cases I used can be found in folder "inputPCScenarios" - they can be loaded directly 
