function [indata] = prepareInput(demandCode)
%prepareInput: Prepare input variables for the simulator 

% ----load settings and initialize network 
parametersetting() %adjusted
load parameters.mat
Initialization(demandCode) %adjusted 
load FinalInput.mat

% ----simulation parameters 
indata.ksi = ksi; 
indata.vehlength = vehlength; 
indata.v_ff = v_ff;
indata.kmax = kmax; 
indata.DT = DT;
indata.defac = defac; %demand scaling and pattern 
indata.WU = round(warm_up/DT); %first x steps to use the Warm-Up matrix 
% indata.linksGenAttract = linksGenAttract;
indata.t_win = t_win; %interval for turn rates use 
indata.gn_R = gn_R; %max greeen time change (MP and PC) 
indata.specialInt = specialIntersections; %list of intersections to exclude 
%from MP application
indata.c_int = c_int; %no of time steps composing one control cycle  

% -----demand 
indata.demandGroup2 = demandGroup2;

% ----- About bus lanes: 
% indata.dec_var = dec_var; 
% indata.f = f;
% indata.greentimes = greentimes; 
% indata.DBLplan = zeros(size(dec_var,1),1); %zero lanes reserved 

% ----- MP
indata.MP = MP; 
indata.applyMPfreq = applyMPfreq;
indata.stageDur = stageDur;
indata.stageVarDur = stageVarDur;
indata.minStageDur = minStageDur;

%Turn rates update and PC
indata.updateTR = updateTR; 

% ---- For turn rates update process(comment)
indata.connCounters = connCounters;  %bins for counting vehs at connections 
indata.OD_links = OD_links;  %
indata.stack_OD = stack_OD;  % stack of remaining trips 
indata.A = A; % adjacency matrix  
indata.upstr2 = upstr2; % upstream links 
indata.LinksP = LinksP; %with virtual plus high node vlinks (all) 
indata.downstrP = downstrP; 
indata.upstrP = upstrP;
indata.empty_connCounters = empty_connCounters;
indata.OD_links_2 = OD_links_2; 

% ----clustering for PC 
indata.no_reg = no_reg; 
indata.no_adjReg = no_adjReg; 
indata.reInd = reInd; %indices of links per region 
indata.nodereg = nodereg; %indices of nodes per region 

% ----network characteristics 
indata.greentimes2 = greentimes2; 
indata.group1 = group1;
indata.group2 = group2; 
indata.group3 = group3; 
indata.junct = junct;
indata.junct2 = junct2; %with virtual
%indata.Links = Links;
indata.Links2 = Links2; %with virtual 

indata.NLinks = NLinks; 
indata.NLinks2 = size(Links2,1); %with virtual 
indata.Nodes = Nodes;
% indata.NPairs = NPairs;
indata.NPairs2 = NPairs2;
indata.satFlow = (indata.Links2(:,2))*1800; %[veh/h]
indata.capacity = capacity; 
indata.max_n = max_n;

%set parameters for PC (queue balancing) - the thetas of the 2nd part of PC
indata.alpha = PCalpha;
indata.beta = PCbeta;
   
%save to load directly if necessary to avoid rerunning this function
save('indata.mat','indata')
end

