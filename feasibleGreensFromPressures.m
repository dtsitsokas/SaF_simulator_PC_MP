function [feasibleGreens] = feasibleGreensFromPressures(realGreens,minGreens,prevGreens,maxFluct)
%Solve IQP problem to find feasible green times for the involved phases by 
%minimizing difference from pressure-calculated decimal green times 
%Constraints of min green and max fluctuation are imposed. Use yalmip to
%solve the QP 

yalmip('clear')
%Define variables: the feasible greens of the phases involved 
x = intvar(1,length(realGreens)); %integer
% x = sdpvar(1,length(realGreens));
%Define Constraints
Const = [ x >= 0, x-minGreens >= 0, -maxFluct <= x-prevGreens <= maxFluct, sum(x)-sum(prevGreens)==0]; 

%Define and objective: 
Obj = (x-realGreens)*(x-realGreens)';

% Set some options for YALMIP and solver
opt = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
% opt = sdpsettings('verbose',1,'solver','quadprog','verbose',0);
% Solve the problem
sol = optimize(Const,Obj,opt);

% Analyze error flags 
if sol.problem == 0
 % Extract and display value
 feasibleGreens = value(x);
else
 disp('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

end

