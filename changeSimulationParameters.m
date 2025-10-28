%Update simulation parameters 

clear 
clc
load indata.mat

%change maximum allowed green time change for Perimeter Control  
indata.gn_R = 5; %(sec)

save('indata.mat','indata') 