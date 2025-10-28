% Create file with PC parameters

PCinput_fname = 'PC_62';
% ----- set PC parameters manually here: 
    
%
PC.K_I = zeros(7,3);
PC.K_P = zeros(7,3);
PC.K_I = [ 15 -8 0;
    0 -5 8;
    -15 8 0;
    0 5 -8;
    -20   0   0;
    0  -15   0;
    0   0  -10];

PC.K_P = PC.K_I*10;
%     PC.K_P = [150 -80 0;
%         0 -50 100;
%         -100 100 0;
%         0 50 -100;
%         -80 0 0;
%         0 -50 0;
%         0 0 -50];

PC.n_set = [5500 6500 5000]';
PC.actCrit = 0.0;
PC.deactCrit = 0.15;



%% ----- Proportional-Integral Gain parameters K_I, K_P (global controller)
%     %set PC parameters here from AFT (after running AFT algorithm / formatted parameters) 
% 
%     parameters = [ 1.8479400e+01  -2.0716658e+00   1.3090768e-03   5.0004709e-03  -3.3426715e+00   6.8018654e+00  -1.3338836e+01   5.5974115e+00   5.6366118e-03   6.0699635e-03   4.6427472e+00  -3.4596572e+00  -1.6433106e+01   1.2015834e-03   1.5000000e-02   9.6923063e-03  -9.7607918e+00   7.1794467e-03   5.2617920e-03   5.9078230e-03  -1.0470632e+01   1.8246177e+01  -6.8710639e+01   1.4563794e-02   7.6551466e-03  -6.9131972e+01   6.2057440e+01  -4.4136224e+01   2.4349013e+01   9.8300350e-03   5.6850834e-03   1.2273656e+00  -4.0402415e+01  -5.4153855e+01   4.1670151e-03   3.8046758e-03   1.0004850e-02  -2.9613973e+01   2.5788204e-03   4.0510642e-03   9.2997594e-03  -5.0909055e+01   4.4787930e+03   9.1919128e+03   7.9785226e+03   9.1642507e-03   6.9452729e-02]';
%     s_setP = 3; %number of regions
%     s_kp_ki = [7 3]; %dimensions of matrices 
%     PC.K_P = zeros(s_kp_ki);
%     PC.K_I = zeros(s_kp_ki);
%     
%     %setpoints 
%     %PC.n_set = transpose(parameters(end-s_setP-1:end-2));
%     PC.n_set = parameters(end-s_setP-1:end-2);
%     %parameters
%     k = 1;
%     for i=1:s_kp_ki(1)
%         PC.K_I(i,:) = parameters(k:k+s_kp_ki(2)-1);
%         PC.K_P(i,:) = parameters(k+s_kp_ki(1)*s_kp_ki(2):k+s_kp_ki(1)*s_kp_ki(2)+s_kp_ki(2)-1);
%         k = k + s_kp_ki(2);
%     end
%     
%     %activation criterion
%     %from 0 to 0.20 (activ when at least one region is less than 20% of the setpoint)
%     PC.actCrit = parameters(end-1);
%     PC.deactCrit = parameters(end);


save(PCinput_fname, 'PC')