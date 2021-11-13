function params = q_parameters()
    % Physical parameters
    params.m = 0.100;  % weight (in kg)
    params.g = 9.81;   % gravitational constant
    params.I = [1.43e-5,   0,          0; % inertial tensor in m^2 kg
         0,         1.43e-5,    0;
         0,         0,          2.89e-5];
    params.kp = [15;15;30];
    params.kd = [12;12;10];
    params.kpm = ones(3,1)*3000;
    params.kdm = ones(3,1)*300;
    params.arm_length = 0.08; % arm length in m
    params.maxangle = 30*pi/180; % you can specify the maximum commanded angle here
    params.clearance = 0.9; % uav clearance check
    params.maxF     = 2.5*params.m*params.g;   % left these untouched from the nano plus
    params.minF     = 0.05*params.m*params.g;  % left these untouched from the nano plus
end


