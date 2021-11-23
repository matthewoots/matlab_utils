function [cost,gradient] = CostFunction(cp)

    example = matfile('Q.mat');
    Q = example.Q;
    n = example.n;
    
    % Experiemental values for gain
    gain_smooth = 0.5;
    gain_swarm = 5.0;
    gain_penalty_kp = 0.01;
    boundary = -15;
    
    [cost0,gradient0] = SmoothnessCost(cp);
    [cost1,gradient1] = SwarmCost(cp, Q, n);
    
    % Increasing alpha shrinks the area
    alpha = 0.3; A = 1;
    wpt = Q.wpt; 
    tmp = Q.ccp;
    t = tmp(1,:);
    [costp,gradientp] = KeypointPenalty(cp, wpt, t, alpha, A, boundary);
    cost = gain_smooth * cost0 + gain_swarm * cost1;
    gradient = gain_smooth * gradient0 + gain_swarm * gradient1;
    cost = cost + gain_penalty_kp * costp;
    gradient = gradient + gain_penalty_kp * gradientp;
end

