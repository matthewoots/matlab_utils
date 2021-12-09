function [cost,gradient] = CostFunction(cp)

    example = matfile('Q.mat');
    Q = example.Q;
    n = example.n;
    
    cost = 0; gradient = zeros(3,width(cp));
    
    % Experimental values for gain
    % Motion/avoidance gains
    gain_smooth = 1;
    gain_swarm = 4; % 4
    gain_terminal = 0.5; % 0.5
    
    % Penalty functions 
    gain_feasibility = 30;
    gain_penalty_kp = 0.1; % 0.1
    
    % Smoothness Cost 
    [costsm,gradientsm] = SmoothnessCost(cp);
    cost = cost + gain_smooth .* costsm;
    gradient = gradient + gain_smooth .* gradientsm;
    
    % Swarm Reciprocal Avoidance Cost
    pzhor = 3; pzver = 3;
    [costsw,gradientsw] = SwarmCost(cp, Q, n, pzhor, pzver);
    cost = cost + gain_swarm .* costsw;
    gradient = gradient +  gain_swarm .* gradientsw;
    
    % Terminal Cost
    ccp = Q(n).ccp;
    cpt = Q(n).cpt0;
    t = ccp(1,:);
    s = 1:1:width(cpt);
    ns = s(cpt(1,:)>=t(1) & cpt(1,:)<=t(end));
    ref = cpt(2:4,ns(1):ns(end));
    [costtm,gradienttm] = TerminalCost(cp, ref);
    cost = cost + gain_terminal * costtm;
    gradient = gradient + gain_terminal .* gradienttm;
    
    % Feasibility Cost
    max_acc = 1;
    dt = Q(n).dt;
    int = Q(n).intv;
    cpdt = dt * int(1);
    [costfs,gradientsf] = FeasibilityCost(cp, cpdt, max_acc);
    cost = cost + gain_feasibility * costfs;
    gradient = gradient + gain_feasibility .* gradientsf;
    
 
    % Keypoint Penalty Cost
    % Increasing alpha shrinks the area
    alpha = 0.01; A = 1;
    wpt = Q(n).wpt; 
    tmp = Q(n).ccp;
    t = tmp(1,:);
    boundary = 0.07;
    [costkp,gradientkp] = KeypointPenalty(cp, wpt, t, alpha, A, boundary);
    cost = cost + gain_penalty_kp * costkp;
    gradient = gradient + gain_penalty_kp .* gradientkp;
end

