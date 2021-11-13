function [cost,gradient] = CostFunction(cp)

    example = matfile('Q.mat');
    Q = example.Q;
    n = example.n;
    
    % Experiemental values for gain
    gain_smooth = 0.0001;
    gain_swarm = 3.0;
    
    [cost0,gradient0] = SmoothnessCost(cp, gain_smooth);
    [cost1,gradient1] = SwarmCost(cp, Q, n, gain_swarm);
    
    cost = cost0 + cost1;
    gradient = gradient0 + gradient1;
end

