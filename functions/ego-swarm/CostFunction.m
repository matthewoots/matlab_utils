function [cost,gradient] = CostFunction(cp)
    gain_smooth = 1.5;
    [cost,gradient] = SmoothnessCost(cp, gain_smooth);
end

