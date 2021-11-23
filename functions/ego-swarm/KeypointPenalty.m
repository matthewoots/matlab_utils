function [cost,gradient] = KeypointPenalty(cp, wpt, t, alpha, A, boundary)
    startp = [t(1) ; cp(:,1)];
    endp = [t(end) ; cp(:,end)];
    waypoint_vector = [startp wpt endp];
    cost = 0;
    prev_cost = 0;
    for i=1:width(cp)
        % find if any of the discrete time corresponds to the time we are
        % evaluating
        idx = nan;
        for k = 1:width(waypoint_vector)
            if i > width(t)
                break;
            end
            if t(i) == waypoint_vector(1,k)
                idx = k;
                break;
            end
        end
        if isnan(idx)
            gradient(:,i) = [0;0;0];
            continue;
        end
        
        d_vector = cp(:,i) - waypoint_vector(2:4,idx);
        d = sqrt(d_vector(1)^2 + d_vector(2)^2 + d_vector(3)^2);
        cost = cost + A * (exp(-alpha * d)) / -d;
        if isinf(cost) || cost < boundary
            cost = boundary;
        end
        diff = cost - prev_cost;
        gradient(:,i) = [diff ; diff ; diff];
        prev_cost = cost;

    end
    
end

