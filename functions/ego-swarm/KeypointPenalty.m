function [cost,gradient] = KeypointPenalty(cp, wpt, t, alpha, A, boundary)
    sp = [t(1) ; cp(:,1)];
    % ep = [t(end) ; cp(:,end)];
    count = 1;
    % waypoint_vector = [sp wpt ep];
    % waypoint_vector = [sp wpt];
    waypoint_vector = wpt;
    cost = 0;
    tmp = zeros(3,width(cp));
    gradient = zeros(3,width(cp));
    for i=1:width(cp)
        % find if any of the discrete time corresponds to the time we are
        % evaluating
        idx(count) = nan;
        for k = 1:width(waypoint_vector)
            if i > width(t)
                break;
            end
            if t(i) == waypoint_vector(1,k)
                idx(count) = k;
                count = count + 1;
                break;
            end
        end
        if isnan(idx(end))
            gradient(:,i) = [0;0;0];
            continue;
        end
        
        d_vector = cp(:,i) - waypoint_vector(2:4,idx(end));
        offset = (A * (exp(-alpha * boundary)) / boundary);
        for j = 1:3
            val = abs(d_vector(j));
            if val < boundary
                val = boundary;
            end
            
            tmp(j,i) = (A * (exp(-alpha * val)) / -val);
            gradient(j,i) = -alpha * tmp(j,i) / val;
            cost = cost + tmp(j,i) + offset;
        end
    end
    
end

