clc
clear all 
close all

funct_path = '../../functions';
addpath(strcat(funct_path,'/common'));
addpath(strcat(funct_path,'/bspline'));

%% Use bspline to get uav states
order = 4;
k = order + 1;

M = getM(order);
u = zeros(1,k); du = zeros(1,k); ddu = zeros(1,k); dddu = zeros(1,k);
p = zeros(k,1);

% p = degree
% n = cp - 1
% knots = m + 1
% basics is m = n + p + 1
% Open B-splines are [knot(p), knot(m-p)]
% * Clamping only happens for the knots not for the control points!
% * Control points are not really attached to time

% User inputs
timespan = [0,20];
pos = []; vel = []; acc = []; jrk = []; tt = [];
kn_seg = 10; % Division of 1 span or 1 segment
for i=1:3
    r = rand(1,4);
    for j = 1:order
        r = [r(1) r r(end)];
    end
    numrange = 12; % range for evaluation
    P(i,:) = numrange * r; % multiply the random [0-1] by the number range
    n = numel(r) - 1; % n = controlpoints - 1

    % How much to segment a knotspan (internal evaluation in a spline segment)
    buffer = zeros(1,order + 1);
    knots_total = n + order + 1 + 1;
    % Range of index to evaluate according
    range = [order:numel(P(i,:))]; 
    % range = [order:knots_total-order]; 
    % Not <---order---><----Accepted Range numel(P)-----><----order+1---->
    % But <---order---><----------Accepted Range numel(P)---------------->

    % I may have made the code abit different by not using knots = m + 1 which
    % means using numel(P) rather than knots_total to segment the time horizon
    % Solution for matching time span
    dt = (timespan(2) - timespan(1)) / (numel(range)-1); % span of 1 knot
    time = linspace(timespan(1), timespan(2), (numel(range))); % knots
    % End of solution for matching time span

    [p0,v0,a0,t0] = getBSpline(order, timespan, P(i,:), kn_seg, false);
    pos(i,:) = p0; 
    vel(i,:) = v0; 
    t = t0;
end

figure(1)
hold on
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot3(P(1,range(i)),P(2,range(i)),P(3,range(i)),'o');
end
plot3(pos(1,:),pos(2,:),pos(3,:),'x');
grid on

figure(2)
subplot(4,1,1)
hold on
plot(t,pos(1,:))
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(time(i+1),P(1,range(i)),'o');
end
title("Position x")
hold off
grid on 

subplot(4,1,2)
hold on
plot(t,pos(2,:))
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(time(i+1),P(2,range(i)),'o');
end
title("Position y")
hold off
grid on

subplot(4,1,3)
hold on
plot(t,pos(3,:))
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(time(i+1),P(3,range(i)),'o');
end
title("Position z")
hold off
grid on

subplot(4,1,4)
hold on
plot(t,vel(1,:))
plot(t,vel(2,:))
plot(t,vel(3,:))
title("Velocity")
hold off
grid on

% I'll test the concept with time (since my waypoints fall at a certain
% point in time)
max_range = numel(range);
vect = [time(2:max_range) ; P(:,range(1):range(max_range-1))];
cost = []; d = [];
prev_cost = 0;
% Increasing alpha shrinks the area
alpha = 0.3;
A = 1;
boundary_val = -20;
count = 0;
interval = t(2) - t(1);
rg_eval_tmp = time(2:max_range);
rg_eval = [];
for k = 1:width(rg_eval_tmp)
    tmp_val = [rg_eval_tmp(k)-interval ...
        rg_eval_tmp(k) ...
        rg_eval_tmp(k)+interval];
    rg_eval = [rg_eval tmp_val];
end
% should we add a (-) to the denominator :
%   the whole curve is flipped
for i=1:width(pos)
    % find if any of the discrete time corresponds to the time we are
    % evaluating
    idx = nan;
    for k = 1:width(rg_eval)
        if t(i) == rg_eval(k)
            idx = k;
            break;
        end
    end
    if isnan(idx)
        cost(i) = 0;
        continue;
    end
    count = count + 1;
    super_idx = ceil(idx/3);
            
    d_vector = pos(:,i) - vect(2:4,super_idx);
    d(i) = sqrt(d_vector(1)^2 + d_vector(2)^2 + d_vector(3)^2);
    cost(i) = A * (exp(-alpha * d(i))) / -d(i);
    if isinf(cost(i)) || cost(i) < boundary_val
        cost(i) = boundary_val;
    end
            
    gradient(i) = cost(i);
end


f = uifigure;
plot(t,cost,'o');
grid on;