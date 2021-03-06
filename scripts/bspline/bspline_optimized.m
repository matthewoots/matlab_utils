%% Bspline Segment

% Following the notation from
% (1) https://link.springer.com/article/10.1007/s003710050206 
% ("General matrix representations for B-splines")

% Following the notation from
% (2) (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)

% According to (2) etc quintic bspline definition of u is
% s_t = (time - t(1)) / delta_t
% u(t) = s(t) ā s(i)
% s(i) lies at the time of the control point(i)

% According to (2) etc quintic bspline uses a knot vector of
% [t(iā2), t(iā1), t(i), t(i+1), t(i+2), t(i+3)]

% p (control points) and t (times allocated to each control point) has a
% range of 0 to n (which is i = 0 to n)

% p(t) lies inside [t(i), t(i+1)), and they depend on k (order+1) control points

clc
close all
clear all

funct_path = '../../functions/bspline_utils';
addpath(funct_path);

order = 5;
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

%% User inputs
timespan = [0,25];
r = rand(1,35);
numrange = 12; % range for evaluation
P = numrange * r; % mutltiple the random [0-1] by the number range
n = numel(r) - 1;

% How much to segment a knotspan (internal evaluation in a spline segment)
buffer = zeros(1,order + 1);
knots_total = n + order + 1 + 1;
% Range of index to evaluate according
range = [order:numel(P)-order]; 


%% Main code
% I may have made the code abit different by not using knots = m + 1 which
% means using numel(P) rather than knots_total to segment the time horizon
t = linspace(timespan(1), timespan(2), numel(P)); % knots
% t = linspace(timespan(1), timespan(2), knots_total); % knots
dt = (timespan(2) - timespan(1)) / (numel(P)-1); % span of 1 knot
kn_seg = 15; % Division of 1 span or 1 segment
EqA = zeros(1,numel(range));
EqJ = zeros(1,numel(range));

for l = 1:numel(range)
    idx = range(l) - order;
    nxt_idx = idx + 1; 
    lpos = zeros(1,kn_seg-1);
    lvel = zeros(1,kn_seg-1);
    lacc = zeros(1,kn_seg-1);
    
    span = linspace(idx, nxt_idx, kn_seg);
    
    if idx < 0
        error('idx is below suggested idx for control points');
    end 
    
    if idx + 1>= numel(P)
        error('idx is out of bounds compared to control points');
    end
    
    tic
    % numel(span)-1 becasue we are considering [0,1) hence not including 1
    for m = 1:numel(span)-1

        % Using convention from (1) to get u
        time = span(m); % current time in index form, of course we dont like to play with conversion
        u_t = (time - idx)/((idx+1) - idx); % using index is the same as using time since u_t is a factor
        % Save the time constant (0 to 1) and the actual time of the point
        timeConst(m) = u_t;
        timeActual(m) = span(m) * dt;

        % p have several conventions according to SO many papers but i
        % would use the convention under (1) and it is the correct one
        % etc if order is 5
        % (1) p = [P(idx-5) P(idx-4) P(idx-3) P(idx-2) P(idx-1) P(idx)]';
        
        % Make the u, du, ddu and p matrix
        for n = 1:k
            u(n) = u_t^(n-1);
            p(n) = P(idx + (n-1) + 1); % we add a +1 here for matlab notation
            if n >= 2
                du(n) = (n-1) * u_t^(n-2);
            end
            if n >= 3
                ddu(n) = (n-1) * (n-2) * u_t^(n-3);
            end
            if n >= 4
                dddu(n) = (n-1) * (n-2) * (n-3) * u_t^(n-4);
            end
        end
        
        inv_dt = 1/dt;
        
        % Matrix multiplication to attain the pos, vel and acc
        lpos(m) = u * M * p;
        lvel(m) = inv_dt * du * M * p;
        lacc(m) = inv_dt^2 * ddu * M * p;
        ljrk(m) = inv_dt^3 * dddu * M * p;
        
        % Cost function
        % For cost (integral over squared time derivative) in close form
        % Depends on the [limit of order] and also whether you want minimize
        % [Acceleration, Jerk, Snap]

        % Depends on what is used to minimize, we can use Jerk or
        % Acceleration
        % [1] For acceleration
        ddu0 = zeros(1,k); ddu1 = zeros(1,k);
        for n = 1:k
            if n >= 3
                ddu0(n) = (n-1) * (n-2) * 0^(n-3);
                ddu1(n) = (n-1) * (n-2) * 1^(n-3);
            end
        end
        Qa = inv_dt^3 * ((ddu0' * ddu0) + (ddu1' * ddu1));
        
        % [2] For Jerk
        dddu0 = zeros(1,k); dddu1 = zeros(1,k);
        for n = 1:k
            if n >= 4
                dddu0(n) = (n-1) * (n-2) * (n-3) * 0^(n-4);
                dddu1(n) = (n-1) * (n-2) * (n-3) * 1^(n-4);
            end
        end
        Qj = inv_dt^4 * ((dddu0' * dddu0) + (dddu1' * dddu1));
        
        % For 1 Segment
        EqA(l) = EqA(l) + p'* M' * Qa * M * p;
        EqJ(l) = EqJ(l) + p'* M' * Qj * M * p;
    end
    
    % Add the segment to the plot array
    if l-1>0
        pos = [pos lpos];
        vel = [vel lvel];
        acc = [acc lacc];
        jrk = [jrk ljrk];
        tt = [tt timeActual];
        % tc = [tc timeConst];
    else
        pos = lpos;
        vel = lvel;
        acc = lacc;
        jrk = ljrk;
        tt = timeActual;
        % tc = timeConst;
    end
end

fprintf("Time Elapsed %f\n",toc);

%% Plot 
% For plotting purposes

figure
subplot(2,1,1)
plot(1:numel(range),EqA)
title("Plot of B-Spline Cost [Acceleration]")
grid on

subplot(2,1,2)
plot(1:numel(range),EqJ)
title("Plot of B-Spline Cost [Jerk]")
grid on