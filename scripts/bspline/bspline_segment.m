%% Debug Solutions
% - Total time span was not generic, but was fixed
% - Useful span is [order:numel(P)-order]
% - Input [start end] and output was not synced since useful span cut the
% many segments of the time span making the spline not accurate

%% Bugs
% - Now the trajectory does not fully ultilize all the control points,
% hence not reaching the final destination
% *** CANNOT use https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/
% convention for Deboor algorithm - since knots are seperated and determine
% the curve but we must use control points to define clamping

%% Bspline Segment

% Following the notation from
% (1) https://link.springer.com/article/10.1007/s003710050206 
% ("General matrix representations for B-splines")

% Following the notation from
% (2) (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)

% According to (2) etc quintic bspline definition of u is
% s_t = (time - t(1)) / delta_t
% u(t) = s(t) − s(i)
% s(i) lies at the time of the control point(i)

% According to (2) etc quintic bspline uses a knot vector of
% [t(i−2), t(i−1), t(i), t(i+1), t(i+2), t(i+3)]

% p (control points) and t (times allocated to each control point) has a
% range of 0 to n (which is i = 0 to n)

% p(t) lies inside [t(i), t(i+1)), and they depend on k (order+1) control points

clc
close all
clear all

funct_path = '../../functions';
addpath(strcat(funct_path,'/bspline'));

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

%% User inputs
timespan = [10,30];
r = rand(1,35);
for i = 1:order
    r = [0 r 0];
end
numrange = 12; % range for evaluation
P = numrange * r; % multiply the random [0-1] by the number range
n = numel(r) - 1; % n = controlpoints - 1

% How much to segment a knotspan (internal evaluation in a spline segment)
buffer = zeros(1,order + 1);
knots_total = n + order + 1 + 1;
% Range of index to evaluate according
range = [order:numel(P)]; 
% range = [order:knots_total-order]; 
% Not <---order---><----Accepted Range numel(P)-----><----order+1---->
% But <---order---><----------Accepted Range numel(P)---------------->


%% Main code

% I may have made the code abit different by not using knots = m + 1 which
% means using numel(P) rather than knots_total to segment the time horizon
% Solution for matching time span
dt = (timespan(2) - timespan(1)) / (numel(range)-1); % span of 1 knot
t = linspace(timespan(1), timespan(2), (numel(range))); % knots
% End of solution for matching time span

kn_seg = 15; % Division of 1 span or 1 segment

for l = 1:numel(range)-1
    idx = range(l) - order;
    nxt_idx = idx + 1; 
    lpos = zeros(1,kn_seg-1);
    lvel = zeros(1,kn_seg-1);
    lacc = zeros(1,kn_seg-1);
    
    span = linspace(idx, nxt_idx, kn_seg); % relative to the start time as 0 regardless of the time 
    actualspan = linspace(t(l), t(l+1), kn_seg); % time in abs (simulation time / given time)
    
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
        timeActual(m) = actualspan(m);

        % p have several conventions according to SO many papers but i
        % would use the convention under (1) and it is the correct one
        % etc if order is 5
        % (1) p = [P(idx-5) P(idx-4) P(idx-3) P(idx-2) P(idx-1) P(idx)]';
        
        % Make the u, du, ddu and p matrix
        for j = 1:k
            u(j) = u_t^(j-1);
            p(j) = P(idx + (j-1) + 1); % we add a +1 here for matlab notation
            if j >= 2
                du(j) = (j-1) * u_t^(j-2);
            end
            if j >= 3
                ddu(j) = (j-1) * (j-2) * u_t^(j-3);
            end
            if j >= 4
                dddu(j) = (j-1) * (j-2) * (j-3) * u_t^(j-4);
            end
        end
        
        inv_dt = 1/dt;
        
        % Matrix multiplication to attain the pos, vel and acc
        lpos(m) = u * M * p;
        lvel(m) = inv_dt * du * M * p;
        lacc(m) = inv_dt^2 * ddu * M * p;
        ljrk(m) = inv_dt^3 * dddu * M * p;
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

tic
[p0,v0,a0,t0] = getBSpline(order, timespan, P, kn_seg, false);
fprintf("Time Elapsed %f\n",toc);


%% Plot 
% For plotting purposes

figure
subplot(4,1,1)
hold on
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(t(i+1),P(range(i)),'o');
end
plot(tt,pos,'x')
title("Position")
hold off
grid on

subplot(4,1,2)
hold on
plot(tt,vel,'x')
title("Velocity")
hold off
grid on

subplot(4,1,3)
hold on
plot(tt,acc,'x')
title("Acceleration")
hold off
grid on

subplot(4,1,4)
hold on
plot(tt,jrk,'x')
title("Jerk")
hold off
grid on

figure
subplot(3,1,1)
hold on
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(t(i+1),P(range(i)),'o');
end
plot(t0,p0,'x')
title("Position")
hold off
grid on

subplot(3,1,2)
hold on
plot(t0,v0,'x')
title("Velocity")
hold off
grid on

subplot(3,1,3)
hold on
plot(t0,a0,'x')
title("Acceleration")
hold off
grid on

% For cost (integral over squared time derivative) in close form
% Depends on the [limit of order] and also whether you want minimize
% [Acceleration, Jerk, Snap]
