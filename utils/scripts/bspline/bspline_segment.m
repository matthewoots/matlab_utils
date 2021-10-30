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

funct_path = '../../functions/bspline_utils';
addpath(funct_path);

order = 3;
k = order + 1;
% 5th order was used in the paper
M = getM(order);
u = zeros(1,k); du = zeros(1,k); ddu = zeros(1,k);
p = zeros(k,1);

%% User inputs
timespan = [0,25];
r = rand(1,30);
range = 12;
P = range * r;

% How much to segment a knotspan (internal evaluation in a spline segment)
kn_seg = 20;
% Range of index to evaluate according
range = [7:15]; 


%% Main code
t = linspace(timespan(1), timespan(2), numel(P));

for l = 1:numel(range)
    % Matlab notation have to be changed since 0 is 1 in Matlab
    % Because idx 1 should represent 0
    idx = range(l)-1;
    nxt_idx = idx +1; 
    lpos = zeros(1,kn_seg-1);
    lvel = zeros(1,kn_seg-1);
    lacc = zeros(1,kn_seg-1);
    
    span = linspace(idx, nxt_idx, kn_seg);
    dt = (timespan(2) - timespan(1))/(numel(P)-1);

    if idx - order <= 0
        error('idx is below suggested idx for control points');
    end 
    
    if idx + 1>= numel(P)
        error('idx is out of bounds compared to control points');
    end
    
    for m = 1:numel(span)-1
        
        % Using convention from (1) to get u
        time = span(m) * dt;
        u_t = (time - t(range(l)))/(t(range(l)+1) - t(range(l)));
        % Save the time constant (0 to 1) and the actual time of the point
        timeConst(m) = u_t;
        timeActual(m) = span(m) * dt;

        % p have several conventions according to SO many papers but i
        % would use the convention under (1)
        % (1) p = [P(idx-5) P(idx-4) P(idx-3) P(idx-2) P(idx-1) P(idx)]';
        % (2) p = [P(idx-2) P(idx-1) P(idx) P(idx+1) P(idx+2) P(idx+3)]';
        
        % Make the u, du, ddu and p matrix
        for n = 1:k
            u(n) = u_t^(n-1);
            p(n) = P(idx - order + (n-1));
            if n >= 2
                du(n) = (n-1) * u_t^(n-2);
            end
            if n >= 3
                ddu(n) = (n-1) * (n-2) * u_t^(n-3);
            end
        end
        
        inv_dt = 1/dt;
        
        % Matrix multiplication to attain the pos, vel and acc
        lpos(m) = u * M * p;
        lvel(m) = inv_dt * du * M * p;
        lacc(m) = inv_dt^2 * ddu * M * p;
    end
    
    % Add the segment to the plot array
    if l-1>0
        pos = [pos lpos];
        vel = [vel lvel];
        acc = [acc lacc];
        tt = [tt timeActual];
        tc = [tc timeConst];
    else
        pos = lpos;
        vel = lvel;
        acc = lacc;
        tt = timeActual;
        tc = timeConst;
    end
end

%% Plot 
% For plotting purposes

figure
subplot(3,1,1)
hold on
for i = 1:numel(range)
    plot(t(range(i)),P(range(i)),'o');
end
plot(tt,pos,'x')
hold off
grid on

subplot(3,1,2)
hold on
plot(tt,vel,'x')
hold off
grid on

subplot(3,1,3)
hold on
plot(tt,acc,'x')
hold off
grid on

% For cost (integral over squared time derivative) in close form
% Depends on the [limit of order] and also whether you want minimize
% [Acceleration, Jerk, Snap]
