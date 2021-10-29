% Bspline Segment

% Following the notation from
% (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)

% Following the notation from
% https://link.springer.com/article/10.1007/s003710050206 
% ("General matrix representations for B-splines")

% p(t) lies inside [t(i), t(i+1)), and they depend on (order+1) control points
% etc quintic
% [t(i−2), t(i−1), t(i), t(i+1), t(i+2), t(i+3)]

% p (control points) and t (times allocated to each control point) has a
% range of 0 to n (which is i = 0 to n)

clc
close all
clear all

order = 6;
% 5th order was used in the paper
M = getM(order);
u = zeros(1,order); du = zeros(1,order); ddu = zeros(1,order);
p = zeros(order,1);

%% User inputs
timespan = [0,10];
P = [5, 2, 1, 4, 5, 3, 5, 8, 2, 4, 8, 9, 10];
% how much to segment a knotspan
kn_seg = 20;
% start index
range = [4 5 6 7]; 


%% Main code
t = linspace(timespan(1), timespan(2), numel(P));

for k = 1:numel(range)
    % because idx 1 should represent 0
    idx = range(k)-1;
    nxt_idx = idx +1; 
    lpos = zeros(1,kn_seg-1);
    lvel = zeros(1,kn_seg-1);
    lacc = zeros(1,kn_seg-1);
    
    span = linspace(idx, nxt_idx, kn_seg);
    dt = (timespan(2) - timespan(1))/(numel(P)-1);

    if idx > numel(P)
        error('i is greater than control points');
    end

    if idx-round(order/2) < 0 || idx+round(order/2) > numel(P)
        error('i+%d and i-%d is out of bounds compared to control points', round(order/2), round(order/2));
    end
    
    for m = 1:numel(span)-1
        time = span(m) * dt - idx * dt;
        s_t = (time - t(1)) / dt;
        s_i = find(t <= time , 1, 'last');
        u_t = s_t - s_i;
        timeConst(m) = time;
        timeActual(m) = span(m) * dt;

        % Make the u, du, ddu matrix
        for n = 1:order
            u(n) = u_t^(n-1);
            if n >= 2
                du(n-1) = (n-1) * u_t^(n-2);
            end
            if n >= 3
                ddu(n-2) = (n-1) * (n-2) * u_t^(n-3);
            end
        end

        % Make the p matrix
        % **hardcode for now
        p = [P(idx-2) P(idx-1) P(idx) P(idx+1) P(idx+2) P(idx+3)]';
        
        inv_dt = 1/dt;
        
        lpos(m) = u * M * p;
        lvel(m) = inv_dt * du * M * p;
        lacc(m) = inv_dt^2 * ddu * M * p;
    end
    if k-1>0
        pos = [pos lpos];
        vel = [vel lvel];
        acc = [acc lacc];
        tt = [tt timeActual];
    else
        pos(1:k*numel(span)-1) = lpos(1:numel(span)-1);
        vel(1:k*numel(span)-1) = lvel(1:numel(span)-1);
        acc(1:k*numel(span)-1) = lacc(1:numel(span)-1);
        tt(1:k*numel(span)-1) = timeActual(1:numel(span)-1);
    end
end

% figure
% plot(timeConst,lpos)
% grid on

figure
subplot(3,1,1)
hold on
for i = 1:numel(range)
    plot(t(range(i)),P(range(i)),'x');
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

% s_t = (time - 0) / s_interval;
% u(t) = s_t - timespan(1);

% https://link.springer.com/article/10.1007/s003710050206 
% u(i) = (time - timespan(1)) / (timespan(2) - timespan(1));

% (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)
% s_t = (time - t(1)) / delta_t
% u(t) = s(t) − s(i)
% s(i) lies at the time of the control point(i)

% For cost (integral over squared time derivative) in close form
% Depends on the [limit of order] and also whether you want minimize
% [Acceleration, Jerk, Snap]

% Q = 1/delta_t^3 * 