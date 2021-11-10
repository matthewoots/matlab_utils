clc
clear all
close all
%% Issues
% - Knot span cannot be the same as division of the knots when calculating Bspline 
%   - Issue is in num = ceil(interval / dt) when making cp

%% Intro and Format

% 1. Setup simulation parameters 
% 2. Setup quadcopters according to nquad
%    - Quadcopters are given waypoints (size_wp) to follow
%    - Uav will initialize a path of CP (uniform distribution)
%    - This will be fed into the bspline to get [p,v] at controller hz
% 3. Each [control point] is assigned a [time]
% 4. Run simulation per iter / controller hz
% 5. Check for the [control points] and the [time] that the uav is at
%    - Trim/Truncate the uav cp
%    - Log it into the quad.cp_h


%% Add Libraries and Paths

funct_path = '../../functions';
addpath(strcat(funct_path,'/common'));
addpath(strcat(funct_path,'/bspline'));
addpath(strcat(funct_path,'/optimizer'));
addpath(strcat(funct_path,'/ego-swarm'));

param_path = '../../params';
addpath(param_path);

class_path = '../../class';
addpath(class_path);

%% Choose Trajectory profile and method

%% Setup [Time]
tt = 15; % total time
dt = 0.1; % time interval / based on control output hz
int = tt/dt; % units

%% Setup [Boundary]
axes = 3; % How many axis are we using
xy_range = 9.0; z_range = 3.0; % range for xy and z
sprd = 10; % Spread of plotting (how much can we see, for visual)
c = [2;2;2]; % Center of boundary sphere/play area
r = [xy_range; xy_range; z_range]; % Radius of boundary sphere
zlim = 0.1; % Z limit cannot be more than 1
[xbnd, ybnd, zbnd] = sphere(r,c,10); % can plot this as this is the sphere

%% Setup [UAV]
nquad = 3; % number of quad
param = q_parameters();
isIdeal = true;
planning_horizon = 5;

%% Setup other variables
stop = false;
xy_rand_pos = 3;

%% BSpline settings
order = 4;
size_wp = 2;
knot_factor = 4; % Changes the dist estimation for 1 knot 
knot_span = dt * knot_factor; % Changes the dist estimation for 1 knot

%% Setup [Start/End/intermediate waypoints]
% Change this to what you need
% Most important is [output] here is x(1:4,:) which is (t, xpos, ypos, zpos)

% p0 = start pose
% pf = end pose
% ss = start state (p, v, a)
% es = end state (p, v, a)
% We can choose the trajectory representation or method
x = [];

% Generating waypoints is not reliant on other uav in syncing up the
% timestamps and uniting the cp segments
for n = 1:nquad
    %% Multi waypoints bspline
    % Keep this the same because this generates a random point
    p0 = (randspherepoint(r,c,zlim))';
    % Set the end point via pf
    pf = [c(1) - (p0(1) - c(1)); ...
              c(2 )- (p0(2) - c(2)); ...
              c(3) - (p0(3) - c(3))]; % To get the opposite point from start
    % Random 3 axes points according to cp size
    cp_tmp = [];
    cp = [];
    x = [];
    seg = [];
    tic
    for j = 1:axes
        % Setup waypoints 
        for q = 1:size_wp
            % Generate intermediate waypoints 
            wp_t0(j,q) = (rand(1)*2-1)*xy_rand_pos + c(j);
        end
        % Generate array of waypoints 
        wp(j,:) = [p0(j) wp_t0(j,:) pf(j)];
    end
    
    % Size of wp segments
    wp_count = numel(wp(j,:)) - 1;
    
    % Get seperation between waypoints
    for q = 1:wp_count
        % get 3d vector difference
        vec_diff = wp(:,q+1) - wp(:,q);
        % get 3d vector distance 
        dist_vec_diff = sqrt(vec_diff(1)^2 + vec_diff(2)^2 + vec_diff(3)^2);
        % upload it into segments to be used
        wp_sep(q) = dist_vec_diff;
    end
    
    seg_total = 0;

    % Estimated velocity is total time needed to complete
    % all the wp seperation over total time given
    est_vel = sum(wp_sep)/tt;

    % Estimation of distance to velocity
    dist = est_vel * knot_span; 
    % Check to see which axis has the biggest segment count
    for q = 1:wp_count
        % ceil helps to push values above 0 to 1 or more
        % or else seg count is 0 and causes an error
        seg(q) = ceil(wp_sep(q)/dist);
        % total number of segments will be according to the distance
        seg_total = seg_total + ceil(wp_sep(q)/dist);
    end
    
    for j = 1:axes
        % Get the total count for the segments for that axis
        t_seg = sum(seg(:));
        cp_t = [];
        for q = 1:wp_count
            % Why x axis would determine the rest of the axis segment
            % length?
            % It is to unsure the unification of the data, decreasing the x
            % axis length to make sure that the y and z axis have a higher
            % number of cp
            if j == 1
                split = floor(seg(q) / t_seg * seg_total);
            else
                split = ceil(seg(q) / t_seg * seg_total);
            end
            cp0 = linspace(wp(j,q),wp(j,q+1),split);
            cp_t = [cp_t cp0(1:end-1)];
        end
        % Create uniform size of segments for all axis with the first
        if j > 1 && numel(cp_tmp(1,:)) ~= numel(cp_t)
            % Over here we crop the cp of the other axis so that it will 
            % fit into our baseline
            cp_tmp(j,:) = cp_t(1:numel(cp_tmp(1,:)));
        else 
            % Over here we use x axis segments as the baseline
            cp_tmp(j,:) = cp_t;
        end
    end
    fprintf("Time elapsed to setup [waypoints & cp] %.4f\n",toc);
    

    t_abs = linspace(0,tt,int);
    ts = [0 tt];

    for j = 1:axes
        % We will clamp the bspline by adding to the first and the last cp
        cp(j,:) = getClampedCP(cp_tmp(j,:),order);
        range = [order:length(cp)]; 
        dtcheck = (ts(2) - ts(1)) / (numel(range)-1); % span of 1 knot
        % try matching with dt set at the beginning
        intv(j) = ceil(dtcheck/dt); 
        fprintf("[%d] range %d dtcheck %.4f dt %.4f\n",j,numel(range),dtcheck,dt);
        if intv(j) == 1
            error("interval matching is 1\n");
        end
        [x(j+1,:), x(j+4,:), accel, x(1,:)] = getBSpline(order, ts, cp(j,:), intv(j), false);
    end
    
    t = linspace(ts(1), ts(2), numel(range));
    
    %% Setup Quadcopter
    ss = zeros(13,1);
    ss(1:3) = p0;
    ss(7:10) = [1;0;0;0];
    Q(n) = quadcopter(nquad, n, param, ss, dt, isIdeal, wp);
    Q(n).setOtherDebugParam(seg,seg_total,dtcheck,dt,intv); 
    Q(n).setPath(x);
    Q(n).setControlPoint(cp, t, planning_horizon, order)
    Q(n).setTerminal([p0 pf]);

end

%% Initialise Figure
fig = figure;

%% Recursive Loop
fprintf('Starting UAV Path and Loop...\n');
time = 0;

for iter = 1:int

    tic
    for n = 1:nquad
        %% Stopping Criteria
        % If the trajectory ended, we should stop the simulation
        if isempty(find(Q(n).path(1,:) >= iter*dt))
            stop = true;
        end
        % Collision check with respective uavs
        for m = setdiff(1:nquad, n)
            col = Q(n).getCollisionCheck(Q(m));
            if (col)
                fprintf('Collision! [%d with %d]\n', n, m); 
                %  stop = true;
            end
        end
    end
    
    if stop == true
        break;
    end
    
    %% State updates
    clf % Clear figure
    %  plot3(xbnd,ybnd,zbnd,'.','MarkerSize',1); % Displays the bounding sphere
    hold on
    for n = 1:nquad
        %% Update 
        Q(n).updateState(iter);
        %% Plotting
        Q(n).plotTerminalPose();
        Q(n).plotCurrentPose();
        % Q(n).plotFuturePath();
        Q(n).plotCurrentPath();
        Q(n).plotWaypoint();
        Q(n).plotControlPoint();
    end

    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
    axis([c(1)-sprd  c(1)+sprd  c(2)-sprd c(1)+sprd c(3)-sprd  c(3)+sprd]);
    grid on
    title(sprintf('Iteration: %d, time: %4.2f', iter, time));
    timer = toc;
    fprintf('Iteration %d, time %.2f\n', iter, timer);    

    %% Set delay to real time 
    if timer < dt
       pause(dt-timer); 
    end
    time = time + dt;

end

fprintf('Completed Simulation\n');

%% State plot
f2 = figure(2);
% Left Bottom Width Height
f2.Position = [100 50 500 500];
for i = 1:nquad 
    Q(i).plotStateHistory(i);
end

f3 = figure(3);
% Left Bottom Width Height
f3.Position = [650 50 500 500];
for i = 1:nquad 
    Q(i).plotTotalVelHistory(i);
end
