clc
clear all
close all

%% Intro and Format
% x format = (time, x, y, z)


%% Add Libraries and Paths

funct_path = '../../functions';
common_path = '../../functions/common';
bspline_path = '../../functions/bspline_utils';
opt_path = '../../functions/optimizer';
addpath(funct_path);
addpath(common_path);
addpath(bspline_path);
addpath(opt_path);

%% Choose Trajectory profile and method
% method = "bvp";
method = "bspline";
size_cp = 3;
% waypointinput = "single";
waypointinput = "multi-random";

%% Set Global Variables
% Shared Variables for UAVs 
global nquad; % how many drones participating
global ctrlpts; % control points for each uav being circulated
global p; % pos of each uav that is being circulated 

%% Setup [Time]
tt = 15; % total time
dt = 0.2; % time interval
int = tt/dt; % units

%% Setup [Boundary]
axes = 3; % How many axis are we using
xy_range = 5.0; z_range = 3.0; % range for xy and z
sprd = 6; % Spread of plotting (how much can we see, for visual)
c = [2;2;2]; % Center of boundary sphere/play area
r = [xy_range; xy_range; z_range]; % Radius of boundary sphere
zlim = 0.1; % Z limit cannot be more than 1
[xbnd, ybnd, zbnd] = sphere(r,c,10); % can plot this as this is the sphere

%% Setup [UAV]
nquad = 5; % number of quad
id = 1:nquad; % give an id to each quad 
clearance = 0.4; % uav clearance check

%% Setup other variables
stop = false;
order = 4;

%% Setup [Start/End/intermediate waypoints]
% Change this to what you need
% Most important is [output] here is x(1:4,:) which is (t, xpos, ypos, zpos)

% p0 = start pose
% pf = end pose
% ss = start state (p, v, a)
% es = end state (p, v, a)
% We can choose the trajectory representation or method

for n = 1:nquad
    % Keep this the same because this generates a random point
    p0 = (randspherepoint(r,c,zlim))';
    % Set the end point via pf
    pf = [c(1) - (p0(1) - c(1)); ...
              c(2 )- (p0(2) - c(2)); ...
              c(3) - (p0(3) - c(3))]; % To get the opposite point from start
    
    %% Single waypoint 
    if method == "bvp"

        % Set [start end] pose
        p_range = [p0 pf];

        % Set [start end] state
        ss = zeros(9,2);
        ss(1:3,1) = p_range(:,1);
        ss(1:3,2) = p_range(:,2);

        % Trajectory
        x{n}(1,:) = linspace(0,tt,int);
        s1 = bvp(dt, tt, ss(:,2), ss(:,1));
        for j = 1:axes
            x{n}(j+1,:) = s1(j,:); % fit the trajectory inside
            v_offset = [x{n}(j+1,1) x{n}(j+1,1:end-1)];
            
            if j > 1
                x{n}(j+1+3,:) = (x{n}(j+1,:)-v_offset)/dt;
            end
        end
        
    %% Multiple waypoints bspline
    elseif method == "bspline"
        t_abs = linspace(0,tt,int);
        ts = [0 tt];
        
        % Random 3 axes points according to cp size
        if waypointinput == "multi-random"
            for q = 1:size_cp
                for j = 1:axes
                    cp_tmp(j,q) = (rand(1)*2-1)*xy_range + c(j);
                    
                end
            end
            for j = 1:axes
                clampstart = []; clampend = [];
                for i = 1:order
                    clampstart = [clampstart p0(j)];
                    clampend = [clampend pf(j)];
                end
                cp(j,:) = [clampstart cp_tmp(j,:) clampend];
            end
        % Single 3 axes points according to cp size
        elseif waypointinput == "single"
            
            for j = 1:axes
                clampstart = []; clampend = [];
                for i = 1:order
                    clampstart = [clampstart p0(j)];
                    clampend = [clampend pf(j)];
                end
                cp(j,:) = [clampstart linspace(p0(j), pf(j), size_cp) clampend];
            end
        end
        
        %  size_cp = int / intv + (2*order-1);
        range = [order:length(cp)]; 
        dtcheck = (ts(2) - ts(1)) / (numel(range)-1); % span of 1 knot
        intv = ceil(dtcheck/dt); % try matching with dt set at the beginning
        
        for j = 1:axes
            [p_0, v_0, a_0, t_0] = getBSpline(order, ts, cp(j,:), intv, false);
            x{n}(1,:) = t_0;
            x{n}(j+1,:) = p_0;
            x{n}(j+4,:) = v_0;
            x{n}(j+7,:) = a_0;
        end
        
    end
    
end

%% Initialise Figure
fig = figure;

%% Recursive Loop
fprintf('Starting UAV Path and Loop...\n');
time = 0;

for iter = 1:int
    clf % Clear figure
    tic

    %% Plotting process
    %  plot3(xbnd,ybnd,zbnd,'.','MarkerSize',1); % Displays the bounding sphere
    hold on
    for n = 1:nquad
        plot3(x{n}(2,1),x{n}(3,1),x{n}(4,1),'o','MarkerSize',10); % Start pose
        plot3(x{n}(2,end),x{n}(3,end),x{n}(4,end),'o','MarkerSize',10); % End pose
        plot3(x{n}(2,iter),x{n}(3,iter),x{n}(4,iter),'x','MarkerSize',10); % Current pose
        plot3(x{n}(2,iter:end),x{n}(3,iter:end),x{n}(4,iter:end),'o','MarkerSize',2); % Displays the path
    end

    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
    axis([c(1)-sprd  c(1)+sprd  c(2)-sprd c(1)+sprd c(3)-sprd  c(3)+sprd]);
    grid on
    title(sprintf('Iteration: %d, time: %4.2f', iter, time));
    timer = toc;

    %% Stopping Criteria
    for i = 1:nquad
        % Collision check with respective uavs
        for m = setdiff(1:nquad, i)
            dv = x{i}(1:3,iter) - x{m}(1:3,iter); % get the distance vector
            a = 1.0; b = 1.0; inv_a2 = 1 / a / a; inv_b2 = 1 / b / b; 
            e_d = sqrt(dv(3)^2 * inv_a2 + (dv(1)^2 + dv(2)^2) * inv_b2); % get the ellipsoidal distance 
            d_e = clearance - e_d; % get the distance error
            if d_e > 0
                %  stop = true; 
                fprintf('Collision! Self(%d) Other(%d)\n', i, m); 
            end
        end
        % If the trajectory ended, we should stop the simulation
        if iter == numel(x{n}(1,:))
            stop = true; 
        end
    end

    if stop
        % error('Stop initiated'); 
        break;
    end

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
    subplot(nquad,2,i+1*(i-1))
    plot(x{i}(1,:),x{i}(2,:),x{i}(1,:),x{i}(3,:),x{i}(1,:),x{i}(4,:));
    xlabel('t [s]'); ylabel('Pos [m]');
    grid on
    title(sprintf('Position UAV %d', i));
    subplot(nquad,2,2*i)
    plot(x{i}(1,:),x{i}(5,:),x{i}(1,:),x{i}(6,:),x{i}(1,:),x{i}(7,:));
    xlabel('t [s]'); ylabel('Vel [m/s]');
    grid on
    title(sprintf('Velocity UAV %d', i));
    
end
