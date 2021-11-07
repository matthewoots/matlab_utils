clc
clear all
close all

%% Intro and Format

%% Add Libraries and Paths

funct_path = '../../functions';
addpath(strcat(funct_path,'/common'));
addpath(strcat(funct_path,'/bspline'));
addpath(strcat(funct_path,'/optimizer'));

param_path = '../../params';
addpath(param_path);

class_path = '../../class';
addpath(class_path);

%% Choose Trajectory profile and method

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
nquad = 2; % number of quad
param = q_parameters();
isIdeal = false;

%% Setup other variables
stop = false;
order = 4;
size_cp = 1;

%% Setup [Start/End/intermediate waypoints]
% Change this to what you need
% Most important is [output] here is x(1:4,:) which is (t, xpos, ypos, zpos)

% p0 = start pose
% pf = end pose
% ss = start state (p, v, a)
% es = end state (p, v, a)
% We can choose the trajectory representation or method

for n = 1:nquad
    %% Multi waypoints bspline
    % Keep this the same because this generates a random point
    p0 = (randspherepoint(r,c,zlim))';
    % Set the end point via pf
    pf = [c(1) - (p0(1) - c(1)); ...
              c(2 )- (p0(2) - c(2)); ...
              c(3) - (p0(3) - c(3))]; % To get the opposite point from start
    % Random 3 axes points according to cp size
    for j = 1:axes
        for q = 1:size_cp
            cp_t0(j,q) = (rand(1)*2-1)*xy_range + c(j);  
        end
        cp_tmp(j,:) = [p0(j) cp_t0(j,:) pf(j)];
    end
          

    t_abs = linspace(0,tt,int);
    ts = [0 tt];

    for j = 1:axes
        cp(j,:) = getClampedCP(cp_tmp(j,:),order);
        range = [order:length(cp)]; 
        dtcheck = (ts(2) - ts(1)) / (numel(range)-1); % span of 1 knot
        intv = ceil(dtcheck/dt); % try matching with dt set at the beginning
        [x(j+1,:), x(j+4,:), accel, x(1,:)] = getBSpline(order, ts, cp(j,:), intv, false);
    end
    
    %% Setup Quadcopter
    ss = zeros(13,1);
    ss(1:3) = p0;
    ss(7:10) = [1;0;0;0];
    Q(n) = quadcopter(nquad, n, param, ss, dt, isIdeal);
    Q(n).setPath(x);
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
        for m = setdiff(1:nquad, i)
            if (Q(n).collisionCheck(Q(m)))
                fprintf('Collision! Self(%d) Other(%d)\n', i, m); 
                %  stop = true;
            end
        end
    end
    
    if stop == true
        break;
    end

    clf % Clear figure
    %  plot3(xbnd,ybnd,zbnd,'.','MarkerSize',1); % Displays the bounding sphere
    hold on
    for n = 1:nquad
        %% Update 
        Q(n).updateState(iter);
        %% Plotting
        Q(n).plotTerminalPose();
        Q(n).plotCurrentPose();
        Q(n).plotFuturePath();
        Q(n).plotCurrentPath();
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
