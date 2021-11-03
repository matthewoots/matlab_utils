clc
clear all
close all

%% Intro and Format
% x format = (time, x, y, z)


%% Add Libraries and Paths

funct_path = '../functions';
common_path = '../functions/common';
bspline_path = '../functions/bspline_utils';
opt_path = '../functions/optimizer';
addpath(funct_path);
addpath(common_path);
addpath(bspline_path);
addpath(opt_path);

%% Setup [Time]

tt = 10; % total time
dt = 0.2; % time interval
int = tt/dt; % units
x(1,:) = linspace(0,tt,int);

%% Setup [Boundary]

axes = 3; % How many axis are we using
sprd = 5; % Spread of plotting (how much can we see, for visual)
c = [2;2;5]; % Center of boundary sphere/play area
r = [5.0;4.0;5.0]; % Radius of boundary sphere
zlim = 0.4; % Z limit cannot be more than 1
[xbnd, ybnd, zbnd] = sphere(r,c,10); % can plot this as this is the sphere

%% Setup [Start/End/intermediate waypoints]
% Change this to what you need
% Most important is [output] here is x(1:4,:) which is (t, xpos, ypos, zpos)

% p0 = start pose
% pf = end pose
% ss = start state (p, v, a)
% es = end state (p, v, a)
% We can choose the trajectory representation or method

p0 = randspherepoint(r,c,zlim);
% Set the end point via pf
pf = [c(1) - (p0(1) - c(1)), ...
      c(2 )- (p0(2) - c(2)), ...
      c(3) - (p0(3) - c(3))]; % To get the opposite point from start

% Set [start] state
ss = zeros(9,1);
ss(1:3) = p0;

% Set [end] state
es = zeros(9,1);
es(1:3) = pf;

% Trajectory
s1 = bvp(dt, tt, es, ss);
for j = 1:axes
    x(j+1,:) = s1(j,:); % Trajectory
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
% plot3(xbnd,ybnd,zbnd,'.','MarkerSize',1); % Displays the bounding sphere
plot3(x(2,iter),x(3,iter),x(4,iter),'x','MarkerSize',10);
hold on
plot3(x(2,iter:end),x(3,iter:end),x(4,iter:end),'--'); % Displays the path

xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
axis([c(1)-sprd  c(1)+sprd  c(2)-sprd c(1)+sprd c(3)-sprd  c(3)+sprd]);
grid on
title(sprintf('Iteration: %d, time: %4.2f', iter, time));
t = toc;

% Set delay to real time 
if t < dt
   pause(dt-t); 
end
time = time + dt;

end

fprintf('Completed Simulation\n');
