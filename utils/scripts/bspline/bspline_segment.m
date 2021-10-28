% Bspline Segment

% Following the notation from
% (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)

% Following the notation from
% https://link.springer.com/article/10.1007/s003710050206 
% ("General matrix representations for B-splines")

% p(t) lies inside [t(i), t(i+1)), and they depend on (order+1) control points
% [t(i−2), t(i−1), t(i), t(i+1), t(i+2), t(i+3)]

clc
close all
clear all

order = 5;
M = getM(order);
u = zeros(1,order + 1);
p = zeros(order + 1,1);

P = [0, 2, 1, 4, 5, 3, 7, 9, 2, 4];

time = 1.5;
timespan = [1,2];

s_interval = (timespan(2) - timespan(1))/(numel(P)-1);
s = linspace(timespan(1), timespan(2), numel(P));

s_t = (time - 0) / s_interval;
u(t) = s_t - timespan(1);

% https://link.springer.com/article/10.1007/s003710050206 
% u(i) = (time - timespan(1)) / (timespan(2) - timespan(1));

% (Real-Time Trajectory Replanning for MAVs using Uniform B-splines 
% and a 3D Circular Buffer)
% s(t) = (t − t(0)) / deltat
% u(t) = s(t) − s(i)
% s(i) lies at the time of the control point(i)

for i = 1:(order+1)
    
end