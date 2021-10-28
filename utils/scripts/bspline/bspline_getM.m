%Following the notation from
%https://link.springer.com/article/10.1007/s003710050206 
%("General matrix representations for B-splines"
% See Theorem 1 of page 182

% Uniform bspline
% There are 2 ways, one is a shortcut, but the other has a real formulation
% View the portion on (section 3.2) on NURBs and (section 4.1)

% t(j) - t(j-1) = constant

% j = start knot point
% order = k-1 
% m(i,j) = (1/(k-1)) * C(k-1-i,k-1) * ...
%   sum(s=j to k-1) {pow(-1,s-j) * C(s-j,k) * pow(k-s-1,k-1-i)} 
% C(i,n) = factorial(n)/(factorial(i) * factorial(n-i));

clc
clear all
close all

funct_path = '../../functions/bspline_utils';
addpath(funct_path);

order = 5;
M = zeros(order,order);

for i = 1:numel(M)
    % notation is different for matlab
    if mod(i,order) == 0
        modv = order;
    else
        modv = mod(i,order);
    end
    idx = [ceil(i/order), modv];
    M(ceil(i/order), modv) = getM(idx,order);
    fprintf("m(%d,%d) = %f\n",ceil(i/order),modv,M(ceil(i/order), modv));
end

M
