function [outputArg1,outputArg2] = TerminalCost(inputArg1,inputArg2)
%TERMINALCOST Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

%% Terminal Cost from C++
% void BsplineOptimizer::calcTerminalCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
% {
% cost = 0.0;
% 
% // zero cost and gradient in hard constraints
% Eigen::Vector3d q_3, q_2, q_1, dq;
% q_3 = q.col(q.cols() - 3);
% q_2 = q.col(q.cols() - 2);
% q_1 = q.col(q.cols() - 1);
% 
% dq = 1 / 6.0 * (q_3 + 4 * q_2 + q_1) - local_target_pt_;
% cost += dq.squaredNorm();
% 
% gradient.col(q.cols() - 3) += 2 * dq * (1 / 6.0);
% gradient.col(q.cols() - 2) += 2 * dq * (4 / 6.0);
% gradient.col(q.cols() - 1) += 2 * dq * (1 / 6.0);
% }