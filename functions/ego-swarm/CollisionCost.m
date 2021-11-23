function [outputArg1,outputArg2] = CollisionCost(inputArg1,inputArg2)
%COLLISIONCOST Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

%% Collsion Cost from C++
% void BsplineOptimizer::calcDistanceCostRebound(const Eigen::MatrixXd &q, double &cost,
%                                              Eigen::MatrixXd &gradient, int iter_num, double smoothness_cost)
% {
% cost = 0.0;
% int end_idx = q.cols() - order_;
% double demarcation = cps_.clearance;
% double a = 3 * demarcation, b = -3 * pow(demarcation, 2), c = pow(demarcation, 3);
% 
% force_stop_type_ = DONT_STOP;
% if (iter_num > 3 && smoothness_cost / (cps_.size - 2 * order_) < 0.1) // 0.1 is an experimental value that indicates the trajectory is smooth enough.
% {
%   check_collision_and_rebound();
% }
% 
% /*** calculate distance cost and gradient ***/
% for (auto i = order_; i < end_idx; ++i)
% {
%   for (size_t j = 0; j < cps_.direction[i].size(); ++j)
%   {
%     double dist = (cps_.points.col(i) - cps_.base_point[i][j]).dot(cps_.direction[i][j]);
%     double dist_err = cps_.clearance - dist;
%     Eigen::Vector3d dist_grad = cps_.direction[i][j];
% 
%     if (dist_err < 0)
%     {
%       /* do nothing */
%     }
%     else if (dist_err < demarcation)
%     {
%       cost += pow(dist_err, 3);
%       gradient.col(i) += -3.0 * dist_err * dist_err * dist_grad;
%     }
%     else
%     {
%       cost += a * dist_err * dist_err + b * dist_err + c;
%       gradient.col(i) += -(2.0 * a * dist_err + b) * dist_grad;
%     }
%   }
% }
% }
