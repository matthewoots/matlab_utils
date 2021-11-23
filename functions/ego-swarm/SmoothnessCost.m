function [cost,gradient] = SmoothnessCost(cp)
    cost = 0; gradient = zeros(3,width(cp));
    cp_tmp = cp;
    for j = 1:width(cp_tmp)-3    
        jerk = cp_tmp(:,j + 3) - ...
            3 * cp_tmp(:,j + 2) + ...
            3 * cp_tmp(:,j + 1) - ...
            cp_tmp(:,j);
        cost = cost + norm(jerk)^2;
        temp_j = 2.0 * jerk;

        gradient(:,j + 0) = gradient(:,j + 0) + (-temp_j);
        gradient(:,j + 1) = gradient(:,j + 1) + (3.0 * temp_j);
        gradient(:,j + 2) = gradient(:,j + 2) + (-3.0 * temp_j);
        gradient(:,j + 3) = gradient(:,j + 3) + (temp_j);
    end
end
    
%% Smoothness Cost from C++
% void BsplineOptimizer::calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost,
%                                             Eigen::MatrixXd &gradient, bool falg_use_jerk /* = true*/)
%   {
% 
%     cost = 0.0;
% 
%     if (falg_use_jerk)
%     {
%       Eigen::Vector3d jerk, temp_j;
% 
%       for (int i = 0; i < q.cols() - 3; i++)
%       {
%         /* evaluate jerk */
%         jerk = q.col(i + 3) - 3 * q.col(i + 2) + 3 * q.col(i + 1) - q.col(i);
%         cost += jerk.squaredNorm();
%         temp_j = 2.0 * jerk;
%         /* jerk gradient */
%         gradient.col(i + 0) += -temp_j;
%         gradient.col(i + 1) += 3.0 * temp_j;
%         gradient.col(i + 2) += -3.0 * temp_j;
%         gradient.col(i + 3) += temp_j;
%       }
%     }
%     else
%     {
%       Eigen::Vector3d acc, temp_acc;
% 
%       for (int i = 0; i < q.cols() - 2; i++)
%       {
%         /* evaluate acc */
%         acc = q.col(i + 2) - 2 * q.col(i + 1) + q.col(i);
%         cost += acc.squaredNorm();
%         temp_acc = 2.0 * acc;
%         /* acc gradient */
%         gradient.col(i + 0) += temp_acc;
%         gradient.col(i + 1) += -2.0 * temp_acc;
%         gradient.col(i + 2) += temp_acc;
%       }
%     }
%   }

