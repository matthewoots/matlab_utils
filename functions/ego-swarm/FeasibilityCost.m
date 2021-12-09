function [cost,gradient] = FeasibilityCost(cp, cpdt, max_acc_)
    
    ts_inv2 = 1 / cpdt / cpdt;
    cost = 0;
    gradient = zeros(3,width(cp));

    for i = 1:width(cp) - 2
      ai = (cp(:,i+2) - 2 * cp(:,i+1) + cp(:,i)) * ts_inv2;
      for j = 1:3
        if (ai(j) > max_acc_)
          cost = cost + (ai(j) - max_acc_)^2;

          gradient(j, i + 0) = gradient(j, i + 0) + 2 * (ai(j) - max_acc_) * ts_inv2;
          gradient(j, i + 1) = gradient(j, i + 1) + (-4 * (ai(j) - max_acc_) * ts_inv2);
          gradient(j, i + 2) = gradient(j, i + 2) + 2 * (ai(j) - max_acc_) * ts_inv2;

        elseif (ai(j) < -max_acc_)
          cost = cost + (ai(j) + max_acc_)^2;

          gradient(j, i + 0) = gradient(j, i + 0) + 2 * (ai(j) + max_acc_) * ts_inv2;
          gradient(j, i + 1) = gradient(j, i + 1) + (-4 * (ai(j) + max_acc_) * ts_inv2);
          gradient(j, i + 2) = gradient(j, i + 2) + 2 * (ai(j) + max_acc_) * ts_inv2;

        else

        end
      end

    end

end

%% Feasibility Cost from C++
% void BsplineOptimizer::calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost,
%                                              Eigen::MatrixXd &gradient)
% {
%     cost = 0.0;
%     /* abbreviation */
%     double ts, /*vm2, am2, */ ts_inv2;
%     // vm2 = max_vel_ * max_vel_;
%     // am2 = max_acc_ * max_acc_;
% 
%     ts = bspline_interval_;
%     ts_inv2 = 1 / ts / ts;
% 
%     /* velocity feasibility */
%     for (int i = 0; i < q.cols() - 1; i++)
%     {
%       Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;
% 
%       //cout << "temp_v * vi=" ;
%       for (int j = 0; j < 3; j++)
%       {
%         if (vi(j) > max_vel_)
%         {
%           // cout << "zx-todo VEL" << endl;
%           // cout << vi(j) << endl;
%           cost += pow(vi(j) - max_vel_, 2) * ts_inv2; // multiply ts_inv3 to make vel and acc has similar magnitude
% 
%           gradient(j, i + 0) += -2 * (vi(j) - max_vel_) / ts * ts_inv2;
%           gradient(j, i + 1) += 2 * (vi(j) - max_vel_) / ts * ts_inv2;
%         }
%         else if (vi(j) < -max_vel_)
%         {
%           cost += pow(vi(j) + max_vel_, 2) * ts_inv2;
% 
%           gradient(j, i + 0) += -2 * (vi(j) + max_vel_) / ts * ts_inv2;
%           gradient(j, i + 1) += 2 * (vi(j) + max_vel_) / ts * ts_inv2;
%         }
%         else
%         {
%           /* code */
%         }
%       }
%     }
% 
%     /* acceleration feasibility */
%     for (int i = 0; i < q.cols() - 2; i++)
%     {
%       Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;
% 
%       //cout << "temp_a * ai=" ;
%       for (int j = 0; j < 3; j++)
%       {
%         if (ai(j) > max_acc_)
%         {
%           // cout << "zx-todo ACC" << endl;
%           // cout << ai(j) << endl;
%           cost += pow(ai(j) - max_acc_, 2);
% 
%           gradient(j, i + 0) += 2 * (ai(j) - max_acc_) * ts_inv2;
%           gradient(j, i + 1) += -4 * (ai(j) - max_acc_) * ts_inv2;
%           gradient(j, i + 2) += 2 * (ai(j) - max_acc_) * ts_inv2;
%         }
%         else if (ai(j) < -max_acc_)
%         {
%           cost += pow(ai(j) + max_acc_, 2);
% 
%           gradient(j, i + 0) += 2 * (ai(j) + max_acc_) * ts_inv2;
%           gradient(j, i + 1) += -4 * (ai(j) + max_acc_) * ts_inv2;
%           gradient(j, i + 2) += 2 * (ai(j) + max_acc_) * ts_inv2;
%         }
%         else
%         {
%           /* code */
%         }
%       }
%       //cout << endl;
%     }
% }