function [cost,gradient] = SwarmCost(cp, Q, n)
    cost = 0; gradient = zeros(3,width(cp));    
    
    a = 2.0; b = 2.0; inv_a2 = 1 / a / a; inv_b2 = 1 / b / b;  clearance = Q(n).c;
    
    for m = setdiff(1:Q(n).nquad, n)
        ccp1 = []; % cp of other uav
        ccp0 = []; % cp of current uav
        % Check for same size array for the evaluated cp
        if width(cp(1:3,:)) < width(Q(m).ccp(2:4,:))
            ccp1 = Q(m).ccp(2:4,1:width(cp(1:3,:)));
            ccp0 = cp(1:3,:);
        elseif width(cp(1:3,:)) > width(Q(m).ccp(2:4,:))
            ccp1 = Q(m).ccp(2:4,:);
            ccp0 = cp(1:3,1:width(Q(m).ccp(2:4,:)));
        else
            ccp1 = Q(m).ccp(2:4,:);
            ccp0 = cp(1:3,:);
        end

        for l = 1:width(ccp0)
            % The count differs for other UAVs since they are
            % sequentially checked and trajectories are updated in
            % sequence
            
            dist_vec = ccp0(1:3,l) - ccp1(1:3,l);
            ellip_dist = sqrt(dist_vec(3)^2 * inv_a2 + (dist_vec(1)^2 + dist_vec(2)^2) * inv_b2);
            dist_err = clearance - ellip_dist;

            gcoeff(1) = -2 * (clearance / ellip_dist - 1) * inv_b2;
            gcoeff(2) = gcoeff(1);
            gcoeff(3) = -2 * (clearance / ellip_dist - 1) * inv_a2;

            if (dist_err < 0)
                % do nothing
                gradient(:,l) = gradient(:,l) + [0;0;0];
            else
                cost = cost + dist_err^2;
                %  We need to add the gradient together for each of the
                %  columns for the different UAVs
                gradient(:,l) = gradient(:,l) + gcoeff' .* dist_vec;
            end
        end
    end

end

%% Swarm Cost from C++
% void BsplineOptimizer::calcSwarmCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
%   {
%     cost = 0.0;
%     int end_idx = q.cols() - order_ - (double)(q.cols() - 2 * order_) * 1.0 / 3.0; // Only check the first 2/3 points
%     const double CLEARANCE = swarm_clearance_ * 2;
%     double t_now = ros::Time::now().toSec();
%     constexpr double a = 2.0, b = 1.0, inv_a2 = 1 / a / a, inv_b2 = 1 / b / b;
% 
%     for (int i = order_; i < end_idx; i++)
%     {
%       double glb_time = t_now + ((double)(order_ - 1) / 2 + (i - order_ + 1)) * bspline_interval_;
% 
%       for (size_t id = 0; id < swarm_trajs_->size(); id++)
%       {
%         if ((swarm_trajs_->at(id).drone_id != (int)id) || swarm_trajs_->at(id).drone_id == drone_id_)
%         {
%           continue;
%         }
% 
%         double traj_i_satrt_time = swarm_trajs_->at(id).start_time_.toSec();
%         if (glb_time < traj_i_satrt_time + swarm_trajs_->at(id).duration_ - 0.1)
%         {
%           /* def cost=(c-sqrt([Q-O]'D[Q-O]))^2, D=[1/b^2,0,0;0,1/b^2,0;0,0,1/a^2] */
%           Eigen::Vector3d swarm_prid = swarm_trajs_->at(id).position_traj_.evaluateDeBoorT(glb_time - traj_i_satrt_time);
%           Eigen::Vector3d dist_vec = cps_.points.col(i) - swarm_prid;
%           double ellip_dist = sqrt(dist_vec(2) * dist_vec(2) * inv_a2 + (dist_vec(0) * dist_vec(0) + dist_vec(1) * dist_vec(1)) * inv_b2);
%           double dist_err = CLEARANCE - ellip_dist;
% 
%           Eigen::Vector3d dist_grad = cps_.points.col(i) - swarm_prid;
%           Eigen::Vector3d Coeff;
%           Coeff(0) = -2 * (CLEARANCE / ellip_dist - 1) * inv_b2;
%           Coeff(1) = Coeff(0);
%           Coeff(2) = -2 * (CLEARANCE / ellip_dist - 1) * inv_a2;
% 
%           if (dist_err < 0)
%           {
%             /* do nothing */
%           }
%           else
%           {
%             cost += pow(dist_err, 2);
%             gradient.col(i) += (Coeff.array() * dist_grad.array()).matrix();
%           }
% 
%           if (min_ellip_dist_ > dist_err)
%           {
%             min_ellip_dist_ = dist_err;
%           }
%         }
%       }
%     }
%   }
