function [cost,gradient] = TerminalCost(cp, ref)
    cost = 0;
    gradient = zeros(3,width(cp));
    for i=1:width(cp)
        diff = cp(:,i) - ref(:,i);
        sqdiff = diff.^2;
        cost = cost + sqrt(sqdiff(1)^2 + sqdiff(2)^2 + sqdiff(3)^2);
        gradient(:,i) = 2 .* diff;
    end
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