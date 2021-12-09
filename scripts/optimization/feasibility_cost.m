clc 
clear all 
close all

max = 2;
q = rand(3,12) .* max - max/2;
max_acc_ = 8;
dt = 0.4;
ts_inv2 = 1 / dt / dt;
cost = zeros(3,width(q) - 2);
gradient = zeros(3,width(q));
pboundary = ones(1,width(q)-2) * max_acc_;
nboundary = ones(1,width(q)-2) * -max_acc_;

for i = 1:width(q) - 2
  ai(:,i) = (q(:,i+2) - 2 * q(:,i+1) + q(:,i)) * ts_inv2;
  for j = 1:3
    if (ai(j) > max_acc_)
      cost(j,i) = (ai(j,i) - max_acc_)^2;

      gradient(j, i + 0) = gradient(j, i + 0) + 2 * (ai(j,i) - max_acc_) * ts_inv2;
      gradient(j, i + 1) = gradient(j, i + 1) + (-4 * (ai(j,i) - max_acc_) * ts_inv2);
      gradient(j, i + 2) = gradient(j, i + 2) + 2 * (ai(j,i) - max_acc_) * ts_inv2;

    elseif (ai(j) < -max_acc_)
      cost(j,i) = (ai(j,i) + max_acc_)^2;

      gradient(j, i + 0) = gradient(j, i + 0) + 2 * (ai(j,i) + max_acc_) * ts_inv2;
      gradient(j, i + 1) = gradient(j, i + 1) + (-4 * (ai(j,i) + max_acc_) * ts_inv2);
      gradient(j, i + 2) = gradient(j, i + 2) + 2 * (ai(j,i) + max_acc_) * ts_inv2;

    else

    end
  end

end

[X,Y] = meshgrid(-max_acc_ - 2:0.5:max_acc_ + 2,-max_acc_ - 2:0.5:max_acc_ + 2);
Z = zeros(width(X),height(X));
for i=1:width(X)
    for j=1:height(X)
        acc = sqrt(X(i,j)^2 + Y(i,j)^2);
        diff = acc-max_acc_;
        if diff >= 0 
            Z(i,j) = (acc-max_acc_)^2;
        end
    end
end


% for i = 1:sample
%     
%     a0(:,i) = (q(1:2,i+2) - 2 * q(1:2,i+1) + q(1:2,i)) * ts_inv2;
%     for j = 1:2
%         if (ai(j) > max_acc_)
%           cost0(j,i) = (a0(j,i) - max_acc_)^2;
% 
%         elseif (ai(j) < -max_acc_)
%           cost0(j,i) = (a0(j,i) + max_acc_)^2;
% 
%         else
%         end
%     end
%     
% end

figure(1)
subplot(5,1,1)
hold on
for j = 1:3
    plot(1:width(gradient),q(j,1:end),'o');
    plot(1:width(gradient),q(j,1:end));
end
title("points")
hold off
grid on

subplot(5,1,2)
hold on
for j=1:3
    plot(1:width(cost),ai(j,:),'o');
    plot(1:width(cost),ai(j,:));
end
plot(1:width(cost),pboundary(1:end));
plot(1:width(cost),nboundary(1:end));
title("acceleration")
hold off
grid on

subplot(5,1,3)
hold on
for i=1:3
    plot(1:width(cost),cost(i,:),'o');
    plot(1:width(cost),cost(i,:));
end
title("cost")
hold off
grid on

subplot(5,1,4)
hold on
for i = 1:3
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot(1:width(gradient),gradient(i,:),'o');
    plot(1:width(gradient),gradient(i,:));
end
title("gradient")
hold off
grid on

subplot(5,1,5)
hold on
for i = 1:width(gradient)
    t_grad(i) = sum(gradient(:,i));
end
plot(1:width(gradient),t_grad,'o');
plot(1:width(gradient),t_grad);
title("gradient total")
hold off
grid on

figure(2)
surf(X,Y,Z)
title("Feasibility Cost")
xlabel("X")
ylabel("Y")
zlabel("Cost")
%% C++ code
% ts = bspline_interval_;
% ts_inv2 = 1 / ts / ts;
%
% /* acceleration feasibility */
% for (int i = 0; i < q.cols() - 2; i++)
% {
%   Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;
% 
%   //cout << "temp_a * ai=" ;
%   for (int j = 0; j < 3; j++)
%   {
%     if (ai(j) > max_acc_)
%     {
%       // cout << "zx-todo ACC" << endl;
%       // cout << ai(j) << endl;
%       cost += pow(ai(j) - max_acc_, 2);
% 
%       gradient(j, i + 0) += 2 * (ai(j) - max_acc_) * ts_inv2;
%       gradient(j, i + 1) += -4 * (ai(j) - max_acc_) * ts_inv2;
%       gradient(j, i + 2) += 2 * (ai(j) - max_acc_) * ts_inv2;
%     }
%     else if (ai(j) < -max_acc_)
%     {
%       cost += pow(ai(j) + max_acc_, 2);
% 
%       gradient(j, i + 0) += 2 * (ai(j) + max_acc_) * ts_inv2;
%       gradient(j, i + 1) += -4 * (ai(j) + max_acc_) * ts_inv2;
%       gradient(j, i + 2) += 2 * (ai(j) + max_acc_) * ts_inv2;
%     }
%     else
%     {
%       /* code */
%     }
%   }