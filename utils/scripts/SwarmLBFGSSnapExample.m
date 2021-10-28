%% Swarm Cost
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

%% Smoothness Cost
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

%% Idealogy
% 1. Use a uniform multisegment spline to define our trajectory
% 2. Find Control Points (ctrlpts) and also the Trajectory (
% 3. 

%% Main Code
clc
clear all
close all 
figt = figure('Name','Optimization','NumberTitle','off');

% Global Variables for UAVs
global nquad; % Number of Quads 
nquad = 5; % Number of Quads 
ID = 1:nquad; % ID of Quads 
global iter;
factor = 2.5;

% Simulation Boundary
ctr = [2;2;5]; % Center of Boundary Sphere
rad = [25.0;25.0;10.0]; % Radius of Boundary Sphere
zlim = 0.1; % z limit cannot be more than 1
sprd = 25; % Spread of plotting

% Global Simulation parameters
totaltime = 15;
dt = 0.1;
int = totaltime/dt; % Interval
timeint = linspace(0,totaltime,int);
traj{nquad} = []; % s1 is the trajectory that is being used and defined

% Parameters for CA / Planning module
ca_int = 0.5;

poi = cell([1 3]); 
poiS = cell([1 3]); 
pathx = cell([1 3]); 
poiO = cell([nquad-1 3]);
col = cell([nquad-1 3]);
colcount = cell([nquad-1 3]);

for i = 1:width(colcount)*height(colcount)
    colcount{i}=zeros(1);
end
rscan = [2;2;1];
uthres = 1.5; sthres = 1.0;
dist = 4;

ellipsoidFunction = @(p,c,r) ((p(1,1) - c(1,1))^2 / r(1,1)^2 + ...
    (p(2,1) - c(2,1))^2 / r(2,1)^2 + ...
    (p(3,1) - c(3,1))^2 / r(3,1)^2);

[xbnd, ybnd, zbnd] = PlotSphere(rad,ctr,10);

global splinet;
global psplinet;
global ctrlpts;
global pctrlpts;
global traj;
global globalint;
global ntraj;

splinet = linspace(0,totaltime,round(int/factor));
globalint = linspace(0,totaltime,int);
for i = 1:nquad
    StartPose(i,:) = RandSphereCoord(rad,ctr,zlim);
    if i > 1
        while norm(StartPose(i,:) - StartPose(i-1,:)) < dist
            StartPose(i,:) = RandSphereCoord(rad,ctr,zlim);
        end
    end
    EndPose(i,:) = [ctr(1)-(StartPose(i,1) - ctr(1)), ...
        ctr(2)-(StartPose(i,2) - ctr(2)), ...
        ctr(3)-(StartPose(i,3) - ctr(3))];
    % xinit format = (time, x, y, z ,xd, yd, zd, xdd, ydd, zdd)
    xinit{i}(1,:) = linspace(0,totaltime,int);
    xStart{i} = [StartPose(i,:)';0;0;0;0;0;0];
    xEnd{i} = [EndPose(i,:)';0;0;0;0;0;0];
    
    slope0 = [0;0;0];
    slopeF = [0;0;0];
    for b = 1:3
        ctrlpts{i}(b,:) = linspace(StartPose(i,b),EndPose(i,b),round(int/factor));
        q(b,:) = spline(splinet,[slope0(b); ctrlpts{i}(b,:)'; slopeF(b)],xinit{i}(1,:));
        for c=1:width(q(b,:))-1
            qd_temp(c) = (q(b,c+1) - q(b,c))/dt;
        end
        qd_new = [slope0(b),qd_temp];
        q(b+3,:) = qd_new;
    end
    traj{i} = q;
    %  [s1{i},maxvel{i}] = bvp(dt, totaltime, xEnd{i}, xStart{i});
    %  xinit{i}(2:10,:) = s1{i};

    poi{i} = zeros(nquad-1,2);
    pathx{i} = zeros(3,1);   
    exitflag{i} = [0;0;0];
end

pctrlpts = ctrlpts;
psplinet = splinet;
ntraj = traj;

options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false);

%% Recursive Loop
fprintf('Starting UAV Path and Loop...\n');
global time;
time = 0;
ca_iter = 0;
iter = 1;


stop = false;
global curr;

while ~stop 
    tic
    
    % Replanning 
    if time - ca_iter * ca_int > ca_int
        clear ctrlpts traj 
        splinet = psplinet(:,psplinet>time);
        splinet = [time,splinet]; 
        globalint = globalint(:,globalint>=time);
        for i = 1:nquad
            clear n1 qd_temp
            curr = i;
            %  [s1{i},maxvel{i}] = bvp(dt, totaltime - time, s1{i}(:,end), s1{i}(:,iter{i}+1));
            %  [f,g] = costfun(s1{i}(1:3,:));
            [x1,fval2,exitflag{i}(ca_iter+1),output{i}] = fminlbfgs(@costfun,pctrlpts{i}(:,psplinet>time),options);
            opttime{i}(ca_iter+1) = output{i}.timeTotal;
            % x1 are the new control points
            qpoint = ntraj{i}(1:6,iter);
            ctrlpts{i} = [qpoint(1:3),x1];

            nslope0 = [qpoint(4);qpoint(5);qpoint(6)];
            nslopeF = [0;0;0];
  
            for b = 1:3
                n1(b,:) = spline(splinet,[nslope0(b); ctrlpts{i}(b,:)'; nslopeF(b)],globalint);
                for c=1:width(n1(b,:))-1
                    qd_temp(c) = (n1(b,c+1) - n1(b,c))/dt;
                end
                qd_new = [slope0(b),qd_temp];
                n1(b+3,:) = qd_new;
            end
            traj{i} = n1;
            qrec{i,ca_iter+1} = ctrlpts{i};
        end
        ca_iter = ca_iter + 1;
    end
    
    pctrlpts = ctrlpts; 
    psplinet = splinet;
    ntraj = traj;
    clf
    %  figure(figt);
    %  plot3(xbnd,ybnd,zbnd,'.','MarkerSize',1);
    hold on
    for i = 1:nquad
        plot3(StartPose(i,1),StartPose(i,2),StartPose(i,3),'x','MarkerSize',5);
        plot3(traj{i}(1,iter),traj{i}(2,iter),traj{i}(3,iter),'x','MarkerSize',15)
        plot3(traj{i}(1,:),traj{i}(2,:),traj{i}(3,:),'.','MarkerSize',1);
        plot3(ctrlpts{i}(1,:),ctrlpts{i}(2,:),ctrlpts{i}(3,:),'o','MarkerSize',3);
    end
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
    axis([ctr(1)-sprd  ctr(1)+sprd  ctr(2)-sprd ctr(1)+sprd ctr(3)-sprd  ctr(3)+sprd]);
    grid on
    title(sprintf('Iteration: %d, time: %4.2f', iter, time));
    t = toc;
    
    for i = 1:nquad
        for num = setdiff(1:nquad, i)
            dist_vec = traj{i}(1:3,iter) - traj{num}(1:3,iter);
            a = 1.0; b = 1.0; inv_a2 = 1 / a / a; inv_b2 = 1 / b / b; clearance = 0.5;
            ellip_dist = sqrt(dist_vec(3)^2 * inv_a2 + (dist_vec(1)^2 + dist_vec(2)^2) * inv_b2);
            dist_err = clearance - ellip_dist;
            if dist_err > 0
                stop = true; 
                fprintf('Collision! Self(%d) Other(%d)\n',i,num);
            end
        end
        % Stopping criteria
        % 10s -> 25 
        % 15s -> 40
        if width(pctrlpts{i}) < 40
            stop = true;
        end
        tmp = ntraj{i}(1:6,iter);
        qplot{i}(:,iter) = [time;tmp];
    end
    % Set delay to real time 
    if t < dt
       pause(dt-t); 
    end
    time = time + dt;
    iter = iter + 1;
end


%% State plot
f2 = figure(2)
% Left Bottom Width Height
f2.Position = [100 50 500 500];
for i = 1:nquad 
    subplot(nquad,2,i+1*(i-1))
    plot(qplot{i}(1,:),qplot{i}(2,:),qplot{i}(1,:),qplot{i}(3,:),qplot{i}(1,:),qplot{i}(4,:));
    xlabel('t [s]'); ylabel('Pos [m]');
    grid on
    title(sprintf('Position UAV %d', i));
    subplot(nquad,2,2*i)
    plot(qplot{i}(1,:),qplot{i}(5,:),qplot{i}(1,:),qplot{i}(6,:),qplot{i}(1,:),qplot{i}(7,:));
    xlabel('t [s]'); ylabel('Vel [m/s]');
    grid on
    title(sprintf('Velocity UAV %d', i));
    
end

f3 = figure(3);
f3.Position = [600 50 400 500];
for i = 1:nquad 
    subplot(nquad,1,i)
    plot(1:ca_iter,opttime{i},1:ca_iter,opttime{i},'o');
    xlabel('idx'); ylabel('t [s]');
    grid on
    title(sprintf('CA Time UAV %d', i));
end

f = 2;
f4 = figure(4);
f4.Position = [1000 50 900 900];
% Maybe we evaluate 3 out of all paths
for c = 1:ca_iter-5
    subplot(round((ca_iter-5)/f),f,c)
    for i = 1:nquad
        hold on
        plot3(qrec{i,c}(1,:),qrec{i,c}(2,:),qrec{i,c}(3,:),'o','MarkerSize',3);
    end
    hold off
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
    axis([ctr(1)-sprd  ctr(1)+sprd  ctr(2)-sprd ctr(1)+sprd ctr(3)-sprd  ctr(3)+sprd]);
    grid on
    title(sprintf('Control points iteration %d/%d', c,ca_iter));
end


%% Functions
function [f,g] = myfun(x)
    f = sum(sin(x) + 3);
    f = sum(f);
    if ( nargout > 1 ), g = cos(x); end
end

function [f,g] = costfun(x)
    % Function returning the value of global variable “x”
    global curr;
    global nquad;
    global time;
    global pctrlpts;
    global psplinet;
    
    i = curr;
    gain_swarm = 0.2; gain_smooth = 0.1;

    offset = sum(psplinet<time);
    
    % According to EGO Swarm (Swarm Cost)
    cost = 0;
    gradient = zeros(3,width(x));
    a = 2.0; b = 1.0; inv_a2 = 1 / a / a; inv_b2 = 1 / b / b; clearance = 1.2;
    for num = setdiff(1:nquad, i)
        for l = 1:width(x)
            % The count differs for other UAVs since they are
            % sequentially checked and trajectories are updated in
            % sequence
            
            dist_vec = x(1:3,l) - pctrlpts{num}(1:3,l+offset);
            ellip_dist = sqrt(dist_vec(3)^2 * inv_a2 + (dist_vec(1)^2 + dist_vec(2)^2) * inv_b2);
            dist_err = clearance - ellip_dist;

            gcoeff(1) = -2 * (clearance / ellip_dist - 1) * inv_b2;
            gcoeff(2) = gcoeff(1);
            gcoeff(3) = -2 * (clearance / ellip_dist - 1) * inv_a2;

            if (dist_err < 0)
                % do nothing
                gradient(:,l) = gradient(:,l) + [0;0;0];
            else
                cost = cost + gain_swarm * dist_err^2;
                %  We need to add the gradient together for each of the
                %  columns for the different UAVs
                gradient(:,l) = gradient(:,l) + gain_swarm * gcoeff' .* dist_vec;
            end
        end
    end
    
    % According to EGO Swarm (Smoothness Cost)
    for j = 1:width(x)-3
        jerk = x(:,j + 3) - 3 * x(:,j + 2) + 3 * x(:,j + 1) - x(:,j);
        cost = cost + gain_smooth * norm(jerk)^2;
        temp_j = 2.0 * jerk;

        gradient(:,j + 0) = gradient(:,j + 0) + gain_smooth * (-temp_j);
        gradient(:,j + 1) = gradient(:,j + 1) + gain_smooth * (3.0 * temp_j);
        gradient(:,j + 2) = gradient(:,j + 2) + gain_smooth * (-3.0 * temp_j);
        gradient(:,j + 3) = gradient(:,j + 3) + gain_smooth * (temp_j);
    end

    f = cost;
    g = gradient;
end

function [s1,maxvel] = bvp(dt, runtime, dest, s0)
    p0 = s0(1:3);
    v0 = s0(4:6);
    a0 = s0(7:9);
    pf = dest(1:3);
    vf = dest(4:6);
    af = dest(7:9);
    timeint = linspace(0,runtime,runtime/dt);

    % s0 = [p0, v0, a0];
    % sf = [pf, vf, af];
    T = runtime;

    delta =  [(pf - p0 -v0 * T - 0.5 * a0 * T^2) ; ...
              (vf - v0 -a0 * T) ; ...
              (af - a0)];
    m = [720, -360*T, 60*T^2 ; -360*T, 168*T^2, -24*T^3 ; 60*T^2, -24*T^3, 3*T^4];
    M = zeros(length(m)*width(m),length(m)*width(m));
    for i = 1:length(m)*width(m)
        M1 = eye(3);
        l = mod(i,length(m)) + length(m)*(~mod(i,length(m)));
        h = ceil(i/3);
        M1(M1==1) = m(l,h);
        %  fprintf('(%d) %d, %d\n',i,1+(l-1)*length(m),1+(h-1)*length(m));
        M(1+(l-1)*length(m):1+(l-1)*length(m)+2,1+(h-1)*length(m):1+(h-1)*length(m)+2) = M1;
    end

    abg = 1/T^5 * M * delta;
    alpha = abg(1:3); beta = abg(4:6); gamma = abg(7:9);

    for ts = 1:width(timeint)
        t = timeint(ts);
        sOpt(:,ts) = [(alpha/120 * t^5 + beta/24 * t^4 + gamma/6 * t^3 + a0/2 * t^2 + v0 * t + p0); ...
                         (alpha/24 * t^4 + beta/6 * t^3 + gamma/2 * t^2 + a0 * t + v0); ...
                         (alpha/6 * t^3 + beta/2 * t^2 + gamma * t + a0)];
    end
    s1 = sOpt;
    maxvel = max(sOpt(4,:)+sOpt(5,:)+sOpt(6,:));
end

function pos = RandSphereCoord(r,c,zlim)
    % rand(1)*2-1 gives -1 to 1
    theta = (rand(1)*2-1)*pi;
    phi = (rand(1)*2*zlim-zlim)*pi/2;
    cosphi = cos(phi);
    sintheta = sin(theta);
    pos(1) = r(1)*cosphi*cos(theta)  + c(1);
    pos(2) = r(2)*cosphi*sintheta  + c(2);
    pos(3) = r(3)*sin(phi)  + c(3);
end

function [x, y, z] = PlotSphere(r,c,n)
    theta = (-n:2:n)/n*pi;
    phi = (-n:2:n)'/n*pi/2;
    cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
    sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
    x = r(1)*cosphi*cos(theta) + c(1);
    y = r(2)*cosphi*sintheta + c(2);
    z = r(3)*sin(phi)*ones(1,n+1) + c(3);
end

%% LBFGS Functions
function [x,fval,exitflag,output,grad]=fminlbfgs(funfcn,x_init,optim)
    %FMINLBFGS finds a local minimum of a function of several variables. 
    %   This optimizer is developed for image registration methods with large 
    %	amounts of unknown variables.
    %
    %   Optimization methods supported:
    %	- Quasi Newton Broyden芳letcher萌oldfarb亡hanno (BFGS)  
    %   - Limited memory BFGS (L-BFGS)
    %   - Steepest Gradient Descent optimization.
    %   
    %   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINLBFGS(FUN,X0,OPTIONS) 
    %
    %   Inputs,
    %		FUN: Function handle or string which is minimized, returning an
    %				error value and optional the error gradient. 
    %		X0: Initial values of unknowns can be a scalar, vector or matrix
    %	 (optional)
    %		OPTIONS: Structure with optimizer options, made by a struct or
    %				optimset. (optimset doesnot support all input options)
    %
    %   Outputs,
    %		X : The found location (values) which minimize the function.
    %		FVAL : The minimum found
    %		EXITFLAG : Gives value, which explain why the minimizer stopt
    %		OUTPUT : Structure with all important ouput values and parameters
    %		GRAD : The gradient at this location 
    %
    %   Extended description of input/ouput variables 
    %   OPTIONS,
    %		OPTIONS.GoalsExactAchieve : If set to 0, a line search method is
    %               used which uses a few function calls to do a good line
    %               search. When set to 1 a normal line search method with Wolfe 
    %				conditions is used (default).
    %		OPTIONS.GradConstr, Set this variable to true if gradient calls are
    %				cpu-expensive (default). If false more gradient calls are 
    %				used and less function calls.
    %	    OPTIONS.HessUpdate : If set to 'bfgs', Broyden芳letcher萌oldfarb亡hanno 
    %				optimization is used (default), when the number of unknowns is 
    %				larger then 3000 the function will switch to Limited memory BFGS, 
    %				or if you set it to 'lbfgs'. When set to 'steepdesc', steepest 
    %				decent optimization is used.
    %		OPTIONS.StoreN : Number of itterations used to approximate the Hessian,
    %			 	in L-BFGS, 20 is default. A lower value may work better with
    %				non smooth functions, because than the Hessian is only valid for
    %				a specific position. A higher value is recommend with quadratic equations. 
    %		OPTIONS.GradObj : Set to 'on' if gradient available otherwise finited difference
    %				is used.
    %     	OPTIONS.Display : Level of display. 'off' displays no output; 'plot' displays
    %				all linesearch results in figures. 'iter' displays output at  each 
    %               iteration; 'final' displays just the final output; 'notify' 
    %				displays output only if the function does not converge; 
    %	    OPTIONS.TolX : Termination tolerance on x, default 1e-6.
    %	    OPTIONS.TolFun : Termination tolerance on the function value, default 1e-6.
    %		OPTIONS.MaxIter : Maximum number of iterations allowed, default 400.
    % 		OPTIONS.MaxFunEvals : Maximum number of function evaluations allowed, 
    %				default 100 times the amount of unknowns.
    %		OPTIONS.DiffMaxChange : Maximum stepsize used for finite difference gradients.
    %		OPTIONS.DiffMinChange : Minimum stepsize used for finite difference gradients.
    %		OPTIONS.OutputFcn : User-defined function that an optimization function calls
    %				at each iteration.
    %		OPTIONS.rho : Wolfe condition on gradient (c1 on wikipedia), default 0.01.
    %		OPTIONS.sigma : Wolfe condition on gradient (c2 on wikipedia), default 0.9. 
    %		OPTIONS.tau1 : Bracket expansion if stepsize becomes larger, default 3.
    %		OPTIONS.tau2 : Left bracket reduction used in section phase,
    %		default 0.1.
    %		OPTIONS.tau3 : Right bracket reduction used in section phase, default 0.5.
    %   FUN,
    %		The speed of this optimizer can be improved by also providing
    %   	the gradient at X. Write the FUN function as follows
    %   	function [f,g]=FUN(X)
    %       	f , value calculation at X;
    %   	if ( nargout > 1 )
    %       	g , gradient calculation at X;
    %   	end
    %	EXITFLAG,
    %		Possible values of EXITFLAG, and the corresponding exit conditions
    %		are
    %  		1, 'Change in the objective function value was less than the specified tolerance TolFun.';
    %  		2, 'Change in x was smaller than the specified tolerance TolX.'; 
    %  		3, 'Magnitude of gradient smaller than the specified tolerance';
    %  		4, 'Boundary fminimum reached.';
    %  		0, 'Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
    %  		-1, 'Algorithm was terminated by the output function.';
    %  		-2, 'Line search cannot find an acceptable point along the current search';
    %
    %   Examples
    %       options = optimset('GradObj','on');
    %       X = fminlbfgs(@myfun,2,options)
    %
    %   	% where myfun is a MATLAB function such as:
    %       function [f,g] = myfun(x)
    %       f = sin(x) + 3;
    %	    if ( nargout > 1 ), g = cos(x); end
    %
    %   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
    %
    %   Function is written by D.Kroon University of Twente (Updated Nov. 2010)
    % Read Optimalisation Parameters
    defaultopt = struct('Display','final','HessUpdate','bfgs','GoalsExactAchieve',1,'GradConstr',true,  ...
                'TolX',1e-6,'TolFun',1e-6,'GradObj','off','MaxIter',400,'MaxFunEvals',100*numel(x_init)-1,  ...
                'DiffMaxChange',1e-1,'DiffMinChange',1e-8,'OutputFcn',[], ...
                'rho',0.0100,'sigma',0.900,'tau1',3,'tau2', 0.1, 'tau3', 0.5,'StoreN',20);
    if (~exist('optim','var')) 
        optim=defaultopt;
    else
        f = fieldnames(defaultopt);
        for i=1:length(f),
            if (~isfield(optim,f{i})||(isempty(optim.(f{i})))), optim.(f{i})=defaultopt.(f{i}); end
        end
    end

    % Initialize the data structure
    data.fval=0;
    data.gradient=0;
    data.fOld=[]; 
    data.xsizes=size(x_init);
    data.numberOfVariables = numel(x_init);
    data.xInitial = x_init(:);
    data.alpha=1;
    data.xOld=data.xInitial; 
    data.iteration=0;
    data.funcCount=0;
    data.gradCount=0;
    data.exitflag=[];
    data.nStored=0;
    data.timeTotal=tic;
    data.timeExtern=0;
    % Switch to L-BFGS in case of more than 3000 unknown variables
    if(optim.HessUpdate(1)=='b') 
        if(data.numberOfVariables<3000), 
            optim.HessUpdate='bfgs';
        else
            optim.HessUpdate='lbfgs';
        end
    end
    if(optim.HessUpdate(1)=='l')
        succes=false;
        while(~succes)
            try
                data.deltaX=zeros(data.numberOfVariables,optim.StoreN);
                data.deltaG=zeros(data.numberOfVariables,optim.StoreN);
                data.saveD=zeros(data.numberOfVariables,optim.StoreN);
                succes=true;
            catch ME
                warning('fminlbfgs:memory','Decreasing StoreN value because out of memory');
                succes=false;
                data.deltaX=[]; data.deltaG=[]; data.saveD=[];
                optim.StoreN=optim.StoreN-1;
                if(optim.StoreN<1)
                    rethrow(ME);
                end
            end
        end
    end
    exitflag=[];
    % Display column headers
    if(strcmp(optim.Display,'iter'))
        disp('     Iteration  Func-count   Grad-count         f(x)         Step-size');
    end
    % Calculate the initial error and gradient
    data.initialStepLength=1;
    [data,fval,grad]=gradient_function(data.xInitial,funfcn, data, optim);
    data.gradient=grad;
    data.dir = -data.gradient;
    data.fInitial = fval;
    data.fPrimeInitial= data.gradient'*data.dir(:);
    data.fOld=data.fInitial;
    data.xOld=data.xInitial;
    data.gOld=data.gradient;

    gNorm = norm(data.gradient,Inf);  % Norm of gradient
    data.initialStepLength = min(1/gNorm,5); 
    % Show the current iteration
    if(strcmp(optim.Display,'iter'))
            s=sprintf('     %5.0f       %5.0f       %5.0f       %13.6g    ',data.iteration,data.funcCount,data.gradCount,data.fInitial); disp(s);
    end

    % Hessian intialization
    if(optim.HessUpdate(1)=='b')
        data.Hessian=eye(data.numberOfVariables);
    end
    % Call output function
    if(call_output_function(data,optim,'init')), exitflag=-1; end

    % Start Minimizing
    while(true)
        % Update number of itterations
        data.iteration=data.iteration+1; 
        % Set current lineSearch parameters
        data.TolFunLnS = eps(max(1,abs(data.fInitial )));
        data.fminimum = data.fInitial - 1e16*(1+abs(data.fInitial));

        % Make arrays to store linesearch results
        data.storefx=[]; data.storepx=[]; data.storex=[]; data.storegx=[];
        % If option display plot, than start new figure
        if(optim.Display(1)=='p'), figure, hold on; end

        % Find a good step size in the direction of the gradient: Linesearch
        if(optim.GoalsExactAchieve==1)
            data=linesearch(funfcn, data,optim);
        else
            data=linesearch_simple(funfcn, data, optim);
        end

        % Make linesearch plot
        if(optim.Display(1)=='p'); 
            plot(data.storex,data.storefx,'r*');
            plot(data.storex,data.storefx,'b');

            alpha_test= linspace(min(data.storex(:))/3, max(data.storex(:))*1.3, 10);
            falpha_test=zeros(1,length(alpha_test));
            for i=1:length(alpha_test)
                [data,falpha_test(i)]=gradient_function(data.xInitial(:)+alpha_test(i)*data.dir(:),funfcn, data, optim);
            end    
            plot(alpha_test,falpha_test,'g');
            plot(data.alpha,data.f_alpha,'go','MarkerSize',8);
        end

        % Check if exitflag is set
        if(~isempty(data.exitflag)),
            exitflag=data.exitflag;
            data.xInitial=data.xOld; 
            data.fInitial=data.fOld;
            data.gradient=data.gOld;
            break, 
        end;

        % Update x with the alpha step
        data.xInitial = data.xInitial + data.alpha*data.dir;

        % Set the current error and gradient
        data.fInitial =  data.f_alpha;
        data.gradient = data.grad;

        % Set initial steplength to 1
        data.initialStepLength = 1;


        gNorm = norm(data.gradient,Inf);  % Norm of gradient

        % Set exit flags 
        if(gNorm <optim.TolFun), exitflag=1; end
        if(max(abs(data.xOld-data.xInitial)) <optim.TolX), exitflag=2; end
        if(data.iteration>=optim.MaxIter), exitflag=0; end

        % Check if exitflag is set
        if(~isempty(exitflag)), break, end;
        % Update the inverse Hessian matrix
        if(optim.HessUpdate(1)~='s')
            % Do the Quasi-Neton Hessian update.
            data = updateQuasiNewtonMatrix_LBFGS(data,optim);
        else
            data.dir = -data.gradient;
        end

        % Derivative of direction
        data.fPrimeInitial= data.gradient'*data.dir(:);
        % Call output function
        if(call_output_function(data,optim,'iter')), exitflag=-1; end

        % Show the current iteration
        if(strcmp(optim.Display(1),'i')||strcmp(optim.Display(1),'p'))
            s=sprintf('     %5.0f       %5.0f       %5.0f       %13.6g   %13.6g',data.iteration,data.funcCount,data.gradCount,data.fInitial,data.alpha); disp(s);
        end

        % Keep the variables for next iteration
        data.fOld=data.fInitial;
        data.xOld=data.xInitial;
        data.gOld=data.gradient;
    end
    % Set output parameters
    fval=data.fInitial;
    grad=data.gradient;
    x = data.xInitial;
    % Reshape x to original shape
    x=reshape(x,data.xsizes);
    % Call output function
    if(call_output_function(data,optim,'done')), exitflag=-1; end
    % Make exist output structure
    if(optim.HessUpdate(1)=='b'), output.algorithm='Broyden芳letcher萌oldfarb亡hanno (BFGS)';
    elseif(optim.HessUpdate(1)=='l'), output.algorithm='limited memory BFGS (L-BFGS)';
    else output.algorithm='Steepest Gradient Descent'; 
    end
    output.message=getexitmessage(exitflag);
    output.iteration = data.iteration;
    output.funccount = data.funcCount;
    output.fval = data.fInitial;
    output.stepsize = data.alpha;
    output.directionalderivative = data.fPrimeInitial;
    output.gradient = reshape(data.gradient, data.xsizes);
    output.searchdirection = data.dir;
    output.timeTotal=toc(data.timeTotal);    
    output.timeExtern=data.timeExtern;
    oupput.timeIntern=output.timeTotal-output.timeExtern;
    % Display final results
    if(~strcmp(optim.Display,'off'))
        disp('    Optimizer Results')
        disp(['        Algorithm Used: ' output.algorithm]);
        disp(['        Exit message : ' output.message]);
        disp(['        iterations : '  int2str(data.iteration)]);
        disp(['        Function Count : ' int2str(data.funcCount)]);
        disp(['        Minimum found : ' num2str(fval)]);
        disp(['        Intern Time : ' num2str(oupput.timeIntern) ' seconds']);
        disp(['        Total Time : ' num2str(output.timeTotal) ' seconds']);
    end
end

function message=getexitmessage(exitflag)
    switch(exitflag)
        case 1, message='Change in the objective function value was less than the specified tolerance TolFun.';
        case 2, message='Change in x was smaller than the specified tolerance TolX.'; 
        case 3, message='Magnitude of gradient smaller than the specified tolerance';
        case 4, message='Boundary fminimum reached.';
        case 0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
        case -1, message='Algorithm was terminated by the output function.';
        case -2, message='Line search cannot find an acceptable point along the current search';
        otherwise, message='Undefined exit code';
    end
end
    
function stopt=call_output_function(data,optim,where)
    stopt=false;
    if(~isempty(optim.OutputFcn))
        output.iteration = data.iteration;
        output.funccount = data.funcCount;
        output.fval = data.fInitial;
        output.stepsize = data.alpha;
        output.directionalderivative = data.fPrimeInitial;
        output.gradient = reshape(data.gradient, data.xsizes);
        output.searchdirection = data.dir;
        stopt=feval(optim.OutputFcn,reshape(data.xInitial,data.xsizes),output,where); 
    end
end
        
	
function data=linesearch_simple(funfcn, data, optim)
    % Find a bracket of acceptable points
    data = bracketingPhase_simple(funfcn, data, optim);
    if (data.bracket_exitflag  == 2)
      % BracketingPhase found a bracket containing acceptable points; 
      % now find acceptable point within bracket
      data = sectioningPhase_simple(funfcn, data, optim);
      data.exitflag = data.section_exitflag; 
    else
      % Already acceptable point found or MaxFunEvals reached
      data.exitflag = data.bracket_exitflag; 
    end
end

function data = bracketingPhase_simple(funfcn, data,optim)
    % Number of itterations
    itw=0; 
    % Point with smaller value, initial
    data.beta=0; 
    data.f_beta=data.fInitial; 
    data.fPrime_beta=data.fPrimeInitial;
    % Initial step is equal to alpha of previous step.
    alpha = data.initialStepLength;
    % Going up hill
    hill=false;
    % Search for brackets
    while(true)
        % Calculate the error registration gradient
        if(optim.GradConstr)
            [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            fPrime_alpha=nan;
            grad=nan;
        else
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            fPrime_alpha = grad'*data.dir(:);
        end

        % Store values linesearch
        data.storefx=[data.storefx f_alpha]; 
        data.storepx=[data.storepx fPrime_alpha]; 
        data.storex=[data.storex alpha]; 
        data.storegx=[data.storegx grad(:)];

        % Update step value
        if(data.f_beta<f_alpha), 
            % Go to smaller stepsize
            alpha=alpha*optim.tau3;

            % Set hill variable
            hill=true;
        else
            % Save current minium point
            data.beta=alpha; data.f_beta=f_alpha; data.fPrime_beta=fPrime_alpha; data.grad=grad;
            if(~hill)
                alpha=alpha*optim.tau1;  
            end
        end

        % Update number of loop iterations
        itw=itw+1; 

        if(itw>(log(optim.TolFun)/log(optim.tau3))),
          % No new optium found, linesearch failed.
          data.bracket_exitflag=-2; break; 
        end

        if(data.beta>0&&hill)
                % Get the brackets around minimum point
                % Pick bracket A from stored trials
                [t,i]=sort(data.storex,'ascend');
                storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
                [t,i]=find(storex>data.beta,1);
                if(isempty(i)), [t,i]=find(storex==data.beta,1); end
                alpha=storex(i); f_alpha=storefx(i); fPrime_alpha=storepx(i);

                % Pick bracket B from stored trials
                [t,i]=sort(data.storex,'descend');
                storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
                [t,i]=find(storex<data.beta,1);
                if(isempty(i)), [t,i]=find(storex==data.beta,1); end
                beta=storex(i); f_beta=storefx(i); fPrime_beta=storepx(i);

                % Calculate derivatives if not already calculated
                if(optim.GradConstr)
                    gstep=data.initialStepLength/1e6; 
                    if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
                    if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
                    [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
                    [data,f_beta2]=gradient_function(data.xInitial(:)+(beta+gstep)*data.dir(:),funfcn, data, optim);
                    fPrime_alpha=(f_alpha2-f_alpha)/gstep;
                    fPrime_beta=(f_beta2-f_beta)/gstep;
                end
                % Set the brackets A and B
                data.a=alpha; data.f_a=f_alpha; data.fPrime_a=fPrime_alpha;
                data.b=beta; data.f_b=f_beta; data.fPrime_b=fPrime_beta;

                % Finished bracketing phase
                data.bracket_exitflag  = 2; return
        end
        % Reached max function evaluations
        if(data.funcCount>=optim.MaxFunEvals), data.bracket_exitflag=0; return; end
    end
end
    
function data = sectioningPhase_simple(funfcn, data, optim)
    % Get the brackets
    brcktEndpntA=data.a; brcktEndpntB=data.b;
    % Calculate minimum between brackets
    [alpha,f_alpha_estimated] = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  
    if(isfield(data,'beta')&&(data.f_beta<f_alpha_estimated)), alpha=data.beta; end
    [t,i]=find(data.storex==alpha,1);
    if((~isempty(i))&&(~isnan(data.storegx(i))))
        f_alpha=data.storefx(i); grad=data.storegx(:,i);
    else
        % Calculate the error and gradient for the next minimizer itteration
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        if(isfield(data,'beta')&&(data.f_beta<f_alpha)), 
            alpha=data.beta; 
            if((~isempty(i))&&(~isnan(data.storegx(i))))
                f_alpha=data.storefx(i); grad=data.storegx(:,i);
            else
                [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            end
        end
    end
    % Store values linesearch
    data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];
    fPrime_alpha = grad'*data.dir(:);
    data.alpha=alpha; 
    data.fPrime_alpha= fPrime_alpha; 
    data.f_alpha= f_alpha;
    data.grad=grad;
    % Set the exit flag to succes   
    data.section_exitflag=[];
end

function data=linesearch(funfcn, data, optim)
    % Find a bracket of acceptable points
    data = bracketingPhase(funfcn, data,optim);
    if (data.bracket_exitflag  == 2)
      % BracketingPhase found a bracket containing acceptable points; 
      % now find acceptable point within bracket
      data = sectioningPhase(funfcn, data, optim);
      data.exitflag = data.section_exitflag; 
    else
      % Already acceptable point found or MaxFunEvals reached
      data.exitflag = data.bracket_exitflag; 
    end
end

function data = sectioningPhase(funfcn, data, optim)
    %
    % sectioningPhase finds an acceptable point alpha within a given bracket [a,b] 
    % containing acceptable points. Notice that funcCount counts the total number of 
    % function evaluations including those of the bracketing phase. 
    while(true)

        % Pick alpha in reduced bracket
        brcktEndpntA = data.a + min(optim.tau2,optim.sigma)*(data.b - data.a); 
        brcktEndpntB = data.b - optim.tau3*(data.b - data.a);

        % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree 
        % polynomial that interpolates f() and f'() at "a" and at "b".
        alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  
        % No acceptable point could be found
        if (abs( (alpha - data.a)*data.fPrime_a ) <= data.TolFunLnS), data.section_exitflag = -2; return; end

        % Calculate value (and gradient if no extra time cost) of current alpha
        if(~optim.GradConstr)
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            fPrime_alpha = grad'*data.dir(:);
        else
            gstep=data.initialStepLength/1e6; 
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
            fPrime_alpha=(f_alpha2-f_alpha)/gstep;
        end
        % Store values linesearch 
        data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 

        % Store current bracket position of A
        aPrev = data.a; 
        f_aPrev = data.f_a; 
        fPrime_aPrev = data.fPrime_a; 
        % Update the current brackets
        if ((f_alpha > data.fInitial + alpha*optim.rho*data.fPrimeInitial) || (f_alpha >= data.f_a))
            % Update bracket B to current alpha
            data.b = alpha; data.f_b = f_alpha; data.fPrime_b = fPrime_alpha;
        else
            % Wolfe conditions, if true then acceptable point found 
            if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial), 
                if(optim.GradConstr)
                    % Gradient was not yet calculated because of time costs
                    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
                    fPrime_alpha = grad'*data.dir(:);
                end
                % Store the found alpha values
                data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
                data.grad=grad;
                data.section_exitflag = []; return, 
            end

            % Update bracket A
            data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;

            if (data.b - data.a)*fPrime_alpha >= 0
                % B becomes old bracket A;
                data.b = aPrev; data.f_b = f_aPrev;  data.fPrime_b = fPrime_aPrev;
            end
        end

        % No acceptable point could be found
        if (abs(data.b-data.a) < eps), data.section_exitflag = -2; return, end
        % maxFunEvals reached
        if(data.funcCount >optim.MaxFunEvals), data.section_exitflag = -1; return, end
    end
end
    
function data = bracketingPhase(funfcn, data, optim)
    % bracketingPhase finds a bracket [a,b] that contains acceptable points; a bracket 
    % is the same as a closed interval, except that a > b is allowed.
    %
    % The outputs f_a and fPrime_a are the values of the function and the derivative 
    % evaluated at the bracket endpoint 'a'. Similar notation applies to the endpoint 
    % 'b'. 
    % Parameters of bracket A
    data.a = []; 
    data.f_a = []; 
    data.fPrime_a = []; 
    % Parameters of bracket B
    data.b = []; 
    data.f_b = []; 
    data.fPrime_b = [];
    % First trial alpha is user-supplied
    % f_alpha will contain f(alpha) for all trial points alpha
    % fPrime_alpha will contain f'(alpha) for all trial points alpha
    alpha = data.initialStepLength;
    f_alpha = data.fInitial;              
    fPrime_alpha = data.fPrimeInitial;    
    % Set maximum value of alpha (determined by fminimum)
    alphaMax = (data.fminimum - data.fInitial)/(optim.rho*data.fPrimeInitial); 
    alphaPrev = 0;
    while(true) 
      % Evaluate f(alpha) and f'(alpha)
      fPrev = f_alpha;
      fPrimePrev = fPrime_alpha;

      % Calculate value (and gradient if no extra time cost) of current alpha
      if(~optim.GradConstr)
          [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
          fPrime_alpha = grad'*data.dir(:);
      else
          gstep=data.initialStepLength/1e6;
          if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
          if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
          [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
          [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
          fPrime_alpha=(f_alpha2-f_alpha)/gstep;
      end

      % Store values linesearch 
      data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 

      % Terminate if f < fminimum
      if (f_alpha <= data.fminimum), data.bracket_exitflag = 4; return; end

      % Bracket located - case 1 (Wolfe conditions)
      if (f_alpha > (data.fInitial + alpha*optim.rho*data.fPrimeInitial)) || (f_alpha >= fPrev)
        % Set the bracket values
        data.a = alphaPrev; data.f_a = fPrev;  data.fPrime_a = fPrimePrev;
        data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
        % Finished bracketing phase
        data.bracket_exitflag  = 2; return 
      end
      % Acceptable steplength found
      if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial), 
          if(optim.GradConstr)
              % Gradient was not yet calculated because of time costs
              [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
              fPrime_alpha = grad'*data.dir(:);
          end
          % Store the found alpha values
          data.alpha=alpha;
          data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha; data.grad=grad;
          % Finished bracketing phase, and no need to call sectioning phase
          data.bracket_exitflag = [];  return 
      end

      % Bracket located - case 2  
      if (fPrime_alpha >= 0)
        % Set the bracket values
        data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
        data.b = alphaPrev; data.f_b = fPrev; data.fPrime_b = fPrimePrev;
        % Finished bracketing phase
        data.bracket_exitflag  = 2; return
      end

      % Update alpha
      if (2*alpha - alphaPrev < alphaMax )
          brcktEndpntA = 2*alpha-alphaPrev; 
          brcktEndpntB = min(alphaMax,alpha+optim.tau1*(alpha-alphaPrev));
          % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial 
          % that interpolates f() and f'() at alphaPrev and at alpha
          alphaNew = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alphaPrev,alpha,fPrev, ...
                                             fPrimePrev,f_alpha,fPrime_alpha,optim);
          alphaPrev = alpha;
          alpha = alphaNew;
      else
          alpha = alphaMax;
      end
      % maxFunEvals reached
      if(data.funcCount >optim.MaxFunEvals), data.bracket_exitflag = -1; return, end
    end
end

function [alpha,f_alpha]= pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2,optim)
    % finds a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial 
    % that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1, 
    % f(alpha2) = f2, f'(alpha2) = fPrime2.
    % determines the coefficients of the cubic polynomial with c(alpha1) = f1, 
    % c'(alpha1) = fPrime1, c(alpha2) = f2, c'(alpha2) = fPrime2.
    coeff = [(fPrime1+fPrime2)*(alpha2-alpha1)-2*(f2-f1) ...
        3*(f2-f1)-(2*fPrime1+fPrime2)*(alpha2-alpha1) (alpha2-alpha1)*fPrime1 f1];
    % Convert bounds to the z-space
    lowerBound = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
    upperBound = (brcktEndpntB - alpha1)/(alpha2 - alpha1);
    % Swap if lowerbound is higher than the upperbound
    if (lowerBound  > upperBound), t=upperBound; upperBound=lowerBound; lowerBound=t; end 
    % Find minima and maxima from the roots of the derivative of the polynomial.
    sPoints = roots([3*coeff(1) 2*coeff(2) coeff(3)]); 
    % Remove imaginaire and points outside range
    sPoints(imag(sPoints)~=0)=[]; 
    sPoints(sPoints<lowerBound)=[]; sPoints(sPoints>upperBound)=[];
    % Make vector with all possible solutions
    sPoints=[lowerBound sPoints(:)' upperBound];
    % Select the global minimum point
    [f_alpha,index]=min(polyval(coeff,sPoints)); z=sPoints(index);
    % Add the offset and scale back from [0..1] to the alpha domain
    alpha = alpha1 + z*(alpha2 - alpha1);
    % Show polynomial search
    if(optim.Display(1)=='p'); 
        vPoints=polyval(coeff,sPoints);
        plot(sPoints*(alpha2 - alpha1)+alpha1,vPoints,'co');
        plot([sPoints(1) sPoints(end)]*(alpha2 - alpha1)+alpha1,[vPoints(1) vPoints(end)],'c*');
        xPoints=linspace(lowerBound/3, upperBound*1.3, 50);
        vPoints=polyval(coeff,xPoints);
        plot(xPoints*(alpha2 - alpha1)+alpha1,vPoints,'c');
    end
end
	
function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
    % Call the error function for error (and gradient)
    if ( nargout <3 )
        timem=tic;   
        fval=funfcn(reshape(x,data.xsizes)); 
        data.timeExtern=data.timeExtern+toc(timem);
        data.funcCount=data.funcCount+1;
    else
        if(strcmp(optim.GradObj,'on'))
            timem=tic;    
            [fval, grad]=feval(funfcn,reshape(x,data.xsizes)); 
            data.timeExtern=data.timeExtern+toc(timem);
            data.funcCount=data.funcCount+1;
            data.gradCount=data.gradCount+1;
        else
            % Calculate gradient with forward difference if not provided by the function
            grad=zeros(length(x),1);
            fval=funfcn(reshape(x,data.xsizes));
            gstep=data.initialStepLength/1e6; 
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            for i=1:length(x),
                x_temp=x; x_temp(i)=x_temp(i)+gstep;
                timem=tic;    
                [fval_g]=feval(funfcn,reshape(x_temp,data.xsizes)); data.funcCount=data.funcCount+1;
                data.timeExtern=data.timeExtern+toc(timem);
                grad(i)=(fval_g-fval)/gstep;
            end
        end
        grad=grad(:);
    end
end
    
function data = updateQuasiNewtonMatrix_LBFGS(data,optim)
    % updates the quasi-Newton matrix that approximates the inverse to the Hessian.
    % Two methods are support BFGS and L-BFGS, in L-BFGS the hessian is not
    % constructed or stored.
    % Calculate position, and gradient diference between the
    % itterations
    deltaX=data.alpha* data.dir;
    deltaG=data.gradient-data.gOld;

    if ((deltaX'*deltaG) >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaG) ))
        if(optim.HessUpdate(1)=='b')
            % Default BFGS as described by Nocedal
            p_k = 1 / (deltaG'*deltaX);
            Vk = eye(data.numberOfVariables) - p_k*deltaG*deltaX';
            % Set Hessian
            data.Hessian = Vk'*data.Hessian *Vk + p_k * deltaX*deltaX';
            % Set new Direction
            data.dir = -data.Hessian*data.gradient;
        else
            % L-BFGS with scaling as described by Nocedal

            % Update a list with the history of deltaX and deltaG
            data.deltaX(:,2:optim.StoreN)=data.deltaX(:,1:optim.StoreN-1); data.deltaX(:,1)=deltaX;
            data.deltaG(:,2:optim.StoreN)=data.deltaG(:,1:optim.StoreN-1); data.deltaG(:,1)=deltaG;

            data.nStored=data.nStored+1; if(data.nStored>optim.StoreN), data.nStored=optim.StoreN; end
            % Initialize variables
            a=zeros(1,data.nStored);
            p=zeros(1,data.nStored);
            q = data.gradient;
            for i=1:data.nStored
                p(i)= 1 / (data.deltaG(:,i)'*data.deltaX(:,i));
                a(i) = p(i)* data.deltaX(:,i)' * q;
                q = q - a(i) * data.deltaG(:,i);
            end
            % Scaling of initial Hessian (identity matrix)
            p_k = data.deltaG(:,1)'*data.deltaX(:,1) / sum(data.deltaG(:,1).^2); 

            % Make r = - Hessian * gradient
            r = p_k * q;
            for i=data.nStored:-1:1,
                b = p(i) * data.deltaG(:,i)' * r;
                r = r + data.deltaX(:,i)*(a(i)-b);
            end

            % Set new direction
            data.dir = -r;
        end
    end
end