classdef quadcopter < handle
    %QUADCOPTER class
    %     s = [x; y; z; xd; yd; zd; qW; qX; qY; qZ; p; q; r];
    %     sdot = [xd; yd; zd; xdd; ydd; zdd; qWd; qXd; qYd; qZd; pd; qd; rd];
    %     state_hist = [t; s]
    %     d_s = [pos; vel; acc; yaw; yawdot] size (11,1)
    
    properties
        nquad           % number of quadcopters in simulation
        id              % id of the quadcopter
        s               % current state of the quadcopter
        d_s             % desired state of the quadcopter
        s_h             % state history of the quadcopter
        terminal        % terminal positions of the quadcopter
        int             % interval that the quadcopter is running on
        t               % current time
        isIdeal         % bool for whether the quadcopter uses dynamic states or ideal states from path
        order           % order of bspline
        
        path            % path in [time position velocity] format
        path0           % original path in [time position velocity] format
        c               % clearance zone for the quadcopter
        cp              % current cp without time and includes clamped cp
        ccp             % current cp that the uav is using
        cp0             % control points at the start
        wp              % waypoints
        wpt             % waypoints with time
        cpt             % control points with time
        cpt0            % control points with time at t0
        ph              % planning horizon
        
        g               % gravity
        m               % mass
        arm_l           % arm length
        I               % inertia
        maxF            % maximum body force trust
        minF            % minimum body force trust
        
        % position controller params
        Kp              % [gainx;gainy;gainz]
        Kd              % [gainx;gainy;gainz]

        % attitude controller params
        KpM             %
        KdM             %
        
        % Other debug parameters
        dt
        intv
        seg_total_xyz
        seg_indi_xyz
        colcount        % Collision count for the uav
        col             % Collision check
        
    end
    
    methods
        
        %% Constructor and Initialize
        function q = quadcopter(nquad, ...
                id, ...
                param, ...
                state, ...
                dt, ...
                isIdeal, ...
                wp)
            if (nargin~=7)
                error("Not correct initiating size");
            end
            
            q.nquad = nquad;
            q.id = id;
            q.s = state;
            q.g = id;
            q.int = dt;
            q.d_s = zeros(11,1);
            q.d_s(1:6) = state(1:6);
            q.s_h = zeros(14,1);
            q.s_h(2:14) = state;
            q.isIdeal = isIdeal;
            q.wp = wp;
            q.col = zeros(1,3);
            q.colcount = 0;
            
            q.m = param.m;
            q.arm_l = param.arm_length;
            q.I = param.I;
            q.c = param.clearance;
            q.Kp = param.kp;
            q.Kd = param.kd;
            q.KpM = param.kpm;
            q.KdM = param.kdm;
            q.maxF = param.maxF;
            q.minF = param.minF;
        end  
        
        %% Set Functions
        function setTerminal(self, terminal)
            self.terminal = terminal;
        end
        
        function setPath(self, path)
            self.path = path;
            self.path0 = path;
        end
        
        function setControlPoint(self, cp, t, time_horizon, order, wp_t)
            self.cp0 = cp;
            self.cpt(1,:) = t;
            self.order = order;
            self.cpt(2:4,:) = cp(:,order:end);
            self.cpt0 = self.cpt;
            current = 0;
            i = find(self.cpt(1,:) <= current + time_horizon);
            self.ccp = self.cpt(:,i);
            self.ph = time_horizon;
            self.cp = self.cp0;
            self.wpt(1,:) = wp_t;
            self.wpt(2:4,:) = self.wp;
        end
        
        function setOtherDebugParam(self,seg,seg_total,dt,intv)
            self.seg_indi_xyz = seg;
            self.seg_total_xyz = seg_total;
            self.dt = dt;
            self.intv = intv;
        end
        
        %% Update Functions
        function updatePath(self, iter)
            % Take one iteration before so that we do not lose the current
            % segment that the uav is in
            time = (iter-1) * self.int;
            iter_arr = find(self.path(1,:) >= time);
            % Truncate the path
            self_tmp = self.path(1:7,iter_arr);
            self.path = [];
            self.path = self_tmp;
        end
        
        % [updateControlPoint] :Update control points so as to sort out data for the replanning
        % module
        function updateControlPoint(self, iter)
            % Use time to get sort out the data
            % Get start and end according to planning horizon
            time = (iter) * self.int;
            iter_arr_start = find(self.cpt(1,:) >= time);
            iter_arr_end = find(self.cpt(1,:) <= time + self.ph);
            % Do a check here for the size of the array
            if width(iter_arr_start(1):iter_arr_end(end)) ~= width(self.ccp)
                fprintf("[updateControlPoint] input %d output %d mismatch\n", ...
                    width(iter_arr_start(1):iter_arr_end(end)), width(self.ccp));
            end
            
            % Get the current cp from cpt (full control points with time)
            self_tmp = self.cpt(1:4,iter_arr_start(1):iter_arr_end(end)); 
            self.ccp = [];
            self.ccp = self_tmp;
        end
        
        function replanControlPoint(self, ccp)
            iter_arr_start = find(self.cpt(1,:) >= ccp(1,1));
            iter_arr_end = find(self.cpt(1,:) <= ccp(1,end));
            % Must initiate a check here to force out any inconsistencies
            % between [ccp size] and the [iter_arr start and end size]
            if width(iter_arr_start(1):iter_arr_end(end)) ~= width(ccp)
                fprintf("[replanControlPoint0] input %d output %d mismatch\n", ...
                    width(iter_arr_start(1):iter_arr_end(end)), width(ccp));
            end
            
            % No saving of variables here, wipe and overwrite
            self.ccp = [];
            self.ccp = ccp;
            self.cpt(:,iter_arr_start(1):iter_arr_end(end)) = ccp;
            
            if width(self.cp0) ~= width(self.cp)
                fprintf("[replanControlPoint1.1] input %d output %d mismatch\n", ...
                    width(self.cp0), width(self.cp));
            end
                
            % [Error] function here sometimes give extra size for cp
            % Change cp to plan ahead
            self.cp(:,iter_arr_start(1)+self.order-1:iter_arr_end(end)+self.order-1) = ...
                self.cpt(2:4,iter_arr_start(1):iter_arr_end(end));
            
            if width(self.cp0) ~= width(self.cp)
                fprintf("[replanControlPoint1.2] input %d output %d mismatch\n", ...
                    width(self.cp0), width(self.cp));
            end
            
        end
        
        function replanPath(self, path, ts)
            iter_arr_start = find(self.path(1,:) >= ts(1));
            iter_arr_end = find(self.path(1,:) <= ts(2));
            self.path(:,iter_arr_start(1):iter_arr_end(end)) = path;
        end
            
        function updateState(self, iter)
            self.t = iter * self.int;
            iter_arr = find(self.path(1,:) >= self.t);
            i = iter_arr(1);
            self.d_s(1:6) = self.path(2:7,i); 
          
            if self.isIdeal
                state = zeros(13,1);
                state(1:6) = self.d_s(1:6);
                self.s = state;
                state_n_t = [self.t ; state];
                self.s_h = [self.s_h state_n_t];
                
            else
                % We give a arbitary range here > 2
                mm = 5; % Because ODE requires a range of t
                tl = linspace(self.t, self.t + self.int, mm);
                [t_tmp, s_tmp] = ode45(@(t,s) self.quadcopterEOM(t,s), tl, self.s);
                self.s = s_tmp(mm-1,:)';
                state_n_t = [self.t ; self.s];
                self.s_h = [self.s_h state_n_t];
            end
        end
        
        %% Get Functions
        function ccp = getControlPointEval(self)
            ccp = self.ccp;
        end
        
        function state = getState(self)
            state = self.s;
        end

        function col = getCollisionCheck(self, other)
            % get the distance vector
            dv = self.s(1:3) - other.s(1:3); 
            a = 1.0; b = 1.0; inv_a2 = 1 / a / a; inv_b2 = 1 / b / b; 
            e_d = sqrt(dv(3)^2 * inv_a2 + (dv(1)^2 + dv(2)^2) * inv_b2); % get the ellipsoidal distance 
            % get the distance error
            d_e = self.c - e_d; 
            if d_e > 0
                col = true;
                self.col(self.colcount + 1,1:3) = self.s(1:3)';
                self.colcount = self.colcount + 1;
            else 
                col = false;
            end
        end
        
        %% Plot Functions
        
        function plotCurrentPose(self)
            plot3(self.s(1), ...
                self.s(2), ...
                self.s(3), ...
                'x','MarkerSize',10); % Current pose
        end
        
        function plotFuturePath(self)
            iter_arr = find(self.path(1,:) >= self.t);
            iter = iter_arr(1);
            plot3(self.path(2,iter:end), ...
                self.path(3,iter:end), ...
                self.path(4,iter:end), ...
                'o','MarkerSize',2); % Displays the path
        end
        
        function plotCurrentPath(self)
            plot3(self.s_h(2,1:end), ...
                self.s_h(3,1:end), ...
                self.s_h(4,1:end), ...
                '-','MarkerSize',2); % Displays the path
        end
        
        function plotTerminalPose(self)
            plot3(self.terminal(1,1), ...
                self.terminal(2,1), ...
                self.terminal(3,1), ...
                'o','MarkerSize',10); % Start pose
            plot3(self.terminal(1,end), ...
                self.terminal(2,end), ...
                self.terminal(3,end), ...
                'o','MarkerSize',10); % End pose
        end
        
        function plotWaypoint(self)
            plot3(self.wp(1,:), ...
                self.wp(2,:), ...
                self.wp(3,:), ...
                'h','MarkerSize',4); % Waypoint pose
        end
        
        function plotControlPoint0(self)
            plot3(self.cp0(1,:), ...
                self.cp0(2,:), ...
                self.cp0(3,:), ...
                '+','MarkerSize',2); % Control points
        end
        
        function plotControlPoint1(self)
            plot3(self.cpt(2,:), ...
                self.cpt(3,:), ...
                self.cpt(4,:), ...
                '+','MarkerSize',3); % Control points
        end
        
        function plotStateHistory(self, idx)
            subplot(self.nquad,2,idx+1*(idx-1))
            plot(self.s_h(1,:),self.s_h(2,:), ...
                self.s_h(1,:),self.s_h(3,:), ...
                self.s_h(1,:),self.s_h(4,:));
            xlabel('t [s]'); ylabel('Pos [m]');
            grid on
            title(sprintf('Position UAV %d', idx));
            subplot(self.nquad,2,2*idx)
            plot(self.s_h(1,:),self.s_h(5,:), ...
                self.s_h(1,:),self.s_h(6,:), ...
                self.s_h(1,:),self.s_h(7,:));
            xlabel('t [s]'); ylabel('Vel [m/s]');
            grid on
            title(sprintf('Velocity UAV %d', idx));
        end
        
        function plotTotalVelHistory(self, idx)
            subplot(self.nquad,2,idx+1*(idx-1))
            vel = abs(self.s_h(5,:)) + ...
                abs(self.s_h(6,:));
            plot(self.s_h(1,:),vel);
            xlabel('t [s]'); ylabel('Hor Vel [m/s]');
            grid on
            title(sprintf('Total Hor Vel UAV %d', idx));
            
            subplot(self.nquad,2,idx+1*(idx-1))
            plot(self.s_h(1,:),abs(self.s_h(7,:)));
            xlabel('t [s]'); ylabel('Ver Vel [m/s]');
            grid on
            title(sprintf('Total Ver Vel UAV %d', idx));
        end
        
        function plotTotalAccHistory(self, idx)
            subplot(self.nquad,2,idx+1*(idx-1))
            vel = abs(self.s_h(5,:)) + ...
                abs(self.s_h(6,:));
            plot(self.s_h(1,:),vel);
            xlabel('t [s]'); ylabel('Hor Vel [m/s]');
            grid on
            title(sprintf('Total Hor Vel UAV %d', idx));
            
            subplot(self.nquad,2,idx+1*(idx-1))
            plot(self.s_h(1,:),abs(self.s_h(7,:)));
            xlabel('t [s]'); ylabel('Ver Vel [m/s]');
            grid on
            title(sprintf('Total Ver Vel UAV %d', idx));
        end
        
        function plotCollision(self)
            for i = 1:self.colcount
                plot3(self.col(i,1), ...
                    self.col(i,2), ...
                    self.col(i,3), ...
                    'o','MarkerSize',7); % collision point
            end
        end
    end
    
    methods(Access = private)    
        function invI = getInvI(I)
            invI = inv(I);
        end
        
        function R = QuatToRot(self, q)
            % QuatToRot Converts a Quaternion to Rotation matrix
            % by Daniel Mellinger
            
            % normalize q
            q = q ./ sqrt(sum(q.^2));
            qahat(1,2) = -q(4);
            qahat(1,3) = q(3);
            qahat(2,3) = -q(2);
            qahat(2,1) = q(4);
            qahat(3,1) = -q(3);
            qahat(3,2) = q(2);

            R = eye(3) + 2 * qahat * qahat + 2 * q(1) * qahat;
        end
        
        function [phi,theta,psi] = RotToRPY_ZXY(self, R)
            %RotToRPY_ZXY Extract Roll, Pitch, Yaw from a world-to-body Rotation Matrix
            %   The rotation matrix in this function is world to body [bRw] you will
            %   need to transpose the matrix if you have a body to world [wRb] such
            %   that [wP] = [wRb] * [bP], where [bP] is a point in the body frame and
            %   [wP] is a point in the world frame
            % by Daniel Mellinger

            phi = asin(R(2,3));
            psi = atan2(-R(2,1)/cos(phi),R(2,2)/cos(phi));
            theta = atan2(-R(1,3)/cos(phi),R(3,3)/cos(phi));

        end
        
        %% EOM
        function sdot = quadcopterEOM(self, t, s)
            pos_des = self.d_s(1:3);
            vel_des = self.d_s(4:6);

            yaw_des = self.d_s(10);
            yawd_des = self.d_s(11);
            
            pos = self.s(1:3);
            vel = self.s(4:6);
            q_r = self.s(7:10);
            omega = self.s(11:13);
            
            rot = self.QuatToRot(q_r');
            [phi, theta, yaw] = self.RotToRPY_ZXY(rot);
            euler = [phi; theta; yaw];
            
            acc_des = self.Kd .* ( ...
                vel_des - vel) + ...
                self.Kp.*(pos_des - pos);
            
            % Desired roll, pitch and yaw
            phi_des = 1/self.g * (acc_des(1)*sin(yaw_des) - acc_des(2)*cos(yaw_des));
            theta_des = 1/self.g * (acc_des(1)*cos(yaw_des) + acc_des(2)*sin(yaw_des));
            psi_des = yaw_des;
            
            euler_des = [phi_des; theta_des; psi_des];
            pqr_des = [0; 0; yawd_des];
            
            % Thurst
            F  = self.m * (self.g + acc_des(3));
            
            % Moment
            M =  self.I * ...
                (self.KdM.*(pqr_des - omega) + ...
                self.KpM.*(euler_des - euler));
            
            % Limit the force and moments due to actuator limits
            A = [0.25,                      0, -0.5/self.arm_l;
                 0.25,  0.5/self.arm_l,                      0;
                 0.25,                      0,  0.5/self.arm_l;
                 0.25, -0.5/self.arm_l,                      0];

            prop_thrusts = A*[F;M(1:2)]; % Not using moment about Z-axis for limits
            prop_thrusts_clamped = max(min(prop_thrusts, self.maxF/4), self.minF/4);

            B = [                 1,                 1,                 1,                  1;
                                  0, self.arm_l,                 0, -self.arm_l;
                 -self.arm_l,                 0, self.arm_l,                 0];
            F = B(1,:)*prop_thrusts_clamped;
            M = [B(2:3,:)*prop_thrusts_clamped; M(3)];
            
            % Assign states
            x = s(1);
            y = s(2);
            z = s(3);
            xdot = s(4);
            ydot = s(5);
            zdot = s(6);
            qW = s(7);
            qX = s(8);
            qY = s(9);
            qZ = s(10);
            p = s(11);
            q = s(12);
            r = s(13);

            quat = [qW; qX; qY; qZ];
            bRw = self.QuatToRot(quat);
            wRb = bRw';

            % Acceleration
            accel = 1 / self.m * (wRb * [0; 0; F] - [0; 0; self.m * self.g]);

            % Angular velocity
            K_quat = 2; %this enforces the magnitude 1 constraint for the quaternion
            quaterror = 1 - (qW^2 + qX^2 + qY^2 + qZ^2);
            qdot = -1/2*[0, -p, -q, -r;...
                         p,  0, -r,  q;...
                         q,  r,  0, -p;...
                         r, -q,  p,  0] * quat + K_quat * quaterror * quat;

            % Angular acceleration
            omega = [p;q;r];
            pqrdot = inv(self.I) * (M - cross(omega, self.I * omega));

            % Assemble sdot
            sdot = zeros(13,1);
            sdot(1)  = xdot;
            sdot(2)  = ydot;
            sdot(3)  = zdot;
            sdot(4)  = accel(1);
            sdot(5)  = accel(2);
            sdot(6)  = accel(3);
            sdot(7)  = qdot(1);
            sdot(8)  = qdot(2);
            sdot(9)  = qdot(3);
            sdot(10) = qdot(4);
            sdot(11) = pqrdot(1);
            sdot(12) = pqrdot(2);
            sdot(13) = pqrdot(3);
            
        end
    end
end

