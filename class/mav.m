classdef mav < handle
    %QUADCOPTER class
    %     s = [x; y; z; xd; yd; zd; qW; qX; qY; qZ; p; q; r];
    %     sdot = [xd; yd; zd; xdd; ydd; zdd; qWd; qXd; qYd; qZd; pd; qd; rd];
    %     state_hist = [t; s]
    %     d_s = [pos; vel; acc; yaw; yawdot] size (11,1)
    
    properties
        nquad           % number of mav in simulation
        id              % id of the mav
        s               % current state of the mav
        d_s             % desired state of the mav
        s_h             % state history of the mav
        terminal        % terminal positions of the mav
        int             % interval that the mav is running on
        t               % current time
        isIdeal         % bool for whether the mav uses dynamic states or ideal states from path
        order           % order of bspline
        
        path            % path in [time position velocity] format
        path0           % original path in [time position velocity] format
        c               % clearance zone for the mav
        cp              % current cp without time and includes clamped cp
        ccp             % current cp that the uav is using
        cp0             % control points at the start
        wp              % waypoints
        wpt             % waypoints with time
        cpt             % control points with time
        cpt0            % control points with time at t0
        ph              % planning horizon
        
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
        function q = mav(nquad, ...
                id, ...
                state, ...
                dt, ...
                isIdeal, ...
                wp, ...
                clearance)
            if (nargin~=7)
                error("Not correct initiating size");
            end
            
            q.nquad = nquad;
            q.id = id;
            q.c = clearance;
            q.s = state;
            q.int = dt;
            q.d_s = state;
            q.s_h = zeros(10,1);
            q.s_h(2:10) = state;
            q.isIdeal = isIdeal;
            q.wp = wp;
            q.col = zeros(1,3);
            q.colcount = 0;
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
            % self.cpt(2:4,:) = cp(:,order:end);
            self.cpt(2:4,:) = cp(:,order+1:end-order+1);
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
            self.d_s(1:9) = self.path(2:10,i); 
          
            if self.isIdeal
                self.s = self.d_s;
                state_n_t = [self.t ; self.d_s];
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
            xlabel('t [s]'); ylabel('HorV [m/s]');
            grid on
            title(sprintf('Hor Vel UAV %d', idx));
            
            subplot(self.nquad,2,2*idx)
            plot(self.s_h(1,:),abs(self.s_h(7,:)));
            xlabel('t [s]'); ylabel('VerV [m/s]');
            grid on
            title(sprintf('Ver Vel UAV %d', idx));
        end
        
        function plotTotalAccHistory(self, idx)
            subplot(self.nquad,1,idx)
            hold on
            plot(self.s_h(1,:),self.s_h(8,:));
            plot(self.s_h(1,:),self.s_h(9,:));
            plot(self.s_h(1,:),self.s_h(10,:));
            xlabel('t [s]'); ylabel('Acc [m/s^2]');
            grid on
            title(sprintf('Acc UAV %d', idx));
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
    
end

