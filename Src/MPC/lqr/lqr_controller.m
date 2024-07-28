classdef lqr_controller < handle
    properties(Access = public)
        u0(1,1) double;
        mu(1,1) double;
        K;
        integrator;
        prevError;
        prevUAV;
        isCaptured;
        returnGoal;
        dt;
        t;
        seenUAV;
    end


    methods(Access = public)
        function obj = lqr_controller(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.mu = quadrotor.mu;
            obj.integrator = [0;0;0];
            obj.prevError = [0;0;0];
            obj.prevUAV = [0;0;0];
            obj.isCaptured = false;
            obj.dt = 0.01; % TODO pass this in the contructor
            obj.t = 0;
            obj.seenUAV = false;

            % Points to linearize around
            zs = zeros(12, 1);
            zs(3) = 1;
            us = obj.u0 * ones([4,1]); % enough to counteract gravity
            [A, B] = linearize(zs, us);
            C = ctrb(A,B);
            % Check if the system is controllable
            if rank(C) < min(size(C))
                error("This linearized system is not controllable!\n")
            end

            x1_gain = 1;
            x2_gain = 1;
            x3_gain = 1;
            a1_gain = 1000;
            a2_gain = 1000;
            a3_gain = 0.001;
            v1_gain = 1;
            v2_gain = 1;
            v3_gain = 1;
            w1_gain = 1;
            w2_gain = 1;
            w3_gain = 0.001;

            
            Q = diag([x1_gain x2_gain x3_gain a1_gain a2_gain a3_gain ...
                      v1_gain v2_gain v3_gain w1_gain w2_gain w3_gain]);
            u_gain = 1/obj.mu^2;
            R = eye(4,4) * u_gain;
            sys = ss(A, B, eye(12), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.K = K;
        end

        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            % State and Control Values to Linearize Around
            obj.t = obj.t + obj.dt;

            if obj.seenUAV 
                % need to have one sample to calculate UAV velocity
                v = uavPos-obj.prevUAV;
            else
                v = [0;0;0];
            end


            zd = zeros(12, 1);
            if isCaptured && ~obj.isCaptured
                obj.isCaptured = true;
                obj.returnGoal = [0;0;z(3)];
                fprintf('UAV Captured!!!\n')
            end

            if ~obj.isCaptured
                uavPos
                predTarget = uavPos + 3*v/obj.dt
                zd(1:3) = predTarget;
                zd(4:6) = v;
                
            else
                zd(1:3) = obj.returnGoal;
                if norm(z(1:3)-obj.returnGoal) < 0.1
                    obj.returnGoal = [0;0;0.1;];
                end
            end 


            u = @(z) obj.u0*ones([4,1])+obj.K*(zd-z);
            u = u(z);
            obj.prevUAV = uavPos;
            obj.seenUAV = true;
        end
   end
            
            
end