classdef lqr_controller < handle
    properties(Access = public)
        u0(1,1) double;
        mu(1,1) double;
        K;
        
    end


    methods(Access = public)
        function obj = lqr_controller(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.mu = quadrotor.mu;

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

            % Cost values for LQR
            x1_gain = 10;
            x2_gain = 10;
            x3_gain = 10;
            a1_gain = 1000;
            a2_gain = 1000;
            a3_gain = 1;
            v1_gain = 1;
            v2_gain = 1;
            v3_gain = 1;
            w1_gain = 1;
            w2_gain = 1;
            w3_gain = 1;
            
            Q = diag([x1_gain x2_gain x3_gain a1_gain a2_gain a3_gain ...
                      v1_gain v2_gain v3_gain w1_gain w2_gain w3_gain]);
            u_gain = 1/obj.mu^2;
            R = eye(4,4) * u_gain;
            [K,~,~] = lqr(A,B,Q,R);
            obj.K = K;
        end
        
        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            zd = zeros(12, 1);
            if ~isCaptured
                zd(1:3) = uavPos;
            end 
            u = @(z) obj.u0*ones([4,1])+obj.K*(zd-z);
            u = u(z);
        end
   end
            
            
end