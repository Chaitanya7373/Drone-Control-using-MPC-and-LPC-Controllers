classdef lqr_controller_relinearize < handle
    properties(Access = public)
        u0(1,1) double;
        mu(1,1) double;
        A;
        B;
        K_prev;
    end


    methods(Access = public)
        function obj = lqr_controller_relinearize(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.mu = quadrotor.mu;
            [A, B] = linearize_symbolic();
            obj.A = A;
            obj.B = B;
            obj.K_prev = zeros([4,12]);

        end
        
        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            % State and Control Values to Linearize Around
            zd = zeros(12, 1);
            if ~isCaptured
                zd(1:3) = uavPos;
            else
                fprintf('UAV Captured!!!\n')
            end

            % Cost values for LQR
            x1_gain = 10;
            x2_gain = 10;
            x3_gain = 20;
            a1_gain = 100;
            a2_gain = 100;
            a3_gain = 200;
            v1_gain = 1;
            v2_gain = 1;
            v3_gain = 2;
            w1_gain = 1;
            w2_gain = 1;
            w3_gain = 2;

            % Use the u from the previous K to linearize?
            us = -obj.K_prev*(z-zd);
            [A,B] = eval_linear(obj.A, obj.B, z, us);
            C = ctrb(A,B);
            % Check if the system is controllable
            if rank(C) < min(size(C))
                error("This linearized system is not controllable!\n")
            end
            
            Q = diag([x1_gain x2_gain x3_gain a1_gain a2_gain a3_gain ...
                      v1_gain v2_gain v3_gain w1_gain w2_gain w3_gain]);
            u_gain = 1/obj.mu^2;
            R = eye(4,4) * u_gain;
            [K,~,~] = lqr(A,B,Q,R);
            obj.K_prev = K;

            u = @(z) -obj.K*(z-zd);
            u = u(z);

        end
   end
            
            
end