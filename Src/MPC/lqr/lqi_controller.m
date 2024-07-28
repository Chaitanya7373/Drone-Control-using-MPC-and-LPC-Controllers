classdef lqi_controller < handle
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
        function obj = lqi_controller(quadrotor)
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

            % https://math.stackexchange.com/questions/3294817/lqr-with-augumented-state-design
            % integrator on only position variables
            Aa = [A, zeros(12, 3); -eye(3), zeros(3, 12)];
            Ba = [B; zeros(3, 4)];
             if rank(ctrb(Aa, Ba)) < 15
                 rank(ctrb(Aa, Ba))
                 error('Augmented linear system for LQI is uncontrollable')
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
            i1_gain = 1;
            i2_gain = 1;
            i3_gain = 1;
            
            Q = diag([x1_gain x2_gain x3_gain a1_gain a2_gain a3_gain ...
                      v1_gain v2_gain v3_gain w1_gain w2_gain w3_gain ...
                      i1_gain i2_gain i3_gain]);
            u_gain = 1/obj.mu^2;
            R = eye(4,4) * u_gain;
            sys = ss(Aa, Ba, eye(15), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.K = K;
        end

        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            % State and Control Values to Linearize Around
            obj.t = obj.t + obj.dt;

            % Error for integrator
            errorE = uavPos-z(1:3);

            % Rotation matrix to convert error to UAV frame

            R_UE = eye(3); % for the case where we dont calculate it
            if obj.seenUAV 
                % need to have one sample to calculate UAV velocity
                v = uavPos-obj.prevUAV;
                b1 = v;
                b2 = cross(b1, errorE);
                b3 = cross(v, b2);
                
                R_UE = [b1 b2 b3];
                error = R_UE\errorE;
            else
                error = errorE;
                v = [0;0;0];
            end

            if ~isCaptured && obj.t > 5
                % trapezoidal approximation of integral of error
                trap = obj.dt*0.5*(obj.prevError + error);
                obj.integrator = obj.integrator + trap;
            end

            zd = zeros(12, 1);
            if isCaptured && ~obj.isCaptured
                obj.isCaptured = true;
                obj.integrator = [0;0;0];
                obj.returnGoal = [0;0;z(3)];
                fprintf('UAV Captured!!!\n')
            end

            if ~obj.isCaptured
                zd(1:3) = uavPos + 2*v/obj.dt;
            else
                zd(1:3) = obj.returnGoal;
                if norm(z(1:3)-obj.returnGoal) < 0.1
                    obj.returnGoal = [0;0;0.1;];
                end
            end 
            integrate = R_UE*obj.integrator;
            % prevent explosion of integral :(
            maxErrInt = 10;
            integrate(integrate>maxErrInt) = maxErrInt;
            integrate(integrate<-maxErrInt) = -maxErrInt;
            u = @(z) obj.u0*ones([4,1])+obj.K*([zd; 0; 0; 0;]-[z; integrate]);
            u = u(z);
            obj.prevError = error;
            obj.prevUAV = uavPos;
            obj.seenUAV = true;
        end
   end
            
            
end