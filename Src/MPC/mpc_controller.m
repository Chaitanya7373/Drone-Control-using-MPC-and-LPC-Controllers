classdef mpc_controller < handle
    properties(Access = public)
        u0(1,1) double;
        mu(1,1) double;
        K;
        mpcobj;
        quadrotor;
        options;
        xc;
    end

    methods(Access = public)
        function obj = mpc_controller(quadrotor)
            obj.quadrotor = quadrotor;  % Store the quadrotor object
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

            mpcverbosity off;
            sys = ss(A, B, eye(12), 0);
           
             % Convert continuous-time model to discrete-time with a default sample time
            Ts = 0.1; % Default sample time (adjust as needed)
            % sysd = c2d(sys, Ts);
            % [K,~,~] = mpc(sys);
            mpcobj = mpc(sys,Ts);
            % obj.K = mpcobj.K;
             % Set mpc object properties
            mpcobj.PredictionHorizon = 10;
            mpcobj.ControlHorizon = 2;
            mpcobj.Weights.ManipulatedVariables = [2 2 2 2];
            mpcobj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0.1];
            mpcobj.Weights.OutputVariables = [2 2 2 20 20 20 0 0 0 0 0 0];
            mpcobj.ManipulatedVariables(1).Min = 0;
            mpcobj.ManipulatedVariables(1).Max = 3;
            mpcobj.ManipulatedVariables(2).Min = 0;
            mpcobj.ManipulatedVariables(2).Max = 3;
            mpcobj.ManipulatedVariables(3).Min = 0;
            mpcobj.ManipulatedVariables(3).Max = 3;
            mpcobj.ManipulatedVariables(4).Min = 0;
            mpcobj.ManipulatedVariables(4).Max = 3;
            
            mpcobj.MV(1).Target = obj.u0;
            mpcobj.MV(2).Target = obj.u0;
            mpcobj.MV(3).Target = obj.u0;
            mpcobj.MV(4).Target = obj.u0;


            obj.mpcobj = mpcobj;
            obj.xc = mpcstate(obj.mpcobj);

        end
        
        % mv = mpcmove(mpcobj,xc,ym,r,v);
        function setOptions(obj, options)
            % Method to set custom options
            if isa(options, 'mpcmoveopt')
                obj.options = options;
            else
                error('Options must be an instance of mpcmoveopt.');
            end
        end

        function u = output(obj, isCaptured, z, uavPos)
             % State and Control Values to Linearize Around
            zd = zeros(12, 1);
            if ~isCaptured
                zd(1:3) = uavPos;
            end 
            obj.xc.Plant = z;
    
            % Pass empty matrices for ym, r, and v
            ym = z; % Measured plant outputs
            r = zeros([12,1]);  % Output references
            r(1:3)=[1;0;1];
            v = [];  % Measured disturbance input
    
            % Call mpcmove with appropriate arguments
            mv = mpcmove(obj.mpcobj, obj.xc, ym, r, v);
            u = obj.u0*ones([4,1]) + mv;

       end
   end
            
end