
classdef mpc_controller < handle
    properties(Access = public)
        u0(1,1) double;
        mu(1,1) double;
        K;
    end

    methods(Access = public)
        function obj = mpc_controller(quadrotor)
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
            % wtf is up with these costs???
            x_gain = 2;
            v_gain = 1;
            Aa = [A, zeros(12, 12); -eye(12), zeros(12, 12)];
            Ba = [B; zeros(12, 4)];
             if rank(ctrb(Aa, Ba)) < 24
                 error('Augmented linear system for LQI is uncontrollable')
             end

            i_gain = 1;
           % Q = diag([x_gain x_gain x_gain 1/(pi/2)^2 1/(pi/2)^2 1/(pi/2)^2 ...
            %          v_gain v_gain v_gain 1/(pi/2)^2 1/(pi/2)^2 1/(pi/2)^2 ...
             %         i_gain i_gain i_gain i_gain i_gain i_gain ...
              %        i_gain i_gain i_gain i_gain i_gain i_gain]);
            %R = eye(4,4) * 1/obj.mu^2;
            sys = ss(A, B, eye(12,12), 0);
            % [K,~,~] = mpc(sys);
            mpcobj = mpc(sys);
            mpcobj.K = K;
        end
        
        mv = mpcmove(mpcobj,xc,ym,r,v);
        

        function u = output(mpcobj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            % State and Control Values to Linearize Around
            zd = zeros(12, 1);
            if ~isCaptured
                zd(1:3) = uavPos;
            end 
            u = @(z) -mpcobj.K*(z-zd);
            u = u(z);
        end
   end
            

% function cntrl = mpc(plant)
%     % Inputs:
%     %   - model: state-space model of the plant
%     %   - ts: sampling time of the plant (seconds)
%     %   - P: prediction horizon
%     %   - M: control horizon
%     %   - W: weights on prediction errors and control moves
%     %   - MV: manipulated variables (input constraints)
%     %   - OV: measured variables (output constraints)
%     %   - DV: manipulated variables (rate of change constraints)
%     if isempty(P)
%         P = 10;  % Default prediction horizon
%     end
%     if isempty(M)
%         M = 2;   % Default control horizon
%     end
%     if isempty(W)
%         W = 1;   % Default weights
%     end
%     if isempty(MV)
%         MV = struct('Min',[], 'Max',[]);  % Default input constraints
%     end
%     if isempty(OV)
%         OV = struct('Min',[], 'Max',[]);  % Default output constraints
%     end
%     if isempty(DV)
%         DV = struct('Min',[], 'Max',[]);  % Default rate of change constraints
%     end
% 
end