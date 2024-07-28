classdef nlmpc_controller < handle
    properties
        quadrotor % Instance of quadrotor
        nlmpcObj % NL-MPC object
        K_prev % Previous control gains
        Q % State cost matrix
        R % Control cost matrix
        Ts % Sampling time
        Np % Prediction horizon
        lastMV
        opt
        u0
        u
        prevUAV;
        isCaptured;
        returnGoal;
    end

    methods
        function obj = nlmpc_controller(quadrotor)
            obj.quadrotor = quadrotor;
            obj.Ts = 0.1; % Sampling time (adjust as needed)
            obj.Np = 18; % Prediction horizon (adjust as needed)
            obj.K_prev = zeros(4, 12); % Initial control gains
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.prevUAV = [0;0;0];
            obj.isCaptured = false;

            % Initialize NL-MPC object
            nx = 12; % Number of states
            ny = 12; % Number of outputs (same as number of states)
            obj.nlmpcObj = nlmpc(nx, ny, 'MV',[1 2 3 4],'MD',[5 6 7 8 9 10]);

            % Set the state function
            obj.nlmpcObj.Model.StateFcn = @QuadrotorStateFcn;
            obj.nlmpcObj.Jacobian.StateFcn = @QuadrotorStateJacobianFcn;

            validateFcns(obj.nlmpcObj,zeros(nx,1),obj.u0*ones(4,1), ones(6, 1)');

            % Set NL-MPC parameters (sampling time, prediction horizon)
            obj.nlmpcObj.Ts = obj.Ts;
            obj.nlmpcObj.PredictionHorizon = obj.Np;

            % Define state and control bounds
            obj.nlmpcObj.MV = struct( ...
                Min={0;0;0;0}, ...
                Max={3;3;3;3}, ...
                RateMin={-0.6;-0.6;-6;-6}, ...
                RateMax={0.6;0.6;0.6;0.6});

            % Set the cost function
            obj.nlmpcObj.Weights.OutputVariables = [1 1 1 1 1 1 0 0 0 0 0 0];
            obj.nlmpcObj.Weights.ManipulatedVariables = [0.1 0.1 0.1 0.1];
            obj.nlmpcObj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0.1];

            nloptions = nlmpcmoveopt;
            nloptions.MVTarget = obj.u0*ones([4,1])'; 
            obj.opt = nloptions;
            obj.lastMV = nloptions.MVTarget;
            obj.u = obj.lastMV';
        end

        function u = output(obj, isCaptured, z, uavPos)

            if isCaptured && ~obj.isCaptured
                obj.isCaptured = true;
                obj.returnGoal = [0;0;z(3)];
                fprintf('UAV Captured!!!\n')
            end

            if ~obj.isCaptured
%                 ref = zeros([12, obj.Np]);
%                 ref(1, :) = linspace(z(1), 0, obj.Np);
%                 ref(2, :) = linspace(z(2), 0, obj.Np);
%                 ref(3, :) = linspace(z(3), 1, obj.Np);
%                 ref = [1;1;1;0;0;0;0;0;0;0;0;0];

                % disp(v);
                newUAVPos = uavPos;
                ref = [newUAVPos' zeros([1,9])];
                % Solve the NL-MPC problem
%                 tic
                % TODO figure out why this is so slow with MD
                mv = nlmpcmove(obj.nlmpcObj, z, obj.lastMV, ref, zeros(1,6), obj.opt);
%                 toc
                obj.lastMV = mv;
                u = mv
%                 fprintf('sus')

            else
                % disp(v);
                newUAVPos = [0;0;1];
                ref = [newUAVPos' zeros([1,9])];
                % Solve the NL-MPC problem

                % TODO actually estimate these
                rEst = [0.1 0.1 0.1]';
                nEst = [0; 0; 0];
                md = [rEst; nEst];
                tic
                % TODO figure out why this is so slow with MD
                mv = nlmpcmove(obj.nlmpcObj, z, obj.lastMV, ref, md', obj.opt);
                toc
                obj.lastMV = mv;
                u = mv
%                 fprintf('amogus')
            end
            obj.u = [obj.u, u];
            obj.prevUAV = uavPos;
        end
    end
end
