classdef lqr_velocity_controller < handle
    properties(Access = public)
        quad; % quadrotor model
        u0; % nominal thrust for hover
        Kapproach; % LQR gains for approaching from a distance
        Kcatch; % LQR gains for close proximity catching
        Kreturn; % LQR gains for returning to base
        Ktransition; % LQR gains for in between approach and catch
        prevUAV; % previous UAV location
        isCaptured; % true if we have captured the UAV
        returnGoal; % current target location for return phase
        dt; % sim timestep size
        t; % current sim time
        seenUAV; % flag for if we have seen UAV at all
        lr; % learning rate for gradient descent
        maxVel; % cap on velocity
        prevV; % previous UAV velocity observation
        samples; % all seen locations of UAV
        zPred; % prediction of our state z in the next/current timestep
        rEsts; % estimates of disturbances (r) in time
        nEsts; % estimates of disturbances (n) in time
        tEsts; % times the disturbance samples were taken
        prevZ;% state z from the previous timestep
        prevU; % input from the previous timestep
        nlmpcObj; % nlmpc controller
        opt; % fmincon options
        nlopt; % nlmpc options
        mpcCount; % counts the number of mpc calls 
        enableMPC; % flag determining if we use MPC to ret
    end


    methods(Access = public)

        function obj = lqr_velocity_controller(quadrotor, dt, useMpc)
             if ~exist('useMpc','var')
                 % third parameter does not exist, so default it to something
                  useMpc = true;
             end
            % Initialize class variables
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.quad = quadrotor;
            obj.prevUAV = [0;0;0];
            obj.isCaptured = false;
            obj.dt = dt;
            obj.t = 0;
            obj.seenUAV = false;
            obj.lr = 1;
            obj.maxVel = 3;
            obj.prevV = [0;0;0];
            obj.samples = [];
            obj.nEsts = [];
            obj.rEsts = [];
            obj.zPred = zeros(12,1);
            obj.prevZ = zeros(12,1);
            obj.prevU = zeros(4,1);
            obj.mpcCount = 0;
            obj.enableMPC = useMpc;

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

            % Gains for approaching phase
            Q = diag([1 1 1 ...     % position gains
                      1 1 1 ...     % orientation gains
                      10 10 10 ...  % velocity gains
                      1 1 1]);      % angular velocity gains
            u_gain = 4;
            R = eye(4,4) * u_gain;
            sys = ss(A, B, eye(12), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.Kapproach = K;

            % Gains for catching phase
            Q = diag([8 8 8 ...     % position gains
                      1 1 1 ...     % orientation gains
                      1 1 1 ...     % velocity gains
                      1 1 1]);      % angular velocity gains
            u_gain = 4;
            R = eye(4,4) * u_gain;
            sys = ss(A, B, eye(12), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.Kcatch = K;

            % Gains for return phase
            Q = diag([1 1 1 ...     % position gains
                      50 50 50 ...     % orientation gains
                      5 5 5 ...     % velocity gains
                      50 50 50]);   % angular velocity gains
            u_gain = 2;
            R = eye(4,4) * u_gain;
            sys = ss(A, B, eye(12), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.Kreturn = K;

            % Initialize NL-MPC object for disturbance rejection on return
            nx = 12; % Number of states
            ny = 12; % Number of outputs (same as number of states)
            obj.nlmpcObj = nlmpc(nx, ny, 'MV',[1 2 3 4],'MD',[5 6 7 8 9 10]);

            % Set the state function
            obj.nlmpcObj.Model.StateFcn = @QuadrotorStateFcn;
            obj.nlmpcObj.Jacobian.StateFcn = @QuadrotorStateJacobianFcn;

            validateFcns(obj.nlmpcObj,zeros(nx,1),obj.u0*ones(4,1), ones(6, 1)');

            % Set NL-MPC parameters (sampling time, prediction horizon)
            obj.nlmpcObj.Ts = 0.1;
            obj.nlmpcObj.PredictionHorizon = 18;

            % Define state and control bounds
            obj.nlmpcObj.MV = struct( ...
                Min={0;0;0;0}, ...
                Max={3;3;3;3}, ...
                RateMin={-0.6;-0.6;-6;-6}, ...
                RateMax={0.6;0.6;0.6;0.6});

            % Set the cost function
            obj.nlmpcObj.Weights.OutputVariables = [1 1 1 1 1 1 1 1 1 1 1 1];
            obj.nlmpcObj.Weights.ManipulatedVariables = [0.1 0.1 0.1 0.1];
            obj.nlmpcObj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0.1];

            nloptions = nlmpcmoveopt;
            nloptions.MVTarget = obj.u0*ones([4,1])'; 
            obj.opt = optimoptions('fmincon','display','off');
            obj.nlopt = nlmpcmoveopt;
        end

        function vd = optimizeApproachVelocity(obj, vPred, e, uavPos, z)
            vd = z(7:9); % optimize from current velocity
            deltaT = obj.dt;
            n = 1000;

            % try and prevent overshoot by limiting velocity
            s = min(obj.maxVel, norm(vPred));

            % Function for error
            function en = e_next(x, y, vhat, vd, s, dt)
                er = y - x;
                eta = norm(er) / s;
                en = (y + vhat*eta) - (x + vd*dt);
            end

            % do gradient descent for fixed iterations
            for i = 1:n
                x = z(1:3);
                y = uavPos;
                en = e_next(x, y, vPred, vd, s, deltaT);
                J = -deltaT*en/norm(en) + deltaT*vd; % +vd for l2 regularization
                vd = vd - obj.lr*J;
            end
            vd = vd/norm(vd) * s;
        end

        function zp = predictNextState(obj, z, u)
            tspan = [obj.t, obj.t+obj.dt];  
            dist = struct("r", @(t,z)[0;0;0], "n", @(t,z)[0;0;0]);
            [t, z] = obj.quad.solve(tspan, z, u, dist);
            zp = z(end, :)';
        end

        function u = constraintInputs(obj, uStar, z)
            % LQR is not guarenteed to produce valid inputs,
            % this will produce the valid inputs that most closely
            % match the desired ones

            % do some constrained optimization
            lb = [0 0 0 0];
            ub = [3 3 3 3];
            Fdes = F(uStar, z, obj.quad);
            fun = @(u) norm(Fdes - F(u', z, obj.quad));
            u0 = (lb + ub)/2;

            A = [];
            b = [];
            Aeq = [];
            beq = [];

            ua = fmincon(fun, u0, A, b, Aeq, beq, lb, ub, [], obj.opt)';
            u = ua;
        end

        function u = nlmpcDisturbanceRejection(obj, z, rEst, nEst)
                % Calculate the reference trajectory
                goal = obj.returnGoal';
                horizon = obj.nlmpcObj.PredictionHorizon;
                ref = repmat([goal zeros([1,9])], ...
                    horizon, 1);
                Ts = obj.nlmpcObj.Ts;
                xrs = linspace(z(1), goal(1), abs(goal(1)-z(1))/Ts);
                yrs = linspace(z(2), goal(2), abs(goal(2)-z(2))/Ts);
                zrs = linspace(z(3), goal(3), abs(goal(3)-z(3))/Ts);
                if min(size(xrs)) > 0
                    ref(1:min(horizon, max(size(xrs))),1) = xrs(1:min(horizon, max(size(xrs))));
                end
                if min(size(yrs)) > 0
                    ref(1:min(horizon, max(size(yrs))),2) = yrs(1:min(horizon, max(size(yrs))));
                end
                if min(size(zrs)) > 0
                    ref(1:min(horizon, max(size(zrs))),3) = zrs(1:min(horizon, max(size(zrs))));
                end

                % Project the disturbance estimation
                dist_fun = obj.disturbancePredictionFunction();
                ts = linspace(obj.t, obj.t+horizon*Ts, horizon);
                md = cell2mat(arrayfun(dist_fun, ts, 'UniformOutput', false));
                md = md(:,1); % select only first to not use whole prediction

                % Solve the NL-MPC problem
                time = mod(obj.mpcCount, int32(1/obj.dt)) == 0; % only show every sim second
                if time
                    tic
                end
                mv = nlmpcmove(obj.nlmpcObj, z, obj.prevU, ref, md');
                if time
                    toc
                end
                obj.mpcCount = obj.mpcCount + 1;
                u = mv;
        end

        function dist_pred_fun = disturbancePredictionFunction(obj)
            % remove the first sample--it is always zero
            rEsts = obj.rEsts(:, 2:end);
            nEsts = obj.nEsts(:, 2:end);
            tEsts = obj.tEsts(:, 2:end);
            samples_size = size(rEsts);
            samples_size = samples_size(2); % how many samples we got?
            if samples_size == 1
                dist_pred_fun = @(t) [rEsts(:, 1); nEsts(:, 1)];
            else
                n = min(samples_size-1, 3);
                rFun = polyPred(rEsts(:, end-n:end), tEsts(:, end-n:end), n);
                nFun = polyPred(nEsts(:, end-n:end), tEsts(:, end-n:end), n);
                dist_pred_fun = @(t) [rFun(t); nFun(t)];
            end
        end

        function Ktransition = calculateKtransition(obj, distance)
            %transitional K matrix for when the distance is less than 3
            %currently not used, and makes it run very slow
            zs = zeros(12, 1);
            zs(3) = 1;
            us = obj.u0 * ones([4,1]); % enough to counteract gravity
            [A, B] = linearize(zs, us);
            % Transition gains
            Q = diag([8-(distance*7/3) 8-(distance*7/3) 8-(distance*7/3) ...     % position gains
                      1 1 1 ...     % orientation gains
                      1+(distance*3) 1+(distance*3) 1+(distance*3) ...     % velocity gains
                      1 1 1]);      % angular velocity gains
            u_gain = 4;
            R = eye(4,4) * u_gain;
            sys = ss(A, B, eye(12), 0);
            [K,~,~] = lqr(sys,Q,R);
            obj.Ktransition = K;
            Ktransition = K;
        end

        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS
            if abs(uavPos(3)) > 10 || uavPos(3) < 0 || abs(uavPos(2)) > 5 || abs(uavPos(1)) > 5
                uavPos = [0;0;0.1];
                obj.prevUAV = uavPos;
            end
            
            % Error for optimization
            e = uavPos-z(1:3);
            
            if isCaptured
                % only look for disturbance when captured
                % Error in prediction (disturbance)
                [nEst, rEst] = estimateDisturbance(z, obj.zPred, obj.quad, obj.dt);
    
                % store disturbance estimates
                obj.rEsts = [obj.rEsts rEst];
                obj.nEsts = [obj.nEsts nEst];
                obj.tEsts = [obj.tEsts obj.t];
            end

            % calculate uav velocity from previous timestep
            if obj.seenUAV && ~obj.isCaptured
                % need to have one sample to calculate UAV velocity
                v = uavPos-obj.prevUAV;
            else
                v = [0;0;0];
            end

            % estimate the velocity of the UAV at the next timestep
            v = v/obj.dt;
            a = (v - obj.prevV) / obj.dt; % estimate acceleration
            vPred = v+a*obj.dt;

            % estimate (poorly) angular velocity.
            % not precise, but can tell straight from not :)
            if norm(v) > 0 && norm(vPred) > 0
                cosTheta = dot(v, vPred) / norm(v) / norm(vPred);
                theta = acos(cosTheta);
                wPred = theta/obj.dt;
            else
                wPred = 0;
            end

            % ought to be at least straight enough
            uavStraight = wPred < 0.0001;
            if uavStraight
                vd = vPred; % just match with a straight enough target
            else
                % need to be more creative with curving targets
                vd = obj.optimizeApproachVelocity(vPred, e, uavPos, z);
            end

            % set flag for furst capture to set return hover point
            if isCaptured && ~obj.isCaptured
                obj.isCaptured = true;
                % go to current altitude, (x,y)=(0,0)
%                 obj.returnGoal = [0;0;z(3)]; 
                obj.returnGoal = [0;0;0.5]; 
                fprintf(['UAV Captured at time %i!!!\n' ...
                    'The runtime of every 1000th nlmpc call will display\n' ...
                    '(once per sim second, only if nlmpc is called)\n'], int32(obj.t*100)/100);
                
                % for the nlmpc solvers sake
                obj.prevU = obj.u0*ones(4,1);
            end

            % consider cases for applying gains, goal state
            zd = zeros(12, 1);
            rad = 3;
            function newVelocity = calculateNewVelocity(velocity, distance)
                if distance > 3
                    newVelocity = velocity*1.2;
                else
                    newVelocity = velocity*(1 + distance*0.2/3);
                end
            end

            if obj.isCaptured % returning to goal   
                if norm(z(1:3)-obj.returnGoal) < 1
                    obj.returnGoal = [0;0;1];
                end
                K = obj.Kreturn;
                if obj.enableMPC && (norm(rEst) > 1e-10 || norm(nEst) > 1e-10)
                    u = obj.nlmpcDisturbanceRejection(z, rEst, nEst);
                else
                    % no disturbance? just use LQR bro
                    zd(1:3) = obj.returnGoal;
                    u = @(z) obj.u0*ones([4,1])+K*(zd-z);
                    u = u(z);
                    u = obj.constraintInputs(u, z);
                end
            
            elseif norm(uavPos-z(1:3)) > rad % far from target
                K = obj.Kapproach;
                yn1 = uavPos+vPred*obj.dt;
                zd(1:3) = yn1;
                zd(7:9) = vd;
            else % close to target, go for the kill!
%                 disp(norm(uavPos-z(1:3)));
                % K = obj.calculateKtransition(norm(uavPos-z(1:3)));
                    %this was for changing K matrix
                K = obj.Kcatch;
                yn1 = uavPos+vPred*obj.dt;
                zd(1:3) = yn1;
                zd(7:9) = calculateNewVelocity(vd, norm(uavPos-z(1:3)));
            end 

            % actually calculate inputs via LQR gains
            if ~obj.isCaptured
                u = @(z) obj.u0*ones([4,1])+K*(zd-z);
                u = u(z);
                % uncomment this for real, but it slows down testing
                u = obj.constraintInputs(u, z);
            end

            % now that we have you, we can predict where we'll be next
            % (for disturbance rejection)
            obj.zPred = obj.predictNextState(z,u);

            % store data for use in next timestep
            obj.prevUAV = uavPos;
            obj.prevV = v;
            obj.prevZ = z;
            obj.prevU = u;
            obj.seenUAV = true;
            obj.t = obj.t + obj.dt;
            obj.samples = [uavPos obj.samples];
        end
   end
            
            
end

function R_CE = rotation(z)
    a1 = z(4);
    a2 = z(5);
    a3 = z(6);
    R3 = [cos(a3) -sin(a3) 0;
          sin(a3) cos(a3)  0;
          0 0 1];
    R2 = [cos(a2) 0 sin(a2);
          0 1 0;
          -sin(a2) 0 cos(a2)];
    R1 = [1 0 0;
          0 cos(a1) -sin(a1);
          0 sin(a1) cos(a1)];
    R_CE = R3*R2*R1;
end

function f = F(u, z, quad)
    R = rotation(z);

    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
    u4 = u(4);

    fv = R*[0;0;u1+u2+u3+u4];
    fw = [quad.l*(u2-u4); 
          quad.l*(u3-u1); 
          quad.sigma*(u1-u2+u3-u4)];

    f = [fv; fw];
end

function J = dFdu(z, quad)
    R = rotation(z);
    r3 = R(:, 3)';
    l = quad.l;
    sigma = quad.sigma;
    J = [r3 0 -l sigma;
         r3 l 0 -sigma;
         r3 0 l sigma;
         r3 -l 0 -sigma];
end

function u = rejectNwithThrust(fnDes, fvDes, nEst, quad)
    l = quad.l;
    sigma = quad.sigma;
    A = [0 l 0 -l;
        -l 0 l 0;
         sigma -sigma sigma -sigma;
         1 1 1 1];
    b = [fnDes-nEst; fvDes];

    u = linsolve(A,b);
end