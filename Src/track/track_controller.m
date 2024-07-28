classdef track_controller < handle
    properties(Access = public)
        g(1,1) double;
        m;
        t;
        r;
        rd;
        rdd;
        ts;
        dt;
        Kp;
        Kv;
        Kr;
        Kw;
        I;
    end


    methods(Access = public)
        function obj = track_controller(quadrotor, q, qd, qdd, ts, dt)
            obj.g = quadrotor.g;
            obj.m = quadrotor.m;
            obj.dt = dt;
            obj.t = 1;
            obj.I = quadrotor.I;
            
            % Store reference trajectory
            obj.r = q;
            obj.rd = qd;
            obj.rdd = qdd;
            obj.ts = ts;

%             % Gains (from 2010 paper)
%             kp = 16*obj.m;
%             kv = 5.6*obj.m;
%             kr = 8.81;
%             kw = 2.54;
%       
%             obj.Kp = diag([kp kp kp]);
%             obj.Kv = diag([kv kv kv]);
%             obj.Kr = diag([kr kr kr]);
%             obj.Kw = diag([kw kw kw]);

            obj.Kp = diag([10 10 10]);
            obj.Kv = diag([100 100 100]);
            obj.Kr = diag([10 10 10]);
            obj.Kw = diag([10 10 10]);
        end

        function u = output(obj, isCaptured, z, uavPos)
            % obj, IS_CAPTURED, Z, UAV_POS

            % Euler Angle Rotation Matrices
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
            R = R_CE;
% 
%             zb = R(:, 3);
            r_des = obj.r(1:3, obj.t);
            rd_des = obj.rd(1:3, obj.t);
            rdd_des = obj.rdd(1:3, obj.t);
            zb = [r_des(3); rd_des(3); rdd_des(3) + obj.g];
            zb = zb/norm(zb);
            
            % Some diff flat stuff
            xc = [cos(obj.r(4, obj.t)), sin(obj.r(4, obj.t)), 0]';
            yb = cross(zb, xc);
            yb = yb/norm(yb,2);
            xb = cross(yb, zb);  
        
            % Linear Error
            ep = z(1:3) - r_des;
            ev = z(4:6) - rdd_des;
            
            Fdes = -obj.Kp*ep - obj.Kv*ev + obj.g*obj.m*[0;0;1] ...
               + obj.m*obj.rdd(1:3, obj.t);
            f = dot(Fdes, zb);

            % Rotational Error
            Rdes = [xb_des yb_des zb_des];
            eR = 0.5*vee(Rdes'*R-R'*Rdes);

            % Get desired angular velocity using flatness
            adot_des = obj.rdd(1:3, obj.t);
            hw = obj.m/f * (adot_des - dot(zb_des, adot_des)*zb_des);

            psidot = obj.rd(4, obj.t);
            zw = [0;0;1];
            
            pd = dot(-hw, yb_des);
            qd = dot(hw, xb_des);
            rd = dot(psidot*zw, zb_des);
            w_des = [pd; qd; rd];
            w = z(10:12);
            ew = w - R'*Rdes*w_des;

            M = -obj.Kr*eR - obj.Kw*ew; % + cross(w, obj.I*w);
            M = [f; M]; % put them together

            l = 0.2;
            sigma = 0.01; % todo save these from quad
            A = [1 1 1 1;
             0 l 0 -l;
            -l 0 l 0;
            sigma -sigma sigma -sigma];

            torques = linsolve(A, M);
            u = torques;

            obj.t = obj.t + 1;
        end
   end
            
            
end


function v = vee(R)
       v = [R(3,2) R(1,3) R(2,1)]';
end