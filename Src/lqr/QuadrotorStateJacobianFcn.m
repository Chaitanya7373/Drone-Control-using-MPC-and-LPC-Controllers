function [A,B] = QuadrotorStateJacobianFcn(in1,in2)
%QuadrotorStateJacobianFcn
%    [A,B] = QuadrotorStateJacobianFcn(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    25-Apr-2024 12:37:08

n1 = in2(8,:);
n2 = in2(9,:);
n3 = in2(10,:);
phi_t = in1(4,:);
phi_dot_t = in1(10,:);
psi_t = in1(6,:);
psi_dot_t = in1(12,:);
r1 = in2(5,:);
r2 = in2(6,:);
r3 = in2(7,:);
theta_t = in1(5,:);
theta_dot_t = in1(11,:);
u1 = in2(1,:);
u2 = in2(2,:);
u3 = in2(3,:);
u4 = in2(4,:);
t2 = cos(phi_t);
t3 = cos(psi_t);
t4 = cos(theta_t);
t5 = sin(phi_t);
t6 = sin(psi_t);
t7 = sin(theta_t);
t8 = phi_t.*2.0;
t9 = psi_dot_t.^2;
t10 = theta_t.*2.0;
t11 = theta_dot_t.^2;
t25 = u1+u2+u3+u4;
t12 = cos(t8);
t13 = t2.^2;
t14 = t4.^2;
t15 = t4.^3;
t16 = sin(t8);
t17 = t5.^2;
t18 = sin(t10);
t19 = t7.*2.0;
t20 = t7.^2;
t21 = t2.*t3;
t22 = t2.*t6;
t23 = t3.*t5;
t24 = t5.*t6;
t26 = t4.*4.0e+1;
t27 = 1.0./t4;
t31 = t7.*4.0e+1;
t32 = t2.*t4.*2.0;
t49 = t2.*t4.*t5.*2.0e+1;
t56 = n1.*t2.*t5.*t7.*2.0e+2;
t57 = t2.*t5.*t7.*u4.*-4.0e+1;
t61 = phi_dot_t.*t2.*t4.*t5.*theta_dot_t.*4.96e+2;
t66 = psi_dot_t.*t2.*t4.*t5.*t7.*theta_dot_t.*4.96e+2;
t28 = 1.0./t14;
t29 = 1.0./t15;
t30 = -t19;
t33 = t23.*2.0;
t34 = t24.*2.0;
t35 = t7.*t24;
t37 = t7.*t13;
t38 = -t13;
t39 = t13.*2.0e+1;
t40 = t7.*t21;
t41 = t7.*t22;
t42 = t7.*t23;
t44 = t19.*t21;
t45 = t19.*t22;
t46 = t13.*t26;
t50 = t2.*t5.*t31;
t59 = t7.*t49;
t60 = t4.*t13.*theta_dot_t.*2.48e+2;
t65 = -t61;
t36 = -t33;
t43 = -t39;
t47 = -t37;
t48 = t37.*2.0e+1;
t51 = -t41;
t52 = -t42;
t53 = t50.*u2;
t54 = t50.*u4;
t58 = t14.*t39;
t62 = t21+t35;
t63 = t24+t40;
t64 = t4.*t37.*theta_dot_t.*2.48e+2;
t69 = t34+t44;
t55 = -t48;
t67 = t22+t52;
t68 = t23+t51;
et1 = t28.*(n1.*t2.*t5.*2.0e+2+n2.*t4.*t37.*1.0e+2-t9.*t13.*t14.*2.48e+2+t11.*t13.*t14.*2.48e+2+t9.*t14.*t17.*2.48e+2-t11.*t14.*t17.*2.48e+2+t2.*t5.*u2.*4.0e+1-t2.*t5.*u4.*4.0e+1-t4.*t37.*u1.*2.0e+1+t4.*t48.*u3+phi_dot_t.*psi_dot_t.*t14.*t37.*2.48e+2+n3.*t2.*t5.*t7.*2.0e+2-n1.*t2.*t5.*t14.*2.0e+2-n2.*t4.*t7.*t17.*1.0e+2-t2.*t5.*t7.*u2.*2.0-t2.*t5.*t7.*u4.*2.0-t2.*t5.*t14.*u2.*4.0e+1+t2.*t5.*t14.*u4.*4.0e+1+t2.*t5.*t19.*u1+t2.*t5.*t19.*u3+t4.*t7.*t17.*u1.*2.0e+1-t4.*t7.*t17.*u3.*2.0e+1-phi_dot_t.*psi_dot_t.*t7.*t14.*t17.*2.48e+2+psi_dot_t.*t2.*t4.*t5.*theta_dot_t.*4.96e+2+psi_dot_t.*t2.*t5.*t15.*theta_dot_t.*4.96e+2-phi_dot_t.*t2.*t4.*t5.*t7.*theta_dot_t.*4.96e+2);
et2 = 1.0./2.48e+2;
et3 = t28.*(n3.*t4.*2.0e+2+t4.*u1.*2.0-t4.*u2.*2.0+t4.*u3.*2.0-t4.*u4.*2.0-n3.*t4.*t13.*1.0e+2-n1.*t4.*t37.*2.0e+2-psi_dot_t.*t7.*theta_dot_t.*4.96e+2+psi_dot_t.*t37.*theta_dot_t.*2.48e+2+t4.*t13.*u2+t4.*t13.*u4-t4.*t37.*u2.*4.0e+1+t4.*t38.*u1+t4.*t38.*u3+t26.*t37.*u4+n2.*t2.*t5.*t14.*1.0e+2-n2.*t2.*t5.*t20.*1.0e+2+phi_dot_t.*t13.*t14.*theta_dot_t.*2.48e+2-phi_dot_t.*t13.*t20.*theta_dot_t.*2.48e+2+psi_dot_t.*t14.*t37.*theta_dot_t.*7.44e+2-t2.*t5.*t14.*u1.*2.0e+1+t2.*t5.*t14.*u3.*2.0e+1+t2.*t5.*t20.*u1.*2.0e+1-t2.*t5.*t20.*u3.*2.0e+1+phi_dot_t.*psi_dot_t.*t2.*t5.*t15.*2.48e+2+t2.*t4.*t5.*t7.*t9.*4.96e+2-t2.*t4.*t5.*t7.*t11.*4.96e+2-phi_dot_t.*psi_dot_t.*t2.*t4.*t5.*t20.*4.96e+2);
et4 = 1.0./2.48e+2;
et5 = t27.*(n2.*t7.*2.0e+2+n2.*t37.*2.0e+2-t7.*u1.*4.0e+1+t31.*u3-t37.*u1.*4.0e+1+t9.*t13.*t15.*4.96e+2+t13.*t31.*u3-phi_dot_t.*psi_dot_t.*t4.*t7.*9.92e+2+phi_dot_t.*psi_dot_t.*t4.*t37.*9.92e+2-n1.*t2.*t4.*t5.*2.0e+2-t4.*t9.*t13.*t20.*9.92e+2-t2.*t4.*t5.*u2.*4.0e+1+t2.*t5.*t26.*u4-phi_dot_t.*t2.*t5.*t7.*theta_dot_t.*4.96e+2-psi_dot_t.*t2.*t5.*t14.*theta_dot_t.*4.96e+2+psi_dot_t.*t2.*t5.*t20.*theta_dot_t.*4.96e+2).*(-1.0./4.96e+2);
et6 = (t7.*t28.*(t53+t56+t57+t65+t66+n2.*t4.*2.0e+2+n3.*t16.*1.0e+2-t4.*u1.*4.0e+1+t16.*u1-t16.*u2+t16.*u3-t16.*u4+t26.*u3+t46.*u3-phi_dot_t.*psi_dot_t.*t14.*4.96e+2+n2.*t4.*t13.*2.0e+2-t9.*t14.*t37.*4.96e+2-t4.*t13.*u1.*4.0e+1+phi_dot_t.*psi_dot_t.*t13.*t14.*4.96e+2))./4.96e+2;
et7 = t28.*(n1.*t4.*-2.0e+2-t4.*u2.*4.0e+1+t26.*u4-psi_dot_t.*theta_dot_t.*cos(t10).*4.96e+2+n1.*t4.*t13.*1.0e+2+phi_dot_t.*t37.*theta_dot_t.*2.48e+2-t4.*t13.*u4.*2.0e+1+t4.*t39.*u2+n2.*t2.*t5.*t7.*1.0e+2+psi_dot_t.*t13.*t14.*theta_dot_t.*2.48e+2-psi_dot_t.*t13.*t20.*theta_dot_t.*2.48e+2+t2.*t5.*t9.*t15.*2.48e+2-t2.*t5.*t7.*u1.*2.0e+1+t2.*t5.*t7.*u3.*2.0e+1-t2.*t4.*t5.*t9.*t20.*4.96e+2+phi_dot_t.*psi_dot_t.*t2.*t4.*t5.*t7.*4.96e+2).*(-1.0./2.48e+2);
et8 = (t7.*t29.*(n3.*2.0e+2+u1.*2.0-u2.*2.0+u3.*2.0-u4.*2.0+n1.*t7.*2.0e+2-n3.*t13.*1.0e+2-n1.*t37.*1.0e+2+phi_dot_t.*t60-t7.*u4.*4.0e+1+t13.*u2+t13.*u4+t31.*u2-t37.*u2.*2.0e+1+t38.*u1+t38.*u3+t48.*u4+t49.*u3+psi_dot_t.*t18.*theta_dot_t.*2.48e+2+n2.*t2.*t4.*t5.*1.0e+2-psi_dot_t.*t4.*t37.*theta_dot_t.*2.48e+2-t2.*t4.*t5.*u1.*2.0e+1+phi_dot_t.*psi_dot_t.*t2.*t5.*t14.*2.48e+2-t2.*t5.*t7.*t9.*t14.*2.48e+2))./1.24e+2;
mt1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,r2.*t63.*2.0+r3.*t67.*2.0+t25.*t67.*2.0,r3.*t62.*-2.0-r2.*t68.*2.0-t25.*t62.*2.0,r2.*t32-r3.*t4.*t5.*2.0-t4.*t5.*t25.*2.0,et1.*et2];
mt2 = [(t27.*(n3.*t12.*2.0e+2+n1.*t37.*2.0e+2+t12.*u1.*2.0-t12.*u2.*2.0+t12.*u3.*2.0-t12.*u4.*2.0-t37.*u4.*4.0e+1-n1.*t7.*t17.*2.0e+2-t7.*t17.*u2.*4.0e+1+t13.*t31.*u2+t17.*t31.*u4-n2.*t2.*t4.*t5.*4.0e+2-phi_dot_t.*t4.*t13.*theta_dot_t.*4.96e+2+phi_dot_t.*t4.*t17.*theta_dot_t.*4.96e+2+psi_dot_t.*t4.*t37.*theta_dot_t.*4.96e+2+t2.*t4.*t5.*u1.*8.0e+1-t2.*t4.*t5.*u3.*8.0e+1-phi_dot_t.*psi_dot_t.*t2.*t5.*t14.*9.92e+2-psi_dot_t.*t4.*t7.*t17.*theta_dot_t.*4.96e+2+t2.*t5.*t7.*t9.*t14.*9.92e+2))./4.96e+2];
mt3 = [(t28.*(t53+t56+t57+t65+t66+n3.*t2.*t5.*2.0e+2+n2.*t4.*t13.*1.0e+2-n2.*t4.*t17.*1.0e+2-t9.*t14.*t37.*2.48e+2+t2.*t5.*u1.*2.0-t2.*t5.*u2.*2.0+t2.*t5.*u3.*2.0-t2.*t5.*u4.*2.0-t4.*t13.*u1.*2.0e+1+t4.*t17.*u1.*2.0e+1-t4.*t17.*u3.*2.0e+1+t4.*t39.*u3+phi_dot_t.*psi_dot_t.*t13.*t14.*2.48e+2-phi_dot_t.*psi_dot_t.*t14.*t17.*2.48e+2+t7.*t9.*t14.*t17.*2.48e+2))./2.48e+2,0.0,0.0,0.0,0.0,0.0,0.0,r1.*t3.*t7.*-2.0+r3.*t4.*t21.*2.0+r2.*t4.*t33+t4.*t21.*t25.*2.0,r1.*t6.*t7.*-2.0+r3.*t4.*t22.*2.0+r2.*t4.*t34+t4.*t22.*t25.*2.0,r1.*t4.*-2.0-r3.*t2.*t7.*2.0-r2.*t5.*t7.*2.0-t2.*t7.*t25.*2.0];
mt4 = [et3.*et4+(t7.*t29.*(n1.*2.0e+2+u2.*4.0e+1-u4.*4.0e+1+n3.*t7.*2.0e+2-n1.*t13.*1.0e+2-n3.*t37.*1.0e+2+phi_dot_t.*t64-t7.*u2.*2.0-t7.*u4.*2.0-t13.*u2.*2.0e+1+t19.*u1+t19.*u3+t37.*u2+t37.*u4+t39.*u4+t47.*u1+t47.*u3+t58.*u2+t59.*u3+n1.*t13.*t14.*1.0e+2+psi_dot_t.*t4.*theta_dot_t.*4.96e+2-t13.*t14.*u4.*2.0e+1-psi_dot_t.*t4.*t13.*theta_dot_t.*2.48e+2-psi_dot_t.*t13.*t15.*theta_dot_t.*2.48e+2-t2.*t5.*t9.*t14.*2.48e+2+t2.*t5.*t11.*t14.*2.48e+2+n2.*t2.*t4.*t5.*t7.*1.0e+2-t2.*t4.*t5.*t7.*u1.*2.0e+1+phi_dot_t.*psi_dot_t.*t2.*t5.*t7.*t14.*2.48e+2))./1.24e+2,et5+et6,et7+et8,0.0,0.0,0.0,0.0,0.0,0.0];
mt5 = [r2.*t62.*-2.0+r3.*t68.*2.0+t25.*t68.*2.0-r1.*t4.*t6.*2.0,r3.*t63.*2.0-r2.*t67.*2.0+t25.*t63.*2.0+r1.*t3.*t4.*2.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,(t28.*(t64+psi_dot_t.*t2.*t5.*t7.*t14.*2.48e+2))./2.48e+2,t27.*(psi_dot_t.*t14.*4.96e+2-psi_dot_t.*t13.*t14.*4.96e+2+t2.*t4.*t5.*theta_dot_t.*4.96e+2).*(-1.0./4.96e+2),(t28.*(t60+psi_dot_t.*t2.*t5.*t14.*2.48e+2))./2.48e+2,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0];
mt6 = [(t28.*(psi_dot_t.*t4.*4.96e+2+phi_dot_t.*t4.*t37.*2.48e+2-psi_dot_t.*t4.*t13.*2.48e+2-psi_dot_t.*t13.*t15.*2.48e+2+t2.*t5.*t14.*theta_dot_t.*4.96e+2))./2.48e+2,t27.*(phi_dot_t.*t2.*t4.*t5.*4.96e+2-psi_dot_t.*t2.*t4.*t5.*t7.*4.96e+2).*(-1.0./4.96e+2),(t28.*(psi_dot_t.*t18.*2.48e+2+phi_dot_t.*t4.*t13.*2.48e+2-psi_dot_t.*t4.*t37.*2.48e+2))./2.48e+2,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t28.*(t60-t4.*theta_dot_t.*4.96e+2+t13.*t15.*theta_dot_t.*2.48e+2+psi_dot_t.*t2.*t5.*t14.*4.96e+2-phi_dot_t.*t2.*t5.*t7.*t14.*2.48e+2).*(-1.0./2.48e+2),t27.*(phi_dot_t.*t14.*4.96e+2-phi_dot_t.*t13.*t14.*4.96e+2+psi_dot_t.*t14.*t37.*9.92e+2-t2.*t4.*t5.*t7.*theta_dot_t.*4.96e+2).*(-1.0./4.96e+2)];
mt7 = [t28.*(t64-t18.*theta_dot_t.*2.48e+2-phi_dot_t.*t2.*t5.*t14.*2.48e+2+psi_dot_t.*t2.*t5.*t7.*t14.*4.96e+2).*(-1.0./2.48e+2)];
A = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],12,12);
if nargout > 1
    t70 = t36+t45;
    B = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t69,t70,t32,t28.*(t30+t37+t59).*(-1.0./2.48e+2),t27.*(-t16+t26+t46).*(-1.0./4.96e+2),t28.*(t13+t49-2.0).*(-1.0./2.48e+2),0.0,0.0,0.0,0.0,0.0,0.0,t69,t70,t32,(t28.*(t30+t37+t43+t58+4.0e+1))./2.48e+2,t27.*(t16-t2.*t5.*t7.*4.0e+1).*(-1.0./4.96e+2),(t28.*(t13+t31+t55-2.0))./2.48e+2,0.0,0.0,0.0,0.0,0.0,0.0,t69,t70,t32,(t28.*(t19+t47+t59))./2.48e+2,(t27.*(t16+t26+t46))./4.96e+2,(t28.*(t38+t49+2.0))./2.48e+2,0.0,0.0,0.0,0.0,0.0,0.0,t69,t70,t32,t28.*(t19+t43+t47+t58+4.0e+1).*(-1.0./2.48e+2),t27.*(t16+t50).*(-1.0./4.96e+2),(t28.*(t13-t31+t48-2.0))./2.48e+2],[12,4]);
end
