clear all;
close all;

%time-related variables
dT=0.001;
tinit=0;
tend=20;
time=[tinit:dT:tend]';

sigma_0 = 1e2;
sigma_1  = sqrt(1e2);
sigma_2  = 0.4;
Fc = 1;
Fs = 1.5;
vs = 0.001;
z_1 = 0;
z_2 = 0;

%define robot parameters
robot_param.m1 = 5;
robot_param.m2 = 3;
robot_param.l1 = 0.5;
robot_param.l2 = 0.3;
robot_param.g_acc = 9.81;

DOF=2;

% dynamics parameters
M=zeros(DOF,DOF);
h=zeros(DOF,1);
g=zeros(DOF,1);


%initialize
q=zeros(DOF,1); %joint angle
q = [0.7854; -pi/2]; %initialize
q_dot=zeros(DOF,1); %joint velocity
tau_c=zeros(DOF,1); %control input

p = zeros(DOF,1); % task variable (=end effector position)
p_dot = zeros(DOF,1); % task velocity
J = zeros(DOF,DOF); % Jacobian
[p,J] = forwardKinPlanar2DOF(robot_param, q); %initialize
p_dot = J*q_dot; % initialize

x=[q;q_dot]; % x = [q1, q2, q1_dot, q2_dot]'
xdot=zeros(2*DOF,1);


%stacking variables
q_stack=[];
q_dot_stack=[];
p_stack=[];
p_dot_stack=[];
tau_c_stack=[];
tau_f_stack=[];
tau_f_est_stack=[];

w = 30*2*pi; 
Q = tf([w],[1 w]);
pre_LPF_1 = 0;
pre_LPF_2 = 0;
tau_f_est = 0;

for i=1:tend/dT+1 
    % update measurements
    q(1)=x(1);
    q(2)=x(2);
    q_dot(1)=x(3);
    q_dot(2)=x(4);

    % update dynamics
    [M, h, g]=getRobotDyn_planar2DOF(robot_param, q,q_dot);
    
    %update kinematics
    [p,J] = forwardKinPlanar2DOF(robot_param, q);
    p_dot = J*q_dot;
    
    
    %apply only damping torque if you want to see if the  robot falls down
    %correctly under the gravity
    %tau_c = -K_d*q_dot;
    
    %joint pd+gravity compensation controller (set point)
    q_des=[0.2;0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= -K_p*(q-q_des) - K_d*q_dot + g;

    
%     %task pd+gravity compensation controller (set point)
%     p_des=[0.5; 0.0];
%     K_p=50*eye(DOF);
%     K_d=5*eye(DOF);
%     tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot ) ;
    
    [tau_f_1, z_1] = lugre(z_1, q_dot(1), Fc, Fs, vs, sigma_0, sigma_1, sigma_2, dT);
    [tau_f_2, z_2] = lugre(z_2, q_dot(2), Fc, Fs, vs, sigma_0, sigma_1, sigma_2, dT);
    tau_f = [tau_f_1 tau_f_2];
    
    % update states 
%     tau = tau_c;
    tau = tau_c + tau_f;
    tau_cf = tau_c - tau_f_est;
%     tau = tau_cf + tau_f;

    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    q_stack=[q_stack; q(1) q(2)];
    q_dot_stack=[q_dot_stack; q_dot(1) q_dot(2)];   
    p_stack=[p_stack; p(1) p(2)];
    p_dot_stack=[p_dot_stack; p_dot(1) p_dot(2)];
    
    J_inv = inv(J);
    LPF_1 = LPF(tau, pre_LPF_1, 0.001, 0.0001);
    LPF_2 = LPF(tau_cf, pre_LPF_2, 0.001, 0.0001);
    tau_f_est = LPF_1 - LPF_2;
    pre_LPF_1 = LPF_1;
    pre_LPF_1 = LPF_2;
    
    tau_f_stack = [tau_f_stack; tau_f];
    tau_f_est_stack = [tau_f_est_stack; tau_f_est(1) tau_f_est(2)];
end

%% plot

linewidth = 1.5;
fontsize = 12;

subplot(2, 2, 1); 
% figure;
plot(time, q_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, q_stack(:,2), 'linewidth',linewidth);
legend('q_1','q_2','Location','Best')
title('Angle');

subplot(2, 2, 2); 
plot(time, p_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, p_stack(:,2), 'linewidth',linewidth);
legend('p_x','p_y','Location','Best')
title('Position');

subplot(2, 2, 3); 
plot(time, tau_f_est_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, tau_f_est_stack(:,2), 'linewidth',linewidth);
legend('estimated tau 1','estimated tau 2','Location','Best')
title('Estimated tau_f');

subplot(2, 2, 4); 
plot(time, tau_f_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, tau_f_stack(:,2), 'linewidth',linewidth);
legend('tau_f1','tau_f2','Location','Best')
title('tau_f');