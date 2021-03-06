clear all;
close all;

%time-related variables
dT=0.001;
tinit=0;
tend=60;
time=[tinit:dT:tend]';

%define robot parameters
robot_param.m1 = 3;
robot_param.m2 = 3;
robot_param.l1 = 1.5;
robot_param.l2 = 1.5;
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

tau_stack_1 = [];
tau_stack_2 = [];
tau_stack_3 = [];
tau_stack_4 = [];
tau_stack_5 = [];

q_stack_2=[];
q_stack_3=[];
q_stack_4=[];
q_stack_5=[];

p_stack_2=[];
p_stack_3=[];
p_stack_4=[];
p_stack_5=[];

%% D = 0
D = 0;
tic;
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
    J_x = J(1, :);
    J_x_pinv = pinv(J_x);
    z = transpose(null(J_x));
    z_pinv = pinv(z);
    p_dot = J*q_dot;
    
    %task pd+gravity compensation controller (set point)
    p_des=[1.5; 0.0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot);
    
    p_des=1.5;
    K_p=50;
    K_d=5;
    W = M;
    tau_n = -D*(q_dot);
    J_w_pinv = inv(W)*J_x'*inv(J_x*inv(W)*J_x');
    f = - K_d*p_dot(1) + K_p*(p_des - p(1));
    projector = eye(2) - J_x' * J_w_pinv';
    tau_c = J_x'*f + projector*tau_n + g;
    
    % update states 
    tau = tau_c;
    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    tau_stack_1 = [tau_stack_1; tau(1) tau(2)];
    q_stack=[q_stack; q(1) q(2)];
    q_dot_stack=[q_dot_stack; q_dot(1) q_dot(2)];   
    p_stack=[p_stack; p(1) p(2)];
    p_dot_stack=[p_dot_stack; p_dot(1) p_dot(2)];   
end
toc
%% D = 0.5
D = 0.5;
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

tic;
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
    J_x = J(1, :);
    J_x_pinv = pinv(J_x);
    z = transpose(null(J_x));
    z_pinv = pinv(z);
    p_dot = J*q_dot;
    
    %task pd+gravity compensation controller (set point)
    p_des=[1.5; 0.0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot);
    
    p_des=1.5;
    K_p=50;
    K_d=5;
    W = M;
    tau_n = -D*(q_dot);
    J_w_pinv = inv(W)*J_x'*inv(J_x*inv(W)*J_x');
    f = - K_d*p_dot(1) + K_p*(p_des - p(1));
    projector = eye(2) - J_x' * J_w_pinv';
    tau_c = J_x'*f + projector*tau_n + g;
    
    % update states 
    tau = tau_c;
    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    tau_stack_2 = [tau_stack_2; tau(1) tau(2)];
    q_stack_2=[q_stack_2; q(1) q(2)];
    p_stack_2=[p_stack_2; p(1) p(2)];
end
toc
%% D = 1.0
D = 1.0;
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

tic;
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
    J_x = J(1, :);
    J_x_pinv = pinv(J_x);
    z = transpose(null(J_x));
    z_pinv = pinv(z);
    p_dot = J*q_dot;
    
    %task pd+gravity compensation controller (set point)
    p_des=[1.5; 0.0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot);
    
    p_des=1.5;
    K_p=50;
    K_d=5;
    W = M;
    tau_n = -D*(q_dot);
    J_w_pinv = inv(W)*J_x'*inv(J_x*inv(W)*J_x');
    f = - K_d*p_dot(1) + K_p*(p_des - p(1));
    projector = eye(2) - J_x' * J_w_pinv';
    tau_c = J_x'*f + projector*tau_n + g;
    
    % update states 
    tau = tau_c;
    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    tau_stack_3 = [tau_stack_3; tau(1) tau(2)];
    q_stack_3=[q_stack_3; q(1) q(2)];
    p_stack_3=[p_stack_3; p(1) p(2)];
end
toc
%% D = 1.5
D = 1.5;
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

tic;
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
    J_x = J(1, :);
    J_x_pinv = pinv(J_x);
    z = transpose(null(J_x));
    z_pinv = pinv(z);
    p_dot = J*q_dot;
    
    %task pd+gravity compensation controller (set point)
    p_des=[1.5; 0.0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot);
    
    p_des=1.5;
    K_p=50;
    K_d=5;
    W = M;
    tau_n = -D*(q_dot);
    J_w_pinv = inv(W)*J_x'*inv(J_x*inv(W)*J_x');
    f = - K_d*p_dot(1) + K_p*(p_des - p(1));
    projector = eye(2) - J_x' * J_w_pinv';
    tau_c = J_x'*f + projector*tau_n + g;
    
    % update states 
    tau = tau_c;
    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    tau_stack_4 = [tau_stack_4; tau(1) tau(2)];
    q_stack_4=[q_stack_4; q(1) q(2)];
    p_stack_4=[p_stack_4; p(1) p(2)];
end
toc
%% D = 10
D = 10;
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

tic;
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
    J_x = J(1, :);
    J_x_pinv = pinv(J_x);
    z = transpose(null(J_x));
    z_pinv = pinv(z);
    p_dot = J*q_dot;
    
    %task pd+gravity compensation controller (set point)
    p_des=[1.5; 0.0];
    K_p=50*eye(DOF);
    K_d=5*eye(DOF);
    tau_c= g+ J'*( -K_p*(p-p_des) - K_d*p_dot);
    
    p_des=1.5;
    K_p=50;
    K_d=5;
    W = M;
    tau_n = -D*(q_dot);
    J_w_pinv = inv(W)*J_x'*inv(J_x*inv(W)*J_x');
    f = - K_d*p_dot(1) + K_p*(p_des - p(1));
    projector = eye(2) - J_x' * J_w_pinv';
    tau_c = J_x'*f + projector*tau_n + g;
    
    % update states 
    tau = tau_c;
    xdot = set_xdot(M,h,g,q_dot,tau);
    x = x + dT*xdot;

    % stack variables (for plot)
    tau_stack_5 = [tau_stack_5; tau(1) tau(2)];
    q_stack_5=[q_stack_5; q(1) q(2)];
    p_stack_5=[p_stack_5; p(1) p(2)];
end
toc
%% plot
close all;

linewidth = 1;
fontsize = 12;

figure();
subplot(2, 1, 1);
plot(time, q_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, q_stack(:,2), 'linewidth',linewidth); 
legend('q_1','q_2');
title('Joint Space');
subplot(2, 1, 2);
plot(time, p_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, p_stack(:,2), 'linewidth',linewidth);
legend('x', 'y');
title('Task Space');
sgtitle('Joint Damping (D=0)')

figure();
subplot(2, 1, 1);
plot(time, q_stack_3(:,1), 'linewidth',linewidth); hold on;
plot(time, q_stack_3(:,2), 'linewidth',linewidth); 
legend('q_1','q_2');
title('Joint Space');
subplot(2, 1, 2);
plot(time, p_stack_3(:,1), 'linewidth',linewidth); hold on;
plot(time, p_stack_3(:,2), 'linewidth',linewidth);
legend('x', 'y');
title('Task Space');
sgtitle('Joint Damping (D=1.0)')

figure();
plot(time, q_stack(:,1), '-c', 'linewidth',linewidth); hold on;
plot(time, q_stack(:,2), ':c', 'linewidth',linewidth); 
plot(time, q_stack_2(:,1), '-m', 'linewidth',linewidth);
plot(time, q_stack_2(:,2), ':m', 'linewidth',linewidth); 
plot(time, q_stack_3(:,1), '-y', 'linewidth',linewidth); 
plot(time, q_stack_3(:,2), ':y', 'linewidth',linewidth); 
plot(time, q_stack_4(:,1), '-k', 'linewidth',linewidth); 
plot(time, q_stack_4(:,2), ':k', 'linewidth',linewidth);
legend('D=0, q_1','D=0, q_2', 'D=0.5, q_1','D=0.5, q_2', ...
    'D=1.0, q_1', 'D=1.0, q_2', 'D=1.5, q_1','D=1.5, q_2');
title('Degree of Joints (Joint Space)');

% figure();
% plot(time, p_stack(:,1), 'linewidth',linewidth); hold on;
% plot(time, p_stack(:,2), 'linewidth',linewidth); 
% plot(time, p_stack_2(:,1), 'linewidth',linewidth); 
% plot(time, p_stack_2(:,2), 'linewidth',linewidth); 
% plot(time, p_stack_3(:,1), 'linewidth',linewidth); 
% plot(time, p_stack_3(:,2), 'linewidth',linewidth); 
% plot(time, p_stack_4(:,1), 'linewidth',linewidth); 
% plot(time, p_stack_4(:,2), 'linewidth',linewidth);
% legend('p_x','p_y','Location','Best')

figure();
subplot(2, 1, 1);
plot(time, p_stack(:,1), 'linewidth',linewidth); hold on;
plot(time, p_stack_2(:,1), 'linewidth',linewidth); 
plot(time, p_stack_3(:,1), 'linewidth',linewidth); 
plot(time, p_stack_4(:,1), 'linewidth',linewidth); 
legend('D=0','D=0.5','D=1.0','D=1.5')
title('Position X');
subplot(2, 1, 2);
plot(time, p_stack(:,2), 'linewidth',linewidth); hold on;
plot(time, p_stack_2(:,2), 'linewidth',linewidth); 
plot(time, p_stack_3(:,2), 'linewidth',linewidth); 
plot(time, p_stack_4(:,2), 'linewidth',linewidth); 
legend('D=0','D=0.5','D=1.0','D=1.5');
title('Position Y');
sgtitle('Task Space')

figure();
plot(time, tau_stack_1(:,1), '-c', 'linewidth',linewidth); hold on;
plot(time, tau_stack_1(:,2), ':c', 'linewidth',linewidth); 
plot(time, tau_stack_2(:,1), '-m', 'linewidth',linewidth);
plot(time, tau_stack_2(:,2), ':m', 'linewidth',linewidth); 
plot(time, tau_stack_3(:,1), '-y', 'linewidth',linewidth); 
plot(time, tau_stack_3(:,2), ':y', 'linewidth',linewidth); 
plot(time, tau_stack_4(:,1), '-k', 'linewidth',linewidth); 
plot(time, tau_stack_4(:,2), ':k', 'linewidth',linewidth);
legend('D=0, \tau_1','D=0, \tau_2', 'D=0.5, \tau_1','D=0.5, \tau_2', ...
    'D=1.0, \tau_1', 'D=1.0, \tau_2', 'D=1.5, \tau_1','D=1.5, \tau_2');
title('Torque of Joints');

figure();
plot(time, tau_stack_5(:,1), '-k', 'linewidth',linewidth); hold on;
plot(time, tau_stack_5(:,2), ':k', 'linewidth',linewidth); 
legend('D=10.0, \tau_1','D=10.0, \tau_2');
title('Torque of Joints (D=10.0)');

figure();
plot(time, q_stack_5(:,1), '-k', 'linewidth',linewidth); hold on;
plot(time, q_stack_5(:,2), ':k', 'linewidth',linewidth); 
legend('D=10.0, q_1','D=10.0, q_2');
title('Degree of Joints (D=10.0)');

figure();
plot(time, p_stack_5(:,1), '-k', 'linewidth',linewidth); hold on;
plot(time, p_stack_5(:,2), ':k', 'linewidth',linewidth); 
legend('D=10.0, p_1','D=10.0, p_2');
title('Position of X (D=10.0)');