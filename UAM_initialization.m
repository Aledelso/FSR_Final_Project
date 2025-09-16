%% INITIALIZATION
p_b_init = [0; 0; 2];
eta_b_init = [0; 0; 0];
q_init = [-pi/2; 0; pi/2];
% p_e_init = [0.05; 0; 1.75];    % Computed using direct kinematics

% Define Final Conditions 
%p_b_final = [1; 1; 2.5]; % Target UAV position
%p_b_final = [0; 0; 2]; % Target UAV position
%eta_b_final = [0; 0; 0]; % Target UAV attitude 
%p_e_final = [1.05; 1; 2.25]; % Target end-effector position
% p_e_final = [0.05; 0; 1.75]; % Target end-effector position
% p_b_e_final = [0.05; 0; -0.25];
% p_e_final = [0.1561; 0; -0.2061]; % Target end-effector position
%p_e_final = [0.1561; 0; 1.7939]; % Target end-effector position

% Final Conditions
p_b_final = [0; 0; 2];                 % UAV in hovering
%p_b_final = [1; 1; 2.5];

eta_b_final = [0; 0; 0];

p_b_e_final = [0.05; 0; -0.25];        % Manipulator still
% p_b_e_final = [0.1561; 0; -0.2061];     % q_final = [-pi/2+pi/4; 0; pi/2-pi/4]



%% PARAMETERS
m_b = 2.0;
H_b = diag([1.24, 1.24, 2.48]);
m_l1 = 0;        % consider that the first two joints are in the same point
m_l2 = 0.049; 
m_l3 = 0.05;
l1 = 0.15;  % length link 1
l2 = 0.05;  % length link 2
offset = 0.1;   % distance along vertical axis between the center of UAV and the manipulator

% constant inertia matrix of link i
H_l1 = diag([0, 0, 0]);
H_l2 = diag([0, 0, 0.0011]); %inertia around the rotational axis (only about the joint-axis)
H_l3 = diag([0, 0, 1.25e-4]); %inertia around the rotational axis (only about the joint-axis)

%% Time 
Ts=0.001;   % Sampling Time

t_iniz = 0;
ttot = 15; % Total duration of the movement
tdead = 5; % Dead time at the end to hold the final position
tot_time = ttot + tdead;

t1 = linspace(t_iniz, ttot, round(ttot/Ts));
t = linspace(t_iniz, tot_time, round(tot_time/Ts));

%% DIRECT KINEMATICS
% For the Centralized Control you need p_e
% For the Decentralized Control you need p_b_e
[p_e_init, p_b_e_init] = Direct_kinematics(p_b_init, eta_b_init, q_init, l1, l2, offset);
% p_b_e_final = p_e_final - p_b_final;
p_e_final = p_b_final + p_b_e_final;

%% Centralized Control PLANNER
% PLANNER for UAM (Start -> Final - based on the professor's planner)
csi_init = [p_b_init; eta_b_init; q_init];
csi_dot_init = zeros(9,1);

x_d_start = [p_b_init; eta_b_init; p_e_init];
x_d = zeros(9,length(t));
dot_x_d = zeros(9,length(t));
ddot_x_d = zeros(9,length(t));

x_d_end = [p_b_final; eta_b_final; p_e_final];

% Define Boundary Conditions
dot_x_d_start = zeros(9, 1);
dot_x_d_end = zeros(9, 1);
ddot_x_d_start = zeros(9, 1);
ddot_x_d_end = zeros(9, 1);
dddot_x_d_start = zeros(9, 1);
dddot_x_d_end = zeros(9, 1);

p_d_full   = zeros(9, length(t1));
dot_p_d_full  = zeros(9, length(t1));
ddot_p_d_full = zeros(9, length(t1));

for j = 1:9 % Loop over all 9 Cartesian state variables
    % The polynomial coefficient calculation matrix A remains the same
    A = [t_iniz^7, t_iniz^6, t_iniz^5, t_iniz^4, t_iniz^3, t_iniz^2, t_iniz, 1;
        ttot^7, ttot^6, ttot^5, ttot^4, ttot^3, ttot^2, ttot, 1;
        7*t_iniz^6, 6*t_iniz^5, 5*t_iniz^4, 4*t_iniz^3, 3*t_iniz^2, 2*t_iniz, 1, 0;
        7*ttot^6, 6*ttot^5, 5*ttot^4, 4*ttot^3, 3*ttot^2, 2*ttot, 1, 0;
        42*t_iniz^5, 30*t_iniz^4, 20*t_iniz^3, 12*t_iniz^2, 6*t_iniz, 2, 0, 0;
        42*ttot^5, 30*ttot^4, 20*ttot^3, 12*ttot^2, 6*ttot, 2, 0, 0;
        210*t_iniz^4, 120*t_iniz^3, 60*t_iniz^2, 24*t_iniz, 6, 0, 0, 0;
        210*ttot^4, 120*ttot^3, 60*ttot^2, 24*ttot, 6, 0, 0, 0];

    b = [x_d_start(j), x_d_end(j), ...
         dot_x_d_start(j), dot_x_d_end(j), ...
         ddot_x_d_start(j), ddot_x_d_end(j), ...
         dddot_x_d_start(j), dddot_x_d_end(j)]';

    % Solve for polynomial coefficients [a7, a6, a5, a4, a3, a2, a1, a0]'
    a_coeffs = A\b;
    a7=a_coeffs(1); a6=a_coeffs(2); a5=a_coeffs(3); a4=a_coeffs(4);
    a3=a_coeffs(5); a2=a_coeffs(6); a1=a_coeffs(7); a0=a_coeffs(8);

    % Calculate trajectories
    p_d_full(j,:) = a7*t1.^7 + a6*t1.^6 + a5*t1.^5 + a4*t1.^4 + a3*t1.^3 + a2*t1.^2 + a1*t1 + a0;

    dot_p_d_full(j,:) = 7*a7*t1.^6 + 6*a6*t1.^5 + 5*a5*t1.^4 + 4*a4*t1.^3 + 3*a3*t1.^2 + 2*a2*t1 + a1;

    ddot_p_d_full(j,:) = 42*a7*t1.^5 + 30*a6*t1.^4 + 20*a5*t1.^3 + 12*a4*t1.^2 + 6*a3*t1 + 2*a2;

    x_d(j,:) = [p_d_full(j,:) zeros(1,round(tdead/Ts))+x_d_end(j)];
    dot_x_d(j,:) = [dot_p_d_full(j,:) zeros(1,round(tdead/Ts))+ dot_x_d_end(j)];
    ddot_x_d(j,:) = [ddot_p_d_full(j,:) zeros(1,round(tdead/Ts))+ ddot_x_d_end(j)];

end

% data for simulink
x_d_dot = dot_x_d;
x_d_ddot = ddot_x_d;

%% Decentralized Control PLANNER UAV
% The passivity based control needs the z-axis pointing down
% Rotation about X-axis
R_pass_imp = [  1,  0,          0;
                0,  cos(pi),    -sin(pi);
                0,  sin(pi),    cos(pi)];

eta_b_init = R_pass_imp * eta_b_init;
eta_b_final = R_pass_imp * eta_b_final;

p= zeros(4,length(t1)); dot_p= zeros(4,length(t1)); ddot_p= zeros(4,length(t1));
p_d= zeros(4,length(t)); dot_p_d= zeros(4,length(t)); ddot_p_d= zeros(4,length(t));
%Initial and final conditions
x0= p_b_init(1); xf=p_b_final(1); dot_x0= 0; dot_xf=0; ddot_x0=0; ddot_xf=0; dddot_x0 = 0; dddot_xf = 0;
y0= -p_b_init(2); yf=-p_b_final(2); dot_y0= 0; dot_yf=0; ddot_y0=0; ddot_yf=0; dddot_y0 = 0; dddot_yf = 0;
z0= -p_b_init(3); zf=-p_b_final(3); dot_z0= 0; dot_zf=0; ddot_z0=0; ddot_zf=0; dddot_z0 = 0; dddot_zf = 0;
psi0 = eta_b_init(3);
psif=eta_b_final(3); dot_psi0= 0; dot_psif=0; ddot_psi0=0; ddot_psif=0; dddot_psi0 = 0; dddot_psif = 0;
p0=[x0,y0,z0,psi0]; dot_p0=[dot_x0,dot_y0,dot_z0,dot_psi0]; ddot_p0=[ddot_x0,ddot_y0,ddot_z0,ddot_psi0]; dddot_p0=[dddot_x0,dddot_y0,dddot_z0,dddot_psi0];
pf=[xf,yf,zf,psif]; dot_pf=[dot_xf,dot_yf,dot_zf,dot_psif]; ddot_pf=[ddot_xf,ddot_yf,ddot_zf,ddot_psif]; dddot_pf=[dddot_xf,dddot_yf,dddot_zf,dddot_psif];
%5-th order polynomial
a0=zeros(1,4); a1=zeros(1,4); a2=zeros(1,4); a3=zeros(1,4); a4=zeros(1,4); a5=zeros(1,4); a6 = zeros(1,4); a7 = zeros(1,4);

for j=1:4
    A = [t_iniz^7, t_iniz^6, t_iniz^5, t_iniz^4, t_iniz^3, t_iniz^2, t_iniz, 1;
        ttot^7, ttot^6, ttot^5, ttot^4, ttot^3, ttot^2, ttot, 1;
        7*t_iniz^6, 6*t_iniz^5, 5*t_iniz^4, 4*t_iniz^3, 3*t_iniz^2, 2*t_iniz, 1, 0;
        7*ttot^6, 6*ttot^5, 5*ttot^4, 4*ttot^3, 3*ttot^2, 2*ttot, 1, 0;
        42*t_iniz^5, 30*t_iniz^4, 20*t_iniz^3, 12*t_iniz^2, 6*t_iniz, 2, 0, 0;
        42*ttot^5, 30*ttot^4, 20*ttot^3, 12*ttot^2, 6*ttot, 2, 0, 0;
        210*t_iniz^4, 120*t_iniz^3, 60*t_iniz^2, 24*t_iniz, 6, 0, 0, 0;
        210*ttot^4, 120*ttot^3, 60*ttot^2, 24*ttot, 6, 0, 0, 0];
    b = [p0(j) pf(j) dot_p0(j) dot_pf(j) ddot_p0(j) ddot_pf(j) dddot_p0(j) dddot_pf(j)]';
    a_temp = A\b;
    a7(j) = a_temp(1);
    a6(j) = a_temp(2);
    a5(j) = a_temp(3);
    a4(j) = a_temp(4);
    a3(j) = a_temp(5);
    a2(j) = a_temp(6);
    a1(j) = a_temp(7);
    a0(j) = a_temp(8);
    
    %trajectories
    p(j,:)=a7(j)*t1.^7 + a6(j)*t1.^6 + a5(j)*t1.^5 +a4(j)*t1.^4 +a3(j)*t1.^3 +a2(j)*t1.^2 +a1(j)*t1 +a0(j);
    dot_p(j,:) = 7*a7(j)*t1.^6 + 6*a6(j)*t1.^5 + 5*a5(j)*t1.^4 +4*a4(j)*t1.^3 +3*a3(j)*t1.^2 +2*a2(j)*t1 +a1(j);
    ddot_p(j,:) = 42*a7(j)*t1.^5 + 30*a6(j)*t1.^4 + 5*4*a5(j)*t1.^3 +4*3*a4(j)*t1.^2 +3*2*a3(j)*t1 +2*a2(j);
    %addition of the stead-state terms
    p_d(j,:)=[p(j,:) zeros(1,round(tdead/Ts))+pf(j)]; 
    dot_p_d(j,:)=[dot_p(j,:) zeros(1,round(tdead/Ts))+dot_pf(j)];
    ddot_p_d(j,:)=[ddot_p(j,:) zeros(1,round(tdead/Ts))+ddot_pf(j)];
end
%data for Simulink
pos_0 = [x0 y0 z0];
lin_vel_0 = [dot_x0 dot_y0 dot_z0];

p_b_d=p_d(1:3,:); dot_p_b_d=dot_p_d(1:3,:); ddot_p_b_d=ddot_p_d(1:3,:);
psi_d=p_d(4,:); dot_psi_d=dot_p_d(4,:); ddot_psi_d=ddot_p_d(4,:);

p_b_dot_init = lin_vel_0';
p_b_init = pos_0';
eta_b_init = [0;0;psi0];
eta_b_dot_init = [0;0;0];

%% Decentralized Control PLANNER MANIPULATOR
% Define Initial Conditions
x_d_start = p_b_e_init;
p_b_e_d = zeros(3,length(t));
dot_x_d = zeros(3,length(t));
ddot_x_d = zeros(3,length(t));
% Define Final Conditions 
x_d_end = p_b_e_final;

% Define Boundary Conditions
dot_x_d_start = zeros(3, 1);
dot_x_d_end = zeros(3, 1);
ddot_x_d_start = zeros(3, 1);
ddot_x_d_end = zeros(3, 1);
dddot_x_d_start = zeros(3, 1);
dddot_x_d_end = zeros(3, 1);

p_d_full   = zeros(3, length(t1));
dot_p_d_full  = zeros(3, length(t1));
ddot_p_d_full = zeros(3, length(t1));

for j = 1:3 % Loop over all 9 Cartesian state variables
    % The polynomial coefficient calculation matrix A remains the same
    A = [t_iniz^7, t_iniz^6, t_iniz^5, t_iniz^4, t_iniz^3, t_iniz^2, t_iniz, 1;
        ttot^7, ttot^6, ttot^5, ttot^4, ttot^3, ttot^2, ttot, 1;
        7*t_iniz^6, 6*t_iniz^5, 5*t_iniz^4, 4*t_iniz^3, 3*t_iniz^2, 2*t_iniz, 1, 0;
        7*ttot^6, 6*ttot^5, 5*ttot^4, 4*ttot^3, 3*ttot^2, 2*ttot, 1, 0;
        42*t_iniz^5, 30*t_iniz^4, 20*t_iniz^3, 12*t_iniz^2, 6*t_iniz, 2, 0, 0;
        42*ttot^5, 30*ttot^4, 20*ttot^3, 12*ttot^2, 6*ttot, 2, 0, 0;
        210*t_iniz^4, 120*t_iniz^3, 60*t_iniz^2, 24*t_iniz, 6, 0, 0, 0;
        210*ttot^4, 120*ttot^3, 60*ttot^2, 24*ttot, 6, 0, 0, 0];
        
    b = [x_d_start(j), x_d_end(j), ...
         dot_x_d_start(j), dot_x_d_end(j), ...
         ddot_x_d_start(j), ddot_x_d_end(j), ...
         dddot_x_d_start(j), dddot_x_d_end(j)]';

    % Solve for polynomial coefficients [a7, a6, a5, a4, a3, a2, a1, a0]'
    a_coeffs = A\b;
    a7=a_coeffs(1); a6=a_coeffs(2); a5=a_coeffs(3); a4=a_coeffs(4);
    a3=a_coeffs(5); a2=a_coeffs(6); a1=a_coeffs(7); a0=a_coeffs(8);
    
    % Calculate trajectories
    p_d_full(j,:) = a7*t1.^7 + a6*t1.^6 + a5*t1.^5 + a4*t1.^4 + a3*t1.^3 + a2*t1.^2 + a1*t1 + a0;
    
    dot_p_d_full(j,:) = 7*a7*t1.^6 + 6*a6*t1.^5 + 5*a5*t1.^4 + 4*a4*t1.^3 + 3*a3*t1.^2 + 2*a2*t1 + a1;
    
    ddot_p_d_full(j,:) = 42*a7*t1.^5 + 30*a6*t1.^4 + 20*a5*t1.^3 + 12*a4*t1.^2 + 6*a3*t1 + 2*a2;
    
    p_b_e_d(j,:) = [p_d_full(j,:) zeros(1,round(tdead/Ts))+x_d_end(j)];
    dot_x_d(j,:) = [dot_p_d_full(j,:) zeros(1,round(tdead/Ts))+ dot_x_d_end(j)];
    ddot_x_d(j,:) = [ddot_p_d_full(j,:) zeros(1,round(tdead/Ts))+ ddot_x_d_end(j)];
    
end

% data for simulink
p_b_e_d_dot = dot_x_d;
p_b_e_d_ddot = ddot_x_d;
q_dot_init = zeros(3,1);

%% Estimator for the UAV in the Decentralized Control
c0 = 1;
r = 2;

c0 = c0 / sqrt(2^(1/r) - 1);
% Calculate coefficients c and gains K
s = tf('s');
G = c0^r/(s+c0)^r;
c = cell2mat(G.Denominator);
c = c(2:end);

K = zeros(r, 1);
product_K = 1;
for j=1:r
    K(j) = c(j)/product_K;
    product_K = product_K*K(j);
end
K = flip(K);


%% FUNCTIONS
function [p_e, p_b_e] = Direct_kinematics(p_b, eta_b, q, l1, l2, offset)

% Euler angles
phi = eta_b(1);    % Yaw angle
theta = eta_b(2);  % Pitch angle
psi = eta_b(3);    % Roll angle

% Individual Rotation Matrices 
% Rotation about X-axis (Roll)
Rx = [1,           0,          0;
      0, cos(phi), -sin(phi);
      0, sin(phi),  cos(phi)];

% Rotation about Y-axis (Pitch)
Ry = [cos(theta),  0, sin(theta);
      0,           1,          0;
     -sin(theta),  0, cos(theta)];

% Rotation about Z-axis (Yaw)
Rz = [cos(psi), -sin(psi), 0;
      sin(psi),  cos(psi), 0;
      0,         0,         1];

% Combined Rotation Matrix R_b (from UAV body frame to Inertial frame)
% Using Z-Y-X Euler angle convention (Yaw-Pitch-Roll)
R_b = Rz * Ry * Rx;

% Ab: Homogeneous transformation matrix from Inertial frame to UAV Body frame
A_b = [R_b,   p_b;
      0 0 0, 1];

% Manipulator Direct Kinematics (Relative to UAV Body Frame)
% Static Transformation to take into account the offset between UAV and Manipulator
R_b_0 = [1,           0,          0;
      0, cos(pi/2), -sin(pi/2);
      0, sin(pi/2),  cos(pi/2)];

T_b_0 = [   R_b_0, [0;0;-offset];
            0,   0,   0,   1];

% D-H parameters for the manipulator
a = [0, l1, l2];
d = [0, 0, 0];
alpha = [pi/2, -pi/2, 0];

% Standard DH transformation matrix A(alpha, a, d, theta) - here theta is q_i
DH_matrix = @(alpha_i, a_i, d_i, q_i) ...
    [cos(q_i), -sin(q_i)*cos(alpha_i),  sin(q_i)*sin(alpha_i),  a_i*cos(q_i);
     sin(q_i),  cos(q_i)*cos(alpha_i), -cos(q_i)*sin(alpha_i),  a_i*sin(q_i);
     0,         sin(alpha_i),          cos(alpha_i),           d_i;
     0,         0,                     0,                      1];

A_0_1 = DH_matrix(alpha(1), a(1), d(1), q(1));
A_1_2 = DH_matrix(alpha(2), a(2), d(2), q(2));
A_2_3 = DH_matrix(alpha(3), a(3), d(3), q(3));

% Transformations from manipulator base frame (Frame 0) to each link frame
T_0_3 = A_0_1 * A_1_2 * A_2_3; % Transformation from Frame 0 to Frame 3 (end-effector)


% Rotation about the y-axis
R_3_e = [cos(pi/2),  0, sin(pi/2);
        0,           1,          0;
        -sin(pi/2),  0, cos(pi/2)];
% To align the z-axis to the approach axis
T_3_e = [R_3_e,          zeros(3,1);
         zeros(1,3),    1];

% Direct kinematic of the Manipulator from UAV body frame to end-effector frame
A_b_e = T_b_0 * T_0_3 * T_3_e;

% Extract position of end-effector frame origin in UAV body frame (p_b_e)
% p_e^b is position vector of the origin of end-effector frame wrt Sigma_b 
p_b_e = A_b_e(1:3, 4);


% Absolute Pose of the End-Effector
% A_e: Absolute homogeneous transformation matrix of the end-effector (in Inertial frame)
A_e = A_b * A_b_e;

% Extract p_e (position of end-effector in inertial frame)
p_e = A_e(1:3, 4); % This is p_e from x = [p_b^T, phi_b^T, p_e^T]^T 

% % Transformation to bring the z-axis pointing down for the PASSIVITY control
% R_pass_b = [1,           0,          0;
%             0, cos(pi), -sin(pi);
%             0, sin(pi),  cos(pi)];      % rotation about x-axis
% T_pass_b = [R_pass_b, zeros(3,1);
%             zeros(1,3), 1];
% T_pass_e = T_pass_b * A_b_e;
% p_pass_e = T_pass_e(1:3, 4);

end