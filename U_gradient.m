%==========================================================================
% Fully Symbolic Calculation of the Gravity Vector g(csi)
%==========================================================================

clear; clc;

%% 1. Define ALL Symbolic Variables and Constants
fprintf('1. Defining all symbolic variables (including pi)...\n');

% -- Generalized Coordinates (csi) --
p_b = sym('p_b', [3 1], 'real');      % Position of UAV [px; py; pz]
phi_b = sym('phi_b', [3 1], 'real');    
q = sym('q', [3 1], 'real');          % Manipulator joint angles [q1; q2; q3]

% Extract individual Euler angles for clarity
phi = phi_b(1);    % Corresponds to phi_b1 in output
theta = phi_b(2);  % Corresponds to phi_b2 in output
psi = phi_b(3);    % Corresponds to phi_b3 in output

% -- Physical Constants and Parameters --
syms g pi real; % <--- CRITICAL CHANGE: 'pi' is now a symbol
syms m_b m_l1 m_l2 m_l3 real;
syms l1 l2 offset real;

%% 2. Kinematics and Jacobians (Symbolic)
fprintf('2. Building symbolic kinematics and Jacobians...\n');

% -- UAV Rotation Matrix (Body to Inertial) --
Rx = [1,           0,          0;
      0, cos(phi), -sin(phi);
      0, sin(phi),  cos(phi)];
Ry = [cos(theta),  0, sin(theta);
      0,           1,          0;
     -sin(theta),  0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0;
      sin(psi),  cos(psi), 0;
      0,         0,         1];
R_b = Rz * Ry * Rx;

% -- Helper Function for Skew-Symmetric Matrices --
skew = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

% -- Manipulator Kinematics (Relative to UAV Body Frame) --
R_b_0 = [1,           0,          0;
        0, cos(pi/2), -sin(pi/2);
        0, sin(pi/2),  cos(pi/2)];  % rotation of pi/2 about x-axis
T_b_0 = [   R_b_0, [0;0;-offset];      % translation along zb
            0,   0,   0,   1];
p_b_0 = T_b_0(1:3, 4);

a = [0, l1, l2];
d = [0, 0, 0];
alpha = [pi/2, -pi/2, 0];

DH_matrix = @(alpha_i, a_i, d_i, q_i) ...
    [cos(q_i), -sin(q_i)*cos(alpha_i),  sin(q_i)*sin(alpha_i),  a_i*cos(q_i);
     sin(q_i),  cos(q_i)*cos(alpha_i), -cos(q_i)*sin(alpha_i),  a_i*sin(q_i);
     0,         sin(alpha_i),          cos(alpha_i),           d_i;
     0,         0,                     0,                      1];

A_0_1 = DH_matrix(alpha(1), a(1), d(1), q(1));
A_1_2 = DH_matrix(alpha(2), a(2), d(2), q(2));
T_0_1 = A_0_1; T_0_2 = A_0_1 * A_1_2;

% -- Link CoM Positions and Jacobians (w.r.t. UAV Body Frame) --
T_b_1 = T_b_0 * T_0_1; R_b_1 = T_b_1(1:3, 1:3); p_b_1 = T_b_1(1:3, 4);
T_b_2 = T_b_0 * T_0_2; R_b_2 = T_b_2(1:3, 1:3); p_b_2 = T_b_2(1:3, 4);

p_b_l1 = p_b_0; % Link 1 CoM at manipulator base
p_b_l2 = p_b_1 + R_b_1 * [l1/2; 0; 0];
p_b_l3 = p_b_2 + R_b_2 * [l2/2; 0; 0];

z0_b = R_b_0(:, 3); z1_b = R_b_1(:, 3); z2_b = R_b_2(:, 3);

J_P_l1 = zeros(3, 3, 'sym');
J_P_l2 = sym(zeros(3,3)); J_P_l3 = sym(zeros(3,3));

J_P_l2(:,1) = skew(z0_b) * (p_b_l2 - p_b_0);
J_P_l2(:,2) = skew(z1_b) * (p_b_l2 - p_b_1);

J_P_l3(:,1) = skew(z0_b) * (p_b_l3 - p_b_0);
J_P_l3(:,2) = skew(z1_b) * (p_b_l3 - p_b_1);
J_P_l3(:,3) = skew(z2_b) * (p_b_l3 - p_b_2);

%% 3. Formulate Total Potential Energy (U)
fprintf('3. Formulating total potential energy U...\n');
e3 = [0; 0; 1];
m_total = m_b + m_l1 + m_l2 + m_l3;

p_weighted_com_b = m_l1*p_b_l1 + m_l2*p_b_l2 + m_l3*p_b_l3;
U_total = m_total * g * e3.' * p_b + g * e3.' * R_b * p_weighted_com_b;

%% 4. Manually Compute the Gravity Vector g(csi) = (dU/dcsi)^T
fprintf('4. Manually computing gradient g(csi)...\n');

g_pb = gradient(m_total * g * p_b(3), p_b);
U_attitude_dependent = g * e3.' * R_b * p_weighted_com_b;
g_phib = gradient(U_attitude_dependent, phi_b);
J_P_weighted = m_l1*J_P_l1 + m_l2*J_P_l2 + m_l3*J_P_l3;
dU_dq_row_vector = g * e3.' * R_b * J_P_weighted;
g_q = dU_dq_row_vector.';
g_csi = [g_pb; g_phib; g_q];

%% 5. Display the Results
fprintf('\n=========================================================\n');
fprintf('Gravity Vector g(csi) = (dU/dcsi)^T:\n');
fprintf('=========================================================\n');
% Use simplify to get the cleanest possible output
% pretty(simplify(g_csi));
% simplify(g_csi)
simplify(g_phib)
