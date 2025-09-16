%==========================================================================
% Fully Symbolic Calculation of the Coriolis Matrix C(csi, csi_dot)
%==========================================================================
% This script computes the Coriolis and Centrifugal matrix C using the
% Christoffel symbols of the first type. It builds upon the previously
% defined symbolic kinematics and the inertia matrix B.

clear; clc;

%% 1. Define ALL Symbolic Variables and Constants
fprintf('1. Defining all symbolic variables...\n');

% -- Generalized Coordinates (csi) and their derivatives (csi_dot) --
csi = sym('csi', [9 1], 'real');
csi_dot = sym('csi_dot', [9 1], 'real');

% Extract individual coordinates for clarity
p_b = csi(1:3);
phi_b = csi(4:6);
q = csi(7:9);

% Extract individual Euler angles
phi = phi_b(1);
theta = phi_b(2);
psi = phi_b(3);

% -- Physical Constants and Parameters --
syms g pi real;
syms m_b m_l1 m_l2 m_l3 real;
syms l1 l2 offset real;
H_b = sym('H_b', [3 3], 'real');
H_l1 = sym('H_l1', [3 3], 'real');
H_l2 = sym('H_l2', [3 3], 'real');
H_l3 = sym('H_l3', [3 3], 'real');

%% 2. Build Inertia Matrix B (Symbolically)
% This section reconstructs the full 9x9 inertia matrix B.
fprintf('2. Building the symbolic Inertia Matrix B...\n');

% -- UAV Kinematics --
Rx = [1,0,0; 0,cos(phi),-sin(phi); 0,sin(phi),cos(phi)];
Ry = [cos(theta),0,sin(theta); 0,1,0; -sin(theta),0,cos(theta)];
Rz = [cos(psi),-sin(psi),0; sin(psi),cos(psi),0; 0,0,1];
R_b = Rz * Ry * Rx;

% Transformation matrix T_b and Q
Q = [   1     0           -sin(theta);
        0     cos(phi)    cos(theta)*sin(phi);
        0     -sin(phi)   cos(theta)*cos(phi)];    % It's easier to compute first Q   
T_b = R_b * Q;        % it's easier this way
%Q = R_b' * T_b;

% -- Helper Function --
skew = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

% -- Manipulator Kinematics and Jacobians --
R_b_0 = [1,           0,          0;
        0, cos(pi/2), -sin(pi/2);
        0, sin(pi/2),  cos(pi/2)];  % rotation of pi/2 about x-axis
T_b_0 = [   R_b_0, [0;0;-offset];      % translation along zb
            0,   0,   0,   1];
p_b_0 = T_b_0(1:3, 4);

a = [0, l1, l2]; d_dh = [0, 0, 0]; alpha = [pi/2, -pi/2, 0];

DH_matrix = @(alpha_i, a_i, d_i, q_i) ...
    [cos(q_i), -sin(q_i)*cos(alpha_i),  sin(q_i)*sin(alpha_i),  a_i*cos(q_i);
     sin(q_i),  cos(q_i)*cos(alpha_i), -cos(q_i)*sin(alpha_i),  a_i*sin(q_i);
     0,         sin(alpha_i),          cos(alpha_i),           d_i;
     0,         0,                     0,                      1];

A_0_1 = DH_matrix(alpha(1), a(1), d_dh(1), q(1));
A_1_2 = DH_matrix(alpha(2), a(2), d_dh(2), q(2));
T_0_1 = A_0_1; T_0_2 = A_0_1 * A_1_2;

T_b_1 = T_b_0 * T_0_1; R_b_1 = T_b_1(1:3, 1:3); p_b_1 = T_b_1(1:3, 4);
T_b_2 = T_b_0 * T_0_2; R_b_2 = T_b_2(1:3, 1:3); p_b_2 = T_b_2(1:3, 4);

p_b_l1 = p_b_0;
p_b_l2 = p_b_1 + R_b_1 * [l1/2; 0; 0];
p_b_l3 = p_b_2 + R_b_2 * [l2/2; 0; 0];
 
R_b_l1 = R_b_0; 
R_b_l2 = R_b_1; 
R_b_l3 = R_b_2; 

z0_b = R_b_0(:, 3); z1_b = R_b_1(:, 3); z2_b = R_b_2(:, 3);

J_P_l1 = zeros(3, 3, 'sym'); J_O_l1 = zeros(3, 3, 'sym');
J_P_l2 = sym(zeros(3,3)); J_O_l2 = sym(zeros(3,3));
J_P_l3 = sym(zeros(3,3)); J_O_l3 = sym(zeros(3,3));

J_O_l1(:,1) = z0_b;
J_P_l2(:,1) = skew(z0_b) * (p_b_l2 - p_b_0); J_P_l2(:,2) = skew(z1_b) * (p_b_l2 - p_b_1);
J_O_l2(:,1) = z0_b; J_O_l2(:,2) = z1_b;
J_P_l3(:,1) = skew(z0_b) * (p_b_l3 - p_b_0); J_P_l3(:,2) = skew(z1_b) * (p_b_l3 - p_b_1); J_P_l3(:,3) = skew(z2_b) * (p_b_l3 - p_b_2);
J_O_l3(:,1) = z0_b; J_O_l3(:,2) = z1_b; J_O_l3(:,3) = z2_b;

% -- Assemble B matrix blocks --
B11 = (m_b + m_l1 + m_l2 + m_l3) * eye(3);
B22_sum = m_l1*T_b'*skew(R_b*p_b_l1)'*skew(R_b*p_b_l1)*T_b + m_l2*T_b'*skew(R_b*p_b_l2)'*skew(R_b*p_b_l2)*T_b + m_l3*T_b'*skew(R_b*p_b_l3)'*skew(R_b*p_b_l3)*T_b + ...
          Q'*(R_b_l1*H_l1*R_b_l1' + R_b_l2*H_l2*R_b_l2' + R_b_l3*H_l3*R_b_l3')*Q;
B22 = Q'*H_b*Q + B22_sum;
B33 = (m_l1*(J_P_l1'*J_P_l1) + J_O_l1'*R_b_l1*H_l1*R_b_l1'*J_O_l1) + ...
      (m_l2*(J_P_l2'*J_P_l2) + J_O_l2'*R_b_l2*H_l2*R_b_l2'*J_O_l2) + ...
      (m_l3*(J_P_l3'*J_P_l3) + J_O_l3'*R_b_l3*H_l3*R_b_l3'*J_O_l3);
B12 = -(m_l1*skew(R_b*p_b_l1)*T_b + m_l2*skew(R_b*p_b_l2)*T_b + m_l3*skew(R_b*p_b_l3)*T_b);
B13 = (m_l1*R_b*J_P_l1 + m_l2*R_b*J_P_l2 + m_l3*R_b*J_P_l3);
B23_sum1 = Q'*(R_b_l1*H_l1*R_b_l1'*J_O_l1 + R_b_l2*H_l2*R_b_l2'*J_O_l2 + R_b_l3*H_l3*R_b_l3'*J_O_l3);
B23_sum2 = -(m_l1*T_b'*skew(R_b*p_b_l1)'*R_b*J_P_l1 + m_l2*T_b'*skew(R_b*p_b_l2)'*R_b*J_P_l2 + m_l3*T_b'*skew(R_b*p_b_l3)'*R_b*J_P_l3);
B23 = B23_sum1 + B23_sum2;

% Full B matrix
B = [B11, B12, B13;
     B12', B22, B23;
     B13', B23', B33];
% B = simplify(B);

%% 3. Compute Christoffel Symbols and Coriolis Matrix C
fprintf('3. Calculating Christoffel symbols (this may take some time)...\n');

n = length(csi);
C = sym(zeros(n, n));

% Loop through each element of the C matrix
for i = 1:n
    for j = 1:n
        sum_term = sym(0);
        % Summation over k
        for k = 1:n
            % Christoffel symbol of the first type: c_ijk
            dB_ij_d_csi_k = diff(B(i,j), csi(k));
            dB_ik_d_csi_j = diff(B(i,k), csi(j));
            dB_jk_d_csi_i = diff(B(j,k), csi(i));
            
            christoffel_ijk = 0.5 * (dB_ik_d_csi_j - dB_jk_d_csi_i + dB_ij_d_csi_k);
            
            sum_term = sum_term + christoffel_ijk * csi_dot(k);
        end
        C(i,j) = sum_term;
    end
end

% C = simplify(C);

%% 4. CODE GENERATION
fprintf('4. Generating the efficient numerical function ''calculate_C_matrix.m''...\n');

% Define the inputs for the function we are about to create.
% It needs the state, its derivative, and all the physical parameters.
input_vars = {csi, csi_dot, ...
              m_b, m_l1, m_l2, m_l3, ...
              l1, l2, offset, ...
              H_b, H_l1, H_l2, H_l3, ...
              pi};

% Generate the .m file
matlabFunction(C, 'File', 'calculate_C_matrix', 'Vars', input_vars, 'Outputs', {'C'});

fprintf('\nSUCCESS!\n');
fprintf('The file ''calculate_C_matrix.m'' has been generated in your current directory.\n');

