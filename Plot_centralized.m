% === Plot: Quadrotor Position Error Norm ===
figure('Renderer', 'painters', 'Position', [10 10 900 350]);
plot(out.tout, out.e_p_b, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\|e_{p_b}\|$ [m]', 'Interpreter', 'latex');
title('Quadrotor Position Error Norm', 'Interpreter', 'latex');
set(gca, 'FontSize', 16); grid on; box on;
axis padded;
xlim([out.tout(1), out.tout(end)]);
% exportgraphics(gcf, 'dec_quad_pos_error_normD.pdf', 'ContentType', 'vector');

% === Plot: Quadrotor Orientation Error Norm ===
figure('Renderer', 'painters', 'Position', [10 10 900 350]);
plot(out.tout, out.e_eta_b, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\|e_{\eta_b}\|$ [rad]', 'Interpreter', 'latex');
title('Quadrotor Orientation Error Norm', 'Interpreter', 'latex');
set(gca, 'FontSize', 16); grid on; box on;
axis padded;
xlim([out.tout(1), out.tout(end)]);
exportgraphics(gcf, 'quad_ori_error_normA1.pdf', 'ContentType', 'vector');

% === Plot: Manipulator End Effector Position Error Norm ===
figure('Renderer', 'painters', 'Position', [10 10 900 350]);
plot(out.tout, out.e_p_e, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\|e_{p_e}\|$ [m]', 'Interpreter', 'latex');
title('Manipulator End Effector Position Error Norm', 'Interpreter', 'latex');
set(gca, 'FontSize', 16); grid on; box on;
axis padded;
xlim([out.tout(1), out.tout(end)]);
% exportgraphics(gcf, 'dec_ee_pos_error_normD.pdf', 'ContentType', 'vector');

% === Plot: Manipulator End Effector Position Error Components ===
figure('Renderer', 'painters', 'Position', [10 10 900 350]);

plot(out.tout, out.e_p_e_full(:,1), 'k-', 'LineWidth', 2); hold on;
plot(out.tout, out.e_p_e_full(:,2), 'b-', 'LineWidth', 2);
plot(out.tout, out.e_p_e_full(:,3), 'r-', 'LineWidth', 2);

xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_{p_e}$ [m]', 'Interpreter', 'latex');
title('Manipulator End Effector Position Error Components', 'Interpreter', 'latex');
legend({'$e_{x}$','$e_{y}$','$e_{z}$'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 16); grid on; box on;

% exportgraphics(gcf, 'dec_ee_pos_error_components.pdf', 'ContentType', 'vector');

%%
% === Plot: Manipulator Joint Torques ===
figure('Renderer', 'painters', 'Position', [10 10 900 350]);
plot(out.tout, out.q(:,7), 'k-', 'LineWidth', 2); hold on;
plot(out.tout, out.q(:,8), 'b-', 'LineWidth', 2);
plot(out.tout, out.q(:,9), 'r-', 'LineWidth', 2);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\tau$ [Nm]', 'Interpreter', 'latex');
title('Manipulator Joint Torques', 'Interpreter', 'latex');
legend({'$q_{1}$','$q_{2}$','$q_{3}$'}, 'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontSize', 16); grid on; box on;
yticks([-90 -0.5 0 0.5 90]);
% exportgraphics(gcf, 'cen_manipulator_joint_torques.pdf', 'ContentType', 'vector');


