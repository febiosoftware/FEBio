% Define box vertices
vertices = [2 2; 2 -2; -2 -2; -2 2; 2 2];
x_vec = [1,0];
y_vec = [0,1];
d1_vec = [1/sqrt(2),1/sqrt(2)];
d2_vec = [-1/sqrt(2),1/sqrt(2)];

% Plot original box
figure(1);
clf
subplot(2, 3, 1);
plot(vertices(:,1), vertices(:,2), 'k-', 'LineWidth', 2); % Making the lines thicker
hold on;
quiver(0, 0, x_vec(1), x_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, y_vec(1), y_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
quiver(0, 0, d1_vec(1), d1_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, d2_vec(1), d2_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color

title('Original');
axis equal;
xlim([-5 5]); % Setting xlim
ylim([-5 5]); % Setting ylim

% Apply pure shear deformation gradient
F = [1 0; 0 2]; % Pure shear deformation gradient
C = F'*F;
lambda_t = sqrt(C(1,1));
lambda_z = sqrt(C(2,2));

% Polar decomposition of F
[U,S,V] = svd(F);
R = U * V'; % Rotation component

% Apply rotation to the original vertices
deformed_vertices = vertices * F';
rotated_x_vec = x_vec * R';
rotated_y_vec = y_vec * R';
rotated_d1_vec = d1_vec * R';
rotated_d2_vec = d2_vec * R';

% Plot deformed box by rotation component
subplot(2, 3, 2);
plot(deformed_vertices(:,1), deformed_vertices(:,2), 'k-', 'LineWidth', 2); % Making the lines thicker
hold on;
quiver(0, 0, rotated_x_vec(1), rotated_x_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, rotated_y_vec(1), rotated_y_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
quiver(0, 0, rotated_d1_vec(1), rotated_d1_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, rotated_d2_vec(1), rotated_d2_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
title('Deposition by Polar Decomposition');
axis equal;
xlim([-5 5]); % Setting xlim
ylim([-5 5]); % Setting ylim

% Apply pure shear deformation gradient
deformed_vertices = vertices * F';
deformed_x_vec = x_vec * F';
deformed_y_vec = y_vec * F';
deformed_x_vec = deformed_x_vec / norm(deformed_x_vec);
deformed_y_vec = deformed_y_vec / norm(deformed_y_vec);
deformed_d1_vec = d1_vec * F';
deformed_d2_vec = d2_vec * F';
deformed_d1_vec = deformed_d1_vec / norm(deformed_d1_vec);
deformed_d2_vec = deformed_d2_vec / norm(deformed_d2_vec);

% Plot deformed box by pure shear deformation
subplot(2, 3, 3);
plot(deformed_vertices(:,1), deformed_vertices(:,2), 'k-', 'LineWidth', 2); % Making the lines thicker
hold on;
quiver(0, 0, deformed_x_vec(1), deformed_x_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, deformed_y_vec(1), deformed_y_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
quiver(0, 0, deformed_d1_vec(1), deformed_d1_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, deformed_d2_vec(1), deformed_d2_vec(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
title('Deposition by Deformation Gradient');
axis equal;
xlim([-5 5]); % Setting xlim
ylim([-5 5]); % Setting ylim



% Apply pure shear deformation gradient
deformed_vertices = vertices * F';
gamma = 1;
a_x = atan(tan(0)*((lambda_t/lambda_z)^gamma));
a_y = atan(tan(pi/2)*((lambda_t/lambda_z)^gamma));
deformed_x_vec = [cos(a_x), sin(a_x)];
deformed_y_vec = [cos(a_y), sin(a_y)];

% Plot deformed box by pure shear deformation
subplot(2, 3, 4);
plot(deformed_vertices(:,1), deformed_vertices(:,2), 'k-', 'LineWidth', 2); % Making the lines thicker
hold on;
quiver(0, 0, deformed_x_vec(1), deformed_x_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, deformed_y_vec(1), deformed_y_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
title('tan(a_{curr}) = (\lambda_t / \lambda_z)^\gamma tan(a_{o}), \gamma = 1');
axis equal;
xlim([-5 5]); % Setting xlim
ylim([-5 5]); % Setting ylim



% Apply pure shear deformation gradient
deformed_vertices = vertices * F';
gamma = 0;
a_x = atan(tan(0)*((lambda_t/lambda_z)^gamma));
a_y = atan(tan(pi/2)*((lambda_t/lambda_z)^gamma));
deformed_x_vec = [cos(a_x), sin(a_x)];
deformed_y_vec = [cos(a_y), sin(a_y)];

% Plot deformed box by pure shear deformation
subplot(2, 3, 5);
plot(deformed_vertices(:,1), deformed_vertices(:,2), 'k-', 'LineWidth', 2); % Making the lines thicker
hold on;
quiver(0, 0, deformed_x_vec(1), deformed_x_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing up, gray color
quiver(0, 0, deformed_y_vec(1), deformed_y_vec(2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 1); % Vector pointing to the side, gray color
title('tan(a_{curr}) = (\lambda_t / \lambda_z)^\gamma tan(a_{o}), \gamma = 0');
axis equal;
xlim([-5 5]); % Setting xlim
ylim([-5 5]); % Setting ylim