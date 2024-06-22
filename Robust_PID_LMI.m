% Define the range of uncertain parameters
a1_min = -2; a1_max = -25;
a2_min = 8; a2_max = 50;
b2_min = 1; b2_max = 40;

% Define the system matrices for the polytopic vertices
A1 = [0 1 0; -a2_min -a1_min 0; 1 0 0];
A2 = [0 1 0; -a2_max -a1_max 0; 1 0 0];
B1 = [0; b2_min; 0];
B2 = [0; b2_max; 0];
B_r = [0; 0; -1];

% Define Q and R matrices for LQR-LMI
Q = [100 0 0;0 1 0;0 0 1];  % Q as a 3x3 identity matrix
R = 1;       % R as a scalar
gamma = 10;  % The constraint on x^T(0) P x(0)

% Use the LMI toolbox to set up the problem
setlmis([]);

% Define LMI variables
P_hat = lmivar(1, [3, 1]); % P_hat is a 3x3 symmetric matrix
Y = lmivar(2, [1, 3]);     % Y is a 1x3 matrix

% Define the LMIs for all vertices
A_vertices = {A1, A2};
B_vertices = {B1, B2};

for i = 1:length(A_vertices)
    A = A_vertices{i};
    B = B_vertices{i};

    % LMI: [A*P_hat + P_hat*A' + B*Y + Y'*B', P_hat, Y'; P_hat, -Q, 0; Y, 0, -R] < 0
    lmiterm([1 1 1 P_hat], A, 1, 's');  % A*P_hat + P_hat*A'
    lmiterm([1 1 1 Y], B, 1, 's');      % B*Y
    lmiterm([1 2 1 P_hat], 1, 1);       % P_hat
    lmiterm([1 3 1 Y], 1, 1);           % Y
    lmiterm([1 2 2 0], -Q);             % -Q
    lmiterm([1 3 3 0], -R);             % -R
    lmiterm([1 3 2 0], zeros(1, 3));    % 1x3 zero matrix
end

% LMI: P_hat > 0
lmiterm([-2 1 1 P_hat], 1, 1);

% Constraint: x^T(0) P x(0) <= gamma
lmiterm([3 1 1 P_hat], 1, 1);
lmiterm([3 1 1 0], -gamma);

% Solve the LMI
lmis = getlmis;
[tmin, xfeas] = feasp(lmis);

% Check if a feasible solution was found
if tmin < 0
    % Extract P_hat and Y from the solution
    P_hat_sol = dec2mat(lmis, xfeas, P_hat);
    Y_sol = dec2mat(lmis, xfeas, Y);

    % Compute the feedback gain K
    K = -Y_sol / P_hat_sol;

    % Display the results
    disp('The feedback (LMI Model) gain K is:');
    disp(K);
else
    disp('No feasible solution found.');
end

%Defining the system for the lmi-controller simulation
% Randomly select parameters within their intervals
a1 = a1_min + (a1_max - a1_min) * rand;
a2 = a2_min + (a2_max - a2_min) * rand;
b2 = b2_min + (b2_max - b2_min) * rand;
A = [0 1 0; -a2 -a1 0; 1 0 0];
B = [0; b2; 0];
B_r = [0; 0; -1];
C = [1 0 0];
D = 0;

%Define the time scale for the simulations
t=0:0.01:60;

% Simulate the step response for the LMI-LQR controller
sys_cl_lmi = ss(A - B * K, B_r, C, D); % Using extracted parameters
[y_lmi, t_lmi, x_lmi] = step(sys_cl_lmi, t);

% Extract step response characteristics from LMI-LQR controller
info_lmi = stepinfo(y_lmi, t_lmi);


% Define the system matrices for the fixed parameter values
A_min = [0 1 0; -a2_min -a1_min 0; 1 0 0];
A_max = [0 1 0; -a2_max -a1_max 0; 1 0 0];
B_min = [0; b2_min; 0];
B_max = [0; b2_max; 0];
B_r = [0; 0; -1];


% Design traditional controllers using step response characteristics from LMI-LQR
lmi_poles = pole(sys_cl_lmi);

% Design state feedback controllers using pole placement
K_min = place(A_min, B_min, lmi_poles);
K_max = place(A_max, B_max, lmi_poles);

% Closed-loop system matrices for the traditional controllers
Ac_min = A - B * K_min;
Ac_max = A - B * K_max;

% Define the state-space systems
sys_cl_min = ss(Ac_min, B_r, C, D);
sys_cl_max = ss(Ac_max, B_r, C, D);



% Simulate the step response for the traditional controllers
[y_min, t_min, x_min] = step(sys_cl_min,t);
[y_max, t_max, x_max] = step(sys_cl_max,t);

% Extract step response characteristics from min and max controller
info_min = stepinfo(y_min, t_min);
info_max = stepinfo(y_max, t_max);

% Plot the step responses for comparison
figure;
plot(t_lmi, y_lmi, 'k', t_min, y_min, 'r', t_max, y_max, 'b');
xlabel('Time (s)');
ylabel('Output');
ylim([0 1.5]);
title('Step Response Comparison');
legend('LMI-LQR Controller', 'Traditional Controller (Min Parameters)', 'Traditional Controller (Max Parameters)');
grid on;

% Display the step response characteristics of the LMI-LQR controller
disp('Step response characteristics of the LMI-LQR controller:');
disp(['Settling Time: ', num2str(info_lmi.SettlingTime)]);
disp(['Rise Time: ', num2str(info_lmi.RiseTime)]);
disp(['Peak Time: ', num2str(info_lmi.PeakTime)]);
disp(['Maximum Overshoot: ', num2str(info_lmi.Overshoot)]);
disp(['Steady-State Error: ', num2str(abs(1 - y_lmi(end)))]); % Assuming a unit step input


% Display the step response characteristics of the minimum controller
disp('Step response characteristics of the minimum controller:');
disp(['Settling Time: ', num2str(info_min.SettlingTime)]);
disp(['Rise Time: ', num2str(info_min.RiseTime)]);
disp(['Peak Time: ', num2str(info_min.PeakTime)]);
disp(['Maximum Overshoot: ', num2str(info_min.Overshoot)]);
disp(['Steady-State Error: ', num2str(abs(1 - y_min(end)))]); % Assuming a unit step input

% Display the step response characteristics of the maximum controller
disp('Step response characteristics of the maximum controller:');
disp(['Settling Time: ', num2str(info_max.SettlingTime)]);
disp(['Rise Time: ', num2str(info_max.RiseTime)]);
disp(['Peak Time: ', num2str(info_max.PeakTime)]);
disp(['Maximum Overshoot: ', num2str(info_max.Overshoot)]);
disp(['Steady-State Error: ', num2str(abs(1 - y_max(end)))]); % Assuming a unit step input

% Display the controller gains for traditional controllers
disp('Controller gains for minimum parameter values:');
disp(K_min);
disp('Controller gains for maximum parameter values:');
disp(K_max);

% Display the selected parameters
disp('Randomly selected parameters:');
disp(['a1 = ', num2str(a1)]);
disp(['a2 = ', num2str(a2)]);
disp(['b2 = ', num2str(b2)]);