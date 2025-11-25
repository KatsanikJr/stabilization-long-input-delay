%% U1_dwell_time.m
% Average predictor-based controller U1 with dwell-time knowledge.
% Implements (1) and (3)--(9) of the paper for the switched system (1)
% with input delay D = 1 and minimum dwell time tau_d = 0.9 s.
%
% Output:
%   - Saves U1_dwell.mat containing X, U for later use.
%   - Plots X(t), U(t) and the switching signal sigma(t).

clc; clear all; close all;
addpath('../Auxiliary Functions');
%% System Dynamics 

A1 = [1 1; 1 2];
A2 = [0.97 1.15 ; 1.06 2.09];
A3 = [1.08 1.2 ; 1.14 2.13];

B1 = [0;1];
B2 = [0;1.05];
B3 = [0;1.1];

% Closed-loop gains Ki such that Ai + Bi*Ki is Hurwitz
K1 = computeK(A1,B1,-3,-2);   % closed-loop matrix: A1 + B1*K1
K2 = computeK(A2,B2,-3,-2);   % closed-loop matrix: A2 + B2*K2
K3 = computeK(A3,B3,-3,-2);   % closed-loop matrix: A3 + B3*K3

% Average system dynamics (via L2 optimization)
A_bar = optimizeL2(A1, A2, A3); %call optimizeL2.m to compute \bar{A} under (8)
B_bar = optimizeL2(B1, B2, B3); %call optimizeL2.m to compute \bar{B} under (8)
K_bar = optimizeL2(K1, K2, K3); %call optimizeL2.m to compute \bar{K} under (9)

%% Simulation parameters

% Switching signal σ(t)
% For exact reproducibility, load the the specific signal used in the paper. 
% The file swsig_values.mat was generated once using generateSwSignal.m
% It must contain a vector "swsig_values" with entries in {1,2,3}.

load('swsig_values.mat')

D               = 1;          % Actual input delay [sec]
dwell_time_min  = 0.9;        % Minimum dwell time [sec]
dwell_time_max  = 3;          % Maximum dwell time [sec]

T  = 10;                      % Total simulation time [sec]
dt = 0.001;                   % Time step [sec]

numSteps   = T/dt;            % Number of integration steps
time       = 0:dt:T;          % Time vector (length = numSteps+1)
delay_steps = round(D/dt);    % Delay in discrete steps

Dcontroller = D;             % Input Delay that controller takes into account (in this case it is the same)
delay_steps_controller = round(Dcontroller/dt);

dwell_steps = round(dwell_time_min / dt);  % dwell time in discrete steps

%%  Pre-allocations

X  = zeros(2, numSteps+1);                  % Plant state
Xt = zeros(2, numSteps+1);                  % state of equation (7) 
U  = zeros(1, delay_steps+numSteps+1+delay_steps_controller); % Control input

% Initial condition
X(:,1) = [1; -1];

% Mode-dependent matrices and exponentials precomputations

Amode = {A1,A2,A3};
Bmode = {B1,B2,B3};
Kmode = {K1,K2,K3};

Pwr_A    = {expm(A1*dt), expm(A2*dt), expm(A3*dt)};
Pwr_Abar = expm(A_bar*dt);

%Inialize index τ_0(t)
switching = swsig_values(1);

%% Main simulation loop

for i=1:numSteps
    mode_now = swsig_values(i);

    A = Amode{mode_now};
    B = Bmode{mode_now};
    
    % ---------- Equation (6): Find last switching instant τ_0(t)  ----------
    if i>2 && swsig_values(i) ~= swsig_values(i-1)
        switching = i;
    end
   
    % ---------- Equation (5): Compute τ(t) ----------
    tr = max(0,switching + dwell_steps - i ); % τ(t)/dt
    tr_t = tr*dt; % with respect to real seconds

    
    % ---------- Equation (7): compute state Xt=X(t+τ(t))  ----------
    integral_term_t = zeros(2,1);
    for k = max(1,i-delay_steps_controller):(i+tr-delay_steps_controller)
             integral_term_t = integral_term_t + Pwr_A{swsig_values(switching)}^(i+tr-delay_steps_controller - k) * Bmode{swsig_values(switching)} * U(k+delay_steps_controller)  ;   % Left-point rule integration
    end
    Xt(:, i) = Pwr_A{swsig_values(switching)}^(tr_t/dt) * X(:, i) + integral_term_t*dt; 
    
     % ---------- Equation (4): compute using average system dynamics ----------
    integral_term = zeros(2,1);
    for j = max(1,i+tr-delay_steps_controller):i
             integral_term = integral_term + Pwr_Abar^(i-j) * B_bar * U(j+delay_steps_controller)  ;   % Left-point rule integration
    end
     % Control input (4)
     U(i+delay_steps) = K_bar *(Pwr_Abar^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term*dt );  

     % ---------- Plant dynamics (1) ----------
     X(:, i+1) = X(:, i) + dt * ( A * X(:, i) + B * U(i) ); % Update state 

end

%% Save data
save('U1_dwell.mat', 'X', 'U');

%% PLOTS 
% Plot state X(t)
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(time, X(1, :), 'r', 'LineWidth',3);
hold on
plot(time, X(2, :), 'b', 'LineWidth',3);
yticks([-6 -3  0  3 6]);
ylim([-7 7]);
grid on;

ax = gca;
ax.FontSize = 20;
ylabel('$X(t)$', 'Interpreter', 'latex','FontSize', 25);
xlabel('$t$ (sec)', 'Interpreter', 'latex','FontSize', 25);
legend('$X_1(t)$','$X_2(t)$', 'Interpreter', 'latex', 'FontSize', 17);

% Plot control input U(t)
subplot(2,1,2)
plot(time, U(delay_steps:delay_steps+numSteps), 'k', 'LineWidth',3);
yticks([-40 -20 0 20 40]);
ylim([-50 50]);
grid on;

ax = gca;
ax.FontSize = 20; 
xlabel('$t$ (sec)', 'Interpreter', 'latex','FontSize', 25);
ylabel('$U(t)$', 'Interpreter', 'latex','FontSize', 25);
hold off


% Plot the switching signal σ(t)
figure('Position', [100, 100, 800, 600]);
plot(time, swsig_values(1:length(time)), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 5);
grid on;

ax = gca;
ax.FontSize = 20; 
xlabel('$t$ (sec)', 'Interpreter', 'latex','FontSize', 25);
ylabel('$\sigma(t)$', 'Interpreter', 'latex','FontSize', 25);


hold on;

% Connect changes between 1 and 2 with lines
for i = 1:length(time)
    if swsig_values(i) ~= swsig_values(i+1)
        plot(time(i:i+1), swsig_values(i:i+1), 'k-', 'LineWidth', 1);
    end
end

yticks([1 2 3]);
ylim([0.5 3.5]);
ax = gca;
ax.FontSize = 20; 
% Improve plot clarity
alpha(0.2); % Set marker transparency

hold off;
