clc
clear all
close all


% Setup of the system 
A1 = [1 1; 1 2];
A2 = [0.97 1.15 ; 1.06 2.09];
A3 = [1.08 1.2 ; 1.14 2.13];

B1 = [0;1];
B2 = [0;1.05];
B3 = [0;1.1];

K1 = computeK(A1,B1,-3,-2); %closed loop gain A1+B1*K1
K2 = computeK(A2,B2,-3,-2); %closed loop gain A2+B2*K2
K3 = computeK(A3,B3,-3,-2); %closed loop gain A2+B2*K2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initializations for simulation
D = 1;  %  delay
dwell_time = 0.9; 
dwell_time_upper = 5;

dt = 0.001; % Time step
T = 10; % Total simulation time
numSteps = T/dt; %Total steps of simulation
time = 0:dt:T; 
delay_steps = round(D/dt);

Dcontroller = 1;
delay_steps_controller = round(Dcontroller/dt);

%kappa_values = generateSwSignal(dwell_time,dwell_time_upper,T,dt);  %random switching signal 
load('swsig_values.mat')

X = zeros(2, numSteps+1); % Initialize state matrix
Xt = zeros(2, numSteps+1); % Initialize state matrix
U = zeros(1, delay_steps+numSteps+1+delay_steps_controller); % Initialize control input matrix
U1 = zeros(1, delay_steps+numSteps+1+delay_steps_controller); % Initialize control input matrix
U2 = zeros(1, delay_steps+numSteps+1+delay_steps_controller); % Initialize control input matrix
U3 = zeros(1, delay_steps+numSteps+1+delay_steps_controller); % Initialize control input matrix
% Assume X0 is given
X(:,1) = [1; -1]; % Initial state

dwell_steps = round(dwell_time / dt); % Number of steps to dwell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amode = {A1,A2,A3};
Bmode = {B1,B2,B3};
Kmode = {K1,K2,K3};

Pwr_A = {expm(A1*dt),expm(A2*dt),expm(A3*dt)};


switching = swsig_values(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:numSteps
   
    kappa = swsig_values(i);

    A = Amode{kappa};
    B = Bmode{kappa};


    
    %Find last switching instant t_0(t)
    if i>2 && swsig_values(i) ~= swsig_values(i-1)
        switching = i;
    end
   
    %Determine tau_(t)
    tr = max(0,switching + dwell_steps - i ); % with repsect to simulation steps index 
    %tr=0;
    tr_t = tr*dt; % with respect to real seconds

   % Compute integral for equation (7)
    integral_term_t = zeros(2,1);
    for k = max(1,i-delay_steps_controller):(i+tr-delay_steps_controller)
             integral_term_t = integral_term_t + Pwr_A{swsig_values(switching)}^(i+tr-delay_steps_controller - k) * Bmode{swsig_values(switching)} * U(k+delay_steps_controller)  ;   % Left-point rule integration
    end

    % Compute integral for equation (4)
    integral_term1 = zeros(2,1);
    integral_term2 = zeros(2,1);
    integral_term3 = zeros(2,1);
    for j =max(1,i+tr-delay_steps_controller):i
             integral_term1 = integral_term1 +  Pwr_A{1}^(i-j) * B1 * U(j+delay_steps_controller)  ;   % Left-point rule integration
             integral_term2 = integral_term2 +  Pwr_A{2}^(i-j) * B2 * U(j+delay_steps_controller)  ;   % Left-point rule integration
             integral_term3 = integral_term3 +  Pwr_A{3}^(i-j) * B3 * U(j+delay_steps_controller)  ;
    end



    X(:, i+1) = X(:, i) + dt * ( A * X(:, i) + B * U(i) ); % Update state / Implementation of (1)
    Xt(:, i) = Pwr_A{swsig_values(switching)}^(tr_t/dt) * X(:, i) + integral_term_t*dt; % Implementation of (7)
    U1(i+1+delay_steps) = K1 * (Pwr_A{1}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term1*dt );  % Update control
    U2(i+1+delay_steps) = K2 * (Pwr_A{2}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term2*dt );
    U3(i+1+delay_steps) = K3 * (Pwr_A{3}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term3*dt );
    U(i+1+delay_steps) = 1/3*(U1(i+1+delay_steps)+U2(i+1+delay_steps)+U3(i+1+delay_steps));  % Update control
    %U(i+1+delay_steps) = K_bar *1/3* ( (Pwr_A{1}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term1*dt )+(Pwr_A{1}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term1*dt )+(Pwr_A{1}^((Dcontroller-tr_t)/dt)* Xt(:, i) + integral_term1*dt ));  % Update control



end
save('U2_nodwell.mat', 'X', 'U');

% PLOTS 
% Plotting the state X
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(time, X(1, :), 'r', 'LineWidth',3);
hold on
plot(time, X(2, :), 'b', 'LineWidth',3);
yticks([-6 -3 0 3 6]);
ylim([-6 6]);
grid on;

ax = gca;
ax.FontSize = 20;
ylabel('$X(t)$', 'Interpreter', 'latex','FontSize', 25);
xlabel('$t$ (sec)', 'Interpreter', 'latex','FontSize', 25);
legend('$X_1(t)$','$X_2(t)$', 'Interpreter', 'latex', 'FontSize', 17); 
% Plotting the control input U
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

%

% Plot the values of kappa with respect to time
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

% Customize the y-axis to only show 1 and 2
yticks([1 2 3]);
ylim([0.5 3.5]);
ax = gca;
ax.FontSize = 20; 
% Improve plot clarity
alpha(0.2); % Set marker transparency
hold off;