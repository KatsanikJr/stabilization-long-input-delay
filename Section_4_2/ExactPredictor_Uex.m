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

%% Simulation parameters

% Switching signal σ(t)
% For correct reproducibility, load the pre-generated (by generateSwSignal.m) signal stored in 'swsig_values.mat'. 
% It must contain a vector "swsig_values" with entries in {1,2,3}.

load('swsig_values.mat')

D               = 1;          % Input delay [sec]
dwell_time_min  = 0.9;        % Minimum dwell time [sec]
dwell_time_max  = 3;          % Maximum dwell time [sec]

T  = 10;                      % Total simulation time [sec]
dt = 0.001;                   % Time step [sec]

numSteps   = T/dt;            % Number of integration steps
time       = 0:dt:T;          % Time vector (length = numSteps+1)
delay_steps = round(D/dt);    % Delay in discrete steps

%%  Pre-allocations


X   = zeros(2, numSteps+1);       % State trajectory
P_T = zeros(2, numSteps+1);       % Predictor state 
U   = zeros(1, delay_steps+numSteps+1); % Control input 
% Initial condition
X(:,1) = [1; -1];

% Mode-dependent matrices and exponentials precomputations

Amode = {A1,A2,A3};
Bmode = {B1,B2,B3};
Kmode = {K1,K2,K3};

% Pre-computed exponentials e^{Ai*dt} for each mode
Pwr_A = {expm(A1*dt), expm(A2*dt), expm(A3*dt)};

%% Main simulation loop
for i=1:numSteps %i=t
    
    mode_now = swsig_values(i);
    A = Amode{mode_now};
    B = Bmode{mode_now};
   
    
    % Build the delay window [i, i + delay_steps]
    if (i + delay_steps) <= length(swsig_values)
        window = swsig_values(i:i+delay_steps);       
    else
        window = swsig_values(i:end);
    end

    % Mode used to choose gain K: σ at t_i + D (or last mode if beyond)
    if (i + delay_steps + 1) <= length(swsig_values)
        K = Kmode{swsig_values(i+delay_steps+1)};
    else
        K = Kmode{swsig_values(end)};
    end

    % Find switches in the delay window  and construct (m_n, s_n)
    switches_D = find(diff(window) ~= 0) + i; % location of new mode
    num_switches = length(switches_D);
    
    % mn: sequence of modes on [t_i, t_i + D]
    % sn: corresponding discrete indices (in steps) where those modes start
    mn = mode_now;           % initial mode in the window
    sn = i;                  % starting index
    if num_switches > 0
       for n = 1:num_switches
           mn = [mn, swsig_values(switches_D(n))]; 
           sn = [sn, (switches_D(n))];
       end
    end   
    % Last index in the prediction horizon (i + delay_steps)
    sn = [sn, i + delay_steps];
    
    % Build the exact predictor P(t) = X(t + D), t=i, equation (18) of the paper
    prodX=1;
    integral_term = zeros(2,1);
    
    for n=1:num_switches+1

        % Product of exponentials for the X-part over the nth subinterval
        prodX = Pwr_A{mn(n)}^(sn(n+1)-sn(n)) * prodX;

        % Product of exponentials for the U-part (from future segments)
        prodU=1;
        if n < num_switches+1
            for j=n:num_switches
                prodU = Pwr_A{mn(j+1)}^(sn(j+2)-sn(j+1))*prodU;
            end
         end
        % Discrete "integration" over the nth subinterval for U
        % (Left-point rule, using stored U with delay buffer)
         for k = max(1,-delay_steps+sn(n)):-delay_steps+sn(n+1)
                integral_term = integral_term + prodU * (Pwr_A{mn(n)}^(-delay_steps+sn(n+1) - k) * Bmode{mn(n)} * U(k+delay_steps) )  ;   % Left-point rule integration
         end
    end
    % Exact predictor
    P_t = prodX * X(:,i) + integral_term * dt;
    
    % Control law and plant update
    
    % Control at time t, equation (140) of the paper
    U(i + delay_steps) = K * P_t;
    
    % Plant dynamics: X_{i+1} = X_i + dt ( A X_i + B U_i )
    X(:,i+1) = X(:,i) + dt * ( A * X(:,i) + B * U(i) ); 

    % Store predictor for plotting
    % P_T(:,i) = P_t;
   

end

%% Save data for later use (e.g., performance index)
save('U_ex.mat', 'X', 'U');

% Plot state X(t) 
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(time, X(1, :), 'r', 'LineWidth',3);
hold on
plot(time, X(2, :), 'b', 'LineWidth',3);
%yticks([-6 -4 -2 0 2 4 6]);
%ylim([-6 6]);
grid on;

ax = gca;
ax.FontSize = 20;
ylabel('$X(t)$', 'Interpreter', 'latex','FontSize', 25);
xlabel('$t$ (sec)', 'Interpreter', 'latex','FontSize', 25);
legend('$X_1(t)$','$X_2(t)$', 'Interpreter', 'latex', 'FontSize', 17); 

% Plot control input U(t)
subplot(2,1,2)
plot(time, U(delay_steps:delay_steps+numSteps), 'k', 'LineWidth',3);
%yticks([-20 -10 0 20 40 60]);
%ylim([-20 60]);
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
xlabel('$t$ (sec)','Interpreter', 'latex','FontSize',25);
ylabel('$\sigma(t)$','Interpreter', 'latex','FontSize', 25);

hold on;
for i = 1:length(time)
    if swsig_values(i) ~= swsig_values(i+1)
        plot(time(i:i+1), swsig_values(i:i+1), 'k-', 'LineWidth', 1);
    end
end
yticks([1 2 3]);
ylim([0.8 3.2]);
ax = gca;
ax.FontSize = 20; 
% Improve plot clarity
alpha(0.2); % Set marker transparency
hold off