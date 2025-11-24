clear all 
clc
%initializations
D = 1;  %  Delay
dwell_time = 0.9;

T = 10; % Total simulation time
dt = 0.001; % Time step

numSteps = T/dt; %Total steps of simulation
time = 0:dt:T; 
delay_steps = round(D/dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('U1_dwell.mat');
X1_d = X;
U1_d =  U(delay_steps:delay_steps+numSteps);

% Compute |X(t)|^2 = |X1|^2 + |X2|^2
X1d_energy = sum(abs(X1_d).^2, 1);   % 1 x N vector

E1d_X = trapz(time, X1d_energy);   % ∫ |X|² dt
E1d_U = trapz(time, U1_d.^2);       % ∫ U² dt
E1d_total = E1d_X + E1d_U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('U1_nodwell.mat');
X1_nd = X;
U1_nd =  U(delay_steps:delay_steps+numSteps);

% Compute |X(t)|^2 = |X1|^2 + |X2|^2
X1nd_energy = sum(abs(X1_nd).^2, 1);   % 1 x N vector

E1nd_X = trapz(time, X1nd_energy);   % ∫ |X|² dt
E1nd_U = trapz(time, U1_nd.^2);       % ∫ U² dt
E1nd_total = E1nd_X + E1nd_U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('U2_dwell.mat');
X2_d = X;
U2_d =  U(delay_steps:delay_steps+numSteps);

% Compute |X(t)|^2 = |X1|^2 + |X2|^2
X2d_energy = sum(abs(X2_d).^2, 1);   % 1 x N vector

E2d_X = trapz(time, X2d_energy);   % ∫ |X|² dt
E2d_U = trapz(time, U2_d.^2);       % ∫ U² dt
E2d_total = E2d_X + E2d_U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('U2_nodwell.mat');
X2_nd = X;
U2_nd =  U(delay_steps:delay_steps+numSteps);

% Compute |X(t)|^2 = |X1|^2 + |X2|^2
X2nd_energy = sum(abs(X2_nd).^2, 1);   % 1 x N vector

E2nd_X = trapz(time, X2nd_energy);   % ∫ |X|² dt
E2nd_U = trapz(time, U2_nd.^2);       % ∫ U² dt
E2nd_total = E2nd_X + E2nd_U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('U_ex.mat');
Xex = X;
Uex =  U(delay_steps:delay_steps+numSteps);

% Compute |X(t)|^2 = |X1|^2 + |X2|^2
Xex_energy = sum(abs(Xex).^2, 1);   % 1 x N vector

Eex_X = trapz(time, Xex_energy);   % ∫ |X|² dt
Eex_U = trapz(time, Uex.^2);       % ∫ U² dt
Eex_total = Eex_X + Eex_U;






