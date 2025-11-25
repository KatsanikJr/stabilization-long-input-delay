%% Performance_Index.m
% Computes the performance index
%
%       J = ∫_0^T ( |X(t)|^2 + |U(t)|^2 ) dt, equation (141) in paper
%
% for the following controllers:
%
%   - U1_dwell     (average predictor with dwell-time knowledge), equations (3)-(4) under (5)--(7) in paper
%   - U1_nodwell   (average predictor without dwell-time knowledge), equations (3)-(4) with τ(t)=0, for all times t
%   - U2_dwell     (averaging predictors with dwell-time knowledge), equations (10)-(12) under (5)--(7) in paper
%   - U2_nodwell   (averaging predictor without dwell-time knowledge), equations (10)-(12) with τ(t)=0, for all times t
%   - U_ex         (exact predictor, non-implementable), equation (140) under (18) in paper
%
% Output:
%   - Prints performance measure J for all controllers defined in equation (141) of the paper
%   - Organizes results in a MATLAB table
%
% Requires that the .mat files are generated beforehand by running:
%   U1_dwell.m, U1_nodwell.m, U2_dwell.m, U2_nodwell.m, ExactPredictor.m

clear all;
clc;
addpath('../Auxiliary Functions');
addpath('../Section_4_1');


%% PARAMETERS
D  = 1;           % Input delay
dt = 0.001;       % Time step
T  = 10;          % Total simulation time

numSteps   = T/dt;
time       = 0:dt:T;
delay_steps = round(D/dt);

%% ========================================================================
%  function to compute the performance index J
% ========================================================================

compute_J = @(X, U, time) trapz(time, sum(abs(X).^2,1)) + trapz(time, U.^2);

%% ========================================================================
%  Load data for each controller & compute performance index
% ========================================================================

%% --- U1_dwell ------------------------------------------------------------
load('U1_dwell.mat');        % loads X, U, time
X1_d = X;
U1_d = U(delay_steps:delay_steps+numSteps);
J_U1_dwell = compute_J(X1_d, U1_d, time);

%% --- U1_nodwell ----------------------------------------------------------
load('U1_nodwell.mat');
X1_nd = X;
U1_nd = U(delay_steps:delay_steps+numSteps);
J_U1_nodwell = compute_J(X1_nd, U1_nd, time);

%% --- U2_dwell ------------------------------------------------------------
load('U2_dwell.mat');
X2_d = X;
U2_d = U(delay_steps:delay_steps+numSteps);
J_U2_dwell = compute_J(X2_d, U2_d, time);

%% --- U2_nodwell ----------------------------------------------------------
load('U2_nodwell.mat');
X2_nd = X;
U2_nd = U(delay_steps:delay_steps+numSteps);
J_U2_nodwell = compute_J(X2_nd, U2_nd, time);

%% --- Exact predictor (non-implementable) --------------------------------
load('U_ex.mat');
Xex = X;
Uex = U(delay_steps:delay_steps+numSteps);
J_U_ex = compute_J(Xex, Uex, time);

%% ========================================================================
%  DISPLAY RESULTS
% ========================================================================

fprintf('\n=================================================\n');
fprintf('        PERFORMANCE INDEX RESULTS (J)\n');
fprintf('=================================================\n');
fprintf(' U1 with dwell-time knowledge      : %g\n', J_U1_dwell);
fprintf(' U1 without dwell-time knowledge   : %g\n', J_U1_nodwell);
fprintf(' U2 with dwell-time knowledge      : %g\n', J_U2_dwell);
fprintf(' U2 without dwell-time knowledge   : %g\n', J_U2_nodwell);
fprintf(' Exact predictor (non-implement.)  : %g\n', J_U_ex);
fprintf('=================================================\n\n');

%% Table output (optional)
Controller = ["U1-dwell"; "U1-no-dwell"; "U2-dwell"; "U2-no-dwell"; "Exact predictor"];
J_value    = [J_U1_dwell; J_U1_nodwell; J_U2_dwell; J_U2_nodwell; J_U_ex];

Results = table(Controller, J_value)



