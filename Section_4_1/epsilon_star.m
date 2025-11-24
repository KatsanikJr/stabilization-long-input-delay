%% epsilon_star.m
% Computes the constant epsilon, an admissible epsilon_star, and the
% corresponding dwell-time bound tau_d_star for the average predictor-based
% controller U1, for the switched system with input delay D = 1.
%
% This script computes exactly the quantities appearing in Theorem 1 of the paper:
%
%   - \epsilon         : equation (14)
%   - \epsilon_\star   : equation (75)    
%   - \tau_d^\star     : equation (76)     
%
% Output:
%   Prints:  epsilon, epsilon_star, tau_d_star


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
% These correspond to the nominal matrices (A_bar, B_bar, K_bar) in the paper.
A_bar = optimizeL2(A1, A2, A3); % \bar{A}
B_bar = optimizeL2(B1, B2, B3); % \bar{B}
K_bar = optimizeL2(K1, K2, K3); % \bar{K}

%% Simulation Parameters

D              = 1;    % Actual input delay [sec]
dwell_time_min = 0.9;  % Minimum dwell time used in simulations [sec]

%% Step 1: ε in equation (14)

eps = max([
    norm(A1 - A_bar, 2), ...
    norm(A2 - A_bar, 2), ...
    norm(A3 - A_bar, 2), ...
    norm(B1 - B_bar, 2), ...
    norm(B2 - B_bar, 2), ...
    norm(B3 - B_bar, 2), ...
    norm(K1 - K_bar, 2), ...
    norm(K2 - K_bar, 2), ...
    norm(K3 - K_bar, 2)]);

%% Step 2: Lyapunov matrices P_i and Q_i
% These come from the Lyapunov equation (13)

Acl1 = A1 + B1*K1; % closed loop dynamics A1+B1*K1       
Acl2 = A2 + B2*K2; % closed loop dynamics A2+B2*K2         
Acl3 = A3 + B3*K3; % closed loop dynamics A3+B3*K3  

Q1 = [1 0 ; 0 1]; % random choice
Q2 = [3 0 ; 0 3]; % random choice
Q3 = [2 0 ; 0 2]; % random choice

lminQ1 = min(eig(Q1)); 
lminQ2 = min(eig(Q2));
lminQ3 = min(eig(Q3));

% Si matrices defined above (13)
S1 = lyap(Acl1', Q1); 
S2 = lyap(Acl2', Q2);
S3 = lyap(Acl3', Q3);

%minimum and maximum eigenvalues
lminS1 = min(eig(S1));
lminS2 = min(eig(S2));
lminS3 = min(eig(S3));

lmaxS1 = max(eig(S1));
lmaxS2 = max(eig(S2));
lmaxS3 = max(eig(S3));

%% Step 3:  norm bounds M_A, M_B, M_K, M_H, equation (37) and \bar{M_A}, \bar{M_B}, equation (57)

MA = max([norm(A1), norm(A2), norm(A3)]);
MB = max([norm(B1), norm(B2), norm(B3)]);
MK = max([norm(K1), norm(K2), norm(K3)]);
MH = max([norm(A1 + B1*K1), norm(A2 + B2*K2), norm(A3 + B3*K3)]); 

MAb = max(MA, norm(A_bar,2));
MBb = max(MB, norm(B_bar,2));

%% Step 4: nu1, equation (41)

nu1 = max([ 4*MK^2*D*exp(2*MH*D) + 1, 4*MK^2*D^2*exp(2*MH*D)*MB^2 + 2 ]);

%% Step 5: compute lambda(epsilon), equations (72)--(74)

delta1 = max(1,MBb) * (MK*D + 1);
delta2 = (2*MK*MBb*D + eps + MK + MBb);
lambda = eps * exp((MAb + eps)*D) * max(delta1, delta2);

%% Step 6: Solve (75) to find an admissible epsilon_star
% We approximate the matrix exponentials depending on epsilon via the third
% order Taylor series

T00 = min( lminQ1 / norm(B1.*S1), lminQ2 / norm(B2.*S2) );
T0  = min( T00, lminQ3 / norm(B3.*S3) ); 

T1 = T0 /(2*sqrt(2*exp(D)*(D*nu1 + 1)));

T3 = 1 / sqrt(2*exp(D)*(D*nu1));

T5 = 2*MK*MBb*D + MK*MBb;

T2 = T1 /(exp(MAb*D) * max(1,MBb) * (MK*D + 1));

T4 = T3 /(exp(MAb*D) * max(1,MBb) * (MK*D + 1));


digits(50);  % set precision 

% First inequality 
p1     = sym([D^2/2, D, 1, -T2]);
r_sym1 = roots(p1);
r_num1 = vpa(r_sym1, 30);     % numeric with 30-digit precision

% Second inequality 
p2     = sym([D^2/2, D, 1, -T4]);
r_sym2 = roots(p2);
r_num2 = vpa(r_sym2, 30);

% Third inequality 
p3     = sym([exp(MAb*D)*(D^2/2), ...
              exp(MAb*D)*D + T5*exp(MAb*D)*(D^2/2), ...
              exp(MAb*D) + T5*exp(MAb*D)*D, ...
              T5*exp(MAb*D), -T1]);
r_sym3 = roots(p3);
r_num3 = vpa(r_sym3, 30);

% Fourth inequality 
p4     = sym([exp(MAb*D)*(D^2/2), ...
              exp(MAb*D)*D + T5*exp(MAb*D)*(D^2/2), ...
              exp(MAb*D) + T5*exp(MAb*D)*D, ...
              T5*exp(MAb*D), -T3]);
r_sym4 = roots(p4);
r_num4 = vpa(r_sym4, 30);

% Collect all roots in one vector
all_roots = [r_num1; r_num2; r_num3; r_num4];
all_roots = double(all_roots);

% Tolerance
tol = 1e-8;

% Keep only real roots
real_roots = all_roots(abs(imag(all_roots)) < tol);
real_roots = real(real_roots);

positive_real_roots = real_roots(real_roots > 0);

% Smallest positive real root
min_positive_real_root = min(positive_real_roots);

%% Step 7: compute \beta in equation (101)

epsf = 0.9 * min_positive_real_root;

delta1f = max(1,MBb) * (MK*D + 1);
delta2f = (2*MK*MBb*D + epsf + MK + MBb);
lambaf  = epsf * exp((MAb + epsf)*D) * max(delta1f, delta2f);

c1 = 2 * norm(B1.*S1)^2 / lminQ1;
c2 = 2 * norm(B2.*S2)^2 / lminQ2;
c3 = 2 * norm(B3.*S3)^2 / lminQ3;

b1 = min( 1 - 2*exp(D)*lambaf^2*D*nu1, ...
          0.5*lminQ1 - 2*c1*exp(D)*lambaf^2*(D*nu1 + 1) ); %equation (93), i=1
b2 = min( 1 - 2*exp(D)*lambaf^2*D*nu1, ... 
          0.5*lminQ2 - 2*c2*exp(D)*lambaf^2*(D*nu1 + 1) ); %equation (94), i=2
b3 = min( 1 - 2*exp(D)*lambaf^2*D*nu1, ...
          0.5*lminQ3 - 2*c3*exp(D)*lambaf^2*(D*nu1 + 1) ); %equation (95), i=3

beta = min([b1, b2, b3]); %equation (101) 

%% Step 8:  Compute mu, k1, k2 and dwell-time condition tau_d_star

k11 = min( lminS1, c1);
k12 = min( lminS2, c2);
k13 = min( lminS3, c3);

k1 = min([k11, k12, k13]); %equation (109)

k21 = max( lmaxS1, c1*exp(D));
k22 = max( lmaxS2, c2*exp(D));
k23 = max( lmaxS3, c3*exp(D));

k2 = max([k21, k22, k23]); %equation (110)

mu = k2 / k1; %equation (108)

tau_d_star = log(mu) / beta; %equation (76)

%% Display key results

fprintf('epsilon      :                             %g\n', eps);
fprintf('epsilon_star :                             %g\n', min_positive_real_root);

fprintf('tau_d_star (theoretical dwell-time bound): %g\n', tau_d_star);
