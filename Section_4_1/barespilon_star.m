%% barepsilon_star.m
% Computes the constant \bar{epsilon}, an admissible \bar{epsilon}_star,
% and the corresponding dwell-time bound \bar{\tau}_d_star for the
% averaging–predictor controller U2, for the switched system with
% input delay D = 1.
%
% This script implements exactly the quantities appearing in Theorem 2 of the paper:
%
%   - \bar{\epsilon}       
%   - \bar{\epsilon}_\star       
%   - \bar{\tau}_d^\star         
%
% Output:
%   Prints:  bar_epsilon, bar_epsilon_star, bar_tau_d_star

clc; clear all; close all;
addpath('../Auxiliary Functions');
%% System Dynamics  
A1 = [1 1; 1 2];
A2 = [0.97 1.15 ; 1.06 2.09];
A3 = [1.08 1.2 ; 1.14 2.13];

B1 = [0;1];
B2 = [0;1.05];
B3 = [0;1.1];

% Closed-loop gains Ki such that Ai + Bi*Ki are Hurwitz (mode-by-mode)
K1 = computeK(A1,B1,-3,-2);
K2 = computeK(A2,B2,-3,-2);
K3 = computeK(A3,B3,-3,-2);

%% Simulation parameters
D = 1;             % Input delay
dwell_time_min = 0.9;

%% Step 1: ε in equation (14)
barepsilon = max([
    norm(A1-A2), norm(A2-A3), norm(A3-A1), ...
    norm(B1-B2), norm(B2-B3), norm(B3-B1), ...
    norm(K1-K2), norm(K2-K3), norm(K3-K1)
]);

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
S1   = lyap(Acl1', Q1);  % solves Acl'*P + P*Acl + Q = 0 ⇒ Acl'*P + P*Acl = -Q
S2  = lyap(Acl2', Q2);  
S3   = lyap(Acl3', Q3); 

%minimum and maximum eigenvalues
lminS1 = min(eig(S1));
lminS2 = min(eig(S2));
lminS3 = min(eig(S3));

lmaxS1 = max(eig(S1));
lmaxS2 = max(eig(S2));
lmaxS3 = max(eig(S3));

%% Step 3:  norm bounds M_A, M_B, M_K, M_H, equation (37) 
MA = max([norm(A1),norm(A2),norm(A3)]);
MB = max([norm(B1),norm(B2),norm(B3)]);
MK = max([norm(K1),norm(K2),norm(K3)]);
MH = max([norm(A1+B1*K1),norm(A2+B2*K2),norm(A3+B3*K3)]); 

%% Step 4: nu1, equation (41)

nu1 = max([ 4*MK^2*D*exp(2*MH*D) + 1, 4*MK^2*D^2*exp(2*MH*D)*MB^2 + 2 ]);

%% Step 5: compute \hat{lambda}(epsilon), equations (123)--(125)

delta1H = max(1,MB)*(MK*D+1);
delta2H = (2*MK*MB*D + barepsilon + MK + MB);
lambdaH = barepsilon * exp((MA+barepsilon)*D)*max(delta1H,delta2H);

%% Step 6: Solve (126) to find an admissible barepsilon_star
% We approximate the matrix exponentials depending on barepsilon via the third
% order Taylor series

T00=min(lminQ1/norm(B1.*S1),lminQ2/norm(B2.*S2));
T0 = min(T00,lminQ3/norm(B3.*S3));

T1=T0/(2*sqrt(2*exp(D)*(D*nu1+1)));

T3=1/sqrt(2*exp(D)*(D*nu1));

T5 =2*MK*MB*D+MK*MB;

T2 =T1/(exp(MA*D)*max(1,MB)*(MK*D+1));

T4 = T3/(exp(MA*D)*max(1,MB)*(MK*D+1));

digits(50);  % set precision

p1 = sym([D^2/2,D, 1,-T2]);   
r_sym1 = roots(p1);          
r_num1 = vpa(r_sym1, 30) ;     % numeric with 30-digit precision

p2 = sym([D^2/2,D, 1,-T4]);  
r_sym2 = roots(p2);          
r_num2 = vpa(r_sym2, 30);      % numeric with 30-digit precision

p3 = sym([exp(MA*D)*(D^2/2),exp(MA*D)*D+T5*exp(MA*D)*(D^2/2),exp(MA*D)+T5*exp(MA*D)*D,T5*exp(MA*D),-T1]);   
r_sym3 = roots(p3);          
r_num3 = vpa(r_sym3, 30) ;     % numeric with 30-digit precision

p4 = sym([exp(MA*D)*(D^2/2),exp(MA*D)*D+T5*exp(MA*D)*(D^2/2),exp(MA*D)+T5*exp(MA*D)*D,T5*exp(MA*D),-T3]);   
r_sym4 = roots(p4);          
r_num4 = vpa(r_sym4, 30) ;     % numeric with 30-digit precision

% Collect all roots in one vector
all_roots = [r_num1; r_num2; r_num3; r_num4];
all_roots = double(all_roots);

% Tolerance
tol = 1e-8;

% Keep only  real roots
real_roots = all_roots(abs(imag(all_roots)) < tol);
real_roots = real(real_roots);

positive_real_roots = real_roots(real_roots > 0);

% Smallest positive real root
min_positive_real_root = min(positive_real_roots);

%% Step 7: compute \bar{\beta} in equation (128)

epsf = 0.9*min_positive_real_root;
delta1f = max(1,MB)*(MK*D+1);
delta2f = (2*MK*MB*D + epsf + MK + MB);
lambaf = epsf * exp((MA+epsf)*D)*max(delta1f,delta2f);

c1 = 2* norm(B1.*S1)^2 / lminQ1 ;
c2 = 2* norm(B2.*S2)^2 / lminQ2 ;
c3 = 2* norm(B3.*S3)^2 / lminQ3 ;

barb1=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ1-2*c1*exp(D)*lambaf^2*(D*nu1+1));  %equation (128), i=1
barb2=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ2-2*c2*exp(D)*lambaf^2*(D*nu1+1));  %equation (128), i=2
barb3=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ3-2*c3*exp(D)*lambaf^2*(D*nu1+1));  %equation (128), i=3

barbeta = min(min(barb1,barb2),barb3); %equation (128) 

%% Step 8:  Compute mu, k1, k2 and dwell-time condition bartau_d_star


k11 = min( lminS1, c1);
k12 = min( lminS2, c2);
k13 = min( lminS3, c3);

k1 = min(min(k11,k12),k13);  %equation (109)

k21 = max (lmaxS1, c1*exp(D));
k22 = max (lmaxS2, c2*exp(D));
k23 = max (lmaxS3, c3*exp(D));

k2 = max(max(k21,k22),k23); %equation (110)

mu = k2/k1; %equation (108)

bartau_d_star = log(mu)/barbeta; %equation (127)

%% Display key results

fprintf('bar_epsilon      :                             %g\n', barepsilon);
fprintf('bar_epsilon_star :                             %g\n', min_positive_real_root);
fprintf('bar_tau_d_star (theoretical dwell-time bound): %g\n', bartau_d_star);