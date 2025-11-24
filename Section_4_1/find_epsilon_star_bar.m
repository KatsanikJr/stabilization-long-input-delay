clc; clear all; close all

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

D = 1;  %  delay
dwell_time = 0.9; 
dwell_time_upper = 3;

%find epsilon, epsilon* and td*
eps = max([norm(A1-A2),norm(A2-A3),norm(A3-A1),norm(B1-B2),norm(B2-B3),norm(B3-B1),norm(K1-K2),norm(K2-K3),norm(K3-K1)]);


Acl1 = A1 + B1*K1;        
Acl2 = A2 + B2*K2;        
Acl3 = A3 + B3*K3;

Q1 = [1 0 ; 0 1];
Q2 = [3 0 ; 0 3];
Q3 = [2 0 ; 0 2];

lminQ1 = min(eig(Q1));
lminQ2 = min(eig(Q2));
lminQ3 = min(eig(Q3));

P1   = lyap(Acl1', Q1);  % solves Acl'*P + P*Acl + Q = 0 ⇒ Acl'*P + P*Acl = -Q
P2  = lyap(Acl2', Q2);  
P3   = lyap(Acl3', Q3); 

lminP1 = min(eig(P1));
lminP2 = min(eig(P2));
lminP3 = min(eig(P3));

lmaxP1 = max(eig(P1));
lmaxP2 = max(eig(P2));
lmaxP3 = max(eig(P3));

MA = max([norm(A1),norm(A2),norm(A3)]);
MB = max([norm(B1),norm(B2),norm(B3)]);
MK = max([norm(K1),norm(K2),norm(K3)]);
MH = max([norm(A1+B1*K1),norm(A2+B2*K2),norm(A3+B3*K3)]); 

nu1 = max([ 4*MK^2*D*exp(2*MH*D)+1, 4*MK^2*D^2*exp(2*MH*D)*MB^2+2 ]);


delta1 = max(1,MB)*(MK*D+1);
delta2 = (2*MK*MB*D + eps + MK + MB);
lamba = eps * exp((MA+eps)*D)*max(delta1,delta2);

T00=min(lminQ1/norm(B1.*P1),lminQ2/norm(B2.*P2));
T0 = min(T00,lminQ3/norm(B3.*P3));

T1=T0/(2*sqrt(2*exp(D)*(D*nu1+1)));

T3=1/sqrt(2*exp(D)*(D*nu1));

T5 =2*MK*MB*D+MK*MB;

T2 =T1/(exp(MA*D)*max(1,MB)*(MK*D+1));

T4 = T3/(exp(MA*D)*max(1,MB)*(MK*D+1));

digits(50);  % set precision (50 significant digits)

p1 = sym([D^2/2,D, 1,-T2]);   % coefficients as symbolic
r_sym1 = roots(p1);           % symbolic roots
r_num1 = vpa(r_sym1, 30)      % numeric with 30-digit precision

p2 = sym([D^2/2,D, 1,-T4]);   % coefficients as symbolic
r_sym2 = roots(p2);           % symbolic roots
r_num2 = vpa(r_sym2, 30)      % numeric with 30-digit precision

p3 = sym([exp(MA*D)*(D^2/2),exp(MA*D)*D+T5*exp(MA*D)*(D^2/2),exp(MA*D)+T5*exp(MA*D)*D,T5*exp(MA*D),-T1]);   % coefficients as symbolic
r_sym3 = roots(p3);           % symbolic roots
r_num3 = vpa(r_sym3, 30)      % numeric with 30-digit precision

p4 = sym([exp(MA*D)*(D^2/2),exp(MA*D)*D+T5*exp(MA*D)*(D^2/2),exp(MA*D)+T5*exp(MA*D)*D,T5*exp(MA*D),-T3]);   % coefficients as symbolic
r_sym4 = roots(p4);           % symbolic roots
r_num4 = vpa(r_sym4, 30)      % numeric with 30-digit precision

% Collect all roots in one vector
all_roots = [r_num1; r_num2; r_num3; r_num4];

% Convert to double for easier numeric processing
all_roots = double(all_roots);

% Tolerance for deciding "real"
tol = 1e-8;

% Keep only (approximately) real roots
real_roots = all_roots(abs(imag(all_roots)) < tol);

% Take their real part
real_roots = real(real_roots);

% Keep only positive real roots
positive_real_roots = real_roots(real_roots > 0);

% Smallest positive real root
min_positive_real_root = min(positive_real_roots)

epsf = 0.9*min_positive_real_root;
delta1f = max(1,MB)*(MK*D+1);
delta2f = (2*MK*MB*D + epsf + MK + MB);
lambaf = epsf * exp((MA+epsf)*D)*max(delta1f,delta2f);

c1 = 2* norm(B1.*P1)^2 / lminQ1 ;
c2 = 2* norm(B2.*P2)^2 / lminQ2 ;
c3 = 2* norm(B3.*P3)^2 / lminQ3 ;

b1=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ1-2*c1*exp(D)*lambaf^2*(D*nu1+1));
b2=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ2-2*c2*exp(D)*lambaf^2*(D*nu1+1));
b3=min(1-2*exp(D)*lambaf^2*D*nu1, 0.5*lminQ3-2*c3*exp(D)*lambaf^2*(D*nu1+1));

beta = min(min(b1,b2),b3);

k11 = min( lminP1, c1);
k12 = min( lminP2, c2);
k13 = min( lminP3, c3);

k1 = min(min(k11,k12),k13);

k21 = max (lmaxP1, c1*exp(D));
k22 = max (lmaxP2, c2*exp(D));
k23 = max (lmaxP3, c3*exp(D));

k2 = max(max(k21,k22),k23);

mu = k2/k1;

tau_d_star_bar = log(mu)/beta;
