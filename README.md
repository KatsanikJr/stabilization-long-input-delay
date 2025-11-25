Stabilization of Linear Switched Systems with Long Constant Input Delay via Predictor Feedbacks

The purpose of these MATLAB codes is to reproduce the simulation results presented in Section 4 of the paper:

"Stabilization of Linear Switched Systems with Long Constant Input Delay via Average or Averaging Predictor Feedbacks."

This repository includes implementations of:

the Average Predictor controller U1

the Averaging Predictors controller U2

the non-causal exact predictor (used only as a reference)

robustness simulations with delay perturbations

scripts that compute the theoretical constants epsilon, epsilon_bar, epsilon_star, epsilon_bar_star, tau_d_star, tau_d_bar_star

A short mathematical summary and documentation are provided in the file "codesdoc.pdf".

Requirements

MATLAB

CVX toolbox (required by optimizeL2.m)

An SDP solver compatible with CVX (for example: SeDuMi, SDPT3, or MOSEK)

Usage

Section 4.1: Controllers with dwell-time knowledge
Run these scripts inside the folder Section_4.1:

U1_dwell
U2_dwell
epsilon_star
barepsilon_star

These scripts simulate the system under U1 and U2, produce plots of X(t), U(t), and the switching signal, and compute the theoretical constants in Theorem 1 and Theorem 2.

Section 4.2: No dwell-time information and exact predictor
Run these scripts inside the folder Section_4.2:

U1_nodwell
U2_nodwell
ExactPredictor
Performance_Index

These scripts simulate the controllers without dwell-time information, simulate the exact predictor, and compute the quadratic performance index
J = integral from 0 to T of (|X(t)|^2 + |U(t)|^2) dt.

Section 4.3: Robustness to delay perturbations
Run these scripts inside the folder Section_4.3:

U1_robust_0_95
U1_robust_1_05
U2_robust_0_95
U2_robust_1_05

These simulate the system when the actual input delay differs from the nominal delay assumed by the controller.

Repository Structure

Auxiliary Functions/
computeK.m
optimizeL2.m
generateSwSignal.m
swsig_values.mat

Section_4.1/
U1_dwell.m
U2_dwell.m
epsilon_star.m
barepsilon_star.m

Section_4.2/
U1_nodwell.m
U2_nodwell.m
ExactPredictor.m
Performance_Index.m

Section_4.3/
U1_robust_0_95.m
U1_robust_1_05.m
U2_robust_0_95.m
U2_robust_1_05.m

README.pdf
LICENSE

License

Copyright 2025.
See LICENSE for detailed license information.

Acknowledgements

This work was supported by the listed funding agency and grant.
The views expressed are those of the authors and do not necessarily reflect the views of the funding body.

How to Cite

Journal version:

@article{Your2025Article,
author = {...},
journal = {Automatica},
title = {Stabilization of Linear Switched Systems with Long Constant Input Delay via Average or Averaging Predictor Feedbacks},
year = {2025}
}

arXiv preprint:

@unpublished{Your2025Arxiv,
author = {...},
note = {arXiv:xxxx.xxxxx, 2025},
title = {Stabilization of Linear Switched Systems with Long Constant Input Delay via Average or Averaging Predictor Feedbacks}
}
