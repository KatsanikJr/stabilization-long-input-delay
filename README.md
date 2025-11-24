# stabilization-long-input-delay
MATLAB code for Section 4 simulations

## Repository structure (Section 4)

- `Section_4.1/`  
  MATLAB scripts for the controller \(U_1\), the controller \(U_2\), and the
  computation of the theoretical bounds \(\varepsilon^\star\) and \(\bar\varepsilon^\star\).

- `Section_4.2/`  
  MATLAB scripts for:
  - the (inapplicable) exact predictor-based controller \(U_{\mathrm{ex}}\),
  - the implementations of \(U_1\) and \(U_2\) **without** dwell-time knowledge,
  - the performance index script that computes
    \[
      J = \int_0^T (\|X(t)\|^2 + |U(t)|^2)\,dt
    \]
    for all controllers.

- `Section_4.3/`  
  MATLAB scripts for the simulations with delay perturbations (to be added).

- `Auxiliary Functions/`  
  Helper functions and data files:
  - `computeK.m` – computes the state-feedback gains \(K_i\).
  - `optimizeL2.m` – computes the "average" matrices \(A_{\mathrm{bar}}, B_{\mathrm{bar}}, K_{\mathrm{bar}}\) via an L2-type CVX optimization.
  - `generateSwSignal.m` – function that generates a random switching signal with a prescribed dwell-time range.
  - `swsig_values.mat` – **the specific switching signal** used in the paper for the simulations of Section 4.  
    This file was generated once using `generateSwSignal.m` and is provided so that the
    figures in the paper can be reproduced exactly.  
    The function `generateSwSignal.m` is included only for reference: if the user calls
    it again, a *different* switching signal will be produced, leading to different
    trajectories than those shown in the article (swsig_values.mat should be redefined accordingly).
