# Stabilization of Linear Switched Systems with Long Constant Input Delay via Average or Averaging Predictor Feedbacks

The purpose of these MATLAB codes is to reproduce the simulation results from Section 4 of the paper "Stabilization of Linear Switched Systems with Long {Constant} Input Delay via Average or Averaging Predictor Feedbacks" by Andreas Katsanikakis and Nikolaos Bekiaris-Liberis which is submitted, as a revised version, in the Control System Letters. The preprint of the paper is available in "[Link for preprint](https://arxiv.org/abs/2506.03908)".

The mathematical background and the documentation for the codes are briefly desribed in the file "READMEcod.pdf.".

## Requirements

The codes require:

- **MATLAB**
- **CVX toolbox** (with a compatible SDP solver such as SeDuMi, SDPT3, or MOSEK)  
- **Symbolic Math Toolbox** (for the computation of theoretical constants in `epsilon_star.m` and `barepsilon_star.m`)
## Usage

The folder names (Section_4_1, Section_4_2, Section_4_3) correspond directly to the section numbering in the paper, to make it easy to locate the scripts associated with each set of simulation results.

The simulation results of Section 4 can be reproduced by running the scripts contained in the folders:

#### Section_4_1/ — Implements controllers U1 and U2 with dwell time knowledge.
Run: 

U1_dwell_time.m

U2_dwell_time.m

epsilon_star.m

barepsilon_star.m


####  Section_4_2/ — Comparison with the exact predictor and no dwell time knowledge cases.
Run:

U1_nodwell_time.m

U2_nodwell_time.m

ExactPredictor_Uex.m

Then compute the performance index. 

Run:

Performance_Index.m


####  Section_4_3/ — Robustness to delay mismatches (delay perturbations).
Run:

U1_robust_0_95.m

U1_robust_1_05.m

U2_robust_0_95.m

U2_robust_1_05.m


## License

Copyright Andreas Katsanikakis 2025. See LICENSE.txt for licensing information.

## Acknowledgements

Funded by the European Union (ERC, C-NORA, 101088147). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Cite this work

arXiv preprint: (https://arxiv.org/abs/2506.03908)
```
@Unpublished{KatsBek25carxiv,
  author = {A. Katsanikakis and N. Bekiaris-Liberis},
  note   = {{arXiv}, 	2506.03908, 2025},
  title  = {Stabilization of Linear Switched Systems with Long Constant Input Delay via Average or Averaging Predictor Feedbacks},
}


