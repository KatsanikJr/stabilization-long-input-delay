function [M_bar] = optimizeL2(M1,M2,M3)
% optimizeL2  Computes the L2-optimal averaged matrix \bar{M} for three matrices.
%
%   M_bar = optimizeL2(M1, M2, M3) solves the convex optimization problem
%
%       minimize    R
%       subject to  ||M_bar - M_i||_2 <= R      for i = 1,2,3
%
%   where R ≥ 0, and the norm is the |·|_2 norm.
%
%  In the paper, this is used to
%   compute:
%
%       - \bar{A} from {A1, A2, A3},  Equation (8)
%       - \bar{B} from {B1, B2, B3},  Equation (8)
%       - \bar{K} from {K1, K2, K3},  Equation (9)
%
%   Inputs:
%       M1, M2, M3   - matrices of identical dimensions
%
%   Output:
%       M_bar        - the L2-optimal averaged matrix
%
%   Requirements:
%       CVX (http://cvxr.com/cvx) with an SDP solver installed.
%
%   Notes:
%       - Uses CVX in quiet mode so that simulations remain clean.
%       - Norm |·|_2 is used for optimization.

    cvx_begin sdp quiet
        variable M_bar(size(M1,1), size(M1,2))  % same dimension as each R_i
        variable R nonnegative

        minimize( R )
        subject to
            norm( M_bar - M1, 2 ) <= R;
            norm( M_bar - M2, 2 ) <= R;
            norm( M_bar - M3, 2 ) <= R;
    cvx_end
   
end