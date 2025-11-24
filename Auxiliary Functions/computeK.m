function gain = computeK(A,B,pole1,pole2)
% computeK  Designs a state-feedback gain K so that A + B*K is Hurwitz.
%
%   gain = computeK(A,B,pole1,pole2) computes a row vector gain such that
%   the closed-loop matrix (A + B*gain) has eigenvalues placed at
%   [pole1, pole2] chosen by the designer.
%
%   This function:
%     1) Checks controllability of the pair (A,B) via the controllability
%        matrix C = [B, A*B].
%     2) If controllable, uses the built-in function PLACE to compute K.
%     3) Returns gain = -K so that:
%            u(t) = gain * x(t)
%        corresponds to a closed-loop A + B*gain.
%
    C = [B, A * B];
    rankC = rank(C);
    if rankC == size(A, 1)
        disp('The system is controllable.');
        desired_poles = [pole1, pole2];
        K = place(A, B, desired_poles);
    else
        disp('The system is not controllable.');
    end
    gain = - K;

end