function gain = computeK(A,B,pole1,pole2)
    % Compute the controllability matrix for the mean pair and check its rank
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