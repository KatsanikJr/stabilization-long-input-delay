function switchingVals = generateSwSignal(tau_min, tau_max, T, dt)
% Create a random switching signal sigma(i) in {1,2,3} 
% with a lower dwell time = tau_min and upper dwell time = tau_max
    numSteps = floor(T/dt);
    switchingVals = zeros(1, numSteps+1);

    % Start with a random mode among {1,2,3}.
    currentMode = randi([1, 3]);
    
    % Time-step counter
    i = 1;  
    
    while i <= numSteps
        
        switchingVals(i) = currentMode;

        % Randomly pick a dwell time in [tau_min, tau_max].
        thisDwell = tau_min + (tau_max - tau_min)*rand();  % in seconds
        dwellSteps = floor(thisDwell / dt);

        % Stay in current mode
        end_i = min(i + dwellSteps - 1, numSteps);

        switchingVals(i : end_i) = currentMode;

        i = end_i + 1;

        % Pick a new mode from {1,2,3}, but not the same as currentMode.
        newMode = currentMode;
        while newMode == currentMode
            newMode = randi([1, 3]);
        end
        
        % Update currentMode
        currentMode = newMode;
    end

    %  If the last slot is unassigned, match it to the second-to-last
    if switchingVals(end) == 0
        switchingVals(end) = switchingVals(end-1);
    end
end