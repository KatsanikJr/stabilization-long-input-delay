function switchingVals = generateSwSignal(tau_min, tau_max, T, dt)
% Create a random switching signal sigma(i) in {1,2,3} 
% with a lower dwell time = tau_min and upper dwell time = tau_max,
% and ensuring that we never remain in the same mode when we switch.

    numSteps = floor(T/dt);
    switchingVals = zeros(1, numSteps+1);

    % 1) Start with a random mode among {1,2,3}.
    currentMode = randi([1, 3]);
    
    % 2) Time-step counter
    i = 1;  
    
    while i <= numSteps
        
        % 2.1) Fill the current mode at step i
        switchingVals(i) = currentMode;

        % 2.2) Randomly pick a dwell time in [tau_min, tau_max].
        thisDwell = tau_min + (tau_max - tau_min)*rand();  % in seconds
        dwellSteps = floor(thisDwell / dt);

        % 2.3) We'll stay in 'currentMode' from index i up to end_i.
        end_i = min(i + dwellSteps - 1, numSteps);

        % 2.4) Fill these steps with 'currentMode'.
        switchingVals(i : end_i) = currentMode;

        % 2.5) Advance i
        i = end_i + 1;

        % 2.6) Pick a new mode from {1,2,3}, but not the same as currentMode.
        newMode = currentMode;
        while newMode == currentMode
            newMode = randi([1, 3]);
        end
        
        % 2.7) Update currentMode
        currentMode = newMode;
    end

    % 3) If the last slot is unassigned, match it to the second-to-last
    if switchingVals(end) == 0
        switchingVals(end) = switchingVals(end-1);
    end
end