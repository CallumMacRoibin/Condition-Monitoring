function stop = pauseTrainingCallback(info)
    % pauseTrainingCallback - Custom function to pause or exit training.
    % Stops training if "pause_training.flag" or "exit_training.flag" is detected.

    stop = false; % Default: do not stop training

    % Check if we're in the iteration state
    if info.State == "iteration"
        % Check for the pause flag file
        if isfile('pause_training.flag')
            disp('Pause flag detected. Stopping training after the current iteration.');
            stop = true; % Signal to stop training
        end

        % Check for the exit flag file
        if isfile('exit_training.flag')
            disp('Exit flag detected. Terminating training.');
            stop = true; % Signal to stop training
            error('Training terminated by user.'); % Force an exit
        end
    end
end