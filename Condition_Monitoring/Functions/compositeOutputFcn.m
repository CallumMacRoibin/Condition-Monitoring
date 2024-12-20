function stop = compositeOutputFcn(info)
% compositeOutputFcn - Combines multiple custom OutputFcn callbacks.
% Integrates stopIfOverfitting and pauseTrainingCallback functionality.

    % Initialize stop flag
    stop = false;

    % Call the first OutputFcn (stopIfOverfitting)
    stop1 = stopIfOverfitting(info);

    % Call the second OutputFcn (pauseTrainingCallback)
    stop2 = pauseTrainingCallback(info);

    % Combine the stop signals (stop if either returns true)
    stop = stop1 || stop2;
end
