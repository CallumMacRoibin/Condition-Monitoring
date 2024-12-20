function stop = stopIfOverfitting(info)
% stopIfOverfitting - Custom early stopping function for training.
% Stops training if the validation loss increases for 3 consecutive
% validation checks.

    % Initialize stop flag
    stop = false;

    % Check current training state
    if info.State == "iteration"
        % Persistent variables to track validation loss
        persistent bestValLoss valLossCounter

        % Initialize variables if empty
        if isempty(bestValLoss) || all(isnan(info.ValidationLoss))
            bestValLoss = info.ValidationLoss;
            valLossCounter = 0;

        % Check if validation loss increased
        elseif info.ValidationLoss > bestValLoss
            valLossCounter = valLossCounter + 1;
            fprintf('Validation loss increased. Count: %d\n', valLossCounter);

        % Reset counter if validation loss improves
        else
            valLossCounter = 0;
            bestValLoss = info.ValidationLoss;
        end

        % Stop training after 3 consecutive increases
        if valLossCounter >= 3
            stop = true;
            disp('Stopping early due to overfitting: Validation loss increased 3 times.');
        end
    end
end
