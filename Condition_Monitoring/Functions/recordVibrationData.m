function recordVibrationData(sampleRate, currentData, scalogramFolder)
    % RECORDVIBRATIONDATA Generate and save a scalogram from captured data.
    %
    % Inputs:
    %   sampleRate     - Sampling rate of the data (Hz).
    %   currentData    - Array of vibrational data.
    %   scalogramFolder - Folder where the scalogram will be saved.

    fprintf('Recording vibrational data for scalogram...\n');

    % Define buffer size for 5 seconds of data
    bufferSize = round(sampleRate * 5); % 5 seconds of data

    % Check if enough data has been captured
    if length(currentData) >= bufferSize
        % Capture the most recent 5 seconds of data
        capturedData = currentData(end-bufferSize+1:end);

        % Generate a timestamped file name
        timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
        scalogramFile = fullfile(scalogramFolder, ['Scalogram_' timestamp '.jpg']);

        % Generate and save the scalogram
        try
            convertSignalToScalogramLive(capturedData, sampleRate, scalogramFile);
            fprintf('Scalogram saved: %s\n', scalogramFile);
        catch ME
            fprintf('Error generating scalogram: %s\n', ME.message);
        end
    else
        fprintf('Not enough data captured yet for a scalogram. Please wait...\n');
    end
end
