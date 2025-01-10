function recordVibrationData(sampleRate, capturedData, scalogramFile)
    % RECORDVIBRATIONDATA Generate and save a scalogram from captured data.
    %
    % Inputs:
    %   sampleRate   - Sampling rate of the data (Hz).
    %   capturedData - Array of captured vibrational data.
    %   scalogramFile - File path to save the scalogram.

    % Call the existing scalogram generation logic
    convertSignalToScalogramLive(capturedData, sampleRate, scalogramFile);
end