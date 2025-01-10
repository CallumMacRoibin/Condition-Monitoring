function filtered_signal = bandpassfiltering(vibrationData, samplingRate)
    % BANDPASSFILTERING Apply a bandpass filter to the given vibration data.
    %
    % Inputs:
    %   vibrationData - The raw vibration data signal (vector).
    %   samplingRate  - The sampling rate of the vibration data (scalar).
    %
    % Outputs:
    %   filtered_signal - The filtered vibration data signal (vector).
    
    % Remove NaN values
    vibrationData = vibrationData(~isnan(vibrationData));
    if isempty(vibrationData)
        error('Input vibration data contains only NaN values.');
    end

    % Parameters for the bandpass filter
    % kurtogram_level = 9;  % Level used in the Kurtogram
    % Ensure the level does not exceed the maximum allowed for the data length
    maxLevel = floor(log2(length(vibrationData))); % Maximum level for current data length
    kurtogram_level = min(9, maxLevel); % Cap level at 9 or the max possible level
    filterOrder = 200;    % Order of the bandpass FIR filter
    
    % Compute the Kurtogram and extract central frequency (FC) and bandwidth (BW)
    [~, ~, ~, FC, ~, BW] = kurtogram(vibrationData, samplingRate, kurtogram_level);

    % Validate FC and BW
    if isempty(FC) || isempty(BW) || FC <= 0 || BW <= 0
        error('Invalid Kurtogram output: FC or BW is non-positive.');
    end
    
    % Design the bandpass filter
    nyquist = samplingRate / 2;
    low_cutoff = max(0, FC - BW / 2);
    high_cutoff = min(nyquist, FC + BW / 2);

    if low_cutoff >= high_cutoff
        error('Invalid filter design parameters: low_cutoff >= high_cutoff.');
    end

    BPF = designfilt('bandpassfir', 'FilterOrder', filterOrder, ...
        'CutoffFrequency1', low_cutoff, 'CutoffFrequency2', high_cutoff, ...
        'SampleRate', samplingRate);

    % Apply the bandpass filter to the vibration data
    filtered_signal = filter(BPF, vibrationData);
end