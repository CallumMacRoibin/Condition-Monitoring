function convertSignalToScalogramLive(signal, samplingRate, outputFile)
    % CONVERTSIGNALTOSCALOGRAM Convert a segment of a signal to a scalogram and save as an image.
    %
    % Inputs:
    %   signal       - 1-D vibrational signal (vector).
    %   samplingRate - Sampling rate of the signal (scalar).
    %   outputFile   - Full path of the output file to save the scalogram (string).

    % Calculate the interval for segmentation (identical to the original logic)
    ratio = 5000 / 97656;            % Ratio for segment size
    interval = round(ratio * samplingRate); % Ensure interval is an integer

    % Check if the signal is long enough for one segment
    if numel(signal) < interval
        error('Signal length is too short for a scalogram segment.');
    end

    % Extract the first segment (identical logic to the original version)
    segment = envelope(signal(1:interval)); % Apply envelope detection

    % Compute Continuous Wavelet Transform (CWT) for the segment
    cfs = cwt(segment, 'amor', seconds(1 / samplingRate));
    cfs = abs(cfs); % Take magnitude

    % Convert scalogram to image
    img = ind2rgb(round(rescale(flip(cfs), 0, 255)), jet(320));

    % Save the scalogram image
    imwrite(imresize(img, [227, 227]), outputFile);
end