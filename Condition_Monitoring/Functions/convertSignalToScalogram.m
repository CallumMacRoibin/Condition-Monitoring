function numScalograms = convertSignalToScalogram(ensemble, folderName)
    % CONVERTSIGNALTOSCALOGRAM Convert 1-D signals to scalograms and save as images.
    %
    % Inputs:
    %   ensemble   - A table with columns:
    %                'gs' (cell array of signals),
    %                'sr' (cell array of sampling rates),
    %                'Label' (string or char array for subfolder labels),
    %                'FileName' (string or char array for file names).
    %   folderName - Root directory to store scalograms (string).
    %
    % Output:
    %   numScalograms - Number of scalograms generated.
    
    % Loop through each row in the ensemble
    numScalograms = 0; % Initialize counter
    for row = 1:height(ensemble)
        % Extract signal, sampling rate, label, and file name
        x = ensemble.gs{row};       % Signal data
        fs = ensemble.sr{row};      % Sampling rate
        label = char(ensemble.Label(row)); % Subfolder name
        fname = char(ensemble.FileName(row)); % File name for saving

        % Calculate interval for segmentation
        ratio = 5000 / 97656;            % Ratio for segment size
        interval = round(ratio * fs);    % Ensure interval is an integer
        N = floor(numel(x) / interval);  % Number of segments

        % Create subfolder for the current label
        path = fullfile(folderName, label);
        if ~exist(path, 'dir')
            mkdir(path);
        end

        % Process each segment
        for idx = 1:N
            % Extract signal segment
            segmentStart = interval * (idx - 1) + 1;
            segmentEnd = interval * idx;
            sig = envelope(x(segmentStart:segmentEnd));

            % Compute Continuous Wavelet Transform (CWT)
            cfs = cwt(sig, 'amor', seconds(1 / fs));
            cfs = abs(cfs); % Take magnitude

            % Convert scalogram to image
            img = ind2rgb(round(rescale(flip(cfs), 0, 255)), jet(320));

            % Generate output file name
            outfname = fullfile(path, [fname '-' num2str(idx) '.jpg']);

            % Save the scalogram image
            imwrite(imresize(img, [227, 227]), outfname);

            % Increment scalogram counter
            numScalograms = numScalograms + 1;
        end
    end
end
