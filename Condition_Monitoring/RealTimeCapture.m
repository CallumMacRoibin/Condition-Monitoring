%% Simulate Live Data Feed
% Ensure access to custom functions
addpath(fullfile(pwd, 'Functions'));

% Define the path to the Experimental_Data folder
basePath = fullfile(pwd, 'Experimental_Data');
dataSubfolder = 'Baseline_Readings'; % Subfolder name

% Define the full path to the .mat file
dataFileName = 'Raw 19_03_2024 20_17_42.mat'; % Replace with actual file name
dataFilePath = fullfile(basePath, dataSubfolder, dataFileName);

% Load the .mat file
if exist(dataFilePath, 'file')
    loadedData = load(dataFilePath);
    if isfield(loadedData, 'bearing')
        bearingData = loadedData.bearing;

        if isfield(bearingData, 'sr') && isfield(bearingData, 'gs')
            sampleRate = bearingData.sr; % Sample rate
            vibrationalData = bearingData.gs; % Vibrational data

            if isvector(vibrationalData)
                % Normalize the vibrational data for visualization
                vibrationalData = vibrationalData / max(abs(vibrationalData));

                % Initialize real-time plotting
                figureHandle = figure('Name', 'Real-Time Vibrational Data', 'NumberTitle', 'off');
                hPlot = plot(nan, nan);
                xlim([0, 10]); % Initial x-axis range (10 seconds)
                ylim([-1.1, 1.1]); % Normalized amplitude
                xlabel('Time (s)');
                ylabel('Amplitude');
                title('Real-Time Vibrational Data');
                grid on;

                % Folder to save scalograms
                scalogramFolder = fullfile(pwd, 'Temp_Scalogram');
                if ~exist(scalogramFolder, 'dir')
                    mkdir(scalogramFolder);
                end

                % Add a Record button to the figure
                btn = uicontrol('Parent', figureHandle, ...
                    'Style', 'pushbutton', ...
                    'String', 'Record', ...
                    'Position', [20, 50, 100, 30]);

                % Store necessary variables in the button's UserData
                btn.UserData.sampleRate = sampleRate;
                btn.UserData.scalogramFolder = scalogramFolder;
                btn.UserData.currentData = [];

                % Define the callback for the button
                btn.Callback = @(src, event)recordVibrationDataCallback(src);

                % Parameters for real-time plotting
                frameSize = round(sampleRate / 10); % 0.1 second per frame
                totalSamples = length(vibrationalData);
                numFrames = ceil(totalSamples / frameSize);
                timeStep = 1 / sampleRate;

                % Preallocate arrays based on totalSamples
                currentTime = nan(1, totalSamples); % Preallocate with NaN
                currentData = nan(1, totalSamples); % Preallocate with NaN

                % Index tracker for appending data
                currentIndex = 0;

                % Simulate real-time plotting
                fprintf('Playing and plotting vibrational data in real-time...\n');
                for frame = 1:numFrames
                    % Determine frame range
                    startIdx = (frame - 1) * frameSize + 1;
                    endIdx = min(frame * frameSize, totalSamples);
                    newTime = ((startIdx:endIdx) - 1) * timeStep;
                    newData = vibrationalData(startIdx:endIdx);

                    % Append new data by assigning directly to preallocated arrays
                    numNewSamples = length(newTime);
                    currentTime(currentIndex + 1:currentIndex + numNewSamples) = newTime;
                    currentData(currentIndex + 1:currentIndex + numNewSamples) = newData;

                    % Update the current index
                    currentIndex = currentIndex + numNewSamples;

                    % Update plot
                    set(hPlot, 'XData', currentTime(1:currentIndex), 'YData', currentData(1:currentIndex));

                    % Update the button's UserData with the latest currentData
                    btn.UserData.currentData = currentData(1:currentIndex);

                    % Adjust x-axis dynamically
                    if currentTime(currentIndex) > 10
                        xlim([currentTime(currentIndex) - 10, currentTime(currentIndex)]);
                    end

                    % Pause to simulate real-time playback
                    pause(length(newData) / sampleRate);
                end
                fprintf('Playback and visualization finished.\n');
            else
                error('Vibrational data is not a vector.');
            end
        else
            error('The struct does not contain required fields "sr" or "gs".');
        end
    else
        error('The loaded file does not contain a "bearing" struct.');
    end
else
    error('The specified .mat file does not exist: %s', dataFilePath);
end
%% Fault Prediction
% Path to the Temp_Scalogram folder
scalogramFolder = fullfile(pwd, 'Temp_Scalogram');

% Load the trained CNN
modelPath = fullfile(pwd, 'trainedCNN.mat');
if exist(modelPath, 'file')
    loadedModel = load(modelPath);
    if isfield(loadedModel, 'net') % 'net' contains the trained CNN
        trainedCNN = loadedModel.net; % Load the CNN into a variable
    else
        error('The .mat file does not contain a valid trained CNN model.');
    end
else
    error('The trained CNN file does not exist at %s.', modelPath);
end

% Get the most recently saved scalogram image
imageFiles = dir(fullfile(scalogramFolder, '*.jpg'));
if isempty(imageFiles)
    error('No scalogram images found in %s.', scalogramFolder);
end

% Sort files by modification date and select the latest
[~, idx] = max([imageFiles.datenum]);
latestScalogram = fullfile(scalogramFolder, imageFiles(idx).name);

% Read and preprocess the image
try
    img = imread(latestScalogram);
    % Resize the image to match the input size of the CNN
    img = imresize(img, [227, 227]); % Example size, adjust to your CNN's input size

    % Convert to single precision if required
    if isa(trainedCNN.Layers(1), 'nnet.cnn.layer.ImageInputLayer')
        img = single(img);
    end

    % Predict the class using the CNN
    [predictedLabel, scores] = classify(trainedCNN, img);

    % Display the results
    fprintf('Predicted Fault Class: %s\n', string(predictedLabel));
    fprintf('Confidence Score: %.2f%%\n', max(scores) * 100);

catch ME
    error('Error processing or classifying the image: %s', ME.message);
end
%% Callback
function recordVibrationDataCallback(button)
    % Wrapper to call the separate record function
    sampleRate = button.UserData.sampleRate;
    scalogramFolder = button.UserData.scalogramFolder;
    vibrationData = button.UserData.currentData; % Use vibrationData to match the function

    % Apply bandpass filter to the captured data
    try
        filteredData = bandpassfiltering(vibrationData, sampleRate); % Pass as vibrationData
        fprintf('Bandpass filtering applied to captured data.\n');
    catch ME
        warning(ME.identifier, 'Error applying bandpass filter: %s', ME.message);
        filteredData = vibrationData; % Use unfiltered data as fallback
    end

    % Generate a timestamped file name with a proper extension
    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    scalogramFile = fullfile(scalogramFolder, ['Scalogram_' timestamp '.jpg']);

    % Call the function in the Functions folder
    try
        recordVibrationData(sampleRate, filteredData, scalogramFile);
    catch ME
        warning(ME.identifier, 'Error generating scalogram: %s', ME.message);
    end
end
