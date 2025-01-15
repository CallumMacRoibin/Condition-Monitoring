%% Real-Time Graph GUI with Fault Prediction
function realTimeGraphGUI()
    % Create the main figure for the GUI
    fig = uifigure('Name', 'Real-Time Graph GUI', 'Position', [100, 100, 600, 400]);

    % File selection components
    selectFileButton = uibutton(fig, ...
        'Text', 'Select MATLAB File', ...
        'Position', [50, 350, 150, 30], ...
        'ButtonPushedFcn', @(btn, event)selectFileCallback());

    filePathLabel = uilabel(fig, ...
        'Text', 'No file selected.', ...
        'Position', [220, 350, 300, 30], ...
        'HorizontalAlignment', 'left');

    % Start real-time plotting button
    startPlottingButton = uibutton(fig, ...
        'Text', 'Start Real-Time Plot', ...
        'Position', [50, 300, 150, 30], ...
        'ButtonPushedFcn', @(btn, event)startPlottingCallback());

    % Add a Record button
    recordButton = uibutton(fig, ...
        'Text', 'Record Data', ...
        'Position', [220, 300, 150, 30], ...
        'ButtonPushedFcn', @(btn, event)recordDataCallback());

    % Axes for real-time plotting
    ax = uiaxes(fig, 'Position', [50, 50, 500, 200]);
    title(ax, 'Real-Time Vibrational Data');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Amplitude');
    grid(ax, 'on');

    % Store the selected file path in the figure's UserData
    fig.UserData.filePath = '';
    fig.UserData.sampleRate = [];
    fig.UserData.vibrationalData = [];
    fig.UserData.scalogramFolder = fullfile(pwd, 'Temp_Scalogram');

    % Callback for the file selection button
    function selectFileCallback()
        % Open a file selection dialog
        [fileName, filePath] = uigetfile('*.mat', 'Select a MATLAB File');
        if fileName ~= 0 % If a file is selected
            fullFilePath = fullfile(filePath, fileName);
            filePathLabel.Text = ['Selected File: ', fileName];
            fig.UserData.filePath = fullFilePath; % Store the file path

            % Load the .mat file and extract vibrational data
            try
                loadedData = load(fullFilePath);
                if isfield(loadedData, 'bearing')
                    bearingData = loadedData.bearing;
                    if isfield(bearingData, 'sr') && isfield(bearingData, 'gs')
                        fig.UserData.sampleRate = bearingData.sr;
                        fig.UserData.vibrationalData = bearingData.gs;
                        uialert(fig, 'File loaded successfully.', 'Success');
                    else
                        uialert(fig, 'File does not contain required fields: "sr" or "gs".', 'Error');
                    end
                else
                    uialert(fig, 'File does not contain the "bearing" struct.', 'Error');
                end
            catch ME
                uialert(fig, sprintf('Error loading file: %s', ME.message), 'Error');
            end
        else
            filePathLabel.Text = 'No file selected.';
            fig.UserData.filePath = ''; % Clear the file path
        end
    end

    % Callback for starting the real-time plot
    function startPlottingCallback()
        if isempty(fig.UserData.filePath)
            uialert(fig, 'No file selected. Please select a file before starting.', 'File Not Selected');
        elseif isempty(fig.UserData.vibrationalData)
            uialert(fig, 'No valid vibrational data found. Please check the file.', 'Data Not Found');
        else
            % Extract data and parameters
            sampleRate = fig.UserData.sampleRate;
            vibrationalData = fig.UserData.vibrationalData;

            if ~isvector(vibrationalData)
                uialert(fig, 'Vibrational data is not a vector. Cannot plot.', 'Data Error');
                return;
            end

            % Normalize the vibrational data
            vibrationalData = vibrationalData / max(abs(vibrationalData));

            % Real-time plotting
            frameSize = round(sampleRate / 10); % 0.1 second per frame
            totalSamples = length(vibrationalData);
            numFrames = ceil(totalSamples / frameSize);
            timeStep = 1 / sampleRate;

            % Preallocate for performance
            currentTime = nan(1, totalSamples);
            currentData = nan(1, totalSamples);

            % Index tracker
            currentIndex = 0;

            % Real-time plotting loop
            for frame = 1:numFrames
                % Determine frame range
                startIdx = (frame - 1) * frameSize + 1;
                endIdx = min(frame * frameSize, totalSamples);
                newTime = ((startIdx:endIdx) - 1) * timeStep;
                newData = vibrationalData(startIdx:endIdx);

                % Append new data
                numNewSamples = length(newTime);
                currentTime(currentIndex + 1:currentIndex + numNewSamples) = newTime;
                currentData(currentIndex + 1:currentIndex + numNewSamples) = newData;

                % Update the current index
                currentIndex = currentIndex + numNewSamples;

                % Update the plot
                plot(ax, currentTime(1:currentIndex), currentData(1:currentIndex));
                drawnow;

                % Pause to simulate real-time playback
                pause(length(newData) / sampleRate);
            end
        end
    end

    % Callback for the record button
    function recordDataCallback()
        if isempty(fig.UserData.vibrationalData)
            uialert(fig, 'No vibrational data available to record.', 'Data Not Found');
            return;
        end

        % Capture recent data for recording
        sampleRate = fig.UserData.sampleRate;
        scalogramFolder = fig.UserData.scalogramFolder;
        vibrationData = fig.UserData.vibrationalData;

        % Apply bandpass filter to the captured data
        try
            filteredData = bandpassfiltering(vibrationData, sampleRate);
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
            uialert(fig, sprintf('Scalogram saved: %s', scalogramFile), 'Success');
        catch ME
            warning(ME.identifier, 'Error generating scalogram: %s', ME.message);
        end
    end

    % Callback to classify scalogram using trained CNN
    function classifyScalogram()
        % Path to the Temp_Scalogram folder
        scalogramFolder = fig.UserData.scalogramFolder;

        % Load the trained CNN
        modelPath = fullfile(pwd, 'trainedCNN.mat');
        if exist(modelPath, 'file')
            loadedModel = load(modelPath);
            if isfield(loadedModel, 'net') % 'net' contains the trained CNN
                trainedCNN = loadedModel.net; % Load the CNN into a variable
            else
                uialert(fig, 'The .mat file does not contain a valid trained CNN model.', 'Error');
                return;
            end
        else
            uialert(fig, 'The trained CNN file does not exist.', 'Error');
            return;
        end

        % Get the most recently saved scalogram image
        imageFiles = dir(fullfile(scalogramFolder, '*.jpg'));
        if isempty(imageFiles)
            uialert(fig, 'No scalogram images found.', 'Error');
            return;
        end

        % Sort files by modification date and select the latest
        [~, idx] = max([imageFiles.datenum]);
        latestScalogram = fullfile(scalogramFolder, imageFiles(idx).name);

        % Read and preprocess the image
        try
            img = imread(latestScalogram);
            % Resize the image to match the input size of the CNN
            img = imresize(img, [227, 227]); % Adjust to your CNN's input size

            % Convert to single precision if required
            if isa(trainedCNN.Layers(1), 'nnet.cnn.layer.ImageInputLayer')
                img = single(img);
            end

            % Predict the class using the CNN
            [predictedLabel, scores] = classify(trainedCNN, img);

            % Display the results
            uialert(fig, sprintf('Predicted Fault Class: %s\nConfidence Score: %.2f%%', ...
                string(predictedLabel), max(scores) * 100), 'Classification Result');

        catch ME
            uialert(fig, sprintf('Error processing or classifying the image: %s', ME.message), 'Error');
        end
    end

    % Add a Classify button
    classifyButton = uibutton(fig, ...
        'Text', 'Classify Scalogram', ...
        'Position', [380, 300, 150, 30], ...
        'ButtonPushedFcn', @(btn, event)classifyScalogram());
end