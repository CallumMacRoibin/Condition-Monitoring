%% Comprehensive Script for Filtering Data, Scalogram Generation and CNN Training
% Add access to folder
% Specify the folder you want to add to the path
folderToAdd = fullfile(pwd, 'Experimental_Data');

% Generate path for the specified folder and all its subfolders
folderAndSubfoldersPath = genpath(folderToAdd);

% Add the generated path to MATLAB path
addpath(folderAndSubfoldersPath);

% Get a list of all subfolders
subFolders = dir(folderToAdd);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

% Ensure access to custom functions
addpath(fullfile(pwd, 'Functions'));

%% Step 1: Apply Bandpass Filter to Vibrational Data
% Initialize Counters
successfulFilters = 0;
failedFilters = 0;

for i = 1:length(subFolders)
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    matFiles = dir(fullfile(currentFolder, '*.mat'));

    for j = 1:length(matFiles)
        % Load the .mat file
        matFileName = fullfile(currentFolder, matFiles(j).name);
        data = load(matFileName);

        % Check if the required struct and fields exist
        if isfield(data, 'bearing')
            bearing = data.bearing;

            if isfield(bearing, 'gs') && isfield(bearing, 'sr')
                vibrational_data = bearing.gs;
                sample_rate = bearing.sr;

                % Clean the vibrational data
                vibrational_data = vibrational_data(~isnan(vibrational_data));
                if isempty(vibrational_data)
                    warning('Vibrational data in %s contains only NaN values.', matFileName);
                    failedFilters = failedFilters + 1;
                    continue;
                end

                % Apply custom bandpass filter function
                try
                    filtered_signal = bandpassfiltering(vibrational_data, sample_rate);
                    successfulFilters = successfulFilters + 1;
                catch ME
                    warning('Error filtering data in %s: %s', matFileName, ME.message);
                    failedFilters = failedFilters + 1;
                    continue;
                end

                % Save filtered data back to .mat file
                bearing.filtered_signal = filtered_signal;
                save(matFileName, 'bearing');
            else
                warning('Required fields missing in struct "bearing" in %s', matFileName);
                failedFilters = failedFilters + 1;
            end
        else
            warning('Struct "bearing" missing in %s', matFileName);
            failedFilters = failedFilters + 1;
        end
    end
end

% Summary Message
disp('Filtering process complete.');
disp(['Number of successful filters: ', num2str(successfulFilters)]);
disp(['Number of failed filters: ', num2str(failedFilters)]);

%% Step 2: Generate Scalogram Images
outputFolder = fullfile(pwd, 'Scalograms');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Initialize scalogram count per fault condition
scalogramCounts = struct();

for i = 1:length(subFolders)
    % Sanitize folder name for valid struct field name
    sanitizedFolderName = matlab.lang.makeValidName(subFolders(i).name);
    
    % Create subfolder in Scalograms folder matching the current subfolder
    scalogramSubFolder = fullfile(outputFolder, subFolders(i).name);
    if ~exist(scalogramSubFolder, 'dir')
        mkdir(scalogramSubFolder);
    end

    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    matFiles = dir(fullfile(currentFolder, '*.mat'));

    scalogramCounts.(sanitizedFolderName) = 0; % Initialize counter for the current fault condition

    for j = 1:length(matFiles)
        % Load the .mat file
        matFileName = fullfile(currentFolder, matFiles(j).name);
        data = load(matFileName);

        if isfield(data, 'bearing')
            bearing = data.bearing;

            if isfield(bearing, 'filtered_signal') && isfield(bearing, 'sr')
                filtered_signal = bearing.filtered_signal;
                sample_rate = bearing.sr;

                % Create a table to mimic an ensemble for compatibility
                signalTable = table({filtered_signal}, {sample_rate}, ...
                    {subFolders(i).name}, {matFiles(j).name}, ...
                    'VariableNames', {'gs', 'sr', 'Label', 'FileName'});

                % Convert the signal to scalograms
                try
                    % Call function and get the number of scalograms generated
                    numScalograms = convertSignalToScalogram(signalTable, outputFolder);
                    scalogramCounts.(sanitizedFolderName) = scalogramCounts.(sanitizedFolderName) + numScalograms;
                catch ME
                    warning('Error generating scalogram for %s: %s', matFileName, ME.message);
                end
            else
                warning('Filtered signal or sampling rate missing in struct "bearing" in %s', matFileName);
            end
        else
            warning('Struct "bearing" missing in %s', matFileName);
        end
    end
end

% Display scalogram counts per fault condition
disp('Scalogram generation summary:');
fields = fieldnames(scalogramCounts);
for i = 1:numel(fields)
    disp([fields{i}, ': ', num2str(scalogramCounts.(fields{i}))]);
end

%% Step 3: Prepare Dataset for CNN Training

% Ask user whether to balance the dataset
disp('Would you like to balance the dataset? (yes/no)');
userInput = input('Enter your choice: ', 's');

% Define the paths
inputFolder = fullfile(pwd, 'Scalograms'); % Use the Scalograms directory
outputFolder = fullfile(pwd, 'Balanced_Scalograms'); % Path to the balanced dataset

if strcmpi(userInput, 'yes')
    % Create a new directory for the balanced dataset
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Load the image datastore
    imds = imageDatastore(inputFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

    % Determine the minimum number of images across all classes
    classes = unique(imds.Labels);
    minImagesPerClass = inf;
    for i = 1:numel(classes)
        classIdx = find(imds.Labels == classes(i));
        minImagesPerClass = min(minImagesPerClass, numel(classIdx));
    end

    % Balance the dataset
    for i = 1:numel(classes)
        % Get all images for the current class
        classIdx = find(imds.Labels == classes(i));

        % Shuffle and select images
        selectedIdx = classIdx(randperm(numel(classIdx), minImagesPerClass));

        % Create class-specific folder in the output directory
        classOutputFolder = fullfile(outputFolder, char(classes(i)));
        if ~exist(classOutputFolder, 'dir')
            mkdir(classOutputFolder);
        end

        % Copy selected images to the new folder
        for j = 1:numel(selectedIdx)
            sourceFile = imds.Files{selectedIdx(j)};
            [~, fileName, fileExt] = fileparts(sourceFile);
            destinationFile = fullfile(classOutputFolder, [fileName, fileExt]);
            copyfile(sourceFile, destinationFile);
        end

        fprintf('Copied %d images to %s', numel(selectedIdx), classOutputFolder);
    end
    fprintf('Dataset balancing complete. Images are in %s', outputFolder);
else
    disp('Dataset balancing skipped.');
    outputFolder = inputFolder; % Use original dataset
end

%% Step 4: Train CNN
% Check if balanced dataset exists
if exist(fullfile(pwd, 'Balanced_Scalograms'), 'dir')
    datasetFolder = fullfile(pwd, 'Balanced_Scalograms');
    disp('Using balanced dataset for training.');
else
    datasetFolder = fullfile(pwd, 'Scalograms');
    disp('Balanced dataset not found. Using original Scalograms for training.');
end

% Load the balanced dataset
imds = imageDatastore(outputFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

% Split dataset into training and validation sets
[imdsTrain, imdsValidation] = splitEachLabel(imds, 0.8, 'randomize');

% Load and modify the pre-trained SqueezeNet model
net = squeezenet;
lgraph = layerGraph(net);

% Get the number of classes
numClasses = numel(categories(imdsTrain.Labels));

% Replace final layers to match the number of classes
newConvLayer = convolution2dLayer([1, 1], numClasses, 'WeightLearnRateFactor', 10, 'BiasLearnRateFactor', 10, 'Name', 'new_conv');
newClassificationLayer = classificationLayer('Name', 'new_classoutput');
lgraph = replaceLayer(lgraph, 'conv10', newConvLayer);
lgraph = replaceLayer(lgraph, 'ClassificationLayer_predictions', newClassificationLayer);

% Define training options
options = trainingOptions('sgdm', ...
    'InitialLearnRate', 0.001, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 2, ...
    'MaxEpochs', 4, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', imdsValidation, ...
    'ValidationFrequency', 30, ...
    'Verbose', false, ...
    'MiniBatchSize', 20, ...
    'Plots', 'training-progress');

% Train the network
net = trainNetwork(imdsTrain, lgraph, options);

% Save the trained network
save('trainedCNN.mat', 'net', 'options', 'imdsTrain', 'imdsValidation');
disp('Model training complete. Network saved as trainedCNN.mat');

%% Step 5: Validate and Test CNN
% Load the test data
testFolder = datasetFolder; % Use the same folder for test images
imdsTest = imageDatastore(testFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

% Classify the test images
YPred = classify(net, imdsTest, 'MiniBatchSize', 20);
YTest = imdsTest.Labels;

% Calculate accuracy
accuracy = sum(YPred == YTest) / numel(YTest);
fprintf('Test Accuracy: %.2f%%\n', accuracy * 100);

% Visualize classification results
figure;
confusionchart(YTest, YPred);
title('Confusion Matrix for Test Data');


