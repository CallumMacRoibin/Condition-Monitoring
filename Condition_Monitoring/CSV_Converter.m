%% CSV to .mat conversion
% Add access to folder
% Specify the folder you want to add to the path
folderToAdd = fullfile(pwd, 'Experimental_Data');

% Generate path for the specified folder and all its subfolders
folderAndSubfoldersPath = genpath(folderToAdd);

% Add the generated path to MATLAB path
addpath(folderAndSubfoldersPath);

%% Step 1: Create .mat files from .csv files
% Get a list of all subfolders
subFolders = dir(folderToAdd);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

for i = 1:length(subFolders)
    % Get the path of the current subfolder
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    csvFiles = dir(fullfile(currentFolder, '*.csv'));

    for j = 1:length(csvFiles)
        % Read the CSV files
        filename = fullfile(currentFolder, csvFiles(j).name);
        data = readmatrix(filename);

        % Extract the data
        sample_rate = data(1, end-1);
        vibrational_data = data(:, end);

        % Create .mat files
        matFileName = fullfile(currentFolder, [csvFiles(j).name(1:end-4), '.mat']);

        % Save the sample rate and vibrational data to a .mat file
        save(matFileName, 'sample_rate', 'vibrational_data');
    end
end

%% Step 2: Delete .csv files
% Iterate over each folder to delete .csv files
for i = 1:length(subFolders)
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    csvFiles = dir(fullfile(currentFolder, '*.csv'));

    for j = 1:length(csvFiles)
        csvFileName = fullfile(currentFolder, csvFiles(j).name);
        delete(csvFileName);
        fprintf('Deleted: %s\n', csvFileName);
    end
end

fprintf('All .csv files have been deleted.\n');

%% Step 3: Add shaft rate variable to .mat files
input_shaft_rate = 20.425; % Define the shaft rate value

for i = 1:length(subFolders)
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    matFiles = dir(fullfile(currentFolder, '*.mat'));

    for j = 1:length(matFiles)
        matFileName = fullfile(currentFolder, matFiles(j).name);
        dataStruct = load(matFileName);

        % Add the new variable (shaft rate) to the structure
        dataStruct.input_shaft_rate = input_shaft_rate;

        % Save the modified structure back to the .mat file
        save(matFileName, '-struct', 'dataStruct');
        fprintf('Updated: %s\n', matFileName);
    end
end

%% Step 4: Modify .mat files to restructure variables
for i = 1:length(subFolders)
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    matFiles = dir(fullfile(currentFolder, '*.mat'));

    for j = 1:length(matFiles)
        filePath = fullfile(currentFolder, matFiles(j).name);
        data = load(filePath);

        % Check if the required variables exist
        if isfield(data, 'input_shaft_rate') && isfield(data, 'sample_rate') && isfield(data, 'vibrational_data')
            % Rename variables
            rate = data.input_shaft_rate;
            sr = data.sample_rate;
            gs = data.vibrational_data;

            % Create the 1x1 struct with renamed variables
            bearingFolder.rate = rate;
            bearingFolder.sr = sr;
            bearingFolder.gs = gs;

            % Save the new struct in the same file, overwriting the original
            save(filePath, 'bearingFolder');
        else
            warning('Required variables not found in %s', filePath);
        end
    end
end

disp('Processing complete.');

%% Step 5: Rename struct from bearingFolder to bearing
for i = 1:length(subFolders)
    currentFolder = fullfile(folderToAdd, subFolders(i).name);
    matFiles = dir(fullfile(currentFolder, '*.mat'));

    for j = 1:length(matFiles)
        filePath = fullfile(currentFolder, matFiles(j).name);
        data = load(filePath);

        % Check if the 'bearingFolder' struct exists
        if isfield(data, 'bearingFolder')
            % Rename 'bearingFolder' to 'bearing'
            bearing = data.bearingFolder;

            % Save the renamed struct back to the .mat file, overwriting it
            save(filePath, 'bearing');
        else
            warning('Struct ''bearingFolder'' not found in %s', filePath);
        end
    end
end

disp('Renaming complete.');
