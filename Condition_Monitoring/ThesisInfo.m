% Define the main folder containing the subfolders
mainFolder = 'Experimental_data';

% Get a list of all subfolders within the main folder
subfolders = dir(mainFolder);
subfolders = subfolders([subfolders.isdir]);  % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));  % Remove '.' and '..'

% Initialize summary variables
summaryTable = table('Size', [numel(subfolders) 4], ...
                     'VariableTypes', {'string', 'int32', 'int64', 'double'}, ...
                     'VariableNames', {'Subfolder', 'NumFiles', 'TotalDataPoints', 'TotalSizeGB'});

% Loop over each subfolder to summarize the data
for i = 1:numel(subfolders)
    % Get the full path to the subfolder
    subfolderPath = fullfile(mainFolder, subfolders(i).name);
    
    % Get a list of all .mat files in the subfolder
    matFiles = dir(fullfile(subfolderPath, '*.mat'));
    
    % Initialize variables to hold the number of data points and total size
    totalDataPoints = 0;
    totalSizeBytes = 0;
    
    % Loop over each file to calculate the total data points and size
    for j = 1:numel(matFiles)
        % Get the full path to the .mat file
        matFilePath = fullfile(subfolderPath, matFiles(j).name);
        
        % Get file information
        fileInfo = dir(matFilePath);
        
        % Accumulate the file size in bytes
        totalSizeBytes = totalSizeBytes + fileInfo.bytes;
        
        % Load the .mat file
        fileData = load(matFilePath);
        
        % Access the 'gs' field of the 'bearing' struct
        if isfield(fileData, 'bearing') && isfield(fileData.bearing, 'gs')
            % Sum the number of data points in the 'gs' field
            totalDataPoints = totalDataPoints + numel(fileData.bearing.gs);
        end
    end
    
    % Convert total size to GB
    totalSizeGB = totalSizeBytes / (1024^3);
    
    % Update the summary table
    summaryTable.Subfolder(i) = subfolders(i).name;
    summaryTable.NumFiles(i) = numel(matFiles);
    summaryTable.TotalDataPoints(i) = totalDataPoints;
    summaryTable.TotalSizeGB(i) = totalSizeGB;
end

% Display the summary table
disp(summaryTable);

% Optionally, save the summary table to a file
writetable(summaryTable, fullfile(mainFolder, 'ExperimentalDataSummary.csv'));
%%
% Define the path to the Baseline subfolder
baselineFolderPath = fullfile('Experimental_data', 'Imbalance');

% Get a list of all MATLAB files in the Baseline subfolder
matFiles = dir(fullfile(baselineFolderPath, '*.mat'));

% Initialize variables to store results
allAmplitudes = [];
rotationalFreq = 10;  % Define the critical rotational frequency (Hz)

% Loop through each MATLAB file in the Baseline subfolder
for j = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(baselineFolderPath, matFiles(j).name);
    dataStruct = load(filePath);
    
    % Extract the vibrational data and sampling rate
    vibrationalData = dataStruct.bearing.gs;
    samplingRate = dataStruct.bearing.sr;
    
    % Apply FFT to the vibrational data
    L = length(vibrationalData);
    f = samplingRate*(0:(L/2))/L;
    Y = fft(vibrationalData);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Find amplitude at the rotational frequency
    [~, idx] = min(abs(f - rotationalFreq));
    amplitudeAtRotFreq = P1(idx);
    
    % Store the amplitude
    allAmplitudes = [allAmplitudes; amplitudeAtRotFreq];
end

% Calculate statistics of the amplitudes
meanAmplitude = mean(allAmplitudes);
varianceAmplitude = var(allAmplitudes);
stdAmplitude = std(allAmplitudes);

% Display the results
disp(['Mean Amplitude at ', num2str(rotationalFreq), ' Hz: ', num2str(meanAmplitude)]);
disp(['Variance of Amplitude: ', num2str(varianceAmplitude)]);
disp(['Standard Deviation of Amplitude: ', num2str(stdAmplitude)]);

% Plot amplitudes for visual inspection
figure;
plot(allAmplitudes, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.4470 0.7410]); % Blue color for clarity
xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude at Rotational Frequency (m/s^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Amplitude at Rotational Frequency (20.425 Hz) Across Imbalance Measurements', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Set axis properties for better presentation
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 length(allAmplitudes)]);
ylim([0 max(allAmplitudes)*1.1]); % Set y-axis limit slightly above max amplitude for better visibility

% Add a trendline
hold on;
p = polyfit(1:length(allAmplitudes), allAmplitudes, 1);
trendline = polyval(p, 1:length(allAmplitudes));
plot(1:length(allAmplitudes), trendline, '--r', 'LineWidth', 2);
legend('Grouped Amplitude', 'Trendline', 'Location', 'best', 'FontSize', 12);
hold off;


% Add a legend if necessary (optional)
legend('Amplitude at 20.425 Hz', 'Location', 'best', 'FontSize', 12);
%%
% Define the path to the Baseline subfolder
% Define the path to the Baseline subfolder
baselineFolderPath = fullfile('Experimental_data', 'Imbalance');

% Get a list of all MATLAB files in the Baseline subfolder
matFiles = dir(fullfile(baselineFolderPath, '*.mat'));

% Initialize variables to store results
allAmplitudes = [];
rotationalFreq = 20.425;  % Define the critical rotational frequency (Hz)

% Loop through each MATLAB file in the Baseline subfolder
for j = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(baselineFolderPath, matFiles(j).name);
    dataStruct = load(filePath);
    
    % Extract the vibrational data and sampling rate
    vibrationalData = dataStruct.bearing.gs;
    samplingRate = dataStruct.bearing.sr;
    
    % Apply FFT to the vibrational data
    L = length(vibrationalData);
    f = samplingRate*(0:(L/2))/L;
    Y = fft(vibrationalData);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Find amplitude at the rotational frequency
    [~, idx] = min(abs(f - rotationalFreq));
    amplitudeAtRotFreq = P1(idx);
    
    % Store the amplitude
    allAmplitudes = [allAmplitudes; amplitudeAtRotFreq];
end
allAmplitudes = sort([allAmplitudes; amplitudeAtRotFreq]);
% Grouping the amplitudes in sets of three and averaging
groupedAmplitudes = [];
groupIndices = [];

for i = 1:3:length(allAmplitudes)-2
    groupedAmplitudes = [groupedAmplitudes; mean(allAmplitudes(i:i+2))];
    groupIndices = [groupIndices; i:i+2];
end

% Calculate statistics of the grouped amplitudes
meanAmplitude = mean(groupedAmplitudes);
varianceAmplitude = var(groupedAmplitudes);
stdAmplitude = std(groupedAmplitudes);

% Display the results
disp(['Mean Grouped Amplitude at ', num2str(rotationalFreq), ' Hz: ', num2str(meanAmplitude)]);
disp(['Variance of Grouped Amplitude: ', num2str(varianceAmplitude)]);
disp(['Standard Deviation of Grouped Amplitude: ', num2str(stdAmplitude)]);

% Tabulate results: Reading Number vs Imbalance Amplitude
readingNumbers = (1:length(allAmplitudes))';
T = table(readingNumbers, allAmplitudes, 'VariableNames', {'ReadingNo', 'ImbalanceAmplitude'});
disp(T);

% Plot individual amplitudes and indicate grouping visually
figure;
hold on;

colors = ['r', 'g', 'b']; % Define colors for the groups
markers = ['o', 's', 'd']; % Define markers for the groups

% Loop over the groups and plot with different markers/colors
for k = 1:length(groupedAmplitudes)
    currentColor = colors(mod(k-1, length(colors)) + 1); % Cycle through colors
    currentMarker = markers(mod(k-1, length(markers)) + 1); % Cycle through markers
    
    % Plot each group of three points with a line connecting them
    plot(groupIndices(k,:), allAmplitudes(groupIndices(k,:)), '-', 'Color', currentColor, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(groupIndices(k,:), allAmplitudes(groupIndices(k,:)), currentMarker, 'MarkerSize', 8, 'LineWidth', 2, 'Color', currentColor, 'HandleVisibility', 'off'); 
end

% Plot the grouped averages with lines connecting them
plot(groupIndices(:,2), groupedAmplitudes, '-k', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Group Average');

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude at Rotational Frequency (m/s^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Amplitude at Rotational Frequency (20.425 Hz) Across Imbalance Measurements', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Set axis properties for better presentation
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 length(allAmplitudes) + 1]);
ylim([0 max(allAmplitudes)*1.1]); % Set y-axis limit slightly above max amplitude for better visibility

% Add a trendline and calculate R-squared value
p = polyfit(1:length(allAmplitudes), allAmplitudes, 1);
trendline = polyval(p, 1:length(allAmplitudes));
plot(1:length(allAmplitudes), trendline, '--r', 'LineWidth', 2, 'DisplayName', 'Trendline');

% Step 4: Calculate R-squared value
% Total sum of squares
SS_tot = sum((allAmplitudes - mean(allAmplitudes)).^2);

% Residual sum of squares
SS_res = sum((allAmplitudes - trendline).^2);

% R-squared value
R2 = 1 - (SS_res / SS_tot);

% Display R-squared value
disp(['R-squared value: ', num2str(R2)]);
% Add a legend
legend('Group Average', 'Trendline', 'Location', 'best', 'FontSize', 12);

hold off;
%%
% Define the name of the CSV file
csvFileName = 'shaft_data.csv'; % Replace with the actual path to your CSV file

% Read the data from the CSV file into a table
data = readtable(csvFileName);

% Display the first few rows of the table to verify the data
disp(head(data));

% Define the name for the MAT file
matFileName = 'shaft_data.mat'; % Desired name for the MAT file

% Save the table to the MAT file
save(matFileName, 'data');

% Confirm that the MAT file has been created
disp(['MAT file created: ' matFileName]);

% Extract columns from the table
readings = data.Reading;
angularMisalignment = data.AngularMisalignment;
parallelMisalignment = data.ParallelMisalignment;

% Initialize shaft length
L = 380; % Set this to the actual length of the shaft

% Convert misalignment to position
% Assume the shaft is a line from (0,0) to (L,0) in the x-direction.
% Positions are calculated based on misalignment.

% Create figure
figure;
hold on;

% Plot each reading
for i = 1:height(data)
    % Current position of the shaft
    startX = 0; % Starting X position
    startY = 0 + parallelMisalignment(i); % Starting Y position
    
    endX = startX+sqrt(((-endY)^2)+(2*startX*endY)+(L^2)-startY^2); % End X position
    endY = parallelMisalignment(i) + angularMisalignment(i)*(L/100); % End Y position

    % Plot the shaft position
    plot([startX, endX], [startY, endY], 'LineWidth', 2, 'DisplayName', ['Reading ' num2str(readings(i))]);
end

% Customize the plot
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
title('Shaft Positions Relative to Zero Starting Position (X-Axis)');
xlim([0 (round(L/100)*100)]);
ylim([0 2]); % Set y-axis limit slightly above max amplitude for better visibility

legend show; % Show legend
grid on; % Add grid for better visualization
hold off;
%%

% Define the path to the Debris subfolder
debrisFolderPath = fullfile('Experimental_data', 'Debris');

% Get a list of all MATLAB files in the Debris subfolder
matFiles = dir(fullfile(debrisFolderPath, '*.mat'));

% Define the shaft rate and harmonics

shaftRate = 20.425;  % Shaft rate (Hz)
% Bearing Values
rollingElements = 8; % Number of rolling elements in the bearing
ballDiameter = 0.235; % Diameter of the rolling elements (in meters)
pitchDiameter = 1.245; % Pitch diameter of the bearing (in meters)
contactAngle = 0; % Contact angle of the bearing (in degrees)
BPFI = (((rollingElements) * (shaftRate)) / 2) * (1 + ((ballDiameter/pitchDiameter) * cosd(contactAngle))); % Ball Pass Frequency Inner race

harmonics = (1:5) * BPFI;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize a figure for the PSD plots
figure;
hold on;

% Loop through each MATLAB file in the Debris subfolder
for j = 1:length(matFiles)
    try
        % Load the .mat file
        filePath = fullfile(debrisFolderPath, matFiles(j).name);
        dataStruct = load(filePath);
        
        % Check if the structure contains the expected fields
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            % Extract the vibrational data and sample rate
            vibrationalData = dataStruct.bearing.gs;  % Vibrational data
            samplingRate = dataStruct.bearing.sr;     % Sample rate (Hz)
            
            % Compute Power Spectral Density (PSD) using pwelch function
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            % Check the range of the PSD values
            disp(['Trial ', num2str(j), ' PSD Value Range: ', num2str(min(psdValues)), ' to ', num2str(max(psdValues))]);

            % Normalize the PSD values to avoid distortion
            psdValues = psdValues / max(psdValues);  % Normalize to maximum PSD
            
            % Plot the PSD (Convert to dB/Hz)
            plot(freq, 10*log10(psdValues), 'DisplayName', ['Trial ', num2str(j)]);  % Label as Trial 1, 2, 3, etc.
        else
            disp(['File ', matFiles(j).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing file ', matFiles(j).name, ': ', ME.message]);
    end
end

% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--r', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 2, 'LabelOrientation', 'horizontal');
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Power Spectral Density Overlaid for Debris Trials with BPFI Harmonics', 'FontSize', 14, 'FontWeight', 'bold');

% Add a legend for the trials
legend('show', 'Location', 'best', 'FontSize', 10);

% Set axis limits and grid
xlim([0 600]);  % Adjust this range as needed
ylim([-80 0]);  % Set Y limits based on your expected range
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

hold off;

%%
% Define the paths to the Baseline and Debris subfolders
baselineFolderPath = fullfile('Experimental_data', 'Baseline');
debrisFolderPath = fullfile('Experimental_data', 'Debris');

% Get a list of all MATLAB files in each subfolder
baselineFiles = dir(fullfile(baselineFolderPath, '*.mat'));
debrisFiles = dir(fullfile(debrisFolderPath, '*.mat'));

% Define the shaft rate and harmonics
shaftRate = 20.425;  % Shaft rate (Hz)
% Bearing Values
rollingElements = 8; % Number of rolling elements in the bearing
ballDiameter = 0.235; % Diameter of the rolling elements (in meters)
pitchDiameter = 1.245; % Pitch diameter of the bearing (in meters)
contactAngle = 0; % Contact angle of the bearing (in degrees)
BPFI = (((rollingElements) * (shaftRate)) / 2) * (1 + ((ballDiameter/pitchDiameter) * cosd(contactAngle))); % Ball Pass Frequency Inner race

harmonics = (1:5) * shaftRate;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize arrays to store PSD values for averaging
baselinePSDs = [];
referenceFreq = [];  % Reference frequency vector

% Process Baseline Files to compute average PSD
for j = 1:length(baselineFiles)
    try
        % Load the .mat file
        filePath = fullfile(baselineFiles(j).folder, baselineFiles(j).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;

            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));
            
            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            if isempty(referenceFreq)
                referenceFreq = freq;
            else
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
            end
            
            baselinePSDs = [baselinePSDs, psdValues];
        else
            disp(['Baseline file ', baselineFiles(j).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing baseline file ', baselineFiles(j).name, ': ', ME.message]);
    end
end

% Average PSD for Baseline
averageBaselinePSD = mean(baselinePSDs, 2);

% Plotting the average Baseline PSD
figure;
hold on;
plot(referenceFreq, 10*log10(averageBaselinePSD), 'b', 'LineWidth', 1.5, 'DisplayName', 'Average Baseline PSD');

% Process Debris Files to compute average PSD of every three readings and plot each
debrisPlotCount = 1;
for k = 1:3:min(floor(length(debrisFiles)/3)*3, 9)  % Only process first 9 files (3 groups of 3)
    try
        % Initialize temporary array for grouping PSDs
        tempDebrisPSDs = [];
        
        for j = 0:2
            % Load the .mat file
            filePath = fullfile(debrisFiles(k+j).folder, debrisFiles(k+j).name);
            dataStruct = load(filePath);
            
            % Extract the vibrational data and sample rate
            if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
                vibrationalData = dataStruct.bearing.gs;
                samplingRate = dataStruct.bearing.sr;
                
                % Remove NaN values
                vibrationalData = vibrationalData(~isnan(vibrationalData));

                % Compute Power Spectral Density (PSD)
                [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
                
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
                
                tempDebrisPSDs = [tempDebrisPSDs, psdValues];
            else
                disp(['Debris file ', debrisFiles(k+j).name, ' does not contain the expected structure. Skipping...']);
            end
        end
        
        % Average PSD for the current group of three debris readings
        if ~isempty(tempDebrisPSDs)
            avgTempPSD = mean(tempDebrisPSDs, 2);
            % Plot each average PSD of every three debris readings
            plot(referenceFreq, 10*log10(avgTempPSD), 'LineWidth', 0.5, 'DisplayName', ['Average Debris PSD (Group ', num2str(debrisPlotCount), ')']);
            debrisPlotCount = debrisPlotCount + 1;
        end
    catch ME
        disp(['Error processing debris files starting with ', debrisFiles(k).name, ': ', ME.message]);
    end
end

% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--k', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 0.5, 'LabelOrientation', 'horizontal');
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Average PSD Comparison: Baseline vs Debris (Grouped by 3)', 'FontSize', 14, 'FontWeight', 'bold');

% Add a legend
legend('show', 'Location', 'best', 'FontSize', 10);

% Set axis limits and grid
xlim([0 100]);  % Adjust this range as needed
ylim([-80 0]);  % Adjust Y limits based on your data
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

hold off;
%%
% Define the paths to the Baseline and Debris subfolders
baselineFolderPath = fullfile('Experimental_data', 'Baseline');
debrisFolderPath = fullfile('Experimental_data', 'Debris');

% Get a list of all MATLAB files in each subfolder
baselineFiles = dir(fullfile(baselineFolderPath, '*.mat'));
debrisFiles = dir(fullfile(debrisFolderPath, '*.mat'));

% Define the shaft rate and harmonics
shaftRate = 20.425;  % Shaft rate (Hz)
% Bearing Values
rollingElements = 8; % Number of rolling elements in the bearing
ballDiameter = 0.235; % Diameter of the rolling elements (in meters)
pitchDiameter = 1.245; % Pitch diameter of the bearing (in meters)
contactAngle = 0; % Contact angle of the bearing (in degrees)
BPFI = (((rollingElements) * (shaftRate)) / 2) * (1 + ((ballDiameter/pitchDiameter) * cosd(contactAngle))); % Ball Pass Frequency Inner race

harmonics = (1:5) * shaftRate;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize arrays to store PSD values for averaging
baselinePSDs = [];
referenceFreq = [];  % Reference frequency vector

% Process Baseline Files to compute average PSD
for j = 1:length(baselineFiles)
    try
        % Load the .mat file
        filePath = fullfile(baselineFiles(j).folder, baselineFiles(j).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;

            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));
            
            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            if isempty(referenceFreq)
                referenceFreq = freq;
            else
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
            end
            
            baselinePSDs = [baselinePSDs, psdValues];
        else
            disp(['Baseline file ', baselineFiles(j).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing baseline file ', baselineFiles(j).name, ': ', ME.message]);
    end
end

% Average PSD for Baseline
averageBaselinePSD = mean(baselinePSDs, 2);

% Plotting the average Baseline PSD
figure;
hold on;
plot(referenceFreq, 10*log10(averageBaselinePSD), 'b', 'LineWidth', 1.5, 'DisplayName', 'Average Baseline PSD');

% Process Debris Files to compute average PSD of every three readings and plot each
debrisPlotCount = 1;
for k = 1:3:floor(length(debrisFiles)/3)*3  % Process all groups of 3, exclude leftover files
    try
        % Initialize temporary array for grouping PSDs
        tempDebrisPSDs = [];
        
        for j = 0:2
            % Load the .mat file
            filePath = fullfile(debrisFiles(k+j).folder, debrisFiles(k+j).name);
            dataStruct = load(filePath);
            
            % Extract the vibrational data and sample rate
            if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
                vibrationalData = dataStruct.bearing.gs;
                samplingRate = dataStruct.bearing.sr;
                
                % Remove NaN values
                vibrationalData = vibrationalData(~isnan(vibrationalData));

                % Compute Power Spectral Density (PSD)
                [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
                
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
                
                tempDebrisPSDs = [tempDebrisPSDs, psdValues];
            else
                disp(['Debris file ', debrisFiles(k+j).name, ' does not contain the expected structure. Skipping...']);
            end
        end
        
        % Average PSD for the current group of three debris readings
        if ~isempty(tempDebrisPSDs)
            avgTempPSD = mean(tempDebrisPSDs, 2);
            % Plot each average PSD of every three debris readings
            plot(referenceFreq, 10*log10(avgTempPSD), 'LineWidth', 0.5, 'DisplayName', ['Average Debris PSD (Group ', num2str(debrisPlotCount), ')']);
            debrisPlotCount = debrisPlotCount + 1;
        end
    catch ME
        disp(['Error processing debris files starting with ', debrisFiles(k).name, ': ', ME.message]);
    end
end

% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--k', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 0.5, 'LabelOrientation', 'horizontal');
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Average PSD Comparison: Baseline vs Debris (Grouped by 3)', 'FontSize', 14, 'FontWeight', 'bold');

% Add a legend
legend('show', 'Location', 'best', 'FontSize', 10);

% Set axis limits and grid
xlim([0 100]);  % Adjust this range as needed
ylim([-80 0]);  % Adjust Y limits based on your data
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

hold off;
%%

% Define the paths to the Baseline and Debris subfolders
baselineFolderPath = fullfile('Experimental_data', 'Baseline');
debrisFolderPath = fullfile('Experimental_data', 'Lubricant');

% Get a list of all MATLAB files in each subfolder
baselineFiles = dir(fullfile(baselineFolderPath, '*.mat'));
debrisFiles = dir(fullfile(debrisFolderPath, '*.mat'));

% Define the shaft rate and harmonics
shaftRate = 20.425;  % Shaft rate (Hz)
% Bearing Values
rollingElements = 8; % Number of rolling elements in the bearing
ballDiameter = 0.235; % Diameter of the rolling elements (in meters)
pitchDiameter = 1.245; % Pitch diameter of the bearing (in meters)
contactAngle = 0; % Contact angle of the bearing (in degrees)
BPFI = (((rollingElements) * (shaftRate)) / 2) * (1 + ((ballDiameter/pitchDiameter) * cosd(contactAngle))); % Ball Pass Frequency Inner race

harmonics = (1:5) * BPFI;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize arrays to store PSD values for averaging
baselinePSDs = [];
referenceFreq = [];  % Reference frequency vector

% Process Baseline Files to compute average PSD
for j = 1:length(baselineFiles)
    try
        % Load the .mat file
        filePath = fullfile(baselineFiles(j).folder, baselineFiles(j).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;

            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));
            
            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            if isempty(referenceFreq)
                referenceFreq = freq;
            else
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
            end
            
            baselinePSDs = [baselinePSDs, psdValues];
        else
            disp(['Baseline file ', baselineFiles(j).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing baseline file ', baselineFiles(j).name, ': ', ME.message]);
    end
end

% Average PSD for Baseline
averageBaselinePSD = mean(baselinePSDs, 2);

% Plotting the average Baseline PSD
figure;
hold on;


% Process each Debris File and plot against the baseline
for k = 1:length(debrisFiles)
    try
        % Load the .mat file
        filePath = fullfile(debrisFiles(k).folder, debrisFiles(k).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;
            
            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));

            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            % Interpolate PSD to match reference frequency vector
            psdValues = interp1(freq, psdValues, referenceFreq);
            
            % Plot each debris PSD individually
            plot(referenceFreq, 10*log10(psdValues), 'LineWidth', 0.01, 'DisplayName', ['Lubrication PSD (File ', num2str(k), ')']);
        else
            disp(['Debris file ', debrisFiles(k).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing debris file ', debrisFiles(k).name, ': ', ME.message]);
    end
end
plot(referenceFreq, 10*log10(averageBaselinePSD), 'b', 'LineWidth', 0.5, 'DisplayName', 'Average Baseline PSD');
% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--k', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('PSD Comparison: Baseline vs Individual Lubrication Files', 'FontSize', 14, 'FontWeight', 'bold');

% Add a legend
legend('show', 'Location', 'best', 'FontSize', 10);

% Set axis limits and grid
xlim([0 600]);  % Adjust this range as needed
ylim([-80 0]);  % Adjust Y limits based on your data
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

hold off;
%%
% Define the paths to the Baseline and Lubricant subfolders
baselineFolderPath = fullfile('Experimental_data', 'Debris');
lubricantFolderPath = fullfile('Experimental_data', 'Lubricant');

% Get a list of all MATLAB files in each subfolder
baselineFiles = dir(fullfile(baselineFolderPath, '*.mat'));
lubricantFiles = dir(fullfile(lubricantFolderPath, '*.mat'));

% Define the shaft rate and harmonics
shaftRate = 20.425;  % Shaft rate (Hz)
% Bearing Values
rollingElements = 8; % Number of rolling elements in the bearing
ballDiameter = 0.235; % Diameter of the rolling elements (in meters)
pitchDiameter = 1.245; % Pitch diameter of the bearing (in meters)
contactAngle = 0; % Contact angle of the bearing (in degrees)
BPFI = (((rollingElements) * (shaftRate)) / 2) * (1 + ((ballDiameter/pitchDiameter) * cosd(contactAngle))); % Ball Pass Frequency Inner race

harmonics = (1:5) * shaftRate;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize arrays to store PSD values for averaging
baselinePSDs = [];
lubricantPSDs = [];
referenceFreq = [];  % Reference frequency vector

% Process Baseline Files to compute average PSD
for j = 1:length(baselineFiles)
    try
        % Load the .mat file
        filePath = fullfile(baselineFiles(j).folder, baselineFiles(j).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;

            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));
            
            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            if isempty(referenceFreq)
                referenceFreq = freq;
            else
                % Interpolate PSD to match reference frequency vector
                psdValues = interp1(freq, psdValues, referenceFreq);
            end
            
            baselinePSDs = [baselinePSDs, psdValues];
        else
            disp(['Baseline file ', baselineFiles(j).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing baseline file ', baselineFiles(j).name, ': ', ME.message]);
    end
end

% Average PSD for Baseline
averageBaselinePSD = mean(baselinePSDs, 2);

% Process Lubricant Files to compute average PSD
for k = 1:length(lubricantFiles)
    try
        % Load the .mat file
        filePath = fullfile(lubricantFiles(k).folder, lubricantFiles(k).name);
        dataStruct = load(filePath);
        
        % Extract the vibrational data and sample rate
        if isfield(dataStruct, 'bearing') && isfield(dataStruct.bearing, 'gs') && isfield(dataStruct.bearing, 'sr')
            vibrationalData = dataStruct.bearing.gs;
            samplingRate = dataStruct.bearing.sr;
            
            % Remove NaN values
            vibrationalData = vibrationalData(~isnan(vibrationalData));

            % Compute Power Spectral Density (PSD)
            [psdValues, freq] = pwelch(vibrationalData, [], [], [], samplingRate);
            
            % Interpolate PSD to match reference frequency vector
            psdValues = interp1(freq, psdValues, referenceFreq);
            
            lubricantPSDs = [lubricantPSDs, psdValues];
        else
            disp(['Lubricant file ', lubricantFiles(k).name, ' does not contain the expected structure. Skipping...']);
        end
    catch ME
        disp(['Error processing lubricant file ', lubricantFiles(k).name, ': ', ME.message]);
    end
end

% Average PSD for Lubricant
averageLubricantPSD = mean(lubricantPSDs, 2);

% Plotting the average PSDs
figure;
hold on;

% Plot the average Baseline PSD
plot(referenceFreq, 10*log10(averageBaselinePSD), 'b', 'LineWidth', 0.5, 'DisplayName', 'Average Debris PSD');

% Plot the average Lubricant PSD
plot(referenceFreq, 10*log10(averageLubricantPSD), 'r', 'LineWidth', 0.5, 'DisplayName', 'Average Lubricant PSD');

% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--k', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Average PSD Comparison: Lubricant Vs Debris', 'FontSize', 14, 'FontWeight', 'bold');

% Add a legend
legend('show', 'Location', 'best', 'FontSize', 10);

% Set axis limits and grid
xlim([0     100]);  % Adjust this range as needed
ylim([-80 0]);  % Adjust Y limits based on your data
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

hold off;
%%
% Define the path to the Baseline and Imbalance subfolders
baselineFolderPath = fullfile('Experimental_data', 'Baseline');
imbalanceFolderPath = fullfile('Experimental_data', 'Imbalance');

% Get a list of all MATLAB files in each subfolder
baselineFiles = dir(fullfile(baselineFolderPath, '*.mat'));
imbalanceFiles = dir(fullfile(imbalanceFolderPath, '*.mat'));

% Define the critical rotational frequency
rotationalFreq = 20.425;  % Critical rotational frequency (Hz)

% Initialize variables to store results for imbalance and baseline
allImbalanceAmplitudes = [];
allBaselineAmplitudes = [];

% Process Baseline Files to compute amplitudes at the rotational frequency
for j = 1:length(baselineFiles)
    try
        % Load the .mat file
        filePath = fullfile(baselineFiles(j).folder, baselineFiles(j).name);
        dataStruct = load(filePath);

        % Extract the vibrational data and sampling rate
        vibrationalData = dataStruct.bearing.gs;
        samplingRate = dataStruct.bearing.sr;

        % Apply FFT to the vibrational data
        L = length(vibrationalData);
        f = samplingRate * (0:(L/2)) / L;
        Y = fft(vibrationalData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Find amplitude at the rotational frequency
        [~, idx] = min(abs(f - rotationalFreq));
        amplitudeAtRotFreq = P1(idx);

        % Store the amplitude
        allBaselineAmplitudes = [allBaselineAmplitudes; amplitudeAtRotFreq];
    catch ME
        disp(['Error processing baseline file ', baselineFiles(j).name, ': ', ME.message]);
    end
end

% Process Imbalance Files to compute amplitudes at the rotational frequency
for j = 1:length(imbalanceFiles)
    try
        % Load the .mat file
        filePath = fullfile(imbalanceFiles(j).folder, imbalanceFiles(j).name);
        dataStruct = load(filePath);

        % Extract the vibrational data and sampling rate
        vibrationalData = dataStruct.bearing.gs;
        samplingRate = dataStruct.bearing.sr;

        % Apply FFT to the vibrational data
        L = length(vibrationalData);
        f = samplingRate * (0:(L/2)) / L;
        Y = fft(vibrationalData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Find amplitude at the rotational frequency
        [~, idx] = min(abs(f - rotationalFreq));
        amplitudeAtRotFreq = P1(idx);

        % Store the amplitude
        allImbalanceAmplitudes = [allImbalanceAmplitudes; amplitudeAtRotFreq];
    catch ME
        disp(['Error processing imbalance file ', imbalanceFiles(j).name, ': ', ME.message]);
    end
end

% Sort and group the imbalance amplitudes in sets of three
allImbalanceAmplitudes = sort(allImbalanceAmplitudes);
groupedAmplitudes = [];
groupIndices = [];

for i = 1:3:length(allImbalanceAmplitudes)-2
    groupedAmplitudes = [groupedAmplitudes; mean(allImbalanceAmplitudes(i:i+2))];
    groupIndices = [groupIndices; i:i+2];
end

% Plot individual amplitudes and indicate grouping visually
figure;
hold on;

% Loop over the groups and plot with different markers/colors
colors = ['r', 'g', 'b']; % Define colors for the groups
markers = ['o', 's', 'd']; % Define markers for the groups

for k = 1:length(groupedAmplitudes)
    currentColor = colors(mod(k-1, length(colors)) + 1); % Cycle through colors
    currentMarker = markers(mod(k-1, length(markers)) + 1); % Cycle through markers
    
    % Plot each group of three points with a line connecting them
    plot(groupIndices(k,:), allImbalanceAmplitudes(groupIndices(k,:)), '-', 'Color', currentColor, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(groupIndices(k,:), allImbalanceAmplitudes(groupIndices(k,:)), currentMarker, 'MarkerSize', 8, 'LineWidth', 2, 'Color', currentColor, 'HandleVisibility', 'off'); 
end

% Plot the grouped averages with lines connecting them
plot(groupIndices(:,2), groupedAmplitudes, '-k', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Imbalance Group Average');

% Extend the baseline amplitudes to match the length of imbalance amplitudes
extendedBaselineAmplitudes = interp1(1:length(allBaselineAmplitudes), allBaselineAmplitudes, linspace(1, length(allBaselineAmplitudes), length(allImbalanceAmplitudes)), 'linear', 'extrap');
extendedBaselineAmplitudes = transpose(extendedBaselineAmplitudes);

% Plot the baseline amplitudes as a reference
plot(1:length(allImbalanceAmplitudes), extendedBaselineAmplitudes, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Baseline Amplitude');

% Add a trendline for the imbalance amplitudes
pImbalance = polyfit(1:length(allImbalanceAmplitudes), allImbalanceAmplitudes, 1);
trendlineImbalance = polyval(pImbalance, 1:length(allImbalanceAmplitudes));
plot(1:length(allImbalanceAmplitudes), trendlineImbalance, '--r', 'LineWidth', 2, 'DisplayName', 'Imbalance Trendline');

% Add a trendline for the baseline amplitudes
%pBaseline = polyfit(1:length(extendedBaselineAmplitudes), extendedBaselineAmplitudes, 1);
%trendlineBaseline = polyval(pBaseline, 1:length(allImbalanceAmplitudes));
%plot(1:length(allImbalanceAmplitudes), trendlineBaseline, '--b', 'LineWidth', 2, 'DisplayName', 'Baseline Trendline');

% Add a legend
legend('Imbalance Group Average', 'Baseline Amplitude', 'Imbalance Trendline', 'Baseline Trendline', 'Location', 'best', 'FontSize', 12);

% Add labels and title
xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude at Rotational Frequency (m/s^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Amplitude at Rotational Frequency (20.425 Hz) Across Imbalance vs Baseline Measurements', 'FontSize', 14, 'FontWeight', 'bold');

% Grid and axis properties
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 length(allImbalanceAmplitudes) + 1]);
ylim([0 max(max(allImbalanceAmplitudes), max(extendedBaselineAmplitudes)) * 1.1]); % Set y-axis limit slightly above max amplitude for better visibility

hold off;

%%

% Define the paths to the Angular Misalignment and Parallel Misalignment subfolders
angularFolderPath = fullfile('Experimental_data', 'Angular Misalignment');
parallelFolderPath = fullfile('Experimental_data', 'Parallel Misalignment');

% Get a list of all MATLAB files in each subfolder
angularFiles = dir(fullfile(angularFolderPath, '*.mat'));
parallelFiles = dir(fullfile(parallelFolderPath, '*.mat'));

% Initialize figure for plotting
figure;
hold on;

% Set line properties for fine lines
lineWidth = 0.5;

% Process Angular Misalignment Files
for j = 1:length(angularFiles)
    try
        % Load the .mat file
        filePath = fullfile(angularFiles(j).folder, angularFiles(j).name);
        dataStruct = load(filePath);

        % Extract the vibrational data and sampling rate
        vibrationalData = dataStruct.bearing.gs;
        samplingRate = dataStruct.bearing.sr;

        % Apply FFT to the vibrational data
        L = length(vibrationalData);
        f = samplingRate * (0:(L/2)) / L;
        Y = fft(vibrationalData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Plot the frequency amplitude plot
        plot(f, P1, 'r-', 'LineWidth', lineWidth, 'DisplayName', 'Angular Misalignment');
    catch ME
        disp(['Error processing angular misalignment file ', angularFiles(j).name, ': ', ME.message]);
    end
end

% Process Parallel Misalignment Files
for j = 1:length(parallelFiles)
    try
        % Load the .mat file
        filePath = fullfile(parallelFiles(j).folder, parallelFiles(j).name);
        dataStruct = load(filePath);

        % Extract the vibrational data and sampling rate
        vibrationalData = dataStruct.bearing.gs;
        samplingRate = dataStruct.bearing.sr;

        % Apply FFT to the vibrational data
        L = length(vibrationalData);
        f = samplingRate * (0:(L/2)) / L;
        Y = fft(vibrationalData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Plot the frequency amplitude plot
        plot(f, P1, 'b-', 'LineWidth', lineWidth, 'DisplayName', 'Parallel Misalignment');
    catch ME
        disp(['Error processing parallel misalignment file ', parallelFiles(j).name, ': ', ME.message]);
    end
end

% Add labels and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (m/s^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Frequency Amplitude Plots for Angular and Parallel Misalignment', 'FontSize', 14, 'FontWeight', 'bold');

% Set axis limits and grid
xlim([0 100]);  % Adjust this range as needed
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Combine the legend entries to avoid repetition
h = findobj(gca, 'Type', 'line');
legend([h(end), h(1)], {'Angular Misalignment', 'Parallel Misalignment'}, 'Location', 'best', 'FontSize', 10);

hold off;

%%
% Define the path to the Angular Misalignment subfolder
angularFolderPath = fullfile('Experimental_data', 'Parallel Misalignment');

% Get a list of all MATLAB files in the Angular Misalignment subfolder
angularFiles = dir(fullfile(angularFolderPath, '*.mat'));

% Define the shaft rate and harmonics
shaftRate = 20.425;  % Shaft rate (Hz)
harmonics = (1:5) * shaftRate;  % First 5 harmonics (1x, 2x, 3x, 4x, 5x)

% Initialize figure for plotting
figure;
hold on;

% Define a colormap and line styles for better differentiation
colors = lines(ceil(length(angularFiles)/3));  % Create a colormap with distinct colors
lineStyles = {'-', '--', ':', '-.'};  % Different line styles

% Set line properties for fine lines
lineWidth = 1;

% Process every third Angular Misalignment File
groupIndex = 1;
for j = 1:3:length(angularFiles)
    try
        % Load the .mat file
        filePath = fullfile(angularFiles(j).folder, angularFiles(j).name);
        dataStruct = load(filePath);

        % Extract the vibrational data and sampling rate
        vibrationalData = dataStruct.bearing.gs;
        samplingRate = dataStruct.bearing.sr;

        % Apply FFT to the vibrational data
        L = length(vibrationalData);
        f = samplingRate * (0:(L/2)) / L;
        Y = fft(vibrationalData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Select color and line style for this trial group
        color = colors(groupIndex, :);
        lineStyle = lineStyles{mod(groupIndex-1, length(lineStyles)) + 1};

        % Plot the frequency amplitude plot for this trial group
        plot(f, P1, 'Color', color, 'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
            'DisplayName', ['Trial Group ', num2str(groupIndex)]);

        % Increment the group index
        groupIndex = groupIndex + 1;
    catch ME
        disp(['Error processing angular misalignment file ', angularFiles(j).name, ': ', ME.message]);
    end
end

% Highlight the shaft rate frequency and harmonics
for i = 1:length(harmonics)
    xline(harmonics(i), '--k', ['Harmonic ', num2str(i), 'x (' num2str(harmonics(i), '%.2f'), ' Hz)'], ...
        'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
end

% Add labels with units and title
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (m/s^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Frequency Amplitude Plots for Parallel Misalignment', 'FontSize', 14, 'FontWeight', 'bold');

% Set axis limits and grid
xlim([0 100]);  % Adjust this range as needed
grid on;

% Set font size for the axes
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Add a legend to differentiate the trial groups
legend('show', 'Location', 'best', 'FontSize', 10);

hold off;





