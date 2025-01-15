clc;
clear;

% Step 1: Connect to the mobile device
m = mobiledev; % Create a mobiledev object
disp('Connecting to MATLAB Mobile...');
disp(m);

% Enable accelerometer sensor
m.AccelerationSensorEnabled = 1;

% Start logging data
m.Logging = 1;
disp('Logging accelerometer data. Move your phone to see live updates.');

% Step 2: Set up the live plot
figure('Name', 'Live Accelerometer Data', 'NumberTitle', 'off');
hPlotZ = plot(nan, nan, 'b', 'DisplayName', 'Z-axis'); % We will only plot Z-axis
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend('show');
title('Live Accelerometer Data');
grid on;

% Define the time window for display (10 seconds)
timeWindow = 10; % seconds

% Initialize empty buffers with fixed size
bufferSize = 100; % Rolling buffer size
timeBuffer = nan(bufferSize, 1);
zBuffer = nan(bufferSize, 1);

% Variable to store initial offset (for zero-centering)
zOffset = NaN;
offsetSamples = 20; % Number of samples to calculate initial offset

% Step 3: Continuous live plotting
while true
    % Get logged data from the mobile device
    [accelData, timeStamps] = accellog(m);

    % Check if data is being retrieved
    if ~isempty(accelData)
        % Extract only Z-axis acceleration
        zDataNew = accelData(:, 3);  % Z-axis data

        % Approximate time values assuming a 10 Hz sampling rate
        timeStamps = linspace(0, length(zDataNew) * 0.1, length(zDataNew))';

        % Calculate the initial offset from the first few samples
        if isnan(zOffset) && length(zDataNew) >= offsetSamples
            zOffset = mean(zDataNew(1:offsetSamples)); % Estimate rest position
        end

        % Center the Z-axis values around zero
        if ~isnan(zOffset)
            zDataNew = zDataNew - zOffset;
        end

        % Add new data to buffers
        numNewPoints = length(zDataNew);
        timeBuffer = [timeBuffer(numNewPoints+1:end); timeStamps];
        zBuffer = [zBuffer(numNewPoints+1:end); zDataNew];

        % Trim the buffer to the last 'timeWindow' seconds
        currentTime = max(timeBuffer(~isnan(timeBuffer))); % Latest valid time
        validIdx = timeBuffer >= currentTime - timeWindow;

        % Ensure enough valid data exists before plotting
        if sum(validIdx) >= 10  % Plot only if at least 10 valid points exist
            % Extract valid data from the rolling buffer
            timeData = timeBuffer(validIdx);
            zData = zBuffer(validIdx);

            % Update the plot
            set(hPlotZ, 'XData', timeData, 'YData', zData);

            % Adjust plot limits dynamically
            maxY = max(abs(zData));
            yRange = maxY*2;
            if yRange == 0
                yRange = 1; % Prevent zero range issues
            end
            ylim([(maxY*-1) - 0.2 * yRange, maxY + 0.2 * yRange]); % Add some padding

            % Adjust X-axis limits to keep time window
            xlim([currentTime - timeWindow, currentTime]);
            drawnow;
        end
    else
        % Skip plotting if no new data is available
        disp('No new data available.');
    end

    % Pause to prevent CPU overuse
    pause(0.1);
end

% Step 4: Stop logging when finished
m.Logging = 0;
disp('Logging stopped.');
