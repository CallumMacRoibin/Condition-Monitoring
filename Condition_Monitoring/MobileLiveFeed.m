% MATLAB Script for Live Accelerometer Data (Ultimate Fix)
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
hPlotX = plot(nan, nan, 'r', 'DisplayName', 'X-axis'); hold on;
hPlotY = plot(nan, nan, 'g', 'DisplayName', 'Y-axis');
hPlotZ = plot(nan, nan, 'b', 'DisplayName', 'Z-axis');
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
xBuffer = nan(bufferSize, 1);
yBuffer = nan(bufferSize, 1);
zBuffer = nan(bufferSize, 1);

% Boolean to track if enough data has been accumulated
dataInitialized = false;

% Step 3: Continuous live plotting
while true
    % Get logged data from the mobile device
    [accelData, timeStamps] = accellog(m);

    % Check if data is being retrieved
    if ~isempty(accelData)
        % Append new data to the rolling buffer
        numNewPoints = length(timeStamps);
        timeBuffer = [timeBuffer(numNewPoints+1:end); timeStamps];
        xBuffer = [xBuffer(numNewPoints+1:end); accelData(:, 1)];
        yBuffer = [yBuffer(numNewPoints+1:end); accelData(:, 2)];
        zBuffer = [zBuffer(numNewPoints+1:end); accelData(:, 3)];

        % Trim the buffer to the last 'timeWindow' seconds
        currentTime = max(timeBuffer(~isnan(timeBuffer))); % Latest valid time
        validIdx = timeBuffer >= currentTime - timeWindow;

        % Ensure enough valid data exists before plotting
        if sum(validIdx) >= 10  % Plot only if at least 10 valid points exist
            % Extract valid data from the rolling buffer
            timeData = timeBuffer(validIdx);
            xData = xBuffer(validIdx);
            yData = yBuffer(validIdx);
            zData = zBuffer(validIdx);

            % Mark data as initialized
            dataInitialized = true;

            % Update the plot
            set(hPlotX, 'XData', timeData, 'YData', xData);
            set(hPlotY, 'XData', timeData, 'YData', yData);
            set(hPlotZ, 'XData', timeData, 'YData', zData);

            % Adjust plot limits
            xlim([currentTime - timeWindow, currentTime]);
            ylim([-20, 20]); % Adjust based on expected range
            drawnow;
        end
    else
        % Skip plotting if no new data is available
        disp('No new data available.');
    end

    % Clear the plot if data is not initialized to avoid lingering lines
    if ~dataInitialized
        set(hPlotX, 'XData', nan, 'YData', nan);
        set(hPlotY, 'XData', nan, 'YData', nan);
        set(hPlotZ, 'XData', nan, 'YData', nan);
    end

    % Pause to prevent CPU overuse
    pause(0.1);

    % Exit condition (press Ctrl+C to stop)
end

% Step 4: Stop logging when finished
m.Logging = 0;
disp('Logging stopped.');
