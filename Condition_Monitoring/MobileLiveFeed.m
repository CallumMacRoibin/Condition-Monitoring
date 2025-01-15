% MATLAB Script for Live Accelerometer Data (Z-Axis Only with Offset and Symmetric Y-Limits)
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
hPlotZ = plot(nan, nan, 'b', 'DisplayName', 'Z-axis'); % Plot for Z-axis
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend('show');
title('Live Accelerometer Data (Z-Axis)');
grid on;

% Define the time window for display (10 seconds)
timeWindow = 10; % seconds

% Initialize buffers and offset
bufferSize = 100; % Rolling buffer size
timeBuffer = []; % Start with an empty buffer for time
zBuffer = []; % Start with an empty buffer for Z-axis data
zOffset = NaN; % Offset for zero-centering
offsetSamples = 20; % Number of samples to calculate initial offset
bufferFull = false; % Flag to indicate if the buffer is full

% Step 3: Continuous live plotting
while true
    % Get logged data from the mobile device
    [accelData, timeStamps] = accellog(m);

    % Check if data is being retrieved
    if ~isempty(accelData)
        % Extract Z-axis data
        zDataNew = accelData(:, 3);

        % Calculate the initial offset if not already done
        if isnan(zOffset) && length(zDataNew) >= offsetSamples
            zOffset = mean(zDataNew(1:offsetSamples)); % Estimate the offset using the first few samples
        end

        % Apply offset to center Z-axis data around zero
        if ~isnan(zOffset)
            zDataNew = zDataNew - zOffset;
        end

        % If the buffer is not full, plot directly without using the buffer
        if ~bufferFull
            timeBuffer = [timeBuffer; timeStamps]; % Append to time buffer
            zBuffer = [zBuffer; zDataNew]; % Append to Z-axis buffer

            % Plot directly if there are valid timestamps
            if ~isempty(timeStamps) && all(~isnan(timeStamps)) && length(timeStamps) > 1
                % Update the plot with the new data
                set(hPlotZ, 'XData', timeStamps, 'YData', zDataNew);

                % Adjust Y-axis limits with padding
                maxY = max(abs(zDataNew));
                if maxY == 0
                    maxY = 1; % Prevent zero range
                end
                padding = 0.1 * maxY; % 10% padding
                ylim([-maxY - padding, maxY + padding]);

                % Adjust X-axis limits if valid
                xlim([min(timeStamps), max(timeStamps)]);
                drawnow;
            end

            % Check if the buffer is now full
            if length(zBuffer) >= bufferSize
                bufferFull = true; % Switch to using the rolling buffer
                timeBuffer = timeBuffer(end-bufferSize+1:end); % Trim to buffer size
                zBuffer = zBuffer(end-bufferSize+1:end); % Trim to buffer size
            end
        else
            % Append new data to the rolling buffer
            numNewPoints = length(timeStamps);
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

                % Adjust plot limits with padding
                maxY = max(abs(zData)); % Determine the maximum absolute value
                if maxY == 0
                    maxY = 1; % Prevent zero range
                end
                padding = 0.1 * maxY; % Add 10% padding to the limits
                ylim([-maxY - padding, maxY + padding]); % Symmetric Y-limits with padding
                xlim([currentTime - timeWindow, currentTime]); % X-axis limits based on time window
                drawnow;
            end
        end
    else
        % Skip plotting if no new data is available
        disp('No new data available.');
    end

    % Pause to prevent CPU overuse
    pause(0.1);

    % Exit condition (press Ctrl+C to stop)
end

% Step 4: Stop logging when finished
m.Logging = 0;
disp('Logging stopped.');
