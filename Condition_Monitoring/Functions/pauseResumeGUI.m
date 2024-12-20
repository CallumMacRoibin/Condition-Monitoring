function pauseResumeGUI()
    % Create the figure
    f = figure('Name', 'Training Control', 'NumberTitle', 'off', ...
        'Position', [500, 500, 300, 150], 'MenuBar', 'none', ...
        'CloseRequestFcn', @deleteFlagFilesAndClose);

    % Add a Pause button
    uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
        'Position', [20, 50, 80, 40], ...
        'Callback', @(~, ~) createPauseFlag());

    % Add a Resume button
    uicontrol('Style', 'pushbutton', 'String', 'Resume', ...
        'Position', [110, 50, 80, 40], ...
        'Callback', @(~, ~) createResumeFlag());

    % Add an Exit button
    uicontrol('Style', 'pushbutton', 'String', 'Exit', ...
        'Position', [200, 50, 80, 40], ...
        'Callback', @(~, ~) createExitFlag());

    function createPauseFlag()
        flagFile = 'pause_training.flag';
        if ~isfile(flagFile)
            fclose(fopen(flagFile, 'w'));
            disp('Pause flag created. Training will pause after the current iteration.');
        else
            disp('Pause flag already exists.');
        end
    end

    function createResumeFlag()
        pauseFile = 'pause_training.flag';
        resumeFile = 'resume_training.flag';

        if isfile(pauseFile)
            delete(pauseFile); % Remove pause flag
        end

        if ~isfile(resumeFile)
            fclose(fopen(resumeFile, 'w'));
            disp('Resume flag created. Training will resume.');
        else
            disp('Resume flag already exists.');
        end
    end

    function createExitFlag()
        flagFile = 'exit_training.flag';
        if ~isfile(flagFile)
            fclose(fopen(flagFile, 'w'));
            disp('Exit flag created. Training will terminate.');
        else
            disp('Exit flag already exists.');
        end
    end

    function deleteFlagFilesAndClose(~, ~)
        if isfile('pause_training.flag')
            delete('pause_training.flag');
        end
        if isfile('resume_training.flag')
            delete('resume_training.flag');
        end
        if isfile('exit_training.flag')
            delete('exit_training.flag');
        end
        disp('GUI closed and flag files cleaned up.');
        delete(f);
    end
end