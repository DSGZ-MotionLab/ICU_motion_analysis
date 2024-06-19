function ICU_motion_analysis()
    % Matlab code for analyzing motion features based on recordings 
    % from a thigh-fixed motion sensor (triaxial accelerometer) in an ICU setting.
    % This function reads accelerometer data, processes it to extract various motion
    % features, and exports the results to a CSV file.

    % Export header for results
    exportHeader = {'ID', 'LegSide', 'RecordingDuration', 'CutDuration', 'IntensityGeneral', ...
                    'ActivityPercentage', 'ActiveBoutsPerHour', 'IntensityActiveLogMean', 'IntensityActiveVariability', ...
                    'DurationActiveLogMean', 'DurationActiveVariability'};
    
    % Axivity variable names and other constants
    axivityVariableNames = {'Time', 'AccZ', 'AccX', 'AccY'};
    legSides = {'left', 'right'};
    samplingRate = 10; % [Hz]
    windowLength = 5; %[s]
    activityThreshold = 0.135; % Signal magnitude area level in [g]    
    
    % Open a dialog box to select the data folder
    dataPath = uigetdir;
    files = dir(fullfile(dataPath, '*.csv')); % ID_legside.csv / ID_events.csv
    fileNames = {files.name};
    
    % Define the regular expression pattern to extract IDs    
    pattern = '^([^_]+)_.*\.csv$';
    tokens = regexp(fileNames, pattern, 'tokens', 'once');
    ids = [tokens{:}];
    ids = unique(ids(~cellfun('isempty', ids)));
    
    % Loop through records
    results = [];
    for i = 1:length(ids)
        id = ids{i};
        for j = 1:length(legSides)
            legSide = legSides{j};
            fileName = fullfile(dataPath, [id, '_', legSide, '.csv']);
            if isfile(fileName)
                eventFile = fullfile(dataPath, [id, '_events.csv']);
                disp(['Importing: ', fileName]);
                try
                    importedData = importRecording(fileName, samplingRate, axivityVariableNames);
                    disp(['Processing: ', fileName]);
                    result = processRecording(importedData, samplingRate, windowLength, activityThreshold, eventFile);
                    results = [results; [id, legSide, num2cell(result)]];
                catch ME
                    warning('Failed to process file %s: %s', fileName, ME.message);
                end
            end
        end
    end
    
    if ~isempty(results)
        resultsTable = cell2table(results, 'VariableNames', exportHeader);
        [file, path] = uiputfile('*.csv', 'Save as');
        if ischar(file)
            writetable(resultsTable, fullfile(path, file));
        else
            disp('Save cancelled');
        end
    end
end
    
%% Import recording function
function [axivityTimetable] = importRecording(fileName, Fs, axivityVariableNames)
    opts = detectImportOptions(fileName);
    opts.DataLines = [1 inf];
    opts.VariableNamingRule = 'preserve';
    axivityTable = readtable(fileName, opts);
    axivityTable.Properties.VariableNames = axivityVariableNames;
    axivityTimetable = table2timetable(axivityTable(:, 2:4), 'RowTimes', axivityTable.Time);
    axivityTimetable = retime(axivityTimetable, 'regular', 'linear', 'SampleRate', Fs);
end

%% Process recording function
function output = processRecording(input, Fs, winLen, activityThr, eventFile)
    data = input;
    % High-pass filter data
    Fc = 0.2;
    [b, a] = butter(4, Fc / (Fs / 2), 'high');
    data{:,:} = filter(b, a, data{:,:});
    
    % Remove periods indicated by events file
    totalDuration = height(data);
    available = 1;
    if isfile(eventFile)
        events = readtable(eventFile);
        startDate = events{1, 2};
        endDate = events{end, 3};
        data = data(timerange(startDate, endDate), :);
        for i = 2:height(events)-1
            data{timerange(events{i, 2}, events{i, 3}), :} = NaN;
        end
        available = 1 - sum(isnan(data{:,1})) / totalDuration;
    end
    cutDuration = sum(isnan(data{:,1}));
    
    % Compute signal magnitude area (SMA) for disjoint segments of length winLen
    timeVector = 0:1/Fs:seconds(data.Time(end) - data.Time(1));
    smaFunction = @(t, x, y, z) trapz(t, abs(x) + abs(y) + abs(z)) / winLen;
    SMA = matlab.tall.movingWindow(smaFunction, round(winLen * Fs), timeVector', data{:,1}, data{:,2}, data{:,3}, 'Stride', round(winLen * Fs), 'EndPoints', 'discard');
    
    recordingDuration = (totalDuration / Fs) / 3600; % [hours]
    cutDuration = (cutDuration / Fs) / 3600; % [hours]

    % Compute motion features
    % Overall motion intensity level
    intensityMean = nanmean(SMA);
    
    % Percentage of active behavior
    active = SMA > activityThr;
    activePercentage = 100 * sum(active) / (round(available * length(active)));

    % Number, Intensity and duration of activity bouts
    boutsArray = getBoutsArray(active, SMA);
    activeBoutsPerHour = length(boutsArray)/(available*recordingDuration);

    activeIntensity = cellfun(@nanmean, boutsArray);
    activeDuration = cellfun(@(x) size(x(:,1),1), boutsArray) * winLen; % [in s]
    
    % Compute mean and SD based on the log-transformed data
    [intensityActiveLogMean, intensityActiveVariability] = transformData(activeIntensity);
    [durationActiveLogMean, durationActiveVariability] = transformData(activeDuration);

    output = [recordingDuration, cutDuration, intensityMean, activePercentage, activeBoutsPerHour, intensityActiveLogMean, intensityActiveVariability, ...
              durationActiveLogMean, durationActiveVariability];
end

%% Helper function to get bouts array
function boutsArray = getBoutsArray(active, SMA)
    boutsArray = [];
    activityBool = [false active' ~= 0 false];
    starts = find(activityBool(2:end) & ~activityBool(1:end-1));
    ends = find(~activityBool(2:end) & activityBool(1:end-1)) - 1;
    boutsArray = [boutsArray, arrayfun(@(s, e) SMA(s:e), starts, ends, 'uniformoutput', false)];
end

%% Transform data function
function [logMean, logVariability] = transformData(input)
    if length(input) > 2
        output = mle(input, 'distribution', 'logn');
        logMean = output(1); % Mean of logarithmic values
        logVariability = output(2); % Standard deviation of logarithmic values
    else
        logMean = NaN;
        logVariability = NaN;
    end
end