function dataprocessing(configfilename)
% DATAPROCESSING - Load TIFF data (movie, time, phprobe), locate particles,
%                  extract/save trajectories.
%
% User input area variables:
%   channelUsed         [cell array(strings)] specifies which emission channel to use, red/green/both
%   date                [strings] date of experiment and beginning of data file name
%   dataRootdir         [strings] root directory where data is located
%   fileNumbers         [cell array(strings)] stores file number
%   frames              [matrix] frame numbers for each data file where particles stop moving
%   localRootDir        [strings] local directory where data is sent
%   numberOfChannels    [scalar] how many emission channels there are
%   numFrames           [matrix] up to which frame number of each data file to analyze

% Code updates:
% -------------
% 
% 2009-10-24 - Fixed output of configuration file (extra .mat was printed)
%            - Changed passing in configfilename to require full filename
%              including .mat extension ('dataconfigfile.mat').
% 2009-08-25 - Fixed bug that detected max slope point even beyond the end
%              of the selected relevant data time window (1:numFrames).
% 2009-07-03 - Inserted help text at top of file.
%            - Added automatic time and ph extraction from TIFF file (saved
%              to trajectories datafile)
%            - Added ph probe detection (saved to trajectories datafile)
%            - Commented out unused variables localdir and rootdir
%            - Added saving of config data, and user-mode input, to allow a
%              single global dataprocessing.m to be called (instead of having
%              copies in each data top-level ("date") folder.


% date = '090615'; 
% fileNumbers = {'010'}; 
% dataRootDir = 'D:\Tijana\posao\steve harrison\results\matlab_data\'; 
% localRootDir = 'D:\Tijana\posao\steve harrison\results\matlab_data\';
% %movie1: frames = [238 248]; % ; 199 209; 240 250; 444 454; 261 271]; for particle detection; choose frames where particles stop moving
% frames = [200 215]; % movie1b; 199 209; 240 250; 444 454; 261 271]; for particle detection; choose frames where particles stop moving
% numFrames =[2500]; %2155];   % 2078 2060 2120 2030
% numberOfChannels = 2; %don't change
% channelUsed = {'red'}; % enter for EACH data set; specify which emission channel to use: 'green' 'red' or [] (=both channels)
%                        % 'green' 'red' or 'red/green'

%% PART 1) Settings to select data set(s) and options for particle detection and trajectory extraction from movie
% [20090703MP] Automate analysis process
if(nargin < 1) configfilename = []; end
%configfilename = 'datasetconfig.mat';
if(exist(configfilename))
    fprintf('Found dataset configuration file %s... Loading settings.\n', configfilename);
    load(configfilename);
else
    [configfilename,pathname] = uigetfile({['*.mat'],'MITHEMI dataset analysis configuration file (*.mat)'}, ...
                                    'Select dataset analysis configuration file to use, or click CANCEL to enter new settings');
    if(~(isequal(configfilename,0) || isequal(pathname,0)))   % User selected file, so load data...
        fprintf('Found dataset configuration file %s... Loading settings.\n', configfilename);
        load(configfilename);
    else                                                      % ... else user clicked cancel, manual entry:
        fprintf('No dataset configuration file found.  Enter settings for data import:\n');
        [~,dirname] = fileparts(pwd);
        date = input(['Date of data/start of filenames (string) [Default is current folder name: ' dirname '] : ']);
        if(isempty(date)) date = dirname; end
        fileNumbers = input('File numbers (as a cell array of three-digit strings, e.g. {{''001part1'',''001part2''},''002''}) : ');
        frames = input('Frame number range for particle detection where particles aren''t moving, [start end; start end; ...] : ');
        numFrames = input('How many frames to consider from start of movie? (allows clipping) : ');
        channelUsed = input('Which channels, for each file? {''red'',...} or {''green'',...} or {[],...}=both : ');

        configfilename = input('Enter name of dataset configuration file (without .mat extension) [datasetconfig] : ','s');
        if(isempty(configfilename)) configfilename = 'datasetconfig'; end
        fprintf('Saving settings to %s.mat.\n', configfilename);
        save(configfilename,'date','fileNumbers','frames','numFrames','channelUsed');  % [TI 090924] deleted redundunt variables
    end
end

%% PART 2) Get time/phprobe, find particles and extract trajectories from movie files
%variables:
%datadir = strings; name of data directory from where data is read
%files = cell array(strings); stores each data file name
%localdir = strings; name of local directory where data is sent

files = cell(1,size(fileNumbers,2)); %size command gives number of filenumbers; cell command makes a cell array containing that many cells
% %for n=1:size(fileNumbers,2)  %this loop generates the file names [TI 090218]
%nummovies = size(fileNumbers,2);             %[TI 090218]
nummovies = length(fileNumbers);             %[TI 090218]
for n=1:nummovies                            %this loop generates the file names [TI 090218]
    if ~iscell(fileNumbers{n})
        files{n}{1} = [date '_' fileNumbers{n}];   %concatanates the strings and puts it into one cell of the cell array 'files'
    else
        for kk = 1:length(fileNumbers{n})
            files{n}{kk} = [date '_' fileNumbers{n}{kk}];
        end
    end
%    localdir = [localRootDir date];         %concatanates to give name of local directory  % [20090703MP] Commented out unused variables
%    datadir = [dataRootDir date];           %concatanates to give name of data directory  % [20090703MP] Commented out unused variables
end
clear date fileNumbers %dataRootDir localRootDir    % [20090703MP] Commented out unused variables
%%
%variables:
%nummovies = scalar; number of data files to read
%ph = matrix; stores the time stamps from the fluorescein pH trajectory
%centroid, avMovie, filtImage, trajectories

for n = 1:nummovies
    tic

    fprintf('Getting time/pH trajectory from movie %d of %d...', n, nummovies);  % [20090703MP] Uncommented, and added auto time extraction from TIFF
% Make a matrix storing the time stamps of the ph trajectory
    %cd(files{n})
    %ph=dlmread('ph.xld','\t',5,1);

    % [20090703MP] Get time and ph data from the TIFF files
%    cd(localdir)           % [20090703MP] Comment out but pay attention to relative paths
%    datafilename = [files{n} '.tif'];  [timems,info,IsHamamatsuMovie] = getandortiffinfo(datafilename);  % Get time from TIFF
    timems = [];
    for kk = 1:length(files{n})
        datafilename{kk} = [files{n}{kk} '.tif'];
        [temptime,info{kk},IsHamamatsuMovie] = getandortiffinfo(datafilename{kk});  % Get time from TIFF
        timems = [timems; temptime];
        movielengths(kk) = length(info{kk});        % Keep track of # of frames in each movie
    end
    if(sum(timems == 0) > 3)    % If getandor gets more than a couple of zeros in time vector then Hamamatsu HCImage didn't save time in the TIF.
        fprintf('getandortiffinfo: There is no valid time data saved in TIF, attempting to read from text file %s.\n', [files{n}{1} '.txt']);
        timems = dlmread([files{n}{1} '.txt'], '\t', 1, 1) * 1000;
    end
    if(IsHamamatsuMovie)
        fprintf('Detected interleaved-channels movie data (Hamamatsu setup).\n');
        FiltNoiseSize = 1.8;
%        if(strcmp(channelsUsed{n},'green'))  IsHamamatsuMovie = 0;  end  % NOTE: fix for green only should not be interleaved
        IsInterleavedMovie = input('Is the Hamamatsu movie interleaved? [Y/n] : ');
        if(isempty(IsInterleavedMovie) | upper(IsInterleavedMovie)=='y')
            [idxgreenframes, idxredframes] = parseinterleavedmovieframes(datafilename, 1, numFrames(n), 0, movielengths); % Get green/red channel indices, don't save separate TIFs
        else
            idxgreenframes = [1:numFrames(n)]; idxredframes = idxgreenframes;   % Split-frame movies use every frame for both red and green.
        end
    else
        fprintf('Detected split-frame channels movie data (Andor setup).\n');
 %       FiltNoiseSize = 1.0; % [100113TI] playing with filter noise size for low power movies; original setting
         FiltNoiseSize = 1.4; % [100113TI] playing with filter noise size for low power movies
        idxgreenframes = [1:numFrames(n)]; idxredframes = idxgreenframes;   % Split-frame movies use every frame for both red and green.
    end

% Find particles in the movies
% %    cd(localdir)           % [20090703MP] Comment out but pay attention to relative paths
%     if isdir(files{n})==0
%         mkdir(files{n});
    if isdir(files{n}{1})==0
        mkdir(files{n}{1});
    end
%     fprintf(['Finding particles in ' files{n} '.tif ... ']);
    fprintf('Finding particles in '); fprintf('%s.tif ', files{n}{:}); fprintf( '...\n');

    if(~IsHamamatsuMovie)
         redbox = [11 246; 5 502]; greenbox = [267 502; 5 502]; % [100113TI] changing box size for two-color analysis
%        redbox = [11 246; 5 502]; greenbox = [256 502; 5 502]; % [100113TI] this was the last setting 
%        redbox = [1 256; 1 512]; greenbox = [257 512; 1 512];
    else
         redbox = [1 512; 1 512]; greenbox = redbox;
    end

%    switch channelUsed{n}
%        case []
    if isempty(channelUsed{n})
            framesets = {idxredframes,idxgreenframes};
            boxsets   = {redbox,greenbox};
            trajname  = {'redTraj','greenTraj'};
    elseif strcmp(channelUsed{n},'red')
%        case 'red'
            framesets = {idxredframes};
            boxsets   = {redbox};
            trajname  = {'redTraj'};
    elseif strcmp(channelUsed{n},'green')
%        case 'green'
            framesets = {idxgreenframes};
            boxsets   = {greenbox};
            trajname  = {'greenTraj'};
    end

    for kk = 1:length(framesets)
        framesetsforsearch{kk} = find((framesets{kk} >= frames(n,1)) & (framesets{kk} <= frames(n,2)));
    end
%     idxredframesforsearch = find((idxredframes >= frames(n,1)) & (idxredframes <= frames(n,2)));       % Find red indices in frame range for particle detection
%     idxgreenframesforsearch = find((idxgreenframes >= frames(n,1)) & (idxgreenframes <= frames(n,2))); % Find green indices in frame range for particle detection

%     [centroid avMovie filtImage] = findparticle([files{n} '.tif'], ...
%                                                 'Frames', framesets, 'Box', boxsets, 'FiltNoiseSize', FiltNoiseSize);
    [centroid avMovie filtImage] = findparticle(movielengths, datafilename, ...
                                                'Frames', framesets, 'Box', boxsets, 'FiltNoiseSize', FiltNoiseSize);
    fprintf('done.\n');

% Extract fluorescence intensity trajectories for each particle   
%    fprintf(['Extracting temporal intensity trajectories of detected particles from ' files{n} '.tif ...']);
    fprintf('Extracting temporal intensity trajectories of detected particles from '); fprintf('%s.tif ', files{n}{:}); fprintf( '...\n');
%    trajectories = extractdata([files{n} '.tif'], centroid, numFrames(n));
%    temptrajectories = extractdata([files{n} '.tif'], centroid, framesets);
    temptrajectories = extractdata(movielengths, datafilename, centroid, framesets);
    for kk = 1:length(temptrajectories)
        eval(['trajectories.' trajname{kk} ' = temptrajectories{kk};']);    % Remember, propertyName is a string
    end
    % NEED to set/deal with trajectories.redTraj, trajectories.greenTraj, etc.
    fprintf('done.\n');
    processTime = toc;

% Extract time and phprobe from TIF, and detect ph drop frame.
% %    phprobe = getphprobefromtiff(datafilename, numFrames(n));   % Get ph from TIFF
%     phprobe = getphprobefromtiff(datafilename, idxgreenframes); %, greenbox);   % Get ph from (first) TIFF
%    idxgreenframesmovie1 = idxgreenframes(find(idxgreenframes < idxgreenframes( min(numFrames(n),movielengths(1)) )));  % Get frame indices (for green channel) from only the first moviepart if this is a multi-part movie (either up to numFrames, or whole movie if numFrames is in 2nd, etc part).
    idxgreenframesmovie1 = idxgreenframes(find(idxgreenframes < min(numFrames(n),movielengths(1)) ));  % Get frame indices (for green channel) from only the first moviepart if this is a multi-part movie (either up to numFrames, or whole movie if numFrames is in 2nd, etc part).
%     phprobe = getphprobefromtiff(datafilename{1}, idxgreenframesmovie1); %, greenbox);   % Get ph from (first) TIFF
phprobe = getphprobefromtiff(datafilename{1}, idxgreenframesmovie1); %, greenbox);   % Get ph from (first) TIFF; second version (getphprobefromtiff2) does not subtract small background boxes or open iris experiments
%    nphframe = detect_ph_drop(timems(idxgreenframesmovie1), phprobe);    % Detect the ph drop frame number;  [20090825] Fixed bug where max slope was detected beyond "end of movie" after numFrames.
fGaussFit = fittype('yofs + h*erfc((x-xo)/w)/2 * exp(-m*(x-xo))','coefficients',{'h','m','w','xo','yofs'},'independent','x');    
timesec = timems/1000;
% yfit = fit(timesec(1:80), phprobe(1:80), fGaussFit, 'StartPoint', [(max(phprobe) - min(phprobe)) 0 2 10 min(phprobe)], ...
%                                             'Lower',[-inf 0 -inf -inf -inf]);
yfit = fit(timesec(1:40), phprobe(1:40), fGaussFit, 'StartPoint', [50 0 2 10 min(phprobe)], ...
                                                    'Lower',[-inf 0 -inf -inf -inf]);       

tdrop = yfit.xo + yfit.w;
wdrop = 2*yfit.w;

y=abs(tdrop-timesec);
nphframe(2) = find(y==min(y));

figure;
plot(timesec(1:40), yfit(timesec(1:40)), 'k');
hold on;
plot(timesec(1:40), phprobe(1:40), 'b');
hold on;
plot(tdrop, yfit(tdrop), 'ro');
xlabel('time (sec)', 'FontSize', 14); ylabel('fluorescence intensity (AU)', 'FontSize', 14);
title(['pHdrop time/frame is ' num2str(tdrop) ' / ' num2str(nphframe(2)) ', pHdrop width is ' num2str(wdrop)]);

% nphframe(2) =20;
% phprobe figure to verify manual ph drop time detection    
%    figure
%    timems_vec = timems(idxgreenframesmovie1)/1000;
%    plot(timems_vec(1:59), phprobe(1:59))
%    xlabel('time (sec)', 'FontSize', 14); ylabel('fluorescence intensity (AU)', 'FontSize', 14);
% Save the trajectories and other data
    cd(files{n}{1})
    save([files{n}{1} 'trajectories'],'-struct','trajectories')
    save([files{n}{1} 'trajectories'],'centroid','avMovie','filtImage','timems','phprobe','nphframe','-append');   % [20090703MP] 
    save([files{n}{1} 'trajectories'],'idxredframes','idxgreenframes','-append');
    cd('..')    % [20090703MP] Return to movie set ("root") folder after saving data
end
clear datafilename info
clear nummovies channelUsed
%clear localdir datadir     % [20090703MP] Not used
fprintf('Done: %d/%d/%d - %02d:%02d:%02d.\n', fix(clock));
