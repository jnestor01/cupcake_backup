function cupcake

% cupcakePilot
% ___________________________________________________________________
% 
% Target estimation and uncertainty task
% ___________________________________________________________________

% HISTORY
% 05/13/21  kjt, wrote, adapted from rd cupcake-aperture. 
% 05/20/21  kjt, add mouse functions
% 06/02/21  kjt, add timing 
% 07/08/21  kjt, contrast manipulation
% 04/12/24  kjt, resurrection
% 06/27/24  kjt, eeg pilot
   
%% Input 
% subject 
subjectID = input('Enter subject ID:  ','s');
testingLocation = input('Enter testing location (desk, DenBehav):  ','s');
p.debug = 0; debug = 0; 
practice = input(['Practice phase: 0/skip-real run \n1-high contrast, free eye movements, neutral cues' ...
    '\n2-high contrast, enforce fixation, neutral cues \n3-high contrast, enforce fixation, both cue types' ...
    '\n4-low contrast, enforce fixation, both cue types:']);

if ~practice
    showInstruct = input('Show instructions? 1-yes 0-no:  ');
    run = input('Enter run #:  ');
    run = sprintf('%02d', run);
else
    run = sprintf('practice%02d', practice);
    showInstruct = 1;
end

showPositionFeedback = 1; % input('Show position feedback? 1-yes 0-no:  ');
useTrackball = 0; % input('Mouse or trackball? 1-trackball 0-mouse:  ');
resumeDataset = 0; % input('Resume old dataset? 1-yes 0-no:  ');

% Save options. 
saveData = 1;

%% Directories 
% add paths
directory = pwd; % get project directory path 
addpath(genpath(directory))
addpath(genpath('../helper-functions'))
% subject directories 
subjectDataDir = sprintf('%s/data/%s',directory,subjectID); 
if ~exist(subjectDataDir,'dir')
    mkdir(subjectDataDir); 
end

%% If resuming dataset, load data
if resumeDataset
    origDir = cd(subjectDataDir);
    
    % open user interface to load incomplete data set
    [filename, pathname] = uigetfile('*.mat', 'Select the data set to resume');
    
    % change back to expt directory
    cd(directory);
    
    % load incomplete data
    load([pathname filename]);
end

%% Setup 
if ~resumeDataset
    % Get parameters
    p = cupcake_params(testingLocation, debug, showInstruct, showPositionFeedback);
    useTrackball = logical(useTrackball);
    
    % Running on PTB-3? Abort otherwise.
    AssertOpenGL;
    
    if strcmp(p.testingLocation, 'desk')
        Screen('Preference', 'SkipSyncTests', 1);
    end
    if strcmp(subjectID,'test')
        p.eyeTracking = 0;
    end
end

if practice == 1
    p.gratingContrasts = 0.8;
    p.enforceFixITI = 0;
    p.eyeTracking = 0;
    p.triggerEEG = 0;
    p.photodiode = 0;
    p.attConds = 0;
elseif practice == 2
    p.gratingContrasts = 0.8;
    p.enforceFixITI = 1;
    p.enforceFixCTI = 1;
    p.eyeTracking = 1;    
    p.triggerEEG = 0;
    p.photodiode = 0;
    p.attConds = 0;
elseif practice == 3
    p.gratingContrasts = 0.8;
    p.enforceFixITI = 1;
    p.enforceFixCTI = 1;
    p.eyeTracking = 1;    
    p.triggerEEG = 0;
    p.photodiode = 0;
elseif practice == 4
    p.gratingContrasts = p.gratingContrasts(1);
    p.enforceFixITI = 1;
    p.enforceFixCTI = 1;
    p.eyeTracking = 1;    
    p.triggerEEG = 0;
    p.photodiode = 0;    
end

if ~practice
%% Check for existing threshold file
thresholdFile = sprintf('data/%s/%s_%s_threshold_%s.mat', subjectID, subjectID, p.expName);
if isfile(thresholdFile)
    threshRun = load(thresholdFile);
    threshContrast = 10^(mean(QuestMean(threshRun.expt.q(1:threshRun.expt.p.nStaircases))));
    useThresh = input(sprintf('\nFound threshold of %3.2f in file \n%s \nUse this threshold? 1-yes, skip/0-no', ...
        thresholdFile, threshContrast));
    if useThresh
        p.gratingContrasts = threshContrast;
    else
        fprintf('Using default contrast level: %3.2f', p.gratingContrasts);
    end
else
    customContrast = input(sprintf(['\nNo threshold file found under %s ' ...
        '\nEnter threshold contrast level (absolute or log10) or skip/0 for default contrast: %3.2f'], thresholdFile, p.gratingContrasts));
    if customContrast
        if customContrast>0
            p.gratingContrasts = customContrast;
        else
            p.gratingContrasts = 10^customContrast;
        end
    end
end
end

%% Initialize stim tracker for MEG
if p.triggersOn
    PTBInitStimTracker;
    global PTBTriggerLength
    PTBTriggerLength = 0.001;
    
    % to send trigger messages to EyeLink
    global PTBEyeTrackerRecording
    if p.eyeTracking
        PTBEyeTrackerRecording = 1;
    else
        PTBEyeTrackerRecording = 0;
    end
end

%% Set up parallel port
if p.triggerEEG
    ppdev_mex('CloseAll');
    ppdev_mex('Open',1);
    ppdev_mex('Write',1,255); % experiment begins
    WaitSecs(0.05);
    ppdev_mex('Write',1,0);
end

%% Eye data i/o
eyeRad = p.eyeRad*p.ppd; 
eyeDataDir = 'data/eyedata';
eyeFile = sprintf('%s_%s', subjectID, run);
if numel(eyeFile)>8
    eyeFile = eyeFile(1:8);
end
eyeFileFull = sprintf('%s/%s_%s_%s.edf', eyeDataDir, subjectID, p.expName, datestr(now, 'yymmdd'));

% Check to see if this eye file already exists
existingEyeFile = dir(sprintf('%s/%s.edf', eyeDataDir, eyeFile));
if ~isempty(existingEyeFile) && p.eyeTracking
    error('eye file already exists! please choose another name.')
end

%% Display key settings to the experimenter
fprintf('\nExperiment settings:\n')
fprintf('subject = %s\n', subjectID)
fprintf('run = %s\n', run)
fprintf('location = %s\n', p.testingLocation)
fprintf('trackball = %d\n', useTrackball)
fprintf('debugging = %d\n', p.debug)
fprintf('aperture = %s\n', p.aperture)

%% Screen
% Set up screen
screenNumber = max(Screen('Screens'));

% Check screen resolution and refresh rate 
scr = Screen('Resolution', screenNumber); % 
if ~all([scr.width scr.height scr.hz] == [p.screenRes p.refRate]) && ~strcmp(subjectID,'test')
    % error('Screen resolution and/or refresh rate has not been set correctly by the experimenter!')
end

% Set up window
multisample = 8; % draws at higher resolution than screen, then downsampled 
winFrac = p.winFrac * [scr.width,scr.height];  

%% Open screen
% Screen('Preference', 'SkipSyncTests', 1);
if strcmp(p.testingLocation, 'desk')
    [window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 winFrac], [],[],[],multisample); 
else
    [window, rect] = Screen('OpenWindow', screenNumber, [], [], [],[],[],multisample); 
    HideCursor(screenNumber); 
end

white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
black = BlackIndex(window);
p.photodiodeColors = [white;black]; 
[cx, cy] = RectCenter(rect); 
p.cx = cx; p.cy = cy;
Screen('TextSize', window, p.fontSize);
Screen('TextColor', window, white);
Screen('TextFont', window, p.font);
Screen('FillRect', window, white*p.backgroundColor);

% Check screen size
[sw, sh] = Screen('WindowSize', window); % height and width of screen (px)
if ~all([sw sh] == p.screenRes) && ~strcmp(subjectID,'test')
    error('Screen resolution is different from requested!')
end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)
if abs(flipInterval - (1/p.refRate)) > 0.001 && ~strcmp(subjectID,'test')
    error('Refresh rate is different from requested!')
end

% Check font
if ~strcmp(p.font, Screen('TextFont', window))
    error('Font was not set to requested: %s', p.font)
end

%% Keyboard
devices = PsychHID('devices');
% set button box device
if p.useKbQueue
    devNum = [];
    for iD=1:numel(devices)
        if strcmp(devices(iD).usageName,'Keyboard') && ...
                strcmp(devices(iD).product,'904')
            devNum = iD;
        end
    end
    if isempty(devNum)
        error('Did not find button box')
    end
else
    devNum = -1;
end

% Set up KbQueue if desired
if p.useKbQueue
    KbQueueCreate(devNum);
    KbQueueStart();
end

% Set up mouse 
m = GetMouseIndices;
m = m(1);
KbQueueCreate(m);

%% Make stimuli
if ~resumeDataset
    % Calculate stimulus dimensions (px) and position
    imPos = round(p.imPos*p.ppd);
    gratingRadius = round(p.gratingDiameter/2*p.ppd);
    % edgeWidth = round(p.apertureEdgeWidth*p.ppd);
    responseRadius = round(p.responseDotDiameter/2*p.ppd); % response probe 
    
    % Make gratings
    image = zeros(round(p.ppd* (p.imSize(1)*2) ));
    for iC = 1:numel(p.gratingContrasts)
        for iSF = 1:numel(p.edgewidth_pix)
            contrast = p.gratingContrasts(iC);
            edgeWidth = p.edgewidth_pix(iSF); % in pixels, controls spatial frequency
            % Make radial grating
            [~, aps] = rd_aperture(image, p.aperture, gratingRadius*2, edgeWidth);
            aps = rescale(aps, 0.5-0.5*contrast, 0.5+0.5*contrast);
            % aps(aps<0.5) = 0.5 - aps(aps<0.5)*contrast;
            % aps(aps>0.5) = 0.5 + aps(aps>0.5)*contrast;
            % Aperture the radial grating
            % aps = rd_aperture(aps, 'cosine', gratingRadius(1), edgeWidth);
            aps = rd_aperture(aps, 'gaussian', p.gaussianSD*p.ppd, edgeWidth); % change this, gratingRadius(1)/2.5
            % save gratings
            grating(:,:,iC,iSF) = aps;
            % Make texture
            tex{iC,iSF} = Screen('MakeTexture', window, aps*white);
        end
    end
    
    % Make the rects for placing the images
    % Target rect
    imSize = size(image);
    imRect = CenterRectOnPoint([0 0 imSize], cx+imPos(1), cy+imPos(2));
    % Response rect
    responseRect = CenterRectOnPoint([0 0 responseRadius*2 responseRadius*2], cx, cy);
    % Oval rect
    ecc = p.ecc * p.ppd; % eccentricity of circle (deg)
    ovalRect = [0 0 2*ecc 2*ecc];
    centeredOvalRect = CenterRectOnPointd(ovalRect, cx, cy);
    
    buffer = 1*p.ppd;
    ovalRectInner = [0 0 2*ecc-imSize(1)-buffer 2*ecc-imSize(2)-buffer];
    ovalRectOuter = [0 0 2*ecc+imSize(1)+buffer 2*ecc+imSize(2)+buffer];
    centeredOvalRectInner = CenterRectOnPointd(ovalRectInner, cx, cy);
    centeredOvalRectOuter = CenterRectOnPointd(ovalRectOuter, cx, cy);
    % Arc rect (oval + half arc stroke width)
    arcRect = [0 0 2*(ecc+p.arcStroke/2) 2*(ecc+p.arcStroke/2)];
    centeredArcRect = CenterRectOnPointd(arcRect, cx, cy);
    % Att cue rect
    attRect = [0 0 p.ppd/(p.cueLength+p.cueBuffer) p.ppd/(p.cueLength+p.cueBuffer)];
    centeredAttRect = CenterRectOnPointd(attRect, cx, cy);
    
    % Calculate fixation diameter in pixels
    fixSize = round(p.fixDiameter*p.ppd);
end

%% Generate trials
if ~resumeDataset
    % Construct trials matrix
    trials_headers = {'angularPos','contrast','att','iti','fixTargetISI'};
    
    % make sure column indices match trials headers
    angularPosIdx = strcmp(trials_headers,'angularPos');
    contrastIdx = strcmp(trials_headers,'contrast');
    attIdx = strcmp(trials_headers,'att');
    itiIdx = strcmp(trials_headers,'iti');
    fixTargetISIIdx = strcmp(trials_headers,'fixTargetISI');
    
    % full factorial design
    % target position, contrast, and attention condition are fully
    % counterbalanced
    % ITI and fixTargetISI are pseudo-counterbalanced
    trials = fullfact([numel(p.theta) ...
        numel(p.gratingContrasts)...
        numel(p.attConds)]);
    
    % repeat trials matrix according to nReps of all conditions
    trials = repmat(trials, p.nReps, 1);
    nTrials = size(trials,1);
    
    % show trials and blocks
    fprintf('\n%s\n\n%d trials, %1.2f blocks, %d trials per block\n\n', datestr(now), nTrials, nTrials/p.nTrialsPerBlock, p.nTrialsPerBlock)
    
    % Generate ITIs
    switch p.itiType
        case 'uniform'
            nITIs = numel(p.itis);
            itis = repmat(1:nITIs,1,ceil(nTrials/nITIs));
            trials(:,itiIdx) = itis(randperm(nTrials));
        case 'hazard' % constant hazard rate
            trials(:,itiIdx) = rd_sampleDiscretePDF(p.itiPDF, nTrials);
        otherwise
            error('p.itiType not recognized')
    end
    
    % Generate fixation target ISIs for prestim alpha analysis
    nfixTargetISIs = numel(p.fixTargetISIs);
    fixTargetISIs = repmat(1:nfixTargetISIs,1,ceil(nTrials/nfixTargetISIs));
    trials(:,fixTargetISIIdx) = fixTargetISIs(randperm( size(fixTargetISIs,2), nTrials));
    
    %% Choose order of trial presentation
    % trialOrder = randperm(nTrials);
    % This randomizes trial order within reps, but not across reps. So we will
    % present every trial in one rep before moving on to the next rep. In this
    % way the experiment is blocked into complete reps (though the number of
    % trials per block may be different then the number of trials in a rep).
    nt = nTrials/p.nReps;
    trialSet = ones(nt,1)*(1:p.nReps);
    repOrders = [];
    for i=1:p.nReps
        repOrders = [repOrders randperm(nt)'];
    end
    trialOrder0 = repOrders + (trialSet-1).*nt;
    trialOrder = trialOrder0(:);
    
    % shuffle trials
    trials = trials(trialOrder,:);
    
    %% Pregenerate cue offsets
    % cueThetaRnd = normrnd(0, 7, [1,nTrials]); 
end

%% Eyetracker
if p.eyeTracking
    % Initialize eye tracker
    [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end
    
    % Write subject ID into the edf file
    Eyelink('message', 'BEGIN DESCRIPTIONS');
    Eyelink('message', 'Subject code: %s', subjectID);
    Eyelink('message', 'Run: %s', run);
    Eyelink('message', 'END DESCRIPTIONS');
    
    % No sounds indicating success of calibration
    el.drift_correction_target_beep = [0 0 0];
    el.drift_correction_failed_beep = [0 0 0];
    el.drift_correction_success_beep = [0 0 0];

    % Accept input from all keyboards
    el.devicenumber = -1; %

    % Update with custom settings
    EyelinkUpdateDefaults(el);
    
    % Calibrate eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
    
    if p.enforceFixITI||p.enforceFixCTI
        eyeRad = round(ang2pix(p.eyeRad, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
        fixRect = [cx-eyeRad, cy-eyeRad, cx+eyeRad, cy+eyeRad]; 
        % Start recording
        rd_eyeLink('startrecording', window, {el, fixRect});
    else
        % Start recording
        rd_eyeLink('startrecording', window, el);
    end
end

%% Instructions and example target
if p.showInstruct 
    

    upperTextLocation = 1.9;
    lowerTextLocation = 1.75;
    instrRatio = p.instrRatio;
    pInstr = p;
    pInstr.ecc = pInstr.ecc*instrRatio;

     %% Practice phase 1 instructions (or normal run with showInstruct on)
    if practice == 1||practice == 0
    %% Practice phase 1 first instruction slide      
    
    instructText{1} = 'Welcome to the experiment.\nIn this task you will be asked to report the location of a target on a ring.';    
    instructText{2} = 'Click to continue';


    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);
    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);

    instrImage = zeros(round(p.ppd* (p.imSize(1)*2*instrRatio) ));    
    [~, instrGrating] = rd_aperture(instrImage, p.aperture, gratingRadius, edgeWidth);
    instrGrating = rescale(instrGrating, 0.5-0.5*p.gratingContrasts(1), 0.5+0.5*p.gratingContrasts(1));
    instrGrating = rd_aperture(instrGrating, 'gaussian', p.gaussianSD*p.ppd, edgeWidth);
    instrGratingTex = Screen('MakeTexture', window, instrGrating*white);
    
    instrTheta = pi/3;
    instrImRect = floor(CenterRectOnPointd([0 0 size(instrImage)], cx+ecc*cos(instrTheta)*instrRatio, cy+ecc*sin(instrTheta)*instrRatio));
    Screen('DrawTexture', window, instrGratingTex, [], instrImRect);
    Screen('Flip', window);
    keyPressed = 0;

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end

    %% Practice phase 1 second instruction slide
    keyPressed = 0;
    instructText{1} = 'To do so, you will move a circle to your\nbest guess of the center of the target.\nWhen you click to submit your guess, the circle will turn black.';    
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);
    
    rectSize = instrRatio * p.responseDotDiameter * p.ppd;
    instrResponseRect = CenterRectOnPointd([0 0 rectSize rectSize], cx+ecc*cos(instrTheta)*instrRatio, cy+ecc*sin(instrTheta)*instrRatio);
    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);
    
    Screen('FillOval', window, p.estimationColorActive, instrResponseRect);
    keyPressed = 0;
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end    

    %% Practice phase 1 third instruction slide
    keyPressed = 0;
    instructText{1} = ['Next, you will draw an arc indicating your\nuncertainty about your position guess.' ...
        '\nMoving the mouse right grows the arc,\n and moving the mouse left shrinks the arc.'];
    instructText{2} = 'Click to continue';

    instrArcAngle = 35;
    instrCenteredArcRect = CenterRectOnPointd(instrRatio*arcRect, cx, cy);
    
    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);
    
    Screen('FillOval', window, p.estimationColorSubmit, instrResponseRect);
    Screen('FrameArc',window, p.arcColorActive, instrCenteredArcRect, ...
        rad2deg(instrTheta)+90-instrArcAngle/2, instrArcAngle, p.arcStroke*instrRatio, p.arcStroke*instrRatio);
    keyPressed = 0;
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end
    
    %% Practice Phase 1 fourth instructions slide
    keyPressed = 0;
    instructText{1} = ['Your goal is to make the arc as small as possible\nwhile keeping the target''s center within the arc.' ...
        '\nYou will earn more points the smaller the arc is, as shown by the graph below,\nbut you will earn no points if the arc misses the target.']  ;
    instructText{2} = 'Click to continue';    

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);
   
    graphHeight = round(8*p.ppd);
    graphWidth = round(12*p.ppd);
    pointsGraph = figure('Visible', 'off', 'Position', [0 0 graphWidth, graphHeight]); clf; set(0, 'CurrentFigure', pointsGraph);
    xs = 0:60;
    ys = 100*exp(-p.points*xs);
     plot(xs, ys ,'w');
    xlabel('arc length (degrees)');
    ylabel('possible points earned');
   set(gcf, 'Color', p.backgroundColor*[1 1 1]);
    set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
    set(gca, 'Color', p.backgroundColor*[1 1 1]);

    pointsGraphFrame = getframe(pointsGraph);
    graphImage = pointsGraphFrame.cdata;
    graphTex = Screen('MakeTexture', window, graphImage);
    Screen('DrawTexture', window, graphTex);
    keyPressed = 0;
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end

    %% Practice phase 2 instructions (enforce fixation)
    elseif practice == 2
    %% Practice phase 2 first instructions slide
    keyPressed = 0;
    instructText{1} = ['It is important that you keep your eyes\nfixed on the dot at the center of the screen\n' ...
        'when the target appears. In this practice\nthe target will not appear until you fixate on the center.'];
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);
    
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end

    elseif practice==3
    %% Practice phase 3 instructions
    %% Practice phase 3 first instruction slide
    keyPressed = 0;
    instructText{1} = ['On previous practices, the target has been preceded by\na circle surrounding the central dot.\n' ...
        'this is a neutral cue meaning the target is equally\nlikely to appear at any location around the circle.'];
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);

    instrAttRect = [0 0 p.ppd/(p.cueLength+p.cueBuffer)*instrRatio p.ppd/(p.cueLength+p.cueBuffer)*instrRatio];
    instrCenteredAttRect = CenterRectOnPointd(instrAttRect, cx, cy);

    Screen('FrameOval', window, white, instrCenteredAttRect, p.strokeWidth); % 
    
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end   

    %% Practice phase 3 second instruction slide

    keyPressed = 0;
    instructText{1} = 'In this practice, on some trials, a directional cue instead of a neutral cue\nwill appear at the center of the screen before the target appears.';
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);

    instrCueTheta = pi/3.8;
    Screen('DrawLine', window, white, cx + p.cueBuffer*p.ppd*cosd(instrCueTheta), cy + p.cueBuffer*p.ppd*sind(instrCueTheta),...
        cx + p.cueLength*p.ppd*cosd(instrCueTheta)*instrRatio, cy + p.cueLength*p.ppd*sind(instrCueTheta)*instrRatio, p.strokeWidth);
    
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end   

    %% Practice phase 3 third instruction slide

    keyPressed = 0;
    instructText{1} = ['The directional cue will roughly point towards\nthe target''s location, and it will usually be off by a small amount,\n' ...
        'but sometimes it will be off by a larger amount.\nThe distance between the cued direction and the target location follows\na distribution shown by the graph below.'];
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    graphHeight = round(8*p.ppd);
    graphWidth = round(12*p.ppd);
    offsetGraph = figure('Visible', 'off', 'Position', [0 0 graphWidth, graphHeight]); clf; set(0, 'CurrentFigure', offsetGraph);
    xs = -20:20;
    ys = normpdf(xs,0,p.sigma);
    plot(xs, ys ,'w');
    xlabel('distance between cued location and target location (degrees)');
    ylabel('probability');
    set(gcf, 'Color', p.backgroundColor*[1 1 1]);
    set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
    set(gca, 'Color', p.backgroundColor*[1 1 1]);

    offsetGraphFrame = getframe(offsetGraph);
    graphImage = offsetGraphFrame.cdata;
    graphTex = Screen('MakeTexture', window, graphImage);
    Screen('DrawTexture', window, graphTex);
    Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end       

    elseif practice == 4
    %% Practice phase 4 instructions
    %% Practice phase 4 first instruction slide

    keyPressed = 0;
    instructText{1} = ['In this practice, in order to make the task difficult, the contrast of\nthe target will be lower than it was in previous practices.' ...
        'You may not always spot the target, in which case you should guess randomly\nor use the directional cue to make a best guess if possible.'];
    instructText{2} = 'Click to continue';

    DrawFormattedText(window, instructText{1}, 'center', cy-round(ecc*instrRatio*upperTextLocation), [1 1 1]*white);
    DrawFormattedText(window, instructText{2}, 'center', cy+round(ecc*instrRatio*lowerTextLocation), [1 1 1]*white);

    drawCircles(window, pInstr);
    drawFixation(window, cx, cy, fixSize*instrRatio, p.fixColor*white);

        instrImage = zeros(round(p.ppd* (p.imSize(1)*2*instrRatio) ));    
    [~, instrGrating] = rd_aperture(instrImage, p.aperture, gratingRadius, edgeWidth);
    instrGrating = rescale(instrGrating, 0.5-0.5*p.gratingContrasts(1), 0.5+0.5*p.gratingContrasts(1));
    instrGrating = rd_aperture(instrGrating, 'gaussian', p.gaussianSD*p.ppd, edgeWidth);
    instrGratingTex = Screen('MakeTexture', window, instrGrating*white);

    instrTheta = pi/3;
    instrImRect = floor(CenterRectOnPointd([0 0 size(instrImage)], cx+ecc*cos(instrTheta)*instrRatio, cy+ecc*sin(instrTheta)*instrRatio));
    Screen('DrawTexture', window, instrGratingTex, [], instrImRect);  

   Screen('Flip', window);

    WaitSecs(0.2);

    while ~keyPressed
        [x,y,buttons] = GetMouse(window);
        if any(buttons)
            keyPressed = keyPressed+1;
            break
        end
    end  


    end
end



% Begin screen
DrawFormattedText(window, 'Have fun!', 'center', cy-round(3*p.ppd), [1 1 1]*white);
DrawFormattedText(window, 'Click to begin', 'center', cy+round(3*p.ppd), [1 1 1]*white); % example image 
% example target (highest contrast) to prevent stim delay on 1st trial 
Screen('DrawTexture', window, tex{numel(tex)}, []); 

if p.photodiode % strcmp(p.testingLocation,'MEG')
    drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
end

Screen('Flip', window);
pause(1)
keyPressed = 0;
while ~keyPressed
%     if p.useKbQueue
%         [keyIsDown, firstPress] = KbQueueCheck();
%         keyCode = logical(firstPress);
%     else
%         [secs, keyCode] = KbWait(devNum);
%     end
%     
%     if strcmp(KbName(keyCode),'1!')
%         keyPressed = 1;
%     end
    [x,y,buttons] = GetMouse(window);
    if any(buttons)
        keyPressed = keyPressed+1;
    end
end

%% TRIALS ** 

if ~resumeDataset
    % Initialize
    block = 1; % block count
    triggerTimes = [];
    correct = []; % 1 correct (target within arc), 0 incorrect (target outside arc)
    mouseTheta = NaN;
    mouseThetaDeg = NaN;
    arcAngle = NaN;
    timeStartArc = NaN;
    timeMouseArc = NaN;
    points = NaN;
    correct = NaN;
    errorMouse = NaN;
    moveFromCenter = 0;
    
    % Save expt start time
    timeExptStart = GetSecs;
    
    trialIdx = 1; % initialize trial counter
else
    trialIdx = iTrial; % pickup where we left off 
end
 
% Trial start
for iTrial = trialIdx:nTrials
    %%%% set conditions for this trial
    thetaDeg = p.theta(trials(iTrial, angularPosIdx)); % target angular position (deg)
    contrast = trials(iTrial, contrastIdx); 
    theta = deg2rad(thetaDeg); % target angular position (rad)
    iti = p.itis(trials(iTrial, itiIdx)); 
    fixTargetISI = p.fixTargetISIs(trials(iTrial,fixTargetISIIdx));  
    
    iStage = 0; % alternates 1 and 0, for photodiode
    
    % ___________________________________________________________________
    % Generate stimuli (if noise) for this trial, otherwise pregenerated
     imRectOval = floor(CenterRectOnPointd(imRect, cx+ecc*cos(theta), cy+ecc*sin(theta)));
%     
%     p.noiseContrast = 0.5;
%     orientation = 0; % filters all orientations so this value is not used 
%     orientBandwidth = 5; 
%     sfBandLow = p.gratingSF/2;
%     sfBandHigh = p.gratingSF*2;
%     maskWithAperture = 0;
%     maxDim = max([sw sh]); 
%     imSize = maxDim/p.ppd; % 10; % dva
    
%     if p.noiseContrast>0
%         % imSize_noise = round([sw/p.ppd sh/p.ppd]);
%         [filteredNoiseIm, numGenerated, pass] = kt_makeFilteredNoise(imSize, p.noiseContrast, ...
%             orientation, orientBandwidth, ...
%             sfBandLow, sfBandHigh, ...
%             p.ppd, maskWithAperture, 'allOrientations');
%         % filteredNoiseIm = rescale(filteredNoiseIm, 0.5-0.5*p.noiseContrast, 0.5+0.5*p.noiseContrast);
%         % apply square aperture to size of oval
%         rad = ecc+buffer*1.5; % degrees
%         [imout, ap] = rd_aperture(filteredNoiseIm, 'square', rad);
%         % make texture
%         % tex_noise = Screen('MakeTexture', window, imout*white);   
%         
%         halfIm = size(filteredNoiseIm)/2;
%         
%         % Add noise to blank
%         blank_noise = ones([maxDim maxDim])*0.5;
%         blank_noise = imout;
%         % Crop image to size of s reen (uneeded once noise generator can take 2D input)
%         sizeBlank = size(blank_noise);
%         if sizeBlank(2)>sh % check height
%             sizeDiff = sizeBlank(2)-sh;
%             blank_noise(end+1-sizeDiff/2:end,:) = []; % crop end
%             blank_noise(1:sizeDiff/2,:) = []; % crop beginning % why is this throwing an error
%         end
%         % Add grating to blank
%         blank_grating = ones([sh sw])*0.5;
%         blank_grating(imRectOval(1):imRectOval(3)-1, imRectOval(1):imRectOval(3)-1) = grating(:,:,contrast);
%         % Add noise and grating
%         blank_noiseGrating = (blank_noise + blank_grating)/2;
%         
%         % Make texture
%         tex_noiseGrating = Screen('MakeTexture', window, blank_noiseGrating*white);
%     end
    
    % ___________________________________________________________________
    % Fixation 
    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(1), 0); % alternate black and white
    end
    timeFix = Screen('Flip', window);

    if p.triggersOn
        PTBSendTrigger(p.triggers.fixation, 0);
        triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    end
    if p.triggerEEG
        ppdev_mex('Write',1,p.triggers.fixation);
        WaitSecs(0.05);
        ppdev_mex('Write',1,0);
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_FIX');
    end

    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
    Screen('Flip', window, timeFix + 0.05);

    % ___________________________________________________________________
    % Cue
    if p.attConds(trials(iTrial,attIdx))==1 % att valid
        % randomly generate cue angle from gaussian
        cueThetaRnd = normrnd(thetaDeg, p.sigma, [1,1]); % from normal distribution (in deg) 
        % rescale angle to 0-360 space
        if cueThetaRnd < 0
            cueTheta = cueThetaRnd + 360;
        elseif cueThetaRnd > 360
            cueTheta = cueThetaRnd - 360;
        else
            cueTheta = cueThetaRnd;
        end
    else
        cueThetaRnd = NaN; 
        cueTheta = NaN;
    end
    
    if p.attConds(trials(iTrial,attIdx))==1 % att valid 
        % draw valid cue 
        Screen('DrawLine', window, white, cx + p.cueBuffer*p.ppd*cosd(cueTheta), cy + p.cueBuffer*p.ppd*sind(cueTheta),...
            cx + p.cueLength*p.ppd*cosd(cueTheta), cy + p.cueLength*p.ppd*sind(cueTheta), p.strokeWidth);
    else % neutral 
        % draw neutral cue
        Screen('FrameOval', window, white, centeredAttRect, p.strokeWidth); % 
    end
    
    drawCircles(window, p)
    drawFixation(window, cx  , cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(1), 0); % alternate black and white
    end
    timeCue = Screen('Flip', window, timeFix + p.fixDur - p.slack);

    if p.triggersOn
        PTBSendTrigger(p.triggers.cue, 0);
        triggerTimes = [triggerTimes; p.triggers.cue GetSecs];
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_CUE');
    end
    if p.triggerEEG
        ppdev_mex('Write',1,p.triggers.cue);
        WaitSecs(0.05);
        ppdev_mex('Write',1,0);
    end

    % ___________________________________________________________________
    % Fixation2
    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
    timeFix2 = Screen('Flip', window, timeCue + p.cueDur - p.slack);

    if p.triggersOn
        PTBSendTrigger(p.triggers.fixation, 0);
        triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_FIX2');
    end

    
    % ___________________________________________________________________
    % TARGET IMAGE 
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.noiseContrast > 0 
        % cRect = CenterRectOnPointd(rect, cx, cy);
        Screen('DrawTexture', window, tex_noiseGrating, []); % rect; noise + bullseye
    else
        % Screen('DrawTexture', window, tex{contrast}, [], imRectOval); % bullseye
        % spatial frequency dependence 
        % testing 2 settings: 
        % --- low contrast, low spatial frequency 
        % --- high contrast, high spatial frequency 
        if contrast==1 % low contrast 
            iSF = 1; 
        elseif contrast==2
            iSF = 2; 
        end
        sf = p.gratingSF(iSF);
        Screen('DrawTexture', window, tex{contrast,iSF}, [], imRectOval); % bullseye
    end
    drawCircles(window, p)  
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(1), 0); % alternate black and white
    end

    %% Check fixation before presenting target
    isFix = 0;
    if p.enforceFixCTI
        eyeSlack_min = p.eyeSlacks(rd_sampleDiscretePDF(p.eyeSlackPDF, 1, 5));
    else
        eyeSlack_min = 0;
    end
    eyeSlack = 0;

    while ~isFix||GetSecs<timeFix2 + p.fix2Dur + fixTargetISI - 2*p.slack
        if p.enforceFixCTI
            isFix = rd_eyeLink('fixcheck', window, {cx, cy, eyeRad});
            if ~isFix
                remain = timeFix2 + p.fix2Dur + fixTargetISI - 2*p.slack - GetSecs + eyeSlack;
                if remain < eyeSlack_min
                    eyeSlack = eyeSlack_min-remain;
                end
            end
        else
            isFix = 1;
        end
    end


    timeTarget = Screen('Flip', window, timeFix2 + p.fix2Dur + fixTargetISI - p.slack);

    if p.triggersOn
        PTBSendTrigger(p.triggers.image, 0);
        triggerTimes = [triggerTimes; p.triggers.image GetSecs]; 
    end
    if p.triggerEEG
        ppdev_mex('Write',1,p.triggers.image);
        WaitSecs(0.05);
        ppdev_mex('Write',1,0);
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_IMAGE');
    end

    
    % ___________________________________________________________________
    % RESPONSE (dot, estimation)
    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
    timeResponseDotProbe = Screen('Flip', window, timeTarget + p.targetDur - p.slack); 

    if p.triggersOn
        PTBSendTrigger(p.triggers.fixation, 0);
        triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_FIX');
    end
    
    % Target response SOA 
    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
    Screen('Flip', window, timeTarget + p.targetRespCueSOA - p.slack);
    
    % Start mouse movement data collection
    KbQueueStart(m);
    SetMouse(cx,cy,window); % Set mouse to center 
    [mx,my] = GetMouse;   
    drawCircles(window, p)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    Screen('FillOval', window, p.estimationColorActive, responseRect); % response at center
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
            
    drawResponseDot(window, mx, my, responseRadius*2, p.estimationColorActive);
    timeResponseDotProbe = Screen('Flip', window);
    
    timeStartDot = GetSecs;
    clicks = 0;
    
    x0 = cx; y0 = cy; % initialize mouse position to center 
    moveFromCenter = 0;
    % while 1 && (GetSecs-timeStartDot) < p.responseDotDur
    firstMove = 0; 
    circ_x = cx;
    circ_y = cy;
    while (GetSecs-timeStartDot) < p.responseDotDur
        WaitSecs(p.mouseSample);
        [x,y,buttons] = GetMouse(window);
        if ~useTrackball
            dx = (x-x0);
            dy = (y-y0);
            if dx~=0 || dy~=0 % wait for mouse movement
                % if ~moveFromCenter, WaitSecs(p.mouseIntialMovementBuffer), end
                moveFromCenter = 1;
                new_x = circ_x+dx-cx;
                new_y = circ_y+dy-cy;
                if new_x >= 0 % right half
                    if new_y <= 0 % bottom half
                        mouseThetaDeg = 360+atand(new_y/new_x); % angle of mouse movement (degrees)
                    else
                        mouseThetaDeg = atand(new_y/new_x);
                    end
                else % left half
                    mouseThetaDeg = atand(new_y/new_x)+180; % angle of mouse movement (degrees + semicircle)
                end
                mouseTheta = deg2rad(mouseThetaDeg); % radians

                circ_x = cx+ecc*cos(mouseTheta);
                circ_y = cy+ecc*sin(mouseTheta);
            end

            drawCircles(window, p)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
            responseRectTheta = CenterRectOnPointd(responseRect, circ_x, circ_y);
            if p.photodiode % strcmp(p.testingLocation,'MEG')
                drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
            end
            Screen('FillOval', window, p.estimationColorActive, responseRectTheta); % response
        elseif useTrackball
            dy = -dy; % move trackball in natural position 
            if dx >= 0 % right half
                if dy >= 0 % bottom half
                    mouseThetaDeg = atand(dy/dx); % angle of mouse movement (degrees)
                else
                    mouseThetaDeg = 360+atand(dy/dx);
                end
            else % left half
                mouseThetaDeg = atand(dy/dx)+180; % angle of mouse movement (degrees + semicircle)
            end
            mouseTheta = deg2rad(mouseThetaDeg); % radians
            if dx~=0 || dy~=0
                moveFromCenter = 1; 
            end
            if moveFromCenter % mx, my % if mouse isn't at center
                drawCircles(window, p)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
                responseRectTheta = CenterRectOnPointd(responseRect, cx+ecc*cos(mouseTheta), cy+ecc*sin(mouseTheta));
                Screen('FillOval', window, p.estimationColorActive, responseRectTheta); % response
                moveFromCenter = 1;
            elseif ~moveFromCenter % if mouse at center
                drawCircles(window, p)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
                responseRectTheta = CenterRectOnPointd(responseRect, cx, cy);
                Screen('FillOval', window, p.estimationColorActive, responseRectTheta); % response
            end
        end
        
        if p.debug
            [mouseText, errMsg] = sprintf('theta = %0.2f, deg = %0.2f \n\n x = %d, y = %d \n\n dx = %0.2f, dy = %0.2f\n\n circ_x = %0.2f, circ_y = %0.2f', mouseTheta, mouseThetaDeg, x, y, dx, dy, circ_x, circ_y);
            DrawFormattedText(window, mouseText, 'center', cy-round(1*p.ppd), [1 1 1]*white);
        end
        if p.photodiode % strcmp(p.testingLocation,'MEG')
            drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
        end
        Screen('Flip', window, timeStartDot+p.mouseSample);
        SetMouse(cx,cy,window); % Set mouse to center
        %WaitSecs(p.mouseSample);
           
        timeMouseDot = GetSecs; % sample until mouse is clicked. if no mouse click then timeout
        if any(buttons) && moveFromCenter
            clicks = clicks+1;
            [~, ~, dotClick] = GetMouse(window);
            WaitSecs(.001); % wait 1 ms
            
            % response dot color change when mouse clicked
            drawCircles(window, p)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
            Screen('FillOval', window, p.estimationColorSubmit, responseRectTheta); % response dot
            if p.photodiode % strcmp(p.testingLocation,'MEG')
                drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
            end
            timeMouseDot = Screen('Flip', window, timeStartDot + p.mouseSample);
            break
        end
    end

    if clicks==0 % no response submitted, end trial now 
        drawCircles(window, p)
        drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
        missText = 'too slow!';
        DrawFormattedText(window, missText, 'center', cy-round(1*p.ppd), [1 1 1]*white);
        if p.photodiode % strcmp(p.testingLocation,'MEG')
            drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
        end
        timeFeedback = Screen('Flip', window, timeStartDot + p.mouseSample);
        pause(1)

    else % otherwise continue
        if p.triggersOn
            PTBSendTrigger(p.triggers.responseDot, 0);
            triggerTimes = [triggerTimes; p.triggers.responseDot GetSecs];
        end
        if p.triggerEEG
            ppdev_mex('Write',1,p.triggers.responseDot);
            WaitSecs(0.05);
            ppdev_mex('Write',1,0);
        end
        
        % ___________________________________________________________________
        % ARC
        if p.wedge
            drawCircles(window, p)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
            Screen('FillOval', window, p.estimationColorSubmit, responseRectTheta); % response dot
            if p.photodiode % strcmp(p.testingLocation,'MEG')
                drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
            end
            timeResponseArcProbe = Screen('Flip', window, timeMouseDot + p.dotArcSOA);

            arcAngle0 = (360*responseRadius*2)/(2*pi*p.ecc*p.ppd); % minimum arc angle (deg) corresponds to response dot diameter
            SetMouse(cx,cy,window); % center mouse
            timeStartArc = GetSecs;
            clicks = 0;
            arcAngle = arcAngle0;
            x0 = cx;
            DXs = [];
            Screen('FrameArc',window, p.arcColorActive, centeredArcRect, mouseThetaDeg+90-arcAngle0/2, arcAngle, p.arcStroke, p.arcStroke);
            while (GetSecs - timeStartArc) < p.responseArcDur
                [x,~,buttons] = GetMouse(window);
                dx = (x-x0)*p.mouseSensitivity;
                if dx > 0 % mouse to right
                    if arcAngle > 360
                        arcAngle = 360;
                    else
                        arcAngle = arcAngle + p.arcIncrement*dx; % arc grows
                    end
                elseif dx < 0 % mouse to left
                    if arcAngle < arcAngle0
                        arcAngle = arcAngle0;
                    else
                        arcAngle = arcAngle + p.arcIncrement*dx; % arc shrinks
                    end
                end
                drawCircles(window, p)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
                if p.photodiode % strcmp(p.testingLocation,'MEG')
                    drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
                end
                Screen('FrameArc',window, p.arcColorActive, centeredArcRect, mouseThetaDeg+90-arcAngle/2, arcAngle, p.arcStroke, p.arcStroke); % arc
                Screen('FillOval', window, p.estimationColorSubmit, responseRectTheta); % response dot
                Screen('Flip', window, timeStartArc + p.mouseSample);

                x0 = x;
                DXs = [DXs dx];

                if any(buttons)
                    clicks = clicks+1;
                    [~, ~, arcClick] = GetMouse(window);
                    WaitSecs(.001); % wait 1 ms

                    % color change when arc clicked
                    drawCircles(window, p)
                    drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
                    Screen('FrameArc',window, p.arcColorSubmit, centeredArcRect, mouseThetaDeg+90-arcAngle/2, arcAngle, p.arcStroke, p.arcStroke); % arc
                    Screen('FillOval', window, p.estimationColorSubmit, responseRectTheta); % response dot
                    if p.photodiode % strcmp(p.testingLocation,'MEG')
                        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
                    end

                    timeMouseArc = Screen('Flip', window, timeStartArc + p.mouseSample);
                    break
                end
            end
            if clicks==0 % no response submitted
                arcAngle = NaN;
                drawCircles(window, p)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
                missText = 'too slow!';
                DrawFormattedText(window, missText, 'center', cy-round(1*p.ppd), [1 1 1]*white);
                if p.photodiode % strcmp(p.testingLocation,'MEG')
                    drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
                end
                timeMouseArc = Screen('Flip', window, timeStartDot + p.mouseSample);
                pause(1)
            end
            if p.triggersOn
                PTBSendTrigger(p.triggers.responseArc, 0);
                triggerTimes = [triggerTimes; p.triggers.responseArc GetSecs];
            end

            if p.eyeTracking
                Eyelink('Message', 'EVENT_RESPONSE');
            end

            % ___________________________________________________________________
            % Reward
            if mouseThetaDeg - thetaDeg < -180
                mouseThetaDeg = 360 + mouseThetaDeg; % rescale for acute angle [-180 to 180]
            elseif mouseThetaDeg - thetaDeg > 180
                mouseThetaDeg = mouseThetaDeg - 360; % this should probably be circular degrees
            end
            % errorMouse = mouseThetaDeg - thetaDeg;
            errorMouse = circ_dist(mouseTheta, theta);
            if thetaDeg > mouseThetaDeg - arcAngle/2 && thetaDeg < mouseThetaDeg + arcAngle/2 % check if arc missed target
                points = round(100*exp(-p.points*arcAngle));
                if arcAngle > 50
                    pointsText = sprintf('+ %d points,\n try to make the arc as small as possible \n while containing the target', points);
                else
                    pointsText = sprintf('+ %d points!', points);
                end
                correct = 1;
            else
                points = 0; correct = 0;
                pointsText = sprintf('+ %d points, missed the target location', points);
            end
            Screen('FrameArc',window, p.arcColorSubmit, centeredArcRect, mouseThetaDeg+90-arcAngle/2, arcAngle, p.arcStroke, p.arcStroke); % arc
            DrawFormattedText(window, pointsText, 'center', cy-round(1*p.ppd), [1 1 1]*white);
        end
        drawCircles(window, p)
        drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
        Screen('FillOval', window, p.estimationColorSubmit, responseRectTheta); % response dot

        if p.showPositionFeedback
            trueRectTheta = CenterRectOnPointd(responseRect, cx+ecc*cos(theta), cy+ecc*sin(theta));
            Screen('FillOval', window, p.feedbackEstimationColor, trueRectTheta); % true position dot
            if p.wedge
                errorText = sprintf('%0.2f%c error', abs(rad2deg(errorMouse)), char(176));
                DrawFormattedText(window, errorText, 'center', cy+round(1*p.ppd), [1 1 1]*white);
            end
        end
        if p.photodiode % strcmp(p.testingLocation,'MEG')
            drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
        end
        timeFeedback = Screen('Flip', window);
        % pause(1) % needed? make p.feedbackDur longer if neede 
        if p.triggersOn
            PTBSendTrigger(p.triggers.feedback, 0);
            triggerTimes = [triggerTimes; p.triggers.feedback GetSecs];
        end
    end
        
    % ___________________________________________________________________
    % ITI
    drawCircles(window, p)
%     drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end
    timeITIstart = Screen('Flip', window, timeFeedback + p.feedbackDur - p.slack);

    if p.triggerEEG
        ppdev_mex('Write',1,p.triggers.ITI);
        WaitSecs(0.05);
        ppdev_mex('Write',1,0);
    end
    
    drawCircles(window, p)
%     drawFixation(window, cx, cy, fixSize, p.fixColor*white); % fixation
    if p.photodiode % strcmp(p.testingLocation,'MEG')
        drawPhotodiode(window, [cx cy]*2, p.photodiodeColors(2), 0); % alternate black and white
    end

    %check fixation during ITI before continuing
    isFix = 0;
    if p.enforceFixITI
        eyeSlack_min = p.eyeSlacks(rd_sampleDiscretePDF(p.eyeSlackPDF, 1, 5));
    else
        eyeSlack_min = 0;
    end
    eyeSlack = 0;

    while ~isFix||GetSecs<timeITIstart+iti-2*p.slack+eyeSlack
        if p.enforceFixITI
            isFix = rd_eyeLink('fixcheck', window, {cx, cy, eyeRad});
            if ~isFix
                remain = timeITIstart + iti - 2*p.slack - GetSecs + eyeSlack;
                if remain < eyeSlack_min
                    eyeSlack = eyeSlack_min-remain;
                end
            end
        else
            isFix = 1;
        end
    end

    timeITIend = Screen('Flip', window, timeITIstart + iti - p.slack);
    
    % ___________________________________________________________________
    % Save trial info
    timing.rtDot(iTrial) = timeMouseDot - timeStartDot; 
    timing.rtArc(iTrial) = timeMouseArc - timeStartArc;
    timing.fixDur(iTrial) = timeCue - timeFix; 
    timing.cueDur(iTrial) = timeTarget - timeCue; 
    timing.targetDur(iTrial) = timeResponseDotProbe - timeTarget; 
    timing.itiDur(iTrial) = timeITIstart - timeITIend;
     
    trialsPresented.itiIdx(iTrial) = trials(iTrial,itiIdx); 
    trialsPresented.iti(iTrial) = p.itis(trials(iTrial,itiIdx)); 
    
    trialsPresented.fixTargetISIIdx(iTrial) = trials(iTrial,fixTargetISIIdx); 
    trialsPresented.fixTargetISI(iTrial) = p.fixTargetISIs(trials(iTrial,fixTargetISIIdx)); 
    
    trialsPresented.att(iTrial) = p.attConds(trials(iTrial,attIdx)); % 1 valid, 0 neutral 
    
    trialsPresented.contrast(iTrial) = p.gratingContrasts(trials(iTrial, contrastIdx)); 
    trialsPresented.sf(iTrial) = sf; 
    trialsPresented.sfIdx(iTrial) = iSF; 
    trialsPresented.theta(iTrial) = theta; % true stimulus position
    trialsPresented.thetaDeg(iTrial) = thetaDeg; % true stimulus position
    trialsPresented.cueThetaRnd(iTrial) = cueThetaRnd; % gaussian cue degree (uncorrected, 0-1 space)
    trialsPresented.cueTheta(iTrial) = cueTheta; % gaussian cue degree (corrected, 0-365 space)  
    trialsPresented.cueOffset(iTrial) = circ_dist( deg2rad(cueTheta), theta ) ; % Save cue offset (rad) 
    trialsPresented.cueOffset_deg(iTrial) = rad2deg(trialsPresented.cueOffset(iTrial)); 
    
    trialsPresented.responseTheta(iTrial) = mouseTheta; % position estimation (rad) 
    trialsPresented.responseThetaDeg(iTrial) = mouseThetaDeg; % position estimation (deg) 
    trialsPresented.arcAngleDeg(iTrial) = arcAngle; % arc angle 
    trialsPresented.arcAngleWedge(iTrial,:) = [mouseThetaDeg - arcAngle/2, mouseThetaDeg + arcAngle/2]; % degrees that arc angle spans
    trialsPresented.points(iTrial) = points; % points earned 
    trialsPresented.correct(iTrial) = correct; % correct 1, incorrect 0 
    trialsPresented.error(iTrial) = errorMouse; % error between reported and true position
    trialsPresented.block(iTrial) = block; % block number 
    trialsPresented.pointsBlock(block) = sum(trialsPresented.points(trialsPresented.block==block)); % points per block 
    
    % Save the workspace every trial 
    save(sprintf('%s/TEMP',subjectDataDir)) % this can be done better 
        
    if mod(iTrial,p.nTrialsPerBlock)==0 || iTrial==nTrials % last trial in a block
        
        blockStartTrial = (iTrial/p.nTrialsPerBlock)*p.nTrialsPerBlock - p.nTrialsPerBlock + 1;
        if blockStartTrial < 0 % we are doing less than one block
            blockStartTrial = 1;
        end
        trialsInBlock = trials(blockStartTrial:iTrial,:);
        blockMessage = sprintf('%s You''ve completed %d of %d blocks.', highpraise, block, ceil(nTrials/p.nTrialsPerBlock));
        if iTrial==nTrials
            keyMessage = ''; % last block
        else
            keyMessage = 'Press 1 to go on.';
        end
        pointsMessages = []; 
        for i = 1:block
            pointsMessage = sprintf('Block %d: %d points', i, trialsPresented.pointsBlock(i));
            pointsMessages = sprintf('%s\n%s', pointsMessages, pointsMessage); 
        end
        breakMessage = sprintf('%s\n\n\n%s\n\n\n%s', blockMessage, pointsMessages, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
        WaitSecs(1);
        if iTrial < nTrials
            keyPressed = 0;
            while ~keyPressed
                if p.useKbQueue
                    [keyIsDown, firstPress] = KbQueueCheck();
                    keyCode = logical(firstPress);
                else
                    [secs, keyCode] = KbWait(devNum);
                end          
                if strcmp(KbName(keyCode),'1!')
                    keyPressed = 1;
                end
            end
        end   
        block = block+1; % keep track of block for block message only 


        %% Save data at blocks 
        if saveData
            expt.date = datestr(now,'yyyymmdd_HHMM'); 
            dataFile = sprintf('data/%s/%s_%s_%s_run%s.mat', subjectID, subjectID, p.expName, expt.date, run);
            save(dataFile, 'expt')
            disp('data saved')
        end

    end
end

% Completion message 
WaitSecs(1); 
DrawFormattedText(window, 'All done! Thanks for your effort','center','center', white); 
Screen('Flip', window); 
WaitSecs(1);  
timeExptEnd = GetSecs; 

% Clean up
ShowCursor;

%% Store expt info
timing.triggers = triggerTimes; 

expt.subjectID = subjectID;
expt.run = run;
expt.p = p;
expt.timing = timing; 
expt.trials_headers = trials_headers;
expt.trials = trials;
expt.trialsPresented = trialsPresented;
expt.whenSaved = datestr(now);    
% expt.time = datestr(now,'HHAM'); 
expt.subjectDataDir = subjectDataDir; 

% Save expt end time
exptDur = (timeExptEnd - timeExptStart)/60; 
expt.exptDur = exptDur; % in min 

%% Analyze and save data
% Is this redundant with the previous save? 
if saveData
    dataFile = sprintf('data/%s/%s_%s_%s_run%s.mat', subjectID, subjectID, p.expName, expt.date, run); 
    save(dataFile, 'expt')
    disp('data saved')
end

%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
    
    % rename eye file
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
    delete(sprintf('%s/%s.edf', eyeDataDir, eyeFile))
end

%% Close screen
Screen('CloseAll')

if p.triggerEEG
    ppdev_mex('Close',1);
end
