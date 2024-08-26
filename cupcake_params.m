function p = cupcake_params(testingLocation, debug, showInstruct, showPositionFeedback) 

%% defaults 
p.testingLocation = testingLocation;
p.debug = logical(debug);
p.showInstruct = logical(showInstruct);
p.showPositionFeedback = logical(showPositionFeedback); 
p.wedge = 0; % if no wedge, the task is only position estimate without confidence report (for piloting)  
p.photodiode = 1; % draws flickering rectangle in upper left corner 
p.enforceFix = 0; % only for thresholding and practice

%% 
switch p.testingLocation
    case 'MEG'
        p.keyNames = {'1!'};
        p.refRate = 60;
        p.screenSize = [23.5 17]; % cm
        p.screenRes = [f1024 768];
        p.viewDist = 42; % cm
        p.eyeTracking = 0;
        p.useKbQueue = 0; % using trackball for all 
        p.soundAmp = 0.1;
        p.triggersOn = 1;
        p.displayPath = '/Users/megadmin/Desktop/Experiments/Rachel/vistadisp/exptTools2/displays/meg_lcd_20180420_brightness-32';
        p.winFrac = 1; 
        p.triggerEEG = 0; 
    case 'desk'
        p.keyNames = {'1!'};
        p.refRate = 60;
        p.screenSize = [13 9]; % (in)
        p.screenRes = [1280 1024];
        p.viewDist = 36; % (in)
        p.eyeTracking = 0;
        p.useKbQueue = 0;
        p.soundAmp = 1;
        p.triggersOn = 0;
        p.winFrac = 0.6; 
        p.triggerEEG = 0; 
    case 'CarrascoL1'
        p.keyNames = {'1!'};
        p.refRate = 60;
        p.screenSize = [40 30];
        p.screenRes = [1280 960];
        p.viewDist = 56;
        p.eyeTracking = 0; 
        p.useKbQueue = 0;
        p.soundAmp = 1;
        p.triggersOn = 0;
        p.winFrac = 1; 
        p.triggerEEG = 0; 
    case {'DenBehav','DenEEG'}
        p.keyNames = {'1!'};
        p.refRate = 120;
        p.screenSize = [53 30];
        p.screenRes = [1920 1080];
        p.viewDist = 75;
        p.eyeTracking = 1; 
        p.useKbQueue = 0;
        p.soundAmp = 1;
        p.triggersOn = 0; % legacy for MEG
        p.winFrac = 1;
        p.triggerEEG = 1;
    otherwise
        error('Testing location not found.')
end

KbName('UnifyKeyNames');
p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.nReps = 1; % 6 reps with 32 trials/rep -> 3 blocks
p.nTrialsPerBlock = 72; % 360/12; % 360/6 % 64 trials -> ~1.5 min/block
p.eyeRad = 1.5; % allowed fixation radius (degrees)   

% Calculate pixels per degree
p.ppd = ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'); % pixels per degree

%% Main 
p.expName = 'cupcake'; % CupcakeSpatial

%% Text
p.font = 'Verdana';
p.fontSize = 24;

%% Fixation
p.fixColor = 1; % [inner; outer] 
p.fixDiameter = .3; % [.35 .7]; % deg

%% Timing (all seconds) 
p.itiType = 'uniform'; % 'uniform','hazard'
switch p.itiType
    case 'uniform'
        p.itis = .5:0.05:.7; 
    case 'hazard'
        p.itis = 1:0.05:2.5;
        p.hazardProb = 0.2; % at every time step, the probability of the event is 0.2
        p.itiPDF = p.hazardProb.*(1-p.hazardProb).^(0:numel(p.itis)-1); % f = p.*(1-p).^x;
    otherwise
        error('p.itiType not recognized')
end

p.eyeSlacks = 0.2:0.02:0.4; % min time between regaining fixation and next trial if p.enforceFix
p.eyeSlackPDF = p.hazardProb.*(1-p.hazardProb).^(0:numel(p.eyeSlacks)-1); % f = p.*(1-p).^x;

% stimuli timing (all seconds)
p.fixDur = 1; % initial fixation 
p.cueDur = 0.05; % attention cue duration 
p.fix2Dur = 0.25; % fixation between cue and target
p.fixTargetISIs = .001:0.001:.2; % fixation target ISI (for prestim alpha)
p.targetDur = 0.05; % 0.08; % target duration 
p.targetRespCueSOA = 0.5; % time between target onset and go onset
p.responseDotDur = 10; % 0.5; % 10; % Inf; % response window timeout 
p.dotArcSOA = 0.5; % buffer to prevent immediately registering the arc response after the position estimate 
p.responseArcDur = 10; % 10; % response arc window timeout 
p.respGoSOA = 0.75; % 0.6 % time between resp cue onset and go onset. set to zero for no go cue.
p.feedbackDur = 1; % how long to display feedback info

% peripherals timing 
p.mouseSample = 1/120; % how often to sample mouse movements
p.mouseIntialMovementBuffer = 0.5; % small buffer after initial movement of mouse from center to periphery
p.mouseSensitivity = 0.5; % <1 to make mouse less sensitive, >1 to make more sensitive
p.mouseSensitivityDot = 0.01; % make dot even less sensitive

p.eyeSlack = 0.12; % cushion between last fixation check and next stimulus presentation
p.slack = (1/p.refRate)*.9; % to prevent lag if screen flip aligns with screen refresh 

%% Images
% grating 
p.imPos = [0 0]; % dva
p.imSize = [2 2]; % dva; this is the size of the image container that holds the stim
p.gratingDiameter = [p.imSize(1) 0]; % [outer inner] 
p.theta = 1:5:360; % target angular locations, angular resolution 
p.gaussianSD = 0.4; % x4 visible SDs ~= 1.6dva p.gratingRadius(1)/2.5; % gaussian SD here 
% p.gratingOrientations = 0; 
p.gratingContrasts = 0.1; % low high contrast[0.05 0.1]; % [0.05 0.1 0.25]; % [0.05 0.1 0.25 1]; % 1; [0.01 0.02 0.05 0.2 0.5]
p.aperture = 'radial-sine-ring'; % disk:'cosine', annulus:'cosine-ring', concentric:'radial-sine-ring', roth:'vignette-ring'
% p.apertureEdgeWidth = 1/6; % 0.05, half of a period, so sf of radial-sine aperture is 1/(2*width)
% p.angularFreq = 6; % for 'vignette' only
p.gratingSF = 4; % cpd; low high SF
p.edgewidth_deg = 1./(2.*p.gratingSF); % edgewidth determines SF in radial grating 
p.edgewidth_pix = round(p.edgewidth_deg .* p.ppd); 
% if strfind(p.aperture,'radial')
%     p.apetureSF = 1/(2*p.apertureEdgeWidth);
% end
% p.goCueColor = 0.75;
rng('shuffle'); % necessary to reset RNG sequence 
p.RNGseed = rng; % save current RNG seed

% Oval (diameter, in deg) imPos
if strfind(p.testingLocation, 'desk')
    p.ecc = 6; 
else 
    p.ecc = 8; % try 5.5ecc with cueing (2) and contrast (0.05 0.1); 10 degrees 
end
p.strokeWidth = 2; % pixels 
p.ovalColor = 255*0.75; % light grey 
p.ovalBuffer = 1; % dva spatial buffer between target and ring 

% Noise 
p.noiseContrast = 0; % 0.2, set to 0 for no noise 

% Cue
p.attConds = [0 1]; % 0 neutral, 1 valid 
p.pValid = 0.5; % 0.1, 0.5, probability valid 
p.pNeutral = 0.5; % 0.9, 0.5, probability invalid 
if p.pValid + p.pNeutral ~= 1
    error('pValid and pNeutral must equal 1')
end
p.cueLength = 0.75; % cue length (degrees) 
p.cueBuffer = 0.05; % degrees from center, fixation buffer % where does the nuetral cue get defined? 
% gaussian prior of cue 
p.sigma = 7; % 15; 30 (too big), 60, gaussian standard deviation (degrees)

% Response Dot and Arc 
p.responseDotDiameter = 0.5; % diameter of response cue dot (degree) 
p.responseColor = 1;
p.estimationColorActive = [0 0 255]; % blue 
p.estimationColorSubmit = [0 0 0]; % black  
p.arcColorActive = [0 0 255]; % blue 
p.arcColorSubmit = [255 255 255]; % white 
p.arcStroke = 10; % arc stroke width
p.arcIncrement = 0.8; % 0.3 degrees of arc angle change with each mouse movement 
% p.dotIncrement = 0.1; % 

% Reward
p.feedbackEstimationColor = [0 255 0]; % green
p.points = 0.08; % constant for scaling points (asymptote ~50 deg)
p.feedbackMessageTreshold = 50; % arcL that shows 'try to make the arc as small as possible' message 

% Thresholding
p.threshDegrees = 5;

%n_trials = n_reps * n_orientations
p.nTrialsThresholding = 72;
p.nStaircases = 2;

p.quest.pThreshold = 0.5;
p.quest.tGuess = log10(p.gratingContrasts(1));
p.quest.tGuessSd = 1;
p.quest.beta = 3.5;
p.quest.delta = 0.05;
p.quest.gamma = (2*p.threshDegrees)/360;
p.quest.grain = 0.01;
p.quest.range = 1;

%% M/EEG triggers
p.triggers.fixation = 2^1; %2
p.triggers.cue = 2^2; %4
p.triggers.image = 2^3; %8
p.triggers.responseDot = 2^4; %16
p.triggers.ITI = 2^5; %32

