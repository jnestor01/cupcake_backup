function [out, exitFlag] = rd_eyeLink(command, window, in)

%slight edits by me 2021

% Possible commands and their ins & outs:
% 'eyestart'
%   in = eyeFile
%   out = el
%
% 'calibrate'
%   in = el
%   out = cal
% 'startrecording'
% in = {el fixrect}
%
% 'trialstart'
%   in = {el, trialNum, cx, cy, rad}
%   out = []
% 'trialend'
% in = {fixRect, trialNum}
% out = []
% 'fixholdcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 'fixcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 'driftcorrect'
%   in = {el, cx, cy}
%   out = driftCorrection
%
% 'trialstop'
%   in = []
%   out = []
%
% 'eyestop'
%   in = {eyeFile, eyeDataDir}
%   out = []

%% Initializations
% assume no output unless some is given
out = [];

% assume everything goes ok (exitFlag=0) until proven otherwise        
exitFlag = 0;
        
%% Do command
switch command
    case 'eyestart'
        %% start eyetracker
        eyeFile = in;
        
        % First check if we can get a connection
        if EyelinkInit()~= 1
            fprintf('Couldn''t initialize connection with eyetracker! Exiting ...\n');
            return
        end
        
        % Set up the eyetracker
        el = EyelinkInitDefaults(window);
        
        Eyelink('Command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA');
        % Eyelink('Command', 'calibration_type = HV9');
        [v vs] = Eyelink('GetTrackerVersion');
        fprintf('\nRunning experiment on a ''%s'' tracker.', vs );
        
        % Start the eye file
        edfFile = sprintf('%s.edf', eyeFile);
        edfFileStatus = Eyelink('OpenFile', edfFile);
        if edfFileStatus==0
            fprintf('\nEye file opened ok.\n\n')
        else
            fprintf('Cannot open .edf file. Exiting ...');
            Screen('CloseAll')
            exitFlag = 1;
            return
        end
        
        out = el; % return the el structure as output
        
    case 'calibrate'
        %% calibrate eyetracker
        el = in;
        
%         calString = sprintf('Eye tracker calibration:\n\nPlease fixate the center of the dot!\n\nPress ''space'' to start or ''q'' to quit!');
%         DrawFormattedText(window, calString, 'center', 'center', 1, []);
%         Screen('Flip', window, 0, 1); 
%         
%         contKey = '';
%         while isempty(find(strcmp(contKey,'space'), 1))
%             keyIsDown = 0;
%             while ~keyIsDown
%                 [keyIsDown, keyTime, keyCode] = KbCheck(-1); %% listen to all keyboards
%             end
%             contKey = KbName(find(keyCode));
%         end
%         if strcmp(contKey,'q')
%             ListenChar(0);
%             ShowCursor;
%             Screen('CloseAll')
%             fclose('all');
%             fprintf('User ended program');
%             exitFlag = 1;
%             return
%         end
%         Screen('Flip', window, 0, 1);
        
        cal = EyelinkDoTrackerSetup(el);
        if cal==el.TERMINATE_KEY
            exitFlag = 1;
            return
        end
        
        out = cal;
        
    case 'startrecording'
        el = in{1};
        fixRect=in{2};
        
         % Must be offline to draw to EyeLink screen
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at fix point
        Eyelink('Command', 'clear_screen 0');
%             Eyelink('command', 'draw_box %d %d %d %d 15', width/2-50, height/2-50, width/2+50, height/2+50);
        Eyelink('command', 'draw_box %d %d %d %d 15',round(fixRect(1)), round(fixRect(2)), round(fixRect(3)), round(fixRect(4)));

        Eyelink('Command', 'set_idle_mode');
        record = 0;
        while ~record
            Eyelink('StartRecording');	% start recording
            % start recording 100 msec before just to be safe
            WaitSecs(0.1);
            key=1;
            while key~=0
                key = EyelinkGetKey(el); % dump any pending local keys
            end
            
            err=Eyelink('CheckRecording'); 	% check recording status
            if err==0
                record = 1;
                Eyelink('Message', 'RECORD_START');
            else
                record = 0;	% results in repetition of fixation check
                Eyelink('Message', 'RECORD_FAILURE');
            end
        end
        
    case 'trialstart'
        %% trial start
        % start only when we are recording and the subject is fixating
        el = in{1};
        trialNum = in{2};
%         cx = in{3};
%         cy = in{4};
%         rad = in{5};
% %         
        driftCorrected = 0;
        
        
        % Start the trial only when 1) eyetracker is recording, 2) subject
        % is fixating
        ready = 0; 
%         while ~ready
%             % Check that we are recording
%             err=Eyelink('CheckRecording');
%             if err~=0
%                 rd_eyeLink('startrecording', window, el);
%             end
%             
%             % Verify that the subject is holding fixation for some set
%             % time before allowing the trial to start. A
%             % timeout period is built into this function.
%             fixation = rd_eyeLink('fixholdcheck', window, {cx, cy, rad});
%             
%             % Drift correct if fixation timed out
%             if ~fixation
%                 rd_eyeLink('driftcorrect', window, {el, cx, cy});
%                 driftCorrected = 1;
%                 ready = 0;
%             else
%                 ready = 1;
%             end
%         end
        
        %out = driftCorrected;
        out=0;
        
        Eyelink('Message', 'TRIALID %d', trialNum);
         % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d"',trialNum);
       % Eyelink('Message', 'SYNCTIME');		% zero-plot time for EDFVIEW
 
        
    case 'trialend'
%fixrect, trialnum
fixRect = in{1};
trialNum= in{2};
        % stop the recording of eye-movements for the current trial
       
        %interest area specific to the trial
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, round(fixRect(1)), round(fixRect(2)), round(fixRect(3)), round(fixRect(4)),'fix');
        %usually a count of the trial
        Eyelink('Message', '!V TRIAL_VAR index %d', trialNum);
        Eyelink('Message', '!V TRIAL_VAR iteration %d', trialNum); % Trial iteration
           
        %%give it info about what group belongs to
        %%these two lines would be updated with the condition names from
        %%setup 
        %Eyelink('Message', '!V TRIAL_VAR type %s', 'condName')

        %%give it info about what group belongs to
        %Eyelink('Message', '!V TRIAL_VAR_DATA %s', 'condName');

        % Sending a 'TRIAL_RESULT' message to mark the end of a trial in 
        % Data Viewer. This is different than the end of recording message 
        % END that is logged when the trial recording ends. The viewer will
        % not parse any messages, events, or samples that exist in the data 
        % file after this message. Can use as the end
        % message for a custom time window in data viewer
        Eyelink('Message', 'TRIAL_RESULT 0');
    case 'fixholdcheck'
        %% check that fixation is held for some amount of time
        cx = in{1}; % x coordinate of screen center
        cy = in{2}; % y coordinate of screen center
        rad = in{3}; % acceptable fixation radius %%% in px?
        
        timeout = 3.00; % 3.00 % maximum fixation check time
        tFixMin = 0.30; % 0.10 % minimum correct fixation time
        
        Eyelink('Message', 'FIX_HOLD_CHECK');
        
        tstart = GetSecs;
        fixation = 0; % is the subject fixating now?
        fixStart = 0; % has a fixation already started?
        tFix = 0; % how long has the current fixation lasted so far?
        
        t = tstart;
        while (((t-tstart) < timeout) && (tFix<=tFixMin))
            % get eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            if numel(domEye)>1 % if tracking binocularly
                domEye = domEye(1);
            end
            x = evt.gx(domEye);
            y = evt.gy(domEye);

            % check for blink
            if isempty(x) || isempty(y)
                fixation = 0;
            else
                % check eye position
                if sqrt((x-cx)^2+(y-cy)^2)<rad
                    fixation = 1;
                else
                    fixation = 0;
                end
            end
            
            % update duration of current fixation
            if fixation==1 && fixStart==0
                tFix = 0;
                tFixStart = GetSecs;
                fixStart = 1;
            elseif fixation==1 && fixStart==1
                tFix = GetSecs-tFixStart;
            else
                tFix = 0;
                fixStart = 0;
            end
            
            t = GetSecs;
        end
        
        out = fixation;
        
    case 'fixcheck'
        %% check fixation at one point in time
        cx = in{1}; % x coordinate of screen center
        cy = in{2}; % y coordinate of screen center
        rad = in{3}; % acceptable fixation radius %%% in px?
        
        % determine recorded eye
        evt = Eyelink('newestfloatsample');
        domEye = find(evt.gx ~= -32768);
        
        % if tracking binocularly, just select one eye to be dominant
        if numel(domEye)>1
            domEye = domEye(1);
        end
        
        Eyelink('Message', 'FIX_CHECK');
        
        % get eye position
        x = evt.gx(domEye);
        y = evt.gy(domEye);
        
        % check for blink
        if isempty(x) || isempty(y)
            fixation = 0;
        else
            % check eye position
            if sqrt((x-cx)^2+(y-cy)^2)<rad
                fixation = 1;
            else
                fixation = 0;
            end
        end
        
        if fixation==0
            Eyelink('Message', sprintf('BROKE_FIXATION'));
        end
        
        out = fixation;
        
    case 'driftcorrect'
        %% do a drift correction
        el = in{1};
        cx = in{2};
        cy = in{3};
        
        Eyelink('Message', 'DRIFT_CORRECTION');
        driftCorrection = EyelinkDoDriftCorrect(el, cx, cy, 1, 1);
        
        out = driftCorrection;
        
    case 'stoprecording'
        %% stop recording
        Eyelink('StopRecording');
        Eyelink('Message','RECORD_STOP');
        
    case 'eyestop'
        %% get the eye file and close down the eye tracker
        eyeFile = in{1};
        eyeDataDir = in{2};
        
        % if still recording, stop recording
        err = Eyelink('CheckRecording');
        if err==0
            rd_eyeLink('stoprecording');
        end
        Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
    Eyelink('Command', 'clear_screen 0'); % Clear Host PC backdrop graphics at the end of the experiment
    WaitSecs(0.5); % Allow some time before closing and transferring file    
    Eyelink('CloseFile'); % Close EDF file on Host PC       
  
    
        fprintf('\n\nSaving file %s/%s ...\n', eyeDataDir, eyeFile)
        
        Eyelink('ReceiveFile', eyeFile, eyeDataDir, 1); 
        Eyelink('CloseFile'); 
        Eyelink('Shutdown');
        
    otherwise
        error('[rd_eyeLink]: ''command'' argument not recognized. See help for available commands.')
end
