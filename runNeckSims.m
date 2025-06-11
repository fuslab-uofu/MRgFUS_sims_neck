% =============================================================================
%               runNeckSims script
% -----------------------------------------------------------------------------
%
% Description:  Driver script for running acoustic and thermal simulations on models 
%               for MRgFUS neck pain study. This script is designed to run simulations 
%               on a neck model, breast and bone transducers, and multiple targets and angles. 
%               Uses the hybrid angular spectrum method (HAS) for acoustic simulation and
%               _________ for thermal simulation. 
%
% Authors:      Michelle Kline & Marta M. Iversen  
%               Department of Radiology and Imaging Sciences  
%               University of Utah  
%
% Inputs:       This script must be updated where indicated with file names, 
%               folder paths, and parameters.
%               1. segmented neck model
%               2. simulated pressure source ('ERFA') files for breast and bone transducers
%               3. desired targets for simulated sonications
%
% Saves files:  1. pressure ('P') - 3D complex single 
%               2. power deposition ('Q') - 3D single
%               3. temperatures ('T') - 4D single
%               
% Output:       Vector of results for each simulation:
%               'modelID'   - name of the segmented neck model
%               'ERFAID'    - name of the simulated transducer (ERFA)
%               'targetID'  - name of the target
%               'targetRCP' - target indices in MATLAB row, col, page
%               'angle'     - rotation of simulated transducer
%               'maxP'      - maximum pressure
%               'maxP_RCP'  - max pressure voxel indices
%               'targetP'   - pressure at target voxel
%               'distP'     - Euclidean distance, target voxel to max pressure voxel
%               'P_fileName'- full path to the saved pressure file
%               'maxQ'      - maximum Q 
%               'maxQ_RCP'  - max Q voxel indices
%               'targetQ'   - Q at target voxel
%               'distQ'     - Euclidean distance, target voxel to max Q voxel
%               'Q_fileName'- full path to the saved Q file
%               'metadata'  - HAS-related metadata
%
% Requirements: MATLAB R2020a or later (toolboxes?)  
%
% Dependencies: 1. runAcousticSim_neck.m 
%               2. plotPQSlice_neck.m
%               2. runThermalSim_neck.m  
%
% =============================================================================


%% SETUP: MODIFY AS NEEDED

studyDir = 'D:\HUMAN\IRB137036_Rieke_Shah\neck_system_paper';

% segmented neck model
% NOTE: HAS expects the model to be a 3D matrix of integers specifically named 'Modl' (not a type-o)
% Each integer corresponds to one medium, values must be contiguous and range from 1 to number of media.
modelDir = fullfile(studyDir, '384957\MichelleCodeTests');
modelID = 'Model_384957'; 
load(fullfile(modelDir, modelID), 'Modl');
modelResolution_rcp = [0.5 0.5 0.5]; % row, col, page

% coordinates of points of interest (target and offtarget locations) for
% all five levels (row, col, page indices in segmented model)
load(fullfile(modelDir, 'trgt_offtrgt_coords.mat'), 'levelPOI');
target_ids = {'target_left', 'target_right'};
otherPOI_ids = {'csf', 'nerve_left', 'nerve_right', 'artery_left', 'artery_right', 'sc'};


% electronic steering, xyz (mm)
steering = [0 0 0]; 

% this is transducer rotation (psi, degrees)
% in this case, negative for left-side targets
% eight sets of three angles (one set for each target) 
angles = [0 15 30];

% simulated pressure input (ERFA file)
ERFA_directory = fullfile(studyDir, 'ERFA');
breast_transducer.fileRoot = 'ERFA_A10.mat';
breast_transducer.ERFA = load(fullfile(ERFA_directory, breast_transducer.fileRoot));
breast_transducer.id = 'breast';
bone_transducer.fileRoot = 'ERFA_S11.mat';
bone_transducer.ERFA = load(fullfile(ERFA_directory, bone_transducer.fileRoot));
bone_transducer.id = 'bone';

% output directories
acousticSaveDir = fullfile(modelDir, 'pressures');
thermalSaveDir = fullfile(modelDir, 'temperatures');


%% ACOUSTIC SIMULATIONS

disp("Running acoustic simulations...");
fprintf("\n");

% compute total number of simulations to run
nTargets = length(levelPOI);
nAngles  = length(angles);
nTrans   = 2;  
nTotal   = nTargets * nAngles * nTrans;

% define an "empty” template struct with exactly the fields that runAcousticSim returns
emptySim = struct( ...
    'modelID',    [], ...
    'modelFile',  [], ...
    'ERFAID',     [], ...
    'ERFAFile',   [], ...
    'targetID',   [], ...
    'targetRCP',  [], ...
    'angle',      [], ...
    'maxP',       [], ...
    'maxP_RCP',   [], ...
    'targetP',    [], ...
    'distP',      [], ...
    'P_fileName', [], ...
    'maxQ',       [], ...
    'maxQ_RCP',   [], ...
    'targetQ',    [], ...
    'distQ',      [], ...
    'Q_fileName', [], ...
    'metadata',   [] ...
);

% preallocate the vector of structs to length nTotal 
acousticSimResults = repmat(emptySim, nTotal, 1);

% run the simulation at specified targets using specified transducers at specified rotations
idx = 1;  % index for results output vector
for level = 1:length(levelPOI)
    for trgt_num = 1:length(target_ids)  % one left and one right at each level
        target_id = target_ids{trgt_num};
        target_RCP = levelPOI(level).(target_id);
        for angleNum = 1:length(angles)
            angle = angles(angleNum);
            if strcmp(target_id, 'target_left')
                angle = -1 * angle;
            end
            for transducer = [breast_transducer, bone_transducer]
                acousticSimResults(idx).modelID = modelID;
                acousticSimResults(idx).modelFile = fullfile(modelDir, modelID);
                acousticSimResults(idx).ERFAID = transducer.id;
                acousticSimResults(idx).targetID = target_id;
                acousticSimResults(idx).targetRCP = target_RCP;
                acousticSimResults(idx).angle = angle;
                % show simulation input info on command line
                displayString = [
                    "Model:      "   + modelID
                    "Transducer: "   + transducer.id
                    "Rotation:   "   + string(angle) + "°"
                    "Target:     "   + target_id + " " + ...
                                     strjoin(string(target_RCP)," ")
                ];
                for k = 1:numel(displayString)
                    fprintf("%s\n", displayString(k));
                end
    
                [P, Q, resultsStruct] = runAcousticSim_neck( ...
                    transducer.ERFA, ...
                    transducer.ERFA.R*1000, ...  % R is the transducer focal length in meters
                    Modl, ...
                    modelResolution_rcp, ...
                    target_RCP, ...
                    angle, ...
                    steering);
    
                acousticSimResults(idx).maxP = resultsStruct.maxP;
                acousticSimResults(idx).maxP_RCP = resultsStruct.maxP_RCP;
                acousticSimResults(idx).targetP = resultsStruct.targetP;
                acousticSimResults(idx).distP = resultsStruct.distP;
                acousticSimResults(idx).maxQ = resultsStruct.maxQ;
                acousticSimResults(idx).maxQ_RCP = resultsStruct.maxQ_RCP;
                acousticSimResults(idx).targetQ = resultsStruct.targetQ;
                acousticSimResults(idx).distQ = resultsStruct.distQ;
                acousticSimResults(idx).metadata = resultsStruct.metadata;
    
                % save pressure, Q, and metadata
                disp("Saving pressure and Q...")
                meta = acousticSimResults(idx).metadata;
                saveString = [modelID, '_', target_id, '_psi_' num2str(angle), ...
                    '_str_' 'x' num2str(steering(1)) 'y' num2str(steering(2)) 'z' num2str(steering(3)) 'mm' , ...
                    '_', transducer.id];
                qFileName = ['Q_',saveString,'.mat'];
                save(fullfile(acousticSaveDir, qFileName), 'Q', 'meta');
                acousticSimResults(idx).Q_fileName = fullfile(acousticSaveDir, qFileName);
                pFileName = ['P_', saveString,'.mat'];
                save(fullfile(acousticSaveDir, pFileName),'P', 'meta');
                acousticSimResults(idx).P_fileName = fullfile(acousticSaveDir, pFileName);
    
                % plot pressure and Q slice at target voxel
                % for better visualization, map images to a clipped colormap range
                % adjust as necessary
                lowerDisplayLimit_Q = 1e+06;
                if acousticSimResults(idx).maxQ <= lowerDisplayLimit_Q
                    lowerDisplayLimit_Q = acousticSimResults(idx).maxQ / 2;
                end
                lowerDisplayLimit_P = 1e+06;
                if acousticSimResults(idx).maxP <= lowerDisplayLimit_P
                    lowerDisplayLimit_P = acousticSimResults(idx).maxP / 2;
                end
    
                % plot Q slice at max and at target ---------- 
                titleStr_max = sprintf('Q %0.5e W/M^3 at \n Max Voxel (%dr,%dc,%dp) \nPsi %d°', ...
                    acousticSimResults(idx).maxQ, ...
                    acousticSimResults(idx).maxQ_RCP(1), ...
                    acousticSimResults(idx).maxQ_RCP(2), ...
                    acousticSimResults(idx).maxQ_RCP(3), ...
                    acousticSimResults(idx).angle);
                titleStr_target = sprintf('Q %0.5e W/M^3 at \n Target Voxel (%dr,%dc,%dp) \nPsi %d°', ...
                    acousticSimResults(idx).targetQ, ...
                    acousticSimResults(idx).targetRCP(1), ...
                    acousticSimResults(idx).targetRCP(2), ...
                    acousticSimResults(idx).targetRCP(3), ...
                    acousticSimResults(idx).angle);
                % coords of target (either left or right) at this level
                target_XY = [acousticSimResults(idx).targetRCP(3), acousticSimResults(idx).targetRCP(1)];
                % coords of max voxel at this level
                maxQ_XY = [acousticSimResults(idx).maxQ_RCP(3), acousticSimResults(idx).maxQ_RCP(1)];
                % coords of off-target points of interest at this level
                POI_XY = ones(length(otherPOI_ids), 2);
                for POI_num = 1:length(otherPOI_ids)
                    POI_XY(POI_num + 1, :) = [levelPOI(level).(otherPOI_ids{POI_num})(3), ...
                        levelPOI(level).(otherPOI_ids{POI_num})(1)];
                end
                
                QFigHandle = plotPQSlice_neck( ...
                    squeeze(Modl(:, acousticSimResults(idx).maxQ_RCP(2), :)), ...   % model slice at max
                    squeeze(Q(:, acousticSimResults(idx).maxQ_RCP(2), :)), ...      % Q slice at max
                    titleStr_max, ...
                    squeeze(Modl(:, acousticSimResults(idx).targetRCP(2), :)), ...  % model slice at target
                    squeeze(Q(:, acousticSimResults(idx).targetRCP(2), :)), ...     % Q at target
                    titleStr_target, ...
                    [lowerDisplayLimit_Q acousticSimResults(idx).maxQ], ...         
                    maxQ_XY, ...
                    target_XY, ...
                    POI_XY);
                QFigHandle.Name = ['Model  ', acousticSimResults(idx).modelID, ...
                    ', Target  ', acousticSimResults(idx).targetID, ...
                    ' Transducer  ', acousticSimResults(idx).ERFAID];
                
                % plot pressure slice at max and target ----------
                titleStr_max = sprintf('Pressure %0.5e Pa at \n Max Voxel (%dr,%dc,%dp) \nPsi %d°', ...
                    acousticSimResults(idx).maxP, ...
                    acousticSimResults(idx).maxP_RCP(1), ...
                    acousticSimResults(idx).maxP_RCP(2), ...
                    acousticSimResults(idx).maxP_RCP(3), ...
                    acousticSimResults(idx).angle);
                titleStr_target = sprintf('Pressure %0.5e Pa at \n Target Voxel (%dr,%dc,%dp) \nPsi %d°', ...
                    acousticSimResults(idx).targetP, ...
                    acousticSimResults(idx).targetRCP(1), ...
                    acousticSimResults(idx).targetRCP(2), ...
                    acousticSimResults(idx).targetRCP(3), ...
                    acousticSimResults(idx).angle);
                 % coords of max voxel at this level
                maxP_XY = [acousticSimResults(idx).maxP_RCP(3), acousticSimResults(idx).maxP_RCP(1)];
                PFigHandle = plotPQSlice_neck( ...
                    squeeze(Modl(:, acousticSimResults(idx).maxP_RCP(2),:)), ...    % model slice at max
                    squeeze(abs(P(:,acousticSimResults(idx).maxP_RCP(2),:))),...    % P slice at max
                    titleStr_max, ...
                    squeeze(Modl(:, acousticSimResults(idx).targetRCP(2),:)), ...   % model slice at target
                    squeeze(abs(P(:,acousticSimResults(idx).targetRCP(2),:))),...   % P slice at target
                    titleStr_target, ...
                    [lowerDisplayLimit_P acousticSimResults(idx).maxP], ...
                    maxP_XY, ...
                    target_XY, ...
                    POI_XY);
                PFigHandle.Name = ['Model ', acousticSimResults(idx).modelID, ...
                    ', Target ', acousticSimResults(idx).targetID, ...
                    ' Transducer  ', acousticSimResults(idx).ERFAID];
    
                fprintf("Complete\n\n");
                idx = idx + 1;
            end
        end
    end
end
   

%% THERMAL simulations

% **assumes acousticSimResults struct is still in memory**
disp("Running thermal simulations...");
fprintf("\n");

templateStruct = struct( ...
    'modelID',          modelID, ...
    'modelFile',        fullfile(modelDir, modelID), ...
    'ERFAID',           [], ...
    'ERFAFile',         [], ...
    'targetID',         [], ...
    'targetRCP',        [], ...
    'angle',            [], ...
    'Q_fileName',       [], ...
    'reqTemp',          [], ...
    'tEnd',             [], ...
    'power',            [], ...
    'timePoints',       [], ...
    'maxT',             [], ...
    'maxT_RCP',         [], ...
    'hotVoxTemps',      [], ...
    'hotVolTemps',      [], ...
    'targetT',          [], ...
    'distT',            [], ...
    'targetVoxTemps',   [], ...
    'targetVolTemps',   [], ...
    'averageTissueT',   [], ...
    'heatingSuccess',   [], ...
    'T_fileName',       [], ...
    'TDose_fileName',   [] ...
);
% preallocate the vector of structs to length nTotal 
thermalSimResults = repmat(templateStruct, length(acousticSimResults), 1);

for simNum = 1:length(acousticSimResults)
    thermalSimResults(simNum).ERFAID = acousticSimResults(simNum).ERFAID;
    thermalSimResults(simNum).ERFAFile = acousticSimResults(simNum).ERFAFile;
    thermalSimResults(simNum).targetID = acousticSimResults(simNum).targetID;
    thermalSimResults(simNum).targetRCP = acousticSimResults(simNum).targetRCP;
    thermalSimResults(simNum).angle = acousticSimResults(simNum).angle;

    % load the Q file we are using for this simulation
    load(acousticSimResults(simNum).Q_fileName, 'Q');
    thermalSimResults(simNum).Q_fileName = acousticSimResults(simNum).Q_fileName;
    disp(acousticSimResults(simNum).Q_fileName)

    % set scale factor for 100W power (e.g., scale factor 1.5 = 150W)
    QScaleFactor = 1;
    thermalSimResults(simNum).power = QScaleFactor * 100;   % FIXME - is 100 a constant?

    % FIXME - description?
    thermalSimResults(simNum).reqTemp = 30;
    thermalSimResults(simNum).tEnd = 20;

    calcThermalDose = true;

    % show simulation input info on command line
    displayString = [
        "Model:         "   + modelID
        "Transducer:    "   + thermalSimResults(simNum).ERFAID
        "Rotation:      "   + string(thermalSimResults(simNum).angle) + "°"
        "Target:        "   + thermalSimResults(simNum).targetID + " " + strjoin(string(thermalSimResults(simNum).targetRCP)," ")
        "Power:         "   + string(thermalSimResults(simNum).power)
        "Required Temp: "   + string(thermalSimResults(simNum).reqTemp)
        "Time:          "   + string(thermalSimResults(simNum).tEnd)
    ];
    for k = 1:numel(displayString)
        fprintf("%s\n", displayString(k));
    end

    % run thermal simulation
    [T, TDose, resultsStruct] = runThermalSim_neck( ...
        Modl, ...
        modelResolution_rcp, ...
        Q, ...
        thermalSimResults(simNum).targetRCP, ...
        QScaleFactor, ...
        thermalSimResults(simNum).reqTemp, ...
        thermalSimResults(simNum).tEnd, ...
        calcThermalDose);  % calculate thermal dose if true

    clear Q;
    thermalSimResults(simNum).timePoints = resultsStruct.timePoints;
    thermalSimResults(simNum).maxT = resultsStruct.maxT;
    thermalSimResults(simNum).maxT_RCP = resultsStruct.maxT_RCP;
    thermalSimResults(simNum).hotVoxTemps = resultsStruct.hotVoxTemps;
    thermalSimResults(simNum).hotVolTemps = resultsStruct.hotVolTemps;
    thermalSimResults(simNum).targetT = resultsStruct.targetT;
    thermalSimResults(simNum).distT = resultsStruct.distT;
    thermalSimResults(simNum).targetVoxTemps = resultsStruct.targetVoxTemps;
    thermalSimResults(simNum).targetVolTemps = resultsStruct.targetVolTemps;
    thermalSimResults(simNum).averateTissueT = resultsStruct.averageMediaT;
    thermalSimResults(simNum).heatingSuccess = resultsStruct.failureToHeat;

    disp("Saving temperatures...")
    % these files can be quite large!
    tempsFile = fullfile(thermalSaveDir, ...
       ['T_', modelID,'_W',num2str(100*QScaleFactor), '_max', num2str(thermalSimResults(simNum).maxT, 4),'.mat']);
%     save(tempsFile, 'T', '-v7.3');   
    thermalSimResults(simNum).T_fileName = tempsFile;

    if (calcThermalDose && ~isempty(TDose))
        disp("Saving thermal dose...")
        TDFile = fullfile(thermalSaveDir, ...
           ['TDose_', modelID,'_W',num2str(100*QScaleFactor), '_max', num2str(thermalSimResults(simNum).maxT, 4),'.mat']);
%         save(TDFile, 'TDose', '-v7.3');   
        thermalSimResults(simNum).TDose_fileName = TDFile;
    end

    titleStr_target = sprintf('T %.2f°C at Time %.1f S \n Target Voxel (%dr,%dc,%dp) \n and Psi %d°', ...
        thermalSimResults(simNum).targetT, ...
        thermalSimResults(simNum).tEnd, ...
        thermalSimResults(simNum).targetRCP(1), ...
        thermalSimResults(simNum).targetRCP(2), ...
        thermalSimResults(simNum).targetRCP(3), ...
        thermalSimResults(simNum).angle);

        % coords of target (either left or right) at this level
            target_XY = [thermalSimResults(simNum).targetRCP(3), thermalSimResults(simNum).targetRCP(1)];
            % coords of off-target points of interest at this level
            POI_XY = ones(length(otherPOI_ids), 2);
            for POI_num = 1:length(otherPOI_ids)
                POI_XY(POI_num + 1, :) = [levelPOI(level).(otherPOI_ids{POI_num})(3), ...
                    levelPOI(level).(otherPOI_ids{POI_num})(1)];
            end
            figHandle = plotTSlice_neck( ...
                squeeze(Modl(:, thermalSimResults(simNum).targetRCP(2), :)), ...    % model slice at target
                squeeze(T(:, thermalSimResults(simNum).targetRCP(2), :, end)), ...  % T slice at target
                titleStr_target, ...
                [2 40], ...                                                         % display limits (°C)         
                target_XY,...
                POI_XY);
            figHandle.Name = ['Model  ', thermalSimResults(simNum).modelID, ...
                ', Target  ', thermalSimResults(simNum).targetID, ...
                ' Transducer  ', thermalSimResults(simNum).ERFAID];

    fprintf("Complete\n\n");
end
