% Created by A. Benitez Andonegui, M. Luehrs
% Modified by A. Vorreuther

%%%%% include these %%%%%
% mainPath: path where all participant data is stored. assumed directory format is: 
%            mainPath > P01 > run 1 
%                           > run 2
%                           > ...

%                     > P02 > run 1 
%                           > run 2
%                           > ...
mainPath = 'C:\B_Yes-No_temporal_reanalysis\matlabanalysis\setup b';

% allSubj_Answers: true answers from participants, dimensions: NrOfAnswers (one per run) x NrOfSubj
allSubj_Answers = [ 2 1 1 1 2 2 1 1 1 1;
                    1 2 1 2 2 1 2 1 2 1; 
                    1 2 1 2 2 2 1 2 1 1;
                    1 2 1 2 1 2 1 1 1 1;
                    1 1 2 1 2 1 2 1 1 1];
                
allSubj_Answers = allSubj_Answers';

%%% dependencies
% this code assumes you have installed:
% - neuroelf: https://github.com/neuroelf/neuroelf-matlab 
% - spm12 : https://www.fil.ion.ucl.ac.uk/spm/
% and uses in-house functions to load nirx files onto Matlab. They can be found here: https://github.com/AmaiaBA/TSI-Input-Ouput
addpath(genpath(fullfile(fileparts(pwd), ...
    "tempEncPaper_matlabcode\")));


%%%% modify these %%%%%

SubjNames = {'P04_2019-02-19','P05_2019-02-19','P09_ 2019-01-15','P101_2019-02-19','P102_2019-02-19'};

NrOfCond            = 2; % # of answer options
TaskDuration        = 2; % in seconds
NrOfRep             = 10; % # of times a trial will be presented (full cycle of answer options)
GLMApproach         = 2; % 1 = compute GLM after each condition; 2 = compute GLM when all conditions have been shown
WeightingMethod     = 2; % possible values: 0-4; 0 = weighted sum of  all 4 values [beta tVal HRFcorrVal rVal]; 1 = betaVAL; 2 = tVAL; 3 = HRFcorrVal; 4 = rVAL;
MultiTrialApproach  = 2; % only  "2" allowed at the moment;  0 none; 1 = multi trial (based on single trial, voted majority) to extract beta, tval, corr values, in seconds; 2 = using whole timecourse
OnOffsetWindow      = [5 10]; % window with respect to task onset (1) and offset (2) to extract beta, tval, corr values, in seconds. Always >= 0 
prt_IncludesRest    = 0; % 0 = no; 1 = yes;
BaselinePeriodEnd   = 35; % baseline calculation period end (for concentration changes), in seconds 
showResults         = 0; % 0 = don´t display any results in command window; 1: display decoded value only; set to 2 to show all BCI ouput

selectedOutputType  = ones(numel(SubjNames),1).*2; % 1= HbO; 2 = HbR; 3 = HbT
useSDCasCP          = 1; % 1 if SDC was used as Confound Predictor; 0 if not; dimensions: NrOfSubj x 1
SDC                 = [9 5]; % [source Nr, detector Nr], can contain multiple SDCs (should go in different rows); [] if no SDC was used
NrOfChan            = 31; % without SDC, if it was used. 
NrOfRuns            = ones(numel(SubjNames),1).*10; % excluding localizer run


%%%%%%%%%%%%


%%%% Run

allSubj_Acc = zeros(numel(SubjNames),NrOfChan,NrOfRep); % the third dimension will contain the accuracies when NrOfRep:-1:1 are used for decoding


for ss=1:numel(SubjNames) % loop over participants

    disp(SubjNames{ss})
    disp(' ')
    OutputHb = selectedOutputType(ss); % 1 = HbO; 2= HbR; 3 = HbT;
    selSD = []; % selected channel | col1 = Sources; col2 = Detectors; if empty, output of all channels will be provided

    % get folder locations in directory
    subjPath = fullfile(mainPath,SubjNames{ss});
    tempDir = dir(subjPath);
    posRuns = [];
    counter = 1;
    for int=1:numel(tempDir)
        if isdir(fullfile(subjPath,tempDir(int).name))
            if ~(strcmp(tempDir(int).name,'.') || strcmp(tempDir(int).name,'..'))
                posRuns = [posRuns; int];
                disp([num2str(counter) ' - ' fullfile(subjPath,tempDir(int).name)])
                counter = counter + 1;
            end
        end
    end
    
    % get indices of interest
    if numel(posRuns) ~= NrOfRuns(ss)+1 % to include localizer run
        disp(' ')
%         disp(' WARNING!!! The # of runs in directory does not correspond with # of runs in experiment')
%         indxRun = input(' select run indices, excluding localizer, e.g., 2:5, [2 5 6 7]: ');
        indxRun = 1:10;
    else
        indxRun = 2:NrOfRuns(ss)+1; % assumes first run is localizer
    end

    allDecisionMat =  cell(numel(indxRun),NrOfRep,NrOfChan);

    for qq=1:numel(indxRun) % loop over runs

        disp(['Run # ' num2str((qq))])
        hdrPath = fullfile(subjPath,tempDir(posRuns(indxRun(qq))).name);% define path where hdr file is saved
        ntcPath = hdrPath;% define path where ntc file is saved
        prtPath = hdrPath;% define path where prt file is saved

        [~,NAME,~] = fileparts(hdrPath); % this assumes that run folder name is the same as the files in the folder
        hdrFileName = ['NIRS-', NAME , '.hdr'];
        ntcFileName = ['NIRS-', NAME , '_filtered_TSI.ntc']; % derived from TSI
        prtFileName = ['NIRS-', NAME , '.prt']; % derived from TSI
        
        %%% load files

        [~, ~, ~, SamplingRate, sd_ind, ~,~, ~, ~, ~,~,~]= HDRFile_extractInfo(hdrPath, hdrFileName); % load hdr file
         if qq == 2 && ss == 5
            fileIDBin = fopen(fullfile(ntcPath,ntcFileName));
            fileversion = fread(fileIDBin,[1], 'uint16');
            NrOfSources = fread(fileIDBin,[1],'uint32');
            NrOfDetectors = fread(fileIDBin,[1],'uint32');
            ntcDataPoints = fread(fileIDBin,[1],'uint32');
            NrOfHbTypes = 3;
            ntcData = fread(fileIDBin,[2106 NrOfDetectors*NrOfSources*NrOfHbTypes], 'double');
        else
            [~, ntcData] = readNTCFile(ntcPath, ntcFileName); % load ntc file
        end
        tmp = fullfile(prtPath,prtFileName);
        prt_file = xff(tmp); % load prt file
%         prt_file = xff('*.prt');
        %%% define values for output

        if OutputHb==1
            outputName = 'OXY';
        elseif OutputHb==2
            outputName = 'DEOXY';
        else
            outputName = 'TOTAL';
        end

        if WeightingMethod == 0
            weightName = 'weighted sum';
        elseif WeightingMethod == 1
            weightName = 'betaVal'; 
        elseif WeightingMethod == 2
            weightName = 'tVal';
        elseif WeightingMethod == 3
            weightName = 'bcorrVal';    
        elseif WeightingMethod == 4    
            weightName = 'rVal';
        end

        if MultiTrialApproach == 0
            weightNameMT = 'none';
            showResults = 0 ; % disabled if not multitrial
        elseif MultiTrialApproach == 1
            weightNameMT = 'Weighted single trial';
        elseif MultiTrialApproach == 2
            weightNameMT = 'Timecourse';
        end

        %%% define position of SDC, if used as confound predictor
        
        % get source-detector indices
        maskSetup = sd_ind; % rows = detectors; columns = sources;
        SelCh = find(maskSetup(:)==1); % find channels used in measurement
        [detInd,sourInd] = ind2sub(size(maskSetup),SelCh); % get indices of source and detectors
        SDmat = [sourInd detInd];
        if useSDCasCP % if a SDC was used
            [~, temp] = ismember(SDC,SDmat,'rows'); % find position of SDC in the SD matrix
            SelSDC = SelCh(temp); % find position of SDCs in ntc file
            SDmat(temp,:) = []; % remove SDC from SDmat 
            SelCh(temp)=[]; % remove SDC index from selected channel indices (for ntc file)
        end
        NrOfSelCh = NrOfChan;
        
        %%% extract rest/task onoffset values from prt file

        task_onoffsets  = [];
        if prt_IncludesRest==0
            startVal = 1;
        else
            startVal = 2; % ignore rest condition
        end

        % extract onset and offset values for each condition
        for int=startVal:2
            task_onoffsets = [task_onoffsets;prt_file.Cond(int).OnOffsets, repmat(int,[size(prt_file.Cond(int).OnOffsets,1),1])];    
        end
        task_onoffsets = sortrows(task_onoffsets);

        if prt_file.ResolutionOfTime == 'Volumes' % this is equivalent to frames
            task_onoffsetsFR = task_onoffsets;
        else
            task_onoffsetsFR = floor(task_onoffsets.*SamplingRate); % put it in frames
        end
        task_onoffsetsFR = task_onoffsets;
        
        task_onoffsetsFR_orig = task_onoffsetsFR;
        task_onset = floor(OnOffsetWindow(1)*SamplingRate); % adapt to OnOffsetWindow, in frames
        task_offset = floor(OnOffsetWindow(2)*SamplingRate); % adapt to OnOffsetWindow, in frames

        if GLMApproach ==2 % only select the onset/offset of first condition
            trialStartArr = task_onoffsetsFR(1:NrOfCond:size(task_onoffsetsFR,1),1)- task_onset;
            trialStopArr = task_onoffsetsFR(NrOfCond:NrOfCond:size(task_onoffsetsFR,1),2)+ task_offset;
            ed_task_onoffsetsFR_orig = [trialStartArr trialStopArr];
        end
        
        %%% include descending # of trials

        for lengthRunTrials=0:NrOfRep-1 % this will reduce the number of trials in a run

            NrOfRep2 = NrOfRep - lengthRunTrials;
            indTrialsMat = getConsecutiveTrials(1:NrOfRep,NrOfRep2); % get consecutive number of trials based on current trial/repetition number
            disp(['# Of trials: ' num2str(NrOfRep2)])
            tempDecision = zeros(size(indTrialsMat,1),NrOfSelCh);

            for comb=1:size(indTrialsMat,1)
                startTrial = indTrialsMat(comb,1); % indx of first trial/repetition
                stopTrial= indTrialsMat(comb,size(indTrialsMat,2)); % indx of last trial/repetition
                ed_task_onoffsetsFR = ed_task_onoffsetsFR_orig(startTrial:stopTrial,:); % onset/offset of trial/repetition of interest
                task_onoffsetsFR = task_onoffsetsFR_orig((NrOfCond*(startTrial-1))+1:NrOfCond*stopTrial,:);
                
                %%% define output structure

                if MultiTrialApproach > 0
                    dim2 = NrOfRep2+1;    
                else
                    dim2 = NrOfRep2;    
                end
                trialNames = cell(1,dim2);
                for int=1:NrOfRep2
                   trialNames{int} = ['Trial' num2str(int)]; 
                end
                if MultiTrialApproach > 0
                    trialNames{end} = 'Multitrial';   
                end
                trial = trialNames;

            
                if NrOfSelCh > 1
                    SelChID = cell(1,NrOfSelCh+1); % inlcude average of all selected channels
                else
                    SelChID = cell(1,NrOfSelCh);
                end

                for int=1:NrOfSelCh
                    SelChID{int} =['S' num2str(SDmat(int,1)) '_D' num2str(SDmat(int,2))];
                end
                if NrOfSelCh >1
                    SelChID{end} = 'Mean';    
                end

                channel = SelChID;

                Results = [];
                for int=1:dim2
                    for chan = 1:size(SelChID,2)
                        Results.(trial{int}).(channel{chan}).BCIOutput = [];
                        Results.(trial{int}).(channel{chan}).Decision = [];
                    end
                end       

                %%% output selected choices

                if qq==1 && lengthRunTrials == 0 && comb==1
                    disp(['Selected output Hb type : ' outputName ])
                    disp(['Selected # of channels : ' num2str(NrOfSelCh)])
                    disp(['Selected channel IDs : ' num2str(SelCh')])
                    disp(['Selected feature : ' weightName ])
                    disp(['Selected multi trial approach : ' weightNameMT])
                    disp(' ')
                end
            
                %%% actual computations %%%

                % initialize variables

                TheoBetaMat = eye(2); % for correlation computations
                currRep = 1;
                if GLMApproach==2
                    NrTrials = NrOfRep2; 
                else
                    NrTrials = size(task_onoffsets,1); % nr of task trials
                end

                if NrOfSelCh>1 % include average values for each hb type
                    selCond = zeros(NrOfRep2,NrOfSelCh,6); % [HbO HbR HbT mean(HbO) mean(HbR) mean(HbT)]
                    selCondMultiTrial = zeros(NrOfSelCh,6); 
                else
                    selCond = zeros(NrOfRep2,NrOfSelCh,3); % [HbO HbR HbT]
                    selCondMultiTrial = zeros(NrOfSelCh,3); 
                end

                % create design matrix
                lastFrame = ed_task_onoffsetsFR(end,2)+ floor(5*SamplingRate); % for design matrix definition
                dm = zeros(lastFrame,NrOfCond); 
                for count=1:NrOfCond
                    tempInd = find(task_onoffsetsFR(:,3)==count);
                    for ll=1:numel(tempInd)
                        dm(task_onoffsetsFR(tempInd(ll),1):task_onoffsetsFR(tempInd(ll),2),count)=1; 
                    end
                end  

                % figure, imagesc(dm), hold on, pause()
                
                
                % initializevariables to save [betaVal, tVal, corrVal,rVal] values for each channel and condition
                
                AllResOxy = zeros(NrOfCond,4,NrOfSelCh);
                AllResDeoxy = zeros(NrOfCond,4,NrOfSelCh);
                AllResTotalHb = zeros(NrOfCond,4,NrOfSelCh);

                AllOxyData = zeros(ed_task_onoffsetsFR(end,2),NrOfSelCh);
                AllDeoxyData = zeros(ed_task_onoffsetsFR(end,2),NrOfSelCh);
                AllTotalHbData = zeros(ed_task_onoffsetsFR(end,2),NrOfSelCh);

                isMultiTrial = 1;
                for cc=1:NrOfSelCh

                    tempOxyData = ntcData(1:ed_task_onoffsetsFR(end,2),SelCh(cc)); % load all oxy data until frame of interest
                    tempDeoxyData = ntcData(1:ed_task_onoffsetsFR(end,2),SelCh(cc)+numel(maskSetup)); % HbR
                    tempTotalHbData = tempOxyData-tempDeoxyData; % HbT

                    AllOxyData(:,cc) = tempOxyData;
                    AllDeoxyData(:,cc) = tempDeoxyData; 
                    AllTotalHbData(:,cc) = tempTotalHbData;							

                end

                if useSDCasCP % extract SDC timecourse
                    SDCtempOxy =  ntcData(1:ed_task_onoffsetsFR(end,2),SelSDC);
                    SDCtempDeoxy =  ntcData(1:ed_task_onoffsetsFR(end,2),SelSDC+numel(maskSetup));
                    SDCtotalHbData = SDCtempOxy-SDCtempDeoxy;
                end

                if NrOfSelCh >1 % define variables to include average values
                    AllMeanResOxy = zeros(NrOfCond,4);
                    AllMeanResDeoxy = zeros(NrOfCond,4);
                    AllMeanResTotalHb = zeros(NrOfCond,4);
                    tempMeanOxy = mean(AllOxyData,2);
                    tempMeanDeoxy = mean(AllDeoxyData,2);
                    tempMeanTotalHb = mean(AllTotalHbData,2);
                end

                % select start/stop frame based on current selected
                % trials/repetition
                if NrOfRep2 == NrOfRep
                    startFrame = ceil(BaselinePeriodEnd*SamplingRate);
                else
                    startFrame = ed_task_onoffsetsFR(1,1);
                end
                
                stopFrame = ed_task_onoffsetsFR(end,2);
                
                % define indices for decoding based on GLM
                if GLMApproach ==2
                    tempCond = NrOfCond;
                else
                    tempCond = 1;
                end
                tempIndx = 1;

                for cond=tempCond:NrOfCond	
                    for cc=1:NrOfSelCh
                        if useSDCasCP % include SDC timecourse as CP
                            tempDM = dm;
                            designMat_HbO = [tempDM [SDCtempOxy;zeros(size(dm,1)-size(SDCtempOxy,1),1)]];
                            designMat_HbR = [tempDM [SDCtempDeoxy;zeros(size(dm,1)-size(SDCtempDeoxy,1),1)]];
                            designMat_HbT = [tempDM [SDCtotalHbData;zeros(size(dm,1)-size(SDCtotalHbData,1),1)]];
                        else
                            tempDM = dm;
                            designMat_HbO = tempDM;
                            designMat_HbR = tempDM;
                            designMat_HbT = tempDM;
                        end
                        
                        % extract beta, tval and corrHRF values HbO
                        % selected timecourse range depends on number of NrOfRep2 and comb variables 
                        [~,stats,RValOxy] =computeGLMandRval(designMat_HbO,AllOxyData(startFrame:stopFrame,cc),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);  
                        AllResOxy(tempIndx:cond,1,cc) = stats.beta(2:end-1,:);
                        AllResOxy(tempIndx:cond,2,cc) = stats.t(2:end-1,:);
                        AllResOxy(tempIndx:cond,4,cc) = RValOxy;

                        % extract beta, tval and corrHRF values HbR
                        [~,stats,RValDeoxy] =computeGLMandRval(designMat_HbR,-1*AllDeoxyData(startFrame:stopFrame,cc),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);
                        AllResDeoxy(tempIndx:cond,1,cc) = stats.beta(2:end-1,:);
                        AllResDeoxy(tempIndx:cond,2,cc) = stats.t(2:end-1,:);
                        AllResDeoxy(tempIndx:cond,4,cc) = RValDeoxy;

                        % extract beta, tval and corrHRF values HbT
                        [~,stats,RValTotalHb] = computeGLMandRval(designMat_HbT,AllTotalHbData(startFrame:stopFrame,cc),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);
                        AllResTotalHb(tempIndx:cond,1,cc) = stats.beta(2:end-1,:);
                        AllResTotalHb(tempIndx:cond,2,cc) = stats.t(2:end-1,:);
                        AllResTotalHb(tempIndx:cond,4,cc) = RValTotalHb;

                    end


                    if NrOfSelCh>1
                        if useSDCasCP
                            tempDM = dm;
                            designMat_HbO = [tempDM [SDCtempOxy;zeros(size(dm,1)-size(SDCtempOxy,1),1)]];
                            designMat_HbR = [tempDM [SDCtempDeoxy;zeros(size(dm,1)-size(SDCtempDeoxy,1),1)]];
                            designMat_HbT = [tempDM [SDCtotalHbData;zeros(size(dm,1)-size(SDCtotalHbData,1),1)]];
                        else
                            tempDM = dm;
                            designMat_HbO = tempDM;
                            designMat_HbR = tempDM;
                            designMat_HbT = tempDM;
                        end
                        
                        % extract MEAN beta, tval and corrHRF values for HbO
                        [~,stats,RValMeanOxy] =computeGLMandRval(designMat_HbO,tempMeanOxy(startFrame:stopFrame),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);
                        AllMeanResOxy(tempIndx:cond,1) = stats.beta(2:end-1);
                        AllMeanResOxy(tempIndx:cond,2) = stats.t(2:end-1);
                        AllMeanResOxy(tempIndx:cond,4) = RValMeanOxy;

                        % extract MEAN beta, tval and corrHRF values for HbR
                        [~,stats,RValMeanDeoxy] =computeGLMandRval(designMat_HbR,-1*tempMeanDeoxy(startFrame:stopFrame),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);
                        AllMeanResDeoxy(tempIndx:cond,1) = stats.beta(2:end-1);
                        AllMeanResDeoxy(tempIndx:cond,2) =  stats.t(2:end-1);
                        AllMeanResDeoxy(tempIndx:cond,4) = RValMeanDeoxy;

                        % extract MEAN beta, tval and corrHRF values for HbT
                        [~,stats,RValMeanTotalHb] =computeGLMandRval(designMat_HbT,tempMeanTotalHb(startFrame:stopFrame),SamplingRate,cond,startFrame,stopFrame,stopFrame-startFrame+1,isMultiTrial, GLMApproach,useSDCasCP);
                        AllMeanResTotalHb(tempIndx:cond,1) = stats.beta(2:end-1);
                        AllMeanResTotalHb(tempIndx:cond,2) =  stats.t(2:end-1);
                        AllMeanResTotalHb(tempIndx:cond,4) = RValMeanTotalHb;
                    end
                    if GLMApproach ==1
                        tempIndx = tempIndx + 1;
                    end
                end

                % Prepare to compute decision based on "whole" timecourse
                for cc=1:NrOfSelCh
                    for cond=1:NrOfCond
                        AllResOxy(cond,3,cc)     = corr(AllResOxy(:,1,cc),TheoBetaMat(:,cond));
                        AllResDeoxy(cond,3,cc)   = corr(AllResDeoxy(:,1,cc),TheoBetaMat(:,cond));
                        AllResTotalHb(cond,3,cc) = corr(AllResTotalHb(:,1,cc),TheoBetaMat(:,cond));
                        if NrOfSelCh>1
                            AllMeanResOxy(cond,3)     = corr(AllMeanResOxy(:,1),TheoBetaMat(:,cond));
                            AllMeanResDeoxy(cond,3)   = corr(AllMeanResDeoxy(:,1),TheoBetaMat(:,cond));
                            AllMeanResTotalHb(cond,3) = corr(AllMeanResTotalHb(:,1),TheoBetaMat(:,cond));
                        end
                    end
                end						

                % save results into variable
                if NrOfSelCh>1
                    hbType = {AllResOxy,AllResDeoxy, AllResTotalHb,AllMeanResOxy,AllMeanResDeoxy, AllMeanResTotalHb};                           
                else
                    hbType = {AllResOxy,AllResDeoxy, AllResTotalHb};
                end

                % TSI BCI output
                for cc=1:NrOfSelCh
                    
                    % define saving variables
                    maxBeta = zeros(numel(hbType),1);
                    maxtVal = zeros(numel(hbType),1);
                    maxTheoBeta = zeros(numel(hbType),1);
                    maxCorr = zeros(numel(hbType),1);  

                    for hh=1:numel(hbType) % extract condition with max value
                        if hh>3
                            [~,maxBeta(hh)] = max(hbType{hh}(:,1));
                            [~,maxtVal(hh)] = max(hbType{hh}(:,2));
                            [~,maxTheoBeta(hh)] = max(hbType{hh}(:,3));
                            [~,maxCorr(hh)] = max(hbType{hh}(:,4));  
                        else
                            [~,maxBeta(hh)] = max(hbType{hh}(:,1,cc));
                            [~,maxtVal(hh)] = max(hbType{hh}(:,2,cc));
                            [~,maxTheoBeta(hh)] = max(hbType{hh}(:,3,cc));
                            [~,maxCorr(hh)] = max(hbType{hh}(:,4,cc));   
                        end

                        maxIndx = [maxBeta(hh) maxtVal(hh) maxTheoBeta(hh) maxCorr(hh)]; % save condition with max value

                        % compute BCI plugin output based on weighting method (0 = weighted sum; 1 = betaVAL; 2 = tVAL; 3 = HRFcorrVal; 4 = rVAL;)
                        tempSelCond = computeDecision(WeightingMethod, maxIndx, NrOfCond,hbType,hh,cc,isMultiTrial);
                        selCondMultiTrial(cc,hh) = tempSelCond;
                    end
                end
                
                % displayResults
                Results=displayAndSaveResults(hbType,OutputHb,selCondMultiTrial,currRep,outputName,...
                                                   NrOfCond,NrOfSelCh,weightNameMT, isMultiTrial,showResults,...
                                                   MultiTrialApproach,selCond,Results,trial,channel);
                
                % extract decision values
                for chan=1:NrOfSelCh
                    tempDecision(comb,chan) = Results.Multitrial.(channel{chan}).Decision;
                end
            end
            
            
            for chan=1:NrOfSelCh
                allDecisionMat{qq,lengthRunTrials+1,chan} = tempDecision(:,chan);
            end

            allSubj_DecisionsB(ss,:,:,:) = allDecisionMat(:,10,BestChannelsID_b(ss));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pause(0.5)
            close all;
        end    
    end
    

    % compute accuracies
    for nn=1:NrOfRep                           
        allSubj_Acc(ss,:, nn) = 100.*(sum((squeeze(cell2mat(allDecisionMat(:,nn,:)))-repelem(allSubj_Answers(:,ss),nn)==0))./(nn*NrOfRuns(ss))); % NrOfSubj x NrOfChan x (NRep:-1:1).
    end
%     end
   
end

squeeze(allSubj_Acc(1,7, :))
