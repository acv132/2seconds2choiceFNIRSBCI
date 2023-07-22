function Results=displayAndSaveResults(hbType,OutputHb,selCond,currRep,outputName,NrOfCond,NrOfSelCh,weightName, isMultiTrial, showResults,MultiTrialApproach, allDec,Results,trial,channel )
% input

% hbType: 3D or 4D matrix containing BCIOutput results [bVal bVal tVal rVal] for all channels and conditions
% OutputHb: selected Hb for output (=1:HbO; =2: HbR; =3: HbT)
% selCond: matrix containing selected BCIOutput results (decision) based on weighting method
% currRep: trial/repetition number (max = NrOfRep)
% outputName: 
% NrOfCond: total # of Conditions
% NrOfSelCh: total # of selected channels
% weightName: selected weigthting approach
% isMultiTrial: flag defining if current trial belongs to the multitrial computation
% showResults: only active if MultitrialApproach >0. 0 = don't show results; 1 = show decision; 2 = show all BCI output window results
% MultiTrialApproach: selected multi trial approach
% allDec: similar to selCond, matrix containing ALL BCIOutput results (decision) based on weighting method. Only relevant when MultitrialApproach > 0
% Results: nested structure that contains results per trial and channel
% trial: cell array for Results structure that contains trial names
% channel: cell array for Results structure that contains channel names

% output

% Results: nested structure that contains results per trial and channel

% Author: A.B., Maastricht University
% Last edited: 22.01.2019

%%

fid = 1;

if isMultiTrial == 0
    disp(' ')
    disp(' ')
    disp(['trial #' num2str(currRep) ':'])
    disp('--------')

    disp(['Results based on method: ' weightName ])
    disp('Display Format--> rows: [betaVal bVal tVal rVal]; columns: conditions')
    for int=1:size(hbType{OutputHb},4)
        
        disp(' ')
        disp(['Channel ' num2str(int)])
        tempMat = [squeeze(hbType{OutputHb}(:,1,currRep,int));
                   squeeze(hbType{OutputHb}(:,3,currRep,int));
                   squeeze(hbType{OutputHb}(:,2,currRep,int));
                   squeeze(hbType{OutputHb}(:,4,currRep,int))];
        fprintf(fid, [repmat(' %7.3f ', 1, NrOfCond) '\n'], tempMat)
        disp(' ')
        Results.(trial{currRep}).(channel{int}).BCIOutput = reshape(tempMat,[NrOfCond 4])';
        Results.(trial{currRep}).(channel{int}).Decision = selCond(currRep,int,OutputHb);
    end
    
    
    disp(['DECISION ' outputName  ' (per channel):  ' num2str(selCond(currRep,:,OutputHb))])
    disp('------------------------------')
    disp(' ')
    if NrOfSelCh>1
        disp('MEAN VAL')
        disp('-----')
        disp('Display Format--> rows: [betaVal bVal tVal rVal]; columns: conditions')
        disp(' ')
        tempMat = [squeeze(hbType{OutputHb+ 3}(:,1,currRep));
                   squeeze(hbType{OutputHb+ 3}(:,3,currRep));
                   squeeze(hbType{OutputHb+ 3}(:,2,currRep));
                   squeeze(hbType{OutputHb+ 3}(:,4,currRep))];
        fprintf(fid, [repmat(' %7.3f ', 1, NrOfCond) '\n'], tempMat)
        disp(' ')
        disp(['decision MEAN '  outputName  ':  ' num2str(selCond(currRep,1,OutputHb+3))])
        disp(' ')
        
        Results.(trial{currRep}).(channel{end}).BCIOutput = reshape(tempMat,[NrOfCond 4])';
        Results.(trial{currRep}).(channel{end}).Decision = selCond(currRep,1,OutputHb+3);
    end    
       
else
        
    if MultiTrialApproach == 1
        weightNameMT = 'Weighted single trial';
    elseif MultiTrialApproach == 2
        weightNameMT = 'Whole timecourse';
    end
    if showResults > 0
        disp(' ')
        disp('Multi Trial Approach')
        disp('--------------------')
        disp(' ')
        disp(['Results based on method: ' weightNameMT ])

    end
    if MultiTrialApproach == 1
        if showResults == 2
            disp(['Selected conditions over trials. Selection method: ' weightName])
            disp('Display Format--> rows: selected cond. in trial; columns: channels')
            disp(' ')           
            fprintf(fid, [repmat(' %7.3f ', 1, NrOfSelCh) '\n'], (allDec(:,:,OutputHb)'))
        end
        for int=1:size(allDec,2)
            Results.(trial{end}).(channel{int}).BCIOutput = (allDec(:,int,OutputHb))';
            Results.(trial{end}).(channel{int}).Decision = selCond(int,OutputHb);
        end

    elseif MultiTrialApproach == 2 
         for int=1:size(hbType{OutputHb},3) 
            tempMat = [squeeze(hbType{OutputHb}(:,1,int));
                       squeeze(hbType{OutputHb}(:,3,int));
                       squeeze(hbType{OutputHb}(:,2,int));
                       squeeze(hbType{OutputHb}(:,4,int))];
            if showResults == 2 
                if int==1 
                    disp('Display Format--> rows: [betaVal bVal tVal rVal]; columns: conditions')
                end
                disp(' ')
                disp(['Channel ' num2str(int)])

                fprintf(fid, [repmat(' %7.3f ', 1, NrOfCond) '\n'], tempMat)
                disp(' ')
            end
            Results.(trial{end}).(channel{int}).BCIOutput = reshape(tempMat,[NrOfCond 4])';
            Results.(trial{end}).(channel{int}).Decision = selCond(int,OutputHb);
        end
    end


    if showResults >0 
        disp(' ')
        disp(['DECISION ' outputName  ' (per channel):  ' num2str(selCond(:,OutputHb)')])
        disp('---------------------')
        disp(' ')
    end
    if NrOfSelCh>1
        if showResults >0
            disp('MEAN VAL')
            disp('-----')
        end
        if MultiTrialApproach == 1

            if showResults ==2
                disp(['Selected conditions over trials. Selection method: ' weightName])
                disp('Display Format--> rows: sel condition per trial')
                disp(' ')
                fprintf(fid, [repmat(' %7.3f ', 1, 1) '\n'], (allDec(:,1,OutputHb+3))')
            end
            Results.(trial{end}).(channel{end}).BCIOutput = (allDec(:,1,OutputHb+3))';
            Results.(trial{end}).(channel{end}).Decision = selCond(1,OutputHb+3);
        elseif MultiTrialApproach == 2 
            tempMat = [squeeze(hbType{OutputHb+ 3}(:,1));
                       squeeze(hbType{OutputHb+ 3}(:,3));
                       squeeze(hbType{OutputHb+ 3}(:,2));
                       squeeze(hbType{OutputHb+ 3}(:,4))];
            if showResults ==2
                disp('Display Format--> rows: [betaVal bVal tVal rVal]; columns: conditions')
                disp(' ')

                fprintf(fid, [repmat(' %7.3f ', 1, NrOfCond) '\n'], tempMat)
                disp(' ')
            end
            Results.(trial{end}).(channel{end}).BCIOutput = reshape(tempMat,[NrOfCond 4])';
            Results.(trial{end}).(channel{end}).Decision = selCond(1,OutputHb+3);
        end

        if showResults >0
            disp(' ')
            disp(['decision MEAN '  outputName  ':  ' num2str(selCond(1,OutputHb+3)')])
            disp(' ')
        end
    end
    
end
