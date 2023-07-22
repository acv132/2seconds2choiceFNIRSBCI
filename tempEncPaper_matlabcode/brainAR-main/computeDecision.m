function selCond=computeDecision(WeightingMethod, maxIndx, NrOfCond,hbType,hbCounter,channelCounter,isMultiTrial)

% this function outputs decision values based on weighting method 

% input

% WeightingMethod: approach to compute decision --> 0 = weighted sum; 1 = max BetaVal; 2 = max t-Val; 3 = max bCorr-Val; 4 = max r-Val
% maxIndx: 1x4 array with [betaVal tVal bCorrVal rVal]
% NrOfCond: # of conditions according to protocol file
% hbType: matrix that stores all [betaVal tVal bCorrVal rVal] values for all conditions, hb Type, channels and trial repetitions
% hbCounter: counter that goes through all hb types (oxy, deoxy, totalHb in case NrOfSelCh = 1; oxy,deoxy,totalHb, meanOxy, meanDeoxy, meanTotalHb if NrOfSelCh >1)
% channelCounter: counter that goes through all selected channels
% isMultiTrial: flag to indicate whether decision computation should be
% based on a multitrial approach or not

% output

% selCond: condition that was selected based on Weighting method

% Author: A.B., Maastricht University
% Last edited: 22.01.2019

%%
switch WeightingMethod

    case 0
    
        [val , ~] = histc(maxIndx,1:NrOfCond); % count frequency of items in maxIndx
        tempValue  = max(val); % select max frequency
        tempInd = find(val == tempValue);
        
        if (numel(tempInd )==1)                                       
           selCond  = tempInd; 
        else % look at the sum of channel-wise {betaVal,tVal,b-corrVal,r-Val}
            if isMultiTrial
               if hbCounter>3 % look at the sum of channel-wise val {betaVal,tVal,b-corrVal,r-Val}
                   temp = sum(hbType{hbCounter}(tempInd,:),2); 
               else
                   temp= sum(hbType{hbCounter}(tempInd,:,channelCounter),2); 
               end
               [~, posValMax ] = max(temp); % select the term with the highest value
               selCond  = tempInd(posValMax); 
            else % single trial approach
                if hbCounter>3
                   temp = sum(hbType{hbCounter}(tempInd,:,:),3); 
                else
                   temp = sum(hbType{hbCounter}(tempInd,:,:,channelCounter),3); 
                end

                % apply the same approach as before                                        
                [~,maxBetaT] = max(temp(:,1));
                [~,maxtValT] = max(temp(:,2));
                [~,maxTheoBetaT] = max(temp(:,3));
                [~,maxCorrT] = max(temp(:,4)); 

                maxIndxT = [maxBetaT maxtValT maxTheoBetaT maxCorrT];

                [val , ~] = histc(maxIndxT,1:NrOfCond);
                if numel(unique(val ))>1
                    tempValue  = max(val);
                    tempInd = find(val==tempValue);
                    if (numel(tempInd )==1)
                        selCond = tempInd; 
                    else % select the one with highest t-Val --> arbitrarily decided
                        selCond = maxtValT;
                    end
                end   
            end

        end

    case 1
        selCond  = maxIndx(1);
    case 2
        selCond  = maxIndx(2);                                                                            
    case 3                                                                    
        selCond  = maxIndx(3);
    case 4
        selCond  = maxIndx(4);
end



