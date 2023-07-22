function [DesignMatHB,stats,RValHB] =computeGLMandRval(dm,HbData,SamplingRate,condCounter,startFrame,stopFrame,lengthArr,isMultiTrial, GLMApproach,useSDCasCP)
% this function convolves the design matrix with the hrf model and creates
% the projection matrix for Ordinary Least Squares

% input

% dm: design matrix, of dimensions ( NrOfCond+1 x totalNrOfFrames )
% Sampling Rate: extracted from hdr file 
% CondCount_Task: condition counter
% startFrame - stopFrame: trial window selection frame range
% lengthArr: minimum trial duration across all trials in order to avoid
% matrix dimension mistmatch

% output


% DeisignMatHB: condition wise design matrix (no constant term)
% stats: glm output --> see help from Matlab glmfit function
% RValHB: correlation between HbData and convolved design matrix


% Author: A.B., Maastricht University
% Last edited: 22.01.2019

%%

[hrfunc,~] = spm_hrf(1/SamplingRate); % define hrf using spm12

if GLMApproach ==2
    
    % only select frames of interest depending on reduced trials. Set
    % remaining frames to zero
    tempDM = dm(startFrame:startFrame+lengthArr-1,:);

    % convolve desing matrix with hrf
    tmpDes = [];
    for int=1:condCounter
       tmpDes = [tmpDes  conv(tempDM(:,int),hrfunc)];
    end
    tmpDesHB = tmpDes(1:size(tempDM,1),:); 

    if useSDCasCP
        tmpDesHB = [tmpDesHB tempDM(:,end)]; % add SDC as confound predictor
    end
    if isMultiTrial
       DesignMatHB = tmpDesHB; 
    else
       DesignMatHB = tmpDesHB(1:lengthArr,:); 
    end
    
    % compute HRFcorr-val: computes correlation between timecourse and
    % design matrix

    if useSDCasCP
        RValHB = zeros(size(tmpDesHB,2)-1, size(HbData,2));
        for int=1:size(HbData,2)
            RValHB(:,int) = (corr(HbData(:,int),(DesignMatHB(:,1:end-1))))';
        end
    else
        RValHB = zeros(size(tmpDesHB,2), size(HbData,2));
        for int=1:size(HbData,2)
            RValHB(:,int) = (corr(HbData(:,int),(DesignMatHB)))';
        end
    end

   

else

    tmpDes = conv(dm(:,condCounter),hrfunc); % first column is constant term
    tmpDesHB = tmpDes(startFrame:stopFrame);
    if useSDCasCP
        tmpDesHB = [tmpDesHB dm(startFrame:stopFrame,end)] ;
    end
    if isMultiTrial
       DesignMatHB = tmpDesHB; 
    else
       DesignMatHB = tmpDesHB(1:lengthArr,:); 
    end
    if useSDCasCP
        RValHB = (corr(HbData,(DesignMatHB(:,1:(end-1)))))';
    else
        RValHB = (corr(HbData,(DesignMatHB)))'; 
    end
    
end

[~,~,stats]=glmfit(DesignMatHB,HbData); % Run GLM
% figure, plot(HbData), hold on, plot(DesignMatHB)

%  if useSDCasCP
%      figure, plot(DesignMatHB*stats.beta(2:end),'--b'), hold on, plot(DesignMatHB(:,1:end-1)), hold on, plot(HbData,'k')
% %      pause(),
%  else
%      figure, plot(DesignMatHB*stats.beta(2:end)), hold on, plot(DesignMatHB(:,1:end)), hold on, plot(HbData,'k')
%  end
