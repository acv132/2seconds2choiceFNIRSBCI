function indTrialsMat = getConsecutiveTrials(arr,nrofVal)
% this function creates arrays of consecutive numbers given the desired length of array
% input:
% arr: array of numbers to choose from, e.g., 1:10
% nrofVal: desired length of array, e.g., 3
% output:
% indTrialsMat: matrix of indices. dimensions: nrofcombinations x nrofVal
% example: getConsecutiveTrials(1:10,3) will give the following output
%     1     2     3
%     2     3     4
%     3     4     5
%     4     5     6
%     5     6     7
%     6     7     8
%     7     8     9
%     8     9    10
%%%%%%%%%%%%%

tempMat = combnk(arr,nrofVal); % get all possible combinations
if nrofVal >1
    indTrialsMat = tempMat(ismember(diff(tempMat,1,2),ones([1 nrofVal-1]),'rows'),:); % select only the rows whose columns are all ones (= consecutive numbers)
else
    indTrialsMat = arr';
end
