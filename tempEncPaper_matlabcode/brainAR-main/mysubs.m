SubjNames_a = {'P03_2019-04-08','P06_2019-01-29','P07_2019-01-15','P08_2019-01-22','P15_2019-02-19'};
SubjNames_b = {'P04_2019-02-19','P05_2019-02-19','P09_ 2019-01-15','P101_2019-02-19','P102_2019-02-19'};

BestChannels_a = {'S2-D6','S1-D2','S2-D6','S8-D7','S3-D4'};
BestChannels_b = {'S5-D6','S3-D3','S2-D3','S8-D8','S2-D3'};

BestChannelsID_a = [7,2,7,22,9];
BestChannelsID_b = [19,11,5,30,5];

AccA = zeros(10,5);

AccB = zeros(10,5);


% for ind = 1:numel(SubjNames_a)
%     AccA(:,ind) = allSubj_Acc(ind,BestChannelsID_a(ind), :);
% end

for ind = 1:numel(SubjNames_b)
    AccB(:,ind) = allSubj_Acc(ind,BestChannelsID_b(ind), :);
end

