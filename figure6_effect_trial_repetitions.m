% Created by A. Benitez Andonegui, M. Luehrs
% Modified by A. Vorreuther

% subjects are divided by setup choice of localizer since channels differ
% per setup
SubjNames_a = {'P03_2019-04-08','P06_2019-01-29','P08_2019-01-22','P15_2019-02-19','P07_2019-01-15'};
SubjNames_b = {'P04_2019-02-19','P05_2019-02-19','P09_ 2019-01-15','P101_2019-02-19','P102_2019-02-19'};

% subjects as sorted during analysis (see above) but with names as used in
% paper
tableSubjNames = {'P03','P06','P08','P10','P07','P04','P05','P09','P01','P02'};

BestChannels_a = {'S2-D6','S1-D2','S8-D7','S3-D4','S2-D6'};
BestChannels_b = {'S5-D6','S3-D3','S2-D3','S8-D8','S2-D3'};

BestChannelsID_a = [7,2,22,9,7];
BestChannelsID_b = [19,11,5,30,5];

% accuracy matrices are predefined (matlab-style)
AccA = zeros(10,5);
AccB = zeros(10,5);

% decisions in COI per subject
DecisionsA = cell(10,5);
DecisionsB = cell(10,5);

% path to functions is added (needed for brainAR_main
addpath(genpath("figure6_effect_trial_repetitions_itr/"));
% edit mainPaths in each brainAR function to match the folder where the
% data is stored

% accuracies for each number of trials (1 trial used, 2 consecutive trials
% used, 3 consecutive trials used, ...) for setup A; for each subject, only
% the best channel's data is used (see also table with best Hb-values
% during localizers
brainAR_main_a
for ind = 1:numel(SubjNames_a) -1
    AccA(:,ind) = allSubj_Acc(ind,BestChannelsID_a(ind), :);
end

brainAR_main_P07
AccA(:,5) = allSubj_Acc(:,BestChannelsID_a(5), :);

% same for setup B
brainAR_main_b
for ind = 1:numel(SubjNames_b)
    AccB(:,ind) = allSubj_Acc(ind,BestChannelsID_b(ind), :);
end

trial_decisions = array2table([allSubj_DecisionsA; allSubj_DecisionsB]');
trial_decisions.Properties.VariableNames(1:10) = tableSubjNames;
oldvariables = trial_decisions.Properties.VariableNames;
newvariables = sort(tableSubjNames);
[~,LOCB] = ismember(newvariables,oldvariables);
trial_decisions = trial_decisions(:,LOCB);
writetable(trial_decisions, "trials_decisions.csv");

% table used for calculation of information-transfer rate is created; each
% row corresponds to the effect of trial number, each column is a subject;
effect_trial_repetition = array2table([AccA AccB]);
effect_trial_repetition.Properties.VariableNames(1:10) = tableSubjNames;

% sort table rows by subjects and columns by "1 trial", "2 trials", "3
% trials", etc.
oldvariables = effect_trial_repetition.Properties.VariableNames;
newvariables = sort(tableSubjNames);
[~,LOCB] = ismember(newvariables,oldvariables);
effect_trial_repetition = effect_trial_repetition(:,LOCB);
tmp_array = table2array(effect_trial_repetition);
itr_calculation_table = array2table(tmp_array.');
itr_calculation_table.Properties.RowNames = effect_trial_repetition.Properties.VariableNames;
itr_calculation_table = flip(itr_calculation_table,2);
itr_calculation_table = sortrows(itr_calculation_table,'RowNames');
itr_calculation_table.Properties.VariableNames(1:10) = ...
["1 trial", "2 trials", "3 trials", "4 trials", "5 trials", ...
"6 trials", "7 trials", "8 trials", "9 trials", "10 trials"];

% add average
M = varfun(@mean, itr_calculation_table, 'InputVariables', @isnumeric);
M.Properties.RowNames = "Average";
M.Properties.VariableNames(1:10) = itr_calculation_table.Properties.VariableNames(1:10);
itr_calculation_table = [itr_calculation_table; M];

% save itr_calculation_table for ITR calculation
writetable(itr_calculation_table, "..\data + material\data\derivative\figure6_effect_trial_repetitions_itr_data.xlsx",'WriteRowNames',true);