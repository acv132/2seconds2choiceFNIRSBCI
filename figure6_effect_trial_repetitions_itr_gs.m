% Created by A. Benitez Andonegui, M. Luehrs
% Modified by A. Vorreuther

%%% define variables for loading data, plotting  & computing ITR values %%%
close all 

fileName = "figure6_effect_trial_repetitions_itr_data.xlsx"; % accuracy matrix: nsubj x ntrials
% pathToViolinPlotFunc = "C:\Users\Anna\Desktop\fNIRS\brainAR-main";

addpath(genpath("..\data + material\data\derivative\"));
addpath(genpath(fullfile(fileparts(pwd), ...
    "analysis & results\tempEncPaper_matlabcode\")));

chanceLev = 50; % chance-level value, in %

% variables for ITR computation
n_trials = 10;
n_choices = 2; %YES/NO
n_subj = 10;
Tau = 65:25:300;
% Tau = 45:25:270; % run duration (in seconds) the different # of trials
                  % TAU should have the same convention as the columns in tempMat 
                  % if columns are in increasing order of # trials, TAU durations should also be increasing
                  % TAU includes initial and final baseline

                  
%%% ------- end of definition ------- %%%

% necessary for ITR computation
P = 0: 0.01:1; % set array for ITR model plotting
P = P(2:end-1); % to avoid -infinity of the logarithm (and the NaNs in ITR computation)

% load accuracy matrix
[tempMat,~,~]=xlsread(fileName); % tempMat is (n_subj+1) x n_trials
av = tempMat(end,:)./100; % average accuracies
tempMat = tempMat(1:n_subj,:); % remove average (if more rows than n_subj)

% remove subjects P08 and P01
% tempMat = tempMat(2:end,:); % P01 removed
% tempMat = [tempMat(1:6,:); tempMat(8:end,:)]; % P08 removed
% n_subj = size(tempMat, 1);

% make it between 0 and 1
allSubjAcc = tempMat./100;
allSubjAcc(allSubjAcc==1) = max(P); % to avoid infinites/NaNs
allSubjAcc(allSubjAcc==0) = min(P); % to avoid infinites/NaNs


%%% plot & compute ITR

% addpath(pathToViolinPlotFunc) % add path to violinplot function

figH = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'color','w')

% subplot 1: compute boxplots of accuracies/n_trials, based on violinplot function
subplot(1,2,1) 
vvH = violinplot(tempMat);

% set axis parameters
ax = gca;
ax.FontSize = 24; 
xTickLab = cell(1,n_trials);
for nn=1:n_trials
    %xTickLab{nn} = [num2str(nn) ' trials']; 
    xTickLab{nn} = [num2str(nn)]; 
end
%xTickLab{1}(end) = [];
ax.XLabel.String = 'Number of trial repetitions';
ax.XTickLabel = xTickLab;
%ax.XTickLabelRotation = 45;
ax.YLim = [0 105];
ax.YLabel.String = 'Classification accuracy (%)';
%ax.Title.String = '(A) Accuracy-based assessment';
ax.Title.String = 'A';
ax.Title.HorizontalAlignment = "left";
ax.Title.Position(1)= ax.XLim(1);
ax.FontName = "Arial";

xh = get(ax,'xlabel'); % handle to the label object
p = get(xh,'position'); % get the current position property
set(xh, "position", [5.1500   -8.9670   -1.0000]);

hold on,
ax = axis;
lineH2 = line([ax(1) ax(2)],repelem(chanceLev,2));
lineH2.Color = "#122d46";%[ 1.0000    0.7812    0.4975];
lineH2.LineWidth = 1.25;
lineH2.LineStyle = '--';

%% choose colors for plots
load("cgr_maps\custom_cgr5.mat");
load("cgr_maps\custom_cgr_gs.mat");
load("cgr_maps\customredgreen.mat");

cm = colormap(turbo(n_trials*5));
cm = gray(n_trials*5); % grayscale colormap
cgr = cgr5; % cgr_gs for grayscale figure
%%
cgr = ["diamond","*","o","^","square"].';


avgAcc= mean(allSubjAcc, 2, "omitnan");
avgAcc_sorted = sort(mean(allSubjAcc, 2, "omitnan"));
bins=discretize(avgAcc,6);
cl = [];%zeros(n_subj,1);

for subj=1:n_subj
    if bins(subj) == 1 
        cl{subj} = cgr{1};
    elseif bins(subj) == 2
        cl{subj}  = cgr{2};
    elseif bins(subj) == 3 || bins(subj) == 4
        cl{subj}  = cgr{3};
    elseif bins(subj) == 5
        cl{subj}  = cgr{4};
    elseif bins(subj) == 6
        cl{subj}  = cgr{5};
    end
end

avgAcc= mean(allSubjAcc, 2, "omitnan");
avgAcc_sorted = sort(mean(allSubjAcc, 2, "omitnan"));
bins=discretize(avgAcc,6);
cl2 = zeros(n_subj,3);

% for subj=1:n_subj
%     if bins(subj) == 1 
%         cl2(subj,:) = cgr5(1,:);
%     elseif bins(subj) == 2
%         cl2(subj,:)  = cgr5(2,:);
%     elseif bins(subj) == 3 || bins(subj) == 4
%         cl2(subj,:)  = cgr5(3,:);
%     elseif bins(subj) == 5
%         cl2(subj,:)  = cgr5(4,:);
%     elseif bins(subj) == 6
%         cl2(subj,:)  = cgr5(5,:);
%     end
% end
% 

alpha = 0.6;
cmbox = zeros(size(cm)); % all black boxplot bars
indx = 1:5:size(cm,1); % select indices within this colormap
for vv=1:numel(vvH)
    vvH(vv).ViolinPlot.Visible = 'off';
    vvH(vv).BoxColor = cmbox(indx(vv),:);
    vvH(vv).BoxWidth = 0.05;
    vvH(vv).ScatterPlot.MarkerFaceColor = "flat";
%     if vv == 10
%         vvH(vv).ScatterPlot.CData = [cl(1:6,:); cl(8:end,:)];
%     else
        vvH(vv).ScatterPlot.CData = cl2;%cl;
%     end

    markerShape = cl;
    points = {vvH(vv).ScatterPlot.XData, vvH(vv).ScatterPlot.YData};
    for subj = 1:n_subj
        scatter(points{1}(subj),points{2}(subj),...
        'Marker', markerShape{subj}, ...
        'MarkerFaceAlpha', alpha, ...
        'MarkerFaceColor', 'k',...
        'MarkerEdgeColor', 'k', ...
        'MarkerEdgeAlpha', alpha);
    end
    %vvH(vv).ScatterPlot.CData = cl{vv}; % Assigning marker colors

    vvH(vv).ScatterPlot.SizeData = 50;
    vvH(vv).ScatterPlot.MarkerFaceAlpha = 0;
    vvH(vv).MeanPlot.XData = (vv+[-1,1].*0.17);
%     vvH(vv).MeanPlot.XData = [vvH(vv).MeanPlot.XData(1) + 0.15 vvH(vv).MeanPlot.XData(2) - 0.15];
    vvH(vv).MeanPlot.Color = cmbox(indx(vv),:);
    vvH(vv).ShowMean = 1;   
    uistack(vvH(vv).ScatterPlot, 'top');
    delete(vvH(vv).WhiskerPlot)
    delete(vvH(vv).MedianPlot)
end
axis square

hold on
h = zeros(5, 1);
legendLab = ["90-100%", "80-90%", "70-80%","50-60% (P01)", "<50% (P08)"];
legendCol = flipud(cgr);
markerCol = flipud(cgr);

for ha=1:5
    h(ha) = scatter(NaN,NaN,...
                "Marker", markerCol(ha),... 
                "SizeData", 100,...
                "MarkerFaceAlpha", alpha,...
                "MarkerEdgeColor","k",...
                'MarkerEdgeAlpha', alpha, ...
                "MarkerFaceColor", "k");
end
legH = legend(h,legendLab);
legH.EdgeColor = [1 1 1];
legH.Title.String = 'Average participant-wise accuracy';
legH.Location = 'southwest';
legH.Position = legH.Position + [0.001 0.001 0 0];
legH.FontSize = 18;
legH.FontName = "Arial";

% subplot 2: ITR computation

subplot(1,2,2) 
ax2 = gca;

itrMat = zeros(n_subj+1, n_trials); % initialize matrix storing itr values. last value = average
for tt=1:n_trials
    hold on,

    % first compute & plot theoretical ITR distribution
    ITR_Model = computeITR(n_choices,P,Tau(tt));
    pH = plot(P,ITR_Model);
    pH.Color = cmbox(indx(tt),:);
    
    % now compute % plot subject specific ITR values
    hold on, 
    for ss=1:n_subj
        hold all,
        P_prime = allSubjAcc(ss,tt);
        itrMat(ss,tt) = computeITR(n_choices,P_prime,Tau(tt));
        scH2 = scatter(P_prime, itrMat(ss,tt));
        scH2.MarkerFaceColor = "#a3a3a3";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
        scH2.MarkerEdgeColor = "#a3a3a3";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
        scH2.MarkerFaceAlpha = 0.25;
        scH2.SizeData = 70;
    end
    % lastly, compute & plot the average of the itrs
    hold on, 
    itrMat(end,tt) = computeITR(n_choices,av(tt),Tau(tt));
    scH(tt) = scatter(av(tt), itrMat(end,tt));
    
    % edit axes
    scH(tt).MarkerFaceColor = cm(indx(tt),:);
    scH(tt).MarkerEdgeColor = [0 0 0];
    scH(tt).SizeData = 180;
    scH(tt).Marker = "o";
end
ax2.XLim = ([min(allSubjAcc(:))-0.01 1]); 
ax2.YLim = ([0 max(itrMat(:))+0.01]); 
ax2.YLabel.String = 'ITR (bits/min)';
ax2.XLabel.String = 'Classification accuracy';
ax2.FontSize = 24;
axis square

xh2 = get(ax2,'xlabel'); % handle to the label object
p2 = get(xh2,'position'); % get the current position property

xTickLab = cell(1,n_trials);
for nn=1:n_trials
    xTickLab{nn} = [num2str(nn) ' trials']; 
end
xTickLab{1}(end) = [];
legendLab = xTickLab;
for nn=1:n_trials
    legendLab{nn} = [legendLab{nn}];
end
legH = legend(scH,legendLab);
legH.EdgeColor = [1 1 1];
legH.Title.String = 'Average ITR values';
legH.Location = 'northwest';
legH.Position = legH.Position + [0.03 0 0 0];
legH.FontSize = 18;
legH.FontName = "Arial";
%ax2.Title.String = '(B) ITR-based assessment';
ax2.Title.String = 'B';
ax2.Title.HorizontalAlignment = "left";
ax2.Title.Position(1)= ax2.XLim(1);
FontName = "Arial";

set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf,'figure5_effect_trial_repetitions_itr_gs.svg','svg')
%('-dtiff','-r300','figure5_effect_trial_repetitions_itr.svg') % save fig
%with specific resolution

disp('correlation mean accuracies & number of trial repetitions')
[RHO,PVAL] = corr([1:10]',av','Type','Spearman')

% itr computation
function itr = computeITR(N,P,tau)
    itr = (log2(N)+ P.*log2(P) + (1-P).*log2((1-P)./(N-1))).*(60./tau);
end
