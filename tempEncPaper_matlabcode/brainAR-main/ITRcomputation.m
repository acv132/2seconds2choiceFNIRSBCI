
%%% define variables for loading data, plotting  & computing ITR values %%%

fileName = "C:\Users\Anna\surfdrive\Marble\publication\testAccuracies.xlsx"; % accuracy matrix: nsubj x ntrials
pathToViolinPlotFunc = "C:\Users\Anna\Downloads\brainAR-main";
chanceLev = 50; % chance-level value, in %

% variables for ITR computation
n_trials = 10;
n_choices = 2; %YES/NO
n_subj = 10;
Tau = 120:20:300; % run duration (in seconds) the different # of trials
                  % TAU should have the same convention as the columns in tempMat 
                  % if columns are in increasing order of # trials, TAU durations should also be increasing
                  % TAU includes initial and final baseline

                  
%%% ------- end of definition ------- %%%

% necessary for ITR computation
P = 0: 0.01:1; % set array for ITR model plotting
P = P(2:end-1); % to avoid -infinity of the logarithm (and the NaNs in ITR computation)

% load accuracy matrix
[tempMat,~,~]=xlsread(fileName); % tempMat for Anna is (n_subj+1) x n_trials
av = tempMat(end,:)./100; % average accuracies
tempMat = tempMat(1:n_subj,:); % remove average

% make it between 0 and 1
allSubjAcc = tempMat./100;
allSubjAcc(allSubjAcc==1) = max(P); % to avoid infinites/NaNs
allSubjAcc(allSubjAcc==0) = min(P); % to avoid infinites/NaNs


%%% plot & compute ITR

addpath(pathToViolinPlotFunc) % add path to violinplot function

figH = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'color','w')

% subplot 1: compute boxplots of accuracies/n_trials, based on violinplot function
subplot(1,2,1) 
vvH = violinplot(tempMat);

% set axis parameters
ax = gca;
ax.FontSize = 13; 
xTickLab = cell(1,n_trials);
for nn=1:n_trials
    xTickLab{nn} = [num2str(nn) ' trials']; 
end
xTickLab{1}(end) = [];
ax.XTickLabel = xTickLab;
ax.XTickLabelRotation = 45;
ax.YLim = [0 105];
ax.YLabel.String = 'Classification accuracy (%)';
ax.Title.String = '(A) Accuracy-based assessment';
hold on,
ax = axis;
lineH2 = line([ax(1) ax(2)],repelem(chanceLev,2));
lineH2.Color = "#122d46";%[ 1.0000    0.7812    0.4975];
lineH2.LineWidth = 1.25;
lineH2.LineStyle = '--';

% edit violinplots so that they look like boxplots
cm = colormap(gray(n_trials*5)); % create colormap
indx = 1:5:size(cm,1); % select indices within this colormap
for vv=1:numel(vvH)
    vvH(vv).ViolinPlot.Visible = 'off';
    vvH(vv).BoxColor = cm(indx(vv),:);
    vvH(vv).BoxWidth = 0.05;
    vvH(vv).ScatterPlot.MarkerFaceColor = "#2e75b5";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
    vvH(vv).ScatterPlot.CData = "#2e75b5";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
    vvH(vv).ScatterPlot.SizeData = 50;
    vvH(vv).MeanPlot.XData = [vvH(vv).MeanPlot.XData(1) + 0.15 vvH(vv).MeanPlot.XData(2) - 0.15];
    vvH(vv).MeanPlot.Color = cm(indx(vv),:);
    vvH(vv).ShowMean = 1;   
end
axis square

% subplot 2: ITR computation

subplot(1,2,2) 
ax2 = gca;

itrMat = zeros(n_subj+1, n_trials); % initialize matrix storing itr values. last value = average
for tt=1:n_trials
    hold on,

    % first compute & plot theoretical ITR distribution
    ITR_Model = computeITR(n_choices,P,Tau(tt));
    pH = plot(P,ITR_Model);
    pH.Color = cm(indx(tt),:);
    
    % now compute % plot subject specific ITR values
    hold on, 
    for ss=1:n_subj
        hold all,
        P_prime = allSubjAcc(ss,tt);
        itrMat(ss,tt) = computeITR(n_choices,P_prime,Tau(tt));
        scH2 = scatter(P_prime, itrMat(ss,tt));
        scH2.MarkerFaceColor = "#2e75b5";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
        scH2.MarkerEdgeColor = "#2e75b5";%[0.357 0.608 0.835];%[0.6350 0.0780 0.1840];
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
    scH(tt).Marker = 'o';
end
ax2.XLim = ([min(allSubjAcc(:))-0.01 1]); 
ax2.YLabel.String = 'ITR (bits/min)';
ax2.XLabel.String = 'Classification accuracy';
ax2.FontSize = 13;
axis square

legendLab = xTickLab;
for nn=1:n_trials
    legendLab{nn} = [legendLab{nn} '/run'];
end
legH = legend(scH,legendLab);
legH.EdgeColor = [1 1 1];
legH.Title.String = 'Average ITR values';
legH.Location = 'northwest';
legH.Position = legH.Position + [0.03 0 0 0];
legH.FontSize = 12;


ax2.Title.String = '(B) ITR-based assessment';
%print('-dtiff','-r300',fullfile('path/to/folder','Figure_XX')) % save fig
%with specific resolution

% itr computation

function itr = computeITR(N,P,tau)
    itr = (log2(N)+ P.*log2(P) + (1-P).*log2((1-P)./(N-1))).*(60./tau);
end
