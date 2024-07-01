%% fit von Mises κ to alignment-stretch data for one sample

% NOTES:
% (1) I found it helpful to reference this paper: Wang... Goulbourne. A mixed Von Mises distribution for modeling soft biological tissues with two distributed fiber properties. Int J of Solids and Structures. 49(21). 2012.
% (2) κ is "b" throughout

data = DMSOexp4_data_3Dplot; % copy the name of the "3Dplot" data file here
b_all = zeros(1,size(data,2)-1);
% extract angles used for the centers of histogram bins:
angles = data(4:end,1) .* (pi/180);
% iterate through timepoints and fit von Mises dist. to orientation histogram:
for t = 2:size(data,2)
    % extract histogram weights at each angle for this timepoint:
    weights = data(4:end,t) .* (180/pi);
    % default starting values for fitting:
    mse_min = 1000;
    b_min = 0;
    % iterate through values of b till error with the data is minimized:
    for b = -5:0.01:5
        % compute von Mises distriution for this concentration parameter (b):
        vonMises = exp(b.*cos(2.*angles)) ./ (pi*besseli(0,b));
        % calculate MSE against the fiber distribution histogram:
        mse = mean((weights - vonMises).^2);
        if mse < mse_min
            mse_min = mse;
            b_min = b;
        end
    end
    % store value of b that gave min MSE:
    b_all(t-1) = b_min;
    % show fit on plot of data:
    figure
    plot(angles, exp(b_min.*cos(2.*angles)) ./ (pi*besseli(0,b_min))); hold on
    plot(angles, weights)
end

% plot resulting b values vs axial stretch:
figure
plot(data(2,2:end), b_all, 'k'); hold on
plot([1,4], [0,0], 'k:', 'LineWidth', 1.5)
xlabel('stretch')
ylabel('concentration parameter, b')


%% fit von Mises κ to alignment-stretch data for pooled samples

% all inputs (DMSO or CD):
data_3Dplot = [DMSOexp4_data_3Dplot, DMSOexp5_data_3Dplot, DMSOexp7_data_3Dplot];
lam_star = [3.16, 2.88, 2.03]; % transition stretch, λ*
%data_3Dplot = [CDexp3_data_3Dplot, CDexp4_data_3Dplot, CDexp6_data_3Dplot];
%lam_star = [1.56, 2.04, 1.93]; % transition stretch, λ*
width_lamBins = 0.23;

% remove excess rows containing calculated transverse stretch (1) and measured stress (3):
data_3Dplot([1 3],:) = [];
% extract (rad) angles used for the centers of histogram bins (also remove these from the data):
angles = data_3Dplot(2:end,1) .* (pi/180);
data_3Dplot = data_3Dplot(:,data_3Dplot(1,:)>0);
% rescale each stretch axis (now in row 1) to be "fraction of λ*":
i_sample = 0;
for t = 1:size(data_3Dplot, 2)
    % ones mark seperations between samples:
    if data_3Dplot(1,t) == 1
        i_sample = i_sample + 1;
    end
    data_3Dplot(1,t) = (data_3Dplot(1,t) - 1) / (lam_star(i_sample) - 1);
end
% sort timepoints by rescaled stretch value:
data_3Dplot = sortrows(data_3Dplot')';
% bin along stretch axis according to width_lamBins:
data_lamBinned = zeros(100, size(data_3Dplot, 1), 2); % assume excess of 100 bins (:,,), each storing mean (,,1) and std (,,2) for rescaled stretch (,1,) and all angle weights (,2:end,)
data_currentBin = [];
lamVal_current = 0;
i_currentBin = 0;
for t = 1:size(data_3Dplot, 2)
    % if neccessary, iterate to next bin:
    if data_3Dplot(1,t) > lamVal_current + width_lamBins
        lamVal_current = lamVal_current + width_lamBins;
        i_currentBin = i_currentBin + 1;
        data_lamBinned(i_currentBin,:,1) = mean(data_currentBin, 2); % calculate bin means
        data_lamBinned(i_currentBin,:,2) = std(data_currentBin, 0, 2); % calculate bin stds
        data_currentBin = []; % reset
    end
    % include data in current bin:
    data_currentBin = [data_currentBin, data_3Dplot(:,t)];
end
data_lamBinned = data_lamBinned(data_lamBinned(:,1,1) > 0,:,:); % remove unused trailing bins
% to check: scatter(angles , data_lamBinned(1,2:end,1));

% initialize storage for von Mises concentration parameters (1 per bin):
b_all = zeros(1, size(data_lamBinned,1));
% iterate through bins and fit von Mises dist. to mean orientation histogram:
for t = 1:size(data_lamBinned,1)
    % extract histogram weights at each angle for this bin:
    weights = data_lamBinned(t,2:end,1)' .* (180/pi);
    % default starting values for fitting:
    mse_min = 1000;
    b_min = 0;
    % iterate through values of b till error with the data is minimized:
    for b = -5:0.01:5
        % compute von Mises distriution for this concentration parameter (b):
        vonMises = exp(b.*cos(2.*angles)) ./ (pi*besseli(0,b));
        % calculate MSE against the fiber distribution histogram:
        mse = mean((weights - vonMises).^2);
        if mse < mse_min
            mse_min = mse;
            b_min = b;
        end
    end
    % store value of b that gave min MSE:
    b_all(t) = b_min;
    % show fit on plot of data:
    figure
    plot(angles, exp(b_min.*cos(2.*angles)) ./ (pi*besseli(0,b_min))); hold on
    plot(angles, weights)
end

% plot resulting b values vs axial stretch:
figure
%scatter(data_lamBinned(:,1,1), b_all, 'k'); hold on
errorbar(data_lamBinned(:,1,1), b_all, data_lamBinned(:,1,2), 'horizontal', 'o'); hold on
plot([0,1.3], [0,0], 'k:', 'LineWidth', 1.0)
plot([1,1], [-0.6,1.6], 'k--', 'LineWidth', 1.0)
xlabel('stretch')
ylabel('concentration parameter, b')
xlim([0 1.3])
ylim([-0.6 1.6])


%% plot resulting b values vs axial stretch

figure
%scatter(data_lamBinned(:,1,1), b_all, 'k'); hold on
errorbar(DMSO_data_lamBinned(:,1,1), DMSO_b_all, DMSO_data_lamBinned(:,1,2), 'horizontal', 'o', 'Color', [0 0 0], "MarkerFaceColor", [0 0 0], "MarkerEdgeColor", [0 0 0], "LineWidth", 1.5); hold on
errorbar(CD_data_lamBinned(:,1,1), CD_b_all, CD_data_lamBinned(:,1,2), 'horizontal', 'o', 'Color', [0.5 0.8 1], "MarkerFaceColor", [0.5 0.8 1], "MarkerEdgeColor", [0.5 0.8 1], "LineWidth", 1.5);
plot([0,1.3], [0,0], 'k:', 'LineWidth', 1.0)
plot([1,1], [-0.6,1.6], 'k--', 'LineWidth', 1.0)
xlabel('fraction of λ*')
ylabel('concentration parameter, b')
xlim([0 1.3])
ylim([-0.6 1.6])


%% fitting DMSO and CD κ data (output from above) vs F_11 (axial stretch):
% (this gives the logistic fit for model input of a κ-stretch function)

f = fit([DMSO_b_forFitting(:,1);CD_b_forFitting(:,1)], [DMSO_b_forFitting(:,3);CD_b_forFitting(:,3)],'b1/(1+exp(-x))+b2','start',[8 -6.5]);
plot(DMSO_b_forFitting(1:9,1), DMSO_b_forFitting(1:9,3),'k', 'LineWidth', 1.5); hold on
plot(DMSO_b_forFitting(10:18,1), DMSO_b_forFitting(10:18,3),'k', 'LineWidth', 1.5)
plot(DMSO_b_forFitting(19:end,1), DMSO_b_forFitting(19:end,3),'k', 'LineWidth', 1.5)
plot(CD_b_forFitting(1:8,1), CD_b_forFitting(1:8,3),'c', 'LineWidth', 1.5); hold on
plot(CD_b_forFitting(9:18,1), CD_b_forFitting(9:18,3),'c', 'LineWidth', 1.5)
plot(CD_b_forFitting(19:end,1), CD_b_forFitting(19:end,3),'c', 'LineWidth', 1.5)
xvals = 1:0.1:5;
logisticFit = f.b1./(1+exp(-xvals))+f.b2;
plot(xvals, f.b1./(1+exp(-xvals))+f.b2)
plot([1 70], [0 0], 'k:')
xlabel("axial stretch"); ylabel("\kappa")
ylim([-1,3])
xlim([1,4])
