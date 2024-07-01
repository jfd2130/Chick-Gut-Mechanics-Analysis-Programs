%% calculate "alignment stretch" λ_a

% "alignment" is defined in terms of the proportion of the area of the
% distibution in +/-10 deg around 0 deg (specifically when the rate of
% change of this value is minimized, i.e. when the second derivative = 0)

% all inputs:
data_3Dplot = CDexp6_data_3Dplot; % DMSOexp4,5,7 CDexp3,4,6
lam_star = 1.93; % 3.16 2.88 2.03 1.56 2.04 1.93

% extract sum of bins in +/-10 around 0 deg (bins are 4 deg each, so middle 5 bins):
[~,i_axialRow] = min(abs(data_3Dplot(4:end,1)));
i_axialRow = i_axialRow + 3;
data_axialSum = sum(data_3Dplot(i_axialRow-2:i_axialRow+2, 2:end), 1);

% calculate the proportion of the distribution in the above angle range:
data_allSum = sum(data_3Dplot(4:end, 2:end), 1);
alignment = data_axialSum ./ data_allSum;

% find axial stretch of alignment:
[f,~] = fit(data_3Dplot(2, 2:end)', alignment', 'b1/(1+exp(-b2*x+b3))+b4', 'start', [0.5 5 7 0.1]);
vals_forFit = linspace(1,4,1000);
fitPoints = f.b1 ./ (1 + exp(-f.b2 .* vals_forFit + f.b3)) + f.b4;
secondDeriv = f.b1 .* ( (2 .* f.b2.^2 .* exp(2 .* f.b3 - 2 .* f.b2.* vals_forFit)) ./ (exp(f.b3 - f.b2.* vals_forFit) + 1).^3  -  (f.b2.^2 .* exp(f.b3 - f.b2.* vals_forFit)) ./ (exp(f.b3 - f.b2.* vals_forFit) + 1).^2 );
[~,i_minSecDer] = min(secondDeriv);
stretch_align = vals_forFit(i_minSecDer);

plot_data = plot(data_3Dplot(2,2:end), alignment, 'LineWidth', 2); hold on
plot_isoLine = plot([1,3.5], [20/180,20/180], 'k', 'LineWidth', 0.5); hold on
plot_fit = plot(vals_forFit, fitPoints,'Color', get(plot_data,'Color'), 'LineWidth', 0.5); hold on
plot_line = plot([stretch_align, stretch_align], [0, max(alignment)], 'Color', get(plot_data,'Color'), 'LineStyle', '--', 'LineWidth', 2);
plot_line = plot([lam_star, lam_star], [0, max(alignment)], 'Color', get(plot_data,'Color'), 'LineStyle', ':', 'LineWidth', 1.5);
xlabel("axial stretch")
ylabel("alignment")


















