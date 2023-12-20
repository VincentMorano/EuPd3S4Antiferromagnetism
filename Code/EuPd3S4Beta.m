% Vincent Morano, 08/01/2022
%% Calculate beta for EuPd3S4 BT7 January (1/2,1/2,1/2) T-scan

clear
close all

set(0, 'defaultaxesfontsize', 7)

estTN = 2.89; % To determine cutoffs on data
fitMin = 0.85; % Minimum percent below transition to fit
fitMax = 0.95; % Maximum percent below transition to fit
fact = [0.5, 1, 0.01, 0];
offset = [0, 0, 0, 0.5];
x0 = [20, 0.5]; % Setting initial beta to mean-field so I don't bias the fitting too much. Prefactor, beta.
saveDirFigBin = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\Betas\';
saveDirFig = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\';
saveDirMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\';

% Import data as tables with variable names
filename1 = 'TscanwarmingLong1845818.bt7';
filename2 = 'TscanwarmingLong2845819.bt7';
filename3 = 'TscanwarmingLong3845820.bt7';

data1 = readtable(filename1, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 41);
fid1 = fopen(filename1); % Open file to take column headings as one line of text and edit
title1 = textscan(fid1, '%s', 100, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 40);
fclose(fid1);
title1 = title1{1};
data1.Properties.VariableNames = title1(2:end); % Remove "#Columns" from variable names

data2 = readtable(filename2, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 41);
fid2 = fopen(filename2); % Open file to take column headings as one line of text and edit
title2 = textscan(fid2, '%s', 100, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 40);
fclose(fid2);
title2 = title2{1};
data2.Properties.VariableNames = title2(2:end); % Remove "#Columns" from variable names

data3 = readtable(filename3, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 41);
fid3 = fopen(filename3); % Open file to take column headings as one line of text and edit
title3 = textscan(fid3, '%s', 100, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 40);
fclose(fid3);
title3 = title3{1};
data3.Properties.VariableNames = title3(2:end); % Remove "#Columns" from variable names

data = [data1; data2; data3]; % Concatenate tables

% Extract parameters for beta fitting
int = data.Detector./data.Time; % Normalize to time
intErr = sqrt(data.Detector)./data.Time;
temp = data.Temp;

% Fit beta and plot
[xFit,redChi2Fit,xErr] = beta(int,intErr,temp,estTN,fitMin,fitMax,fact,offset,x0,saveDirFigBin);

figure('Units', 'normalized', 'Position', [0, 0.3, 0.35, 0.6])
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
hold on
xlabel('\it{T}\rm{_{min} (K)}')
ylabel('\beta')
errorbar(temp(temp>fitMin*estTN & temp<fitMax*estTN), xFit(4,:), xErr(4,:), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
box on
pbaspect([16 9 1])
set(gca, 'TickLength', [0.02, 0.01])
hold off
%exportgraphics(gcf, [saveDirFig, 'beta.png'])
%save([saveDirMat, 'beta.mat']) % Save the workspace

% Pick a fitting range and make a colorplot of beta versus Tc

tMin = 2.6; % Minimum Tc for fitting based on the above beta vs Tc plot
[TcGrid, betaGrid] = meshgrid(linspace(2.85, 2.94, 200),linspace(0.25, 0.9, 200)); % Create grid of Tc and beta values

% Set initial values
TcPre = 2.89; % Preliminary Tc for selecting critical data range.
bgr0 = int(end); % Initial value for background fitting.
int0 = int(1) - int(end); % Initial value for intensity fitting.

% Restrict range of data for fitting
fitRange = temp>tMin;
intFit = int(fitRange);
intErrFit = intErr(fitRange);
tempFit = temp(fitRange);

% Plot my redChi^2(Tc ,beta)
numParamGrid = 2; % Now only fitting for two variables i.e. bgr and I0
redChi2Grid = zeros(size(TcGrid)); % Get my reduced chi-squared for every point on this grid
bgrGrid = zeros(size(TcGrid));
I0Grid = zeros(size(TcGrid));
for i = 1:length(TcGrid(:, 1))
    for j = 1:length(TcGrid(1, :))
        probGrid = optimproblem('ObjectiveSense', 'min'); % Stores optimization problem information. Default is minimization of the Objective function.
        xGrid = optimvar('xGrid', 4); % My optimization variables. bgr, I0, Tc, beta
        yGrid = @(bgrVar, I0Var, TcVar, betaVar) bgrVar.*ones(size(intFit)) + (I0Var.*(1 - tempFit./TcVar).^(2.*betaVar)).*(tempFit<=TcVar); % Objective function form for fit. A flat background with the intensity included for temperatures below Tc. Use logical factors to make anonymous piecewise function handles.
        intCalcGrid = fcn2optimexpr(yGrid, xGrid(1), xGrid(2), xGrid(3), xGrid(4)); % Normally optimproblem won't work with expressions of the form ...^x, so use this conversion function to fix that issue.
        objGrid = sum((intFit - intCalcGrid).^2./intErrFit.^2, 'all')/(length(intFit) - numParamGrid); % Reduced chi-squared.
        probGrid.Objective = objGrid; % Chi^2 is the Objective function to be minimized
        x0Grid.xGrid = [bgr0, int0, TcGrid(i, j), betaGrid(i, j)];
        probGrid.Constraints.cons1 = xGrid(3)==TcGrid(i, j); % Fit to specific Tc
        probGrid.Constraints.cons2 = xGrid(4)==betaGrid(i, j); % Fit to specific beta
        [solGrid, fvalGrid] = solve(probGrid, x0Grid); % Solve optimization problem
        bgrGrid(i, j) = solGrid.xGrid(1); % Save grid background
        I0Grid(i, j) = solGrid.xGrid(2); % Save grid I0
        redChi2Grid(i, j) = fvalGrid; % Save grid chi
    end
end

chiBest = min(redChi2Grid, [], 'All'); % Get best chi
indBest = find(redChi2Grid==chiBest); % Use index to get corresponding values
bgrBest = bgrGrid(indBest);
I0Best = I0Grid(indBest);
TcBest = TcGrid(indBest);
betaBest = betaGrid(indBest);

% Make a surface plot of the result
figure(3)
nexttile
hold on
s = surf(TcGrid, betaGrid, redChi2Grid);
s.EdgeColor = 'none'; % Don't show individual pixels
c = colorbar('north', 'Ticks', 1:0.1:1.5, 'TickLength', 0.02, 'Color', 'w');
clim([1, 1.5])
ylabel(c, '\chi_r^2', 'fontsize', 9)
xlabel('\it{T_C}\rm{ (K)}')
ylabel('\beta')
colormap(flipud(jet)) % So most of the plot's background is blue and my best fit is green
view(2) % Set plot along c-axis
deltaChi = chiBest*(1 + 1/(length(intFit) - numParamGrid));
%deltaChi = chiBest*1.2;
deltaChiGrid = deltaChi.*ones(size(TcGrid));
chiDiff = redChi2Grid - deltaChiGrid; % Use when this is zero to draw the intersect between a plane of height deltaChiGrid and my surface plot. From MATLAB Answers site.
C = contour(TcGrid, betaGrid, chiDiff, [0 0]);
TcL = C(1, 2:end);
betaL = C(2, 2:end);
chiL = interp2(TcGrid, betaGrid, redChi2Grid, TcL, betaL);
line(TcL, betaL, chiL, 'Color', 'k', 'LineWidth', 2);
ylim([0.25 0.9])
pbaspect([16 9 1])
set(gca, 'Layer', 'top', 'TickLength', [0.02, 0.01])
hold off

nexttile
hold on
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{I}\rm{ (det. / sec.)}')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
tempPlotPrelim = linspace(temp(1), temp(end), 400);
yPlotPrelim = bgrBest.*ones(size(tempPlotPrelim)) + (I0Best.*(1-tempPlotPrelim./TcBest).^(2.*betaBest)).*(tempPlotPrelim<=TcBest);
tempPlot = linspace(tempFit(1), tempFit(end), 400);
yPlot = bgrBest.*ones(size(tempPlot)) + (I0Best.*(1 - tempPlot./TcBest).^(2.*betaBest)).*(tempPlot<=TcBest);
ind3 = tempPlotPrelim>=min(tempFit);
plot(tempPlotPrelim(ind3), yPlotPrelim(ind3), 'r', 'LineWidth', 1.0)
plot(tempPlotPrelim(~ind3), yPlotPrelim(~ind3), '--r', 'LineWidth', 1.0)
%xline(min(tempFit), '--', {'Lower Fitting Limit'})
xline(min(TcBest), '--', {'\it{T_C}'}, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle', 'fontsize', 7)
pbaspect([16 9 1])
set(gca, 'TickLength', [0.02, 0.01])
box on
hold off
%exportgraphics(gcf, [saveDirFig, 'betaColorPlot.png'])
%save([saveDirMat, 'betaColorPlot.mat']) % Save the workspace

disp(' ')
disp('Optimal Parameters')
disp(['Background: ', num2str(bgrBest)])
disp(['Low-T Intensity: ', num2str(I0Best)])
disp(['T_c: ', num2str(TcBest,5), '+-', num2str((max(TcL) - min(TcL))/2, 2)])
disp(['Beta: ', num2str(betaBest,3), '+-', num2str((max(betaL) - min(betaL))/2, 3)])
disp(['Chi_r^2: ', num2str(chiBest)])