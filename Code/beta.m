function [xFit,redChi2Fit,xErr] = beta(int,intErr,temp,estTN,fitMin,fitMax,fact,offset,x0In,saveDir)
%BETA Fit for beta over a range of cutoffs
%   Send intensity, intensity error, temperature, estimated Neel
%   temperature, a fit min, and a fit max. Returns fitted parameters and
%   errorbars for fits to data above each cutoff value ranging from the fit
%   min to fit max of the estimated Neel temperature. Fact and offset vary
%   each parameter to determine an errorbar. x0In is a column vector of 2
%   initial values, not including the background or TN i.e. prefactor and
%   beta. Save figures and parameters to saveDir.

% Fit for beta
cutoffs = temp(temp>fitMin*estTN & temp<fitMax*estTN) - 0.001*ones(size(temp(temp>fitMin*estTN & temp<fitMax*estTN))); % Minimum temperature cutoffs to test
tempRed = @(TC, T) 1-T./TC; % Reduced temperature
model = @(x, T) x(1) + x(2).*tempRed(x(3), T).^(2.*x(4)).*(T<x(3)); % bg, I0, TC, beta. Flat background and critical scattering, see PAN manual from DAVE except for factor of 2.

% Fit to all points above one of the cutoff temperatures in each iteration
for i=1:length(cutoffs)
    tempFit = temp(temp>cutoffs(i));
    intFit = int(temp>cutoffs(i));
    intErrFit = intErr(temp>cutoffs(i));
    
    modelInput = @(x) model(x, tempFit);
    x0 = [mean(intFit(tempFit>estTN)), x0In(1), estTN, x0In(2)];
    errPts = 5e2;
    [xFit(:,i), redChi2Fit(:,i), xErr(:,i), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(intFit, intErrFit, modelInput, x0, errPts, fact, offset);
    
    % Plot goodness of fit
    tempCalc = linspace(tempFit(1), tempFit(end), 5e2); % Solid line over fitted region
    tempCalcFull = linspace(tempFit(1)*.8, tempFit(1), 5e2); % Dashed line beyond fitted region
    intCalc = model(xFit(:,i), tempCalc);
    intCalcFull = model(xFit(:,i), tempCalcFull);
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['\beta = ', num2str(xFit(4,i), 3), '\pm', num2str(xErr(4,i), 3)])
    xlabel('\it{T}\rm{ (K)}')
    ylabel('\it{I}\rm{ (cts / sec.)}')
    e1 = errorbar(temp, int, intErr, 'LineStyle', 'none', 'Marker', 'o');
    p1 = plot(tempCalc, intCalc, 'LineWidth', 1, 'Color', 'r');
    p2 = plot(tempCalcFull, intCalcFull, '--', 'LineWidth', 1, 'Color', 'r');
    x1 = xline(min(tempFit), '--k');
    legend([e1, p1, p2, x1],{'Data', 'Fit', 'Extrapolation', 'Lower Bound'})
    ylim([0, max(int + intErr)]);
    axis square
    box on
    hold off
    if (saveDir ~= "")
        exportgraphics(gcf, [saveDir, 'beta', strrep(num2str(tempFit(1)), '.', 'p'), '.png'])
    end

    disp(['Beta: ', num2str(xFit(4,i))])
    disp(['Beta error: ', num2str(xErr(4,i))])
    pause(0.1)
end

% Plot beta's versus lowest T point
figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
hold on
xlabel('\it{T}\rm{_{min} (K)}')
ylabel('\beta')
errorbar(temp(temp>fitMin*estTN & temp<fitMax*estTN), xFit(4,:), xErr(4,:), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
box on
pbaspect([16 9 1])
hold off

end