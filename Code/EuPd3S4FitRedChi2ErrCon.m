function [xFit,redChi2Fit,xErr,chiUpper,chiTrial,paramTrial,interpPts,slopes,intercepts,paramLower,paramUpper] = EuPd3S4FitRedChi2ErrCon(data,dataErr,model,x0,errPts,fact,offset,lb,ub)
%fitRedChi2ErrCon Calculate errorbars on parameters fitted with a reduced
%chi-squared and bounds
%   Give data, errorbars for the data, a fitting model, initial fitting
%   values, the number of sampled error points per run, the number of
%   interpolation runs, a factor, and an offset for trial values. Return
%   the fitted parameters, the reduced chi-squared, the errorbars on each
%   fitted parameter, the threshold gof, the errorbar gof's, the errorbar
%   fitted parameters, the indices of the four interpolated points, the
%   slopes and intercepts connecting them, and the upper and lower
%   parameters used for the 2*errorbar calculation. lb and ub are bounds
%   for fmincon().

redChi2 = @(fnObs, fnErr, fnCal, nParam) sum((fnObs-fnCal).^2./fnErr.^2, 'all')/(length(fnObs)-nParam);
obj = @(x) redChi2(data, dataErr, model(x), length(x0));
options = optimoptions(@fmincon, 'Display', 'off');
[xFit, redChi2Fit] = fmincon(obj, x0', [], [], [], [], lb, ub, [], options);

% Use the fitted parameters above to generate trial parameters for the
% errorbar calculation with the input factors and terms in trialParam.
chiUpper = redChi2Fit*(1+1/(length(data)-length(xFit))); % Upper bound of the gof for determining the errorbar. Check this since dof changes and may be undefined.

paramTrial = nan(errPts, length(xFit));
for i = 1:length(xFit)
    paramTrial(:, i) = linspace(xFit(i)*(1-fact(i))-offset(i), xFit(i)*(1+fact(i))+offset(i), errPts);
end
xTrialFit = nan(errPts, length(xFit)-1, length(xFit)); % Each column is the parameter value after fitting, depth is going from one iterated parameter to the next while holding other constant.
chiTrial = nan(errPts, length(xFit)); % Chi2 for every fit for the errorbar
xErr = nan(size(xFit'));
slopes = nan(2, length(xFit)); % First element is for the lower-parameter, second upper
intercepts = nan(2, length(xFit));
interpPts = nan(4, length(xFit));
paramLower = nan(1, length(xFit));
paramUpper = nan(1, length(xFit));

for j = 1:length(xFit)
    initParam = xFit';  % Fitting to one fewer parameter, only send the initial values for the fitted parameters. Can make this more efficient.
    initParam(j) = [];
    if j==1
        modelCalcErr = @(x, iterParam) model([iterParam, x(j:length(xFit)-1)]); % Notation is fun(var, x(1), x(2), x(3), x(4)), fun(x(1), var, x(2), x(3), x(4)), etc.
    elseif j<length(xFit)
        modelCalcErr = @(x, iterParam) model([x(1:j-1), iterParam, x(j:length(xFit)-1)]);
    else
        modelCalcErr = @(x, iterParam) model([x(1:j-1), iterParam]);
    end
    
    % Resize constraints for fitting to length(xFit)-1 variables
    lbErr=lb;
    ubErr=ub;
    lbErr(j)=[];
    ubErr(j)=[];
    
    % Fit to the peak
    for k = 1:errPts
        modelInputErr = @(x) modelCalcErr(x, paramTrial(k, j));
        obj = @(x) redChi2(data, dataErr, modelInputErr(x), length(x0)); % Using nParam=length(x0) rather than length(initParam) so gof's match when iterating parameter
        [xTrialFit(k,:,j), chiTrial(k,j)] = fmincon(obj, initParam, [], [], [], [], lbErr, ubErr, [], options);
    end
    int2 = find(chiTrial(:, j) < chiUpper & paramTrial(:, j) < xFit(j), 1, 'first'); % First point within chiUpper
    int1Pre = find(chiTrial(:, j) > chiUpper & paramTrial(:, j) < xFit(j)); % List of points before first within chiUpper.
    int1 = int1Pre(find(int1Pre < int2, 1, 'last')); % Last point before the first below chiUpper. Make sure to get the index, not just the index of the index (if only use find).
    int3 = find(chiTrial(:, j) < chiUpper & paramTrial(:, j) > xFit(j), 1, 'last'); % Last point within chiUpper
    int4Pre = find(chiTrial(:, j) > chiUpper & paramTrial(:, j) > xFit(j)); % List of points after last within chiUpper.
    int4 = int4Pre(find(int4Pre > int3, 1, 'first')); % First point after the last below chiUpper.
    if isempty(int1) || isempty(int2) || isempty(int3) || isempty(int4)
        disp('Error too large.') % Could not find the intercepts required to calculate an error bar, likely because the parameter wasn't changed enough to pass the gof threshold.
    else
        interpPts(:, j) = [int1; int2; int3; int4]; % Largest parameter with unsuitable chi that's less than the optimal parameter, smallest parameter with suitable chi that's less than the optimal parameter, largest parameter with suitable chi that's greater than the optimal parameter, smallest parameter with unsuitable chi that's greater than the optimal parameter
        slopes(1, j) = (chiTrial(interpPts(2, j), j)-chiTrial(interpPts(1, j), j))/(paramTrial(interpPts(2, j), j)-paramTrial(interpPts(1, j), j));
        slopes(2, j) = (chiTrial(interpPts(4, j), j)-chiTrial(interpPts(3, j), j))/(paramTrial(interpPts(4, j), j)-paramTrial(interpPts(3, j), j));
        intercepts(1, j) = chiTrial(interpPts(1, j), j)-slopes(1, j)*paramTrial(interpPts(1, j), j);
        intercepts(2, j) = chiTrial(interpPts(3, j), j)-slopes(2, j)*paramTrial(interpPts(3, j), j);
        paramLower(j) = (chiUpper-intercepts(1, j))/slopes(1, j); % Greatest integrated intensity with unsuitable chi that's smaller than the optimal intensity
        paramUpper(j) = (chiUpper-intercepts(2, j))/slopes(2, j); % Smallest integrated intensity with unsuitable chi that's greater than the optimal intensity
        xErr(j) = abs(paramUpper(j)-paramLower(j))/2; % Size of the errorbar for integrated intensity
    end
end
end