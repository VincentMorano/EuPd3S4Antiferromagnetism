%% Vincent Morano, 01/25/2021

% Notes on filenames:
% 845748:845750, 845768:845798, 845800:845803 are 5 K, Exp 1
% 845614:845640, 845644:845651, 845654:845656, 845658:845660 are 1.5 K, Exp 1
% 847555:847561 are 4.5 K, Exp 2
% 847443:847472, 847475, 847478:847487 are 0.5 T, Exp 2 AFM. Some must be combined.
% 847424:847442 are 0.5 T, Exp 2 FM.
% 847493:847511 are 3 T, Exp 2

% The monitor had a problem at 845721 so to compare scans before and after
% this (expt 1) I need to normalize to time rather than monitor.

%% Fit for integrated intensities of rocking scans, calculate structure factors.

clear
close all

set(0, 'defaulttextinterpreter', 'latex')
set(0, 'defaultlegendinterpreter', 'latex')

factor1 = 1; % Scale factor to apply to the experiment 1 squared structure factors. Useful if they're too low for FullProf.
factor2 = 1e3; % Scale factor to apply to the experiment 2 squared structure factors. Useful if they're too low for FullProf.
bgThresh = inf; % Don't save points fit with bg higher than this.
dropRange = 1.5; % Angular range about powder lines to drop if within. +-
errPts = 1e3; % Number of parameter iterations when calculating the errorbar
factSig = [0.5, 1.0, 0.0, 0.1]; % Determines range of iterated values for errorbar. [bg, area, cent, sig]
offsetSig = [0.0, 0.0, 0.3, 0.1]; % Determines range of iterated values for errorbar. [bg, area, cent, sig]
fact1 = [0.5, 1.0, 0.0]; % Determines range of iterated values for errorbar. [bg, area, cent] Expt 1.
offset1 = [3.0, 5, 0.1]; % Determines range of iterated values for errorbar. [bg, area, cent]
fact2 = [0.5, 1.0, 0.0]; % Determines range of iterated values for errorbar. [bg, area, cent] Expt 2.
offset2 = [3e-3, 5e-3, 0.1]; % Determines range of iterated values for errorbar. [bg, area, cent]

cuLines = [42.96887, 50.03655, 73.46521, 89.06480, 94.19227, 115.51865, 134.35330, 142.04612]; % A4 values of copper lines. ICSD 52256
alLines = [38.16995, 44.36433, 64.54410, 77.52596, 81.67868, 98.06911, 110.74359, 115.17975, 135.27979, 157.57358]; % ICSD 52611
expt1a = [845748:845750, 845768:845798, 845800:845803]; % 5K, Exp 1
expt1b = [845614:845640, 845644:845651, 845654:845656, 845658:845660]; % 1.5 K, Exp 1
expt2a = 847555:847561; % 4.5 K, Exp 2
expt2b = [847443:847472, 847475, 847478:847487]; % 0.5 T, Exp 2 AFM
expt2c = 847424:847442; % 0.5 T, Exp 2 FM
expt2d = 847493:847511; % 3 T, Exp 2
sumFiles = [847445, 847465; 847446, 847466; 847447, 847467; 847448, 847468; 847449, 847469; 847450, 847470; 847454, 847471; 847455, 847472; 847456, 847475; 847457, 847478; 847458, 847479; 847463, 847480; 847486, 847487];

% Import data
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\EuPd3S4\All\';
fileNum = [expt1a, expt1b, expt2a, expt2b, expt2c, expt2d];
fileName = strcat('fpx', string(fileNum), '.bt7');
headerLines = 41;
scans = importDataBT7(fileLoc, fileName, headerLines);

% Keep track of scan types and desired normalization
for i = 1:length(scans)
    scans(i).fileNum = fileNum(i);
    if ismember(fileNum(i), expt1a)
        scans(i).scanType = 1;
        scans(i).norm = 0; % 0 Time, 1 Monitor
    elseif ismember(fileNum(i), expt1b)
        scans(i).scanType = 2;
        scans(i).norm = 0;
    elseif ismember(fileNum(i), expt2a)
        scans(i).scanType = 3;
        scans(i).norm = 1;
    elseif ismember(fileNum(i), expt2b)
        scans(i).scanType = 4;
        scans(i).norm = 1;
    elseif ismember(fileNum(i), expt2c)
        scans(i).scanType = 5;
        scans(i).norm = 1;
    elseif ismember(fileNum(i), expt2d)
        scans(i).scanType = 6;
        scans(i).norm = 1;
    end
end

% Sum identical scans, including detector and monitor. Drop the old files and append the new ones.
for i = 1:length(sumFiles) % Iterate through each set of files that needs to be summed
    sumInd1 = find(arrayfun(@(x) x.fileNum, scans(1:end-(i-1)))==sumFiles(i,1)); % Find the index of the structure array element i.e. scan to be summed. Be careful that the list of filenames is going to be updated, so don't use the initial copy. arrayfun reccommended in the note here: https://www.mathworks.com/help/matlab/matlab_prog/create-a-structure-array.html. This fails for empty structure element, hence the index.
    sumInd2 = find(arrayfun(@(x) x.fileNum, scans(1:end-(i-1)))==sumFiles(i,2));
    scans(end+1).det = scans(sumInd1).det + scans(sumInd2).det; % Adding new file on to end of shrinking array.
    scans(end).detErr = sqrt(scans(sumInd1).detErr.^2 + scans(sumInd2).detErr.^2);
    scans(end).mon = scans(sumInd1).mon + scans(sumInd2).mon;
    scans(end).monErr = sqrt(scans(sumInd1).monErr.^2 + scans(sumInd2).monErr.^2);
    scans(end).time = scans(sumInd1).time + scans(sumInd2).time;
    scans(end).a3 = (scans(sumInd1).a3 + scans(sumInd2).a3)./2; % These should be the same between scans, be sure.
    scans(end).a4 = (scans(sumInd1).a4 + scans(sumInd2).a4)./2;
    scans(end).T = (scans(sumInd1).T + scans(sumInd2).T)./2;
    scans(end).B = (scans(sumInd1).B + scans(sumInd2).B)./2;
    scans(end).E = (scans(sumInd1).E + scans(sumInd2).E)./2;
    scans(end).H = (scans(sumInd1).H + scans(sumInd2).H)./2;
    scans(end).K = (scans(sumInd1).K + scans(sumInd2).K)./2;
    scans(end).L = (scans(sumInd1).L + scans(sumInd2).L)./2;
    scans(end).HKL = (scans(sumInd1).HKL + scans(sumInd2).HKL)./2;
    scans(end).meanT = (scans(sumInd1).meanT + scans(sumInd2).meanT)./2;
    scans(end).meanB = (scans(sumInd1).meanB + scans(sumInd2).meanB)./2;
    scans(end).meanE = (scans(sumInd1).meanE + scans(sumInd2).meanE)./2;
    scans(end).scanType = (scans(sumInd1).scanType + scans(sumInd2).scanType)./2;
    scans(end).norm = (scans(sumInd1).norm + scans(sumInd2).norm)./2;
    scans(end).meanH = round(mean(scans(end).H), 3);
    scans(end).meanK = round(mean(scans(end).K), 3);
    scans(end).meanL = round(mean(scans(end).L), 3);
    scans(end).meanHKL = round(mean(scans(end).HKL), 3);
    scans(end).meanA3 = round(mean(scans(end).a3), 3);
    scans(end).meanA4 = round(mean(scans(end).a4), 3);
    scans([sumInd1, sumInd2]) = []; % This is to remove the elements simultaneously rather than one after the other with an uncorrected index.
end
for i=1:length(scans)
    scans(i).W=0; % Energy transfer
    scans(i).Exp=EuPd3S4Exp1; % ResLib file. Confirm same for both experiments
    lattice=GetLattice(scans(i).Exp);
    scans(i).a=lattice.a; % Cubic
    scans(i).Q=2*pi/scans(i).a*sqrt(scans(i).meanH.^2 + scans(i).meanK.^2 + scans(i).meanL.^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
    if scans(i).norm==0 % Normalize to time
        scans(i).int = scans(i).det./scans(i).time;
        scans(i).intErr = scans(i).detErr./scans(i).time;
    elseif scans(i).norm==1 % Normalize to monitor
        scans(i).int = scans(i).det./scans(i).mon;
        scans(i).intErr = (scans(i).det./scans(i).mon).*sqrt((scans(i).detErr./scans(i).det).^2 + (scans(i).monErr./scans(i).mon).^2);
    end
end

% Fit for the widths of a few prominent peaks in each phase
% 5 K, Exp 1| 845748:845750, 845768:845798, 845800:845803
% Prominent: 845748:845750, 845774, 845800 so i = 1:3, 10, 35
sigInd1 = [1:3 10 35];
% 1.5 K, Exp 1| 845614:845640, 845644:845651, 845654:845656, 845658:845660
% Prominent: 845614, 845615 so i = 39, 40
sigInd2 = [39 40];
% 4.5 K, Exp 2| 847555:847561
% Prominent: 847558:847560 so i = 83:85
sigInd3 = 83:85;
% 0.5 T, Exp 2 AFM (some combined)| 847443:847472, 847475, 847478:847487
% Prominent: 847443, 847444, 847465 so (1/2,1/2,1/2), (1/2,1/2,3/2), (1/2,1/2,5/2) so i = 87, 88, 140
sigInd4 = [87 88 140];
% 0.5 T, Exp 2 FM| 847424:847442
% Prominent: 847426, 847429, 847432, 847433, 847439, 847441 so i = 123, 126, 129, 130, 136, 138
sigInd5 = [123 126 129 130 136 138];
% 3 T, Exp 2| 847493:847511
% Prominent: 847493, 847495, 847497, 847498, 847501, 847502, 847508, 847510 so i = 102, 104, 106, 107, 110, 111, 117, 119
sigInd6 = [102 104 106 107 110 111 117 119];

model = @(x, a3) x(1) + x(2)./(sqrt(2*pi).*x(4)).*exp(-(a3-x(3)).^2./(2.*x(4).^2)); % bg, area, center, std. Gaussian peak in terms of area and standard deviation.
for j = 1:6
    if j == 1
        sigInd = sigInd1;
    elseif j == 2
        sigInd = sigInd2;
    elseif j == 3
        sigInd = sigInd3;
    elseif j == 4
        sigInd = sigInd4;
    elseif j == 5
        sigInd = sigInd5;
    elseif j == 6
        sigInd = sigInd6;
    end
    tmpSig = zeros(size(sigInd));

    for i = 1:length(sigInd)
        % Set initial values, fit my minimizing reduced chi-squared
        scans(sigInd(i)).farFromPowder = sum( abs(scans(sigInd(i)).meanA4.*ones(size(cuLines)) - cuLines) > dropRange ) == length(cuLines) & sum(abs(mean(scans(sigInd(i)).a4).*ones(size(alLines)) - alLines) > dropRange) == length(alLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
        modelInput = @(x) model(x, scans(sigInd(i)).a3);
        tmp = sort(scans(sigInd(i)).int); % Sorted intensities
        bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
        cent0 = sum(scans(sigInd(i)).a3.*(scans(sigInd(i)).int-bg0))./sum(scans(sigInd(i)).int-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
        sig0 = 0.1;
        area0 = (max(scans(sigInd(i)).int)-bg0)*sqrt(2*pi)*sig0; % The prefactor of the Gaussian is the peak height so multiply by sqrt(2*pi)*sig0 to get the area. Guessing the height is the difference between the max observed intensity and estimated background.
        scans(sigInd(i)).x0Sig = [bg0, area0, cent0, sig0]; % Starting point for fitting.
        [scans(sigInd(i)).xFitSig, scans(sigInd(i)).redChi2FitSig, scans(sigInd(i)).xErrSig, scans(sigInd(i)).chiUpperSig, scans(sigInd(i)).chiTrialSig, scans(sigInd(i)).paramTrialSig, scans(sigInd(i)).interpPtsSig, scans(sigInd(i)).slopesSig, scans(sigInd(i)).interceptsSig, scans(sigInd(i)).paramLowerSig, scans(sigInd(i)).paramUpperSig] = EuPd3S4FitRedChi2Err(scans(sigInd(i)).int, scans(sigInd(i)).intErr, modelInput, scans(sigInd(i)).x0Sig, errPts, factSig, offsetSig);
        scans(sigInd(i)).include = (scans(sigInd(i)).xFitSig(1)<bgThresh) & scans(sigInd(i)).farFromPowder;
    
        % Plot fits
        a3Cal = linspace(min(scans(sigInd(i)).a3), max(scans(sigInd(i)).a3), 200)';
        intCal = model(scans(sigInd(i)).xFitSig, a3Cal);
        
        close all
        figure('Units', 'normalized', 'Position', [0, 0.2, 0.5, 0.7])
        clf
        hold on
        title(['Rocking Scan Fit: $\left(', num2str(scans(sigInd(i)).meanHKL(1)), ', ', num2str(scans(sigInd(i)).meanHKL(2)), ', ', num2str(scans(sigInd(i)).meanHKL(3)), '\right)$, ', num2str(scans(sigInd(i)).meanT, 5), ' K, ', num2str(scans(sigInd(i)).meanE, 3), ' meV, ', '$\left( \frac{', num2str(sigInd(i)), '}{', num2str(length(scans)), '} \right)$'])
        xlabel('A3 (deg.)')
        if scans(sigInd(i)).norm==0 % Normalized to time
            ylabel('Intensity $\left( \frac{\mbox{cts}}{\mbox{sec}} \right)$')
        elseif scans(sigInd(i)).norm==1 % Normalized to monitor
            ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$')
        end
        errorbar(scans(sigInd(i)).a3, scans(sigInd(i)).int, scans(sigInd(i)).intErr, 'o')
        plot(a3Cal, intCal)
        set(gca,'FontSize',12)
        xlim([min(scans(sigInd(i)).a3), max(scans(sigInd(i)).a3)])
        box on
        pbaspect([1 1.618 1])
        hold off
        saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\BT7Reduction\fpx', num2str(scans(sigInd(i)).fileNum), 'Sig', num2str(j), 'Curve.png'])
        
        EuPd3S4PlotBGGaussianErrSig(scans, sigInd(i));
        saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\BT7Reduction\fpx', num2str(scans(sigInd(i)).fileNum), 'Sig', num2str(j), 'Error.png'])
    
        tmpSig(i) = scans(sigInd(i)).xFitSig(4);
    end

    if j == 1
        sig1 = mean(tmpSig);
    elseif j == 2
        sig2 = mean(tmpSig);
    elseif j == 3
        sig3 = mean(tmpSig);
    elseif j == 4
        sig4 = mean(tmpSig);
    elseif j == 5
        sig5 = mean(tmpSig);
    elseif j == 6
        sig6 = mean(tmpSig);
    end
end

% Fit for integrated intensities
for i=1:length(scans)    
    % Set initial values, fit my minimizing reduced chi-squared, apply a scale factor to the intensities if desired
    if ismember(fileNum(i), expt1a)
        sig = sig1;
        fact = fact1;
        offset = offset1;
    elseif ismember(fileNum(i), expt1b)
        sig = sig2;
        fact = fact1;
        offset = offset1;
    elseif ismember(fileNum(i), expt2a)
        sig = sig3;
        fact = fact2;
        offset = offset2;
    elseif ismember(fileNum(i), expt2b)
        sig = sig4;
        fact = fact2;
        offset = offset2;
    elseif ismember(fileNum(i), expt2c)
        sig = sig5;
        fact = fact2;
        offset = offset2;
    elseif ismember(fileNum(i), expt2d)
        sig = sig6;
        fact = fact2;
        offset = offset2;
    end
    scans(i).farFromPowder = sum( abs(scans(i).meanA4.*ones(size(cuLines)) - cuLines) > dropRange ) == length(cuLines) & sum(abs(mean(scans(i).a4).*ones(size(alLines)) - alLines) > dropRange) == length(alLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    model = @(x, a3) x(1) + x(2)./(sqrt(2*pi).*sig).*exp(-(a3-x(3)).^2./(2.*sig.^2)); % bg, area, center. Gaussian peak in terms of area and standard deviation.
    modelInput = @(x) model(x, scans(i).a3);
    tmp = sort(scans(i).int); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    area0 = (max(scans(i).int)-bg0)*sqrt(2*pi)*sig; % The prefactor of the Gaussian is the peak height so multiply by sqrt(2*pi)*sig0 to get the area. Guessing the height is the difference between the max observed intensity and estimated background.
    cent0 = sum(scans(i).a3.*(scans(i).int-bg0))./sum(scans(i).int-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    scans(i).x0 = [bg0, area0, cent0]; % Starting point for fitting.
    lb = [-inf, -inf, min(scans(i).a3)]; % Lower bound for peak fitting.
    ub = [inf, inf, max(scans(i).a3)]; % Upper bound for peak fitting.
    [scans(i).xFit, scans(i).redChi2Fit, scans(i).xErr, scans(i).chiUpper, scans(i).chiTrial, scans(i).paramTrial, scans(i).interpPts, scans(i).slopes, scans(i).intercepts, scans(i).paramLower, scans(i).paramUpper] = EuPd3S4FitRedChi2ErrCon(scans(i).int, scans(i).intErr, modelInput, scans(i).x0, errPts, fact, offset, lb, ub);
    scans(i).include = (scans(i).xFit(1)<bgThresh) & scans(i).farFromPowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.2, 0.5, 0.7])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).meanHKL(1)), ', ', num2str(scans(i).meanHKL(2)), ', ', num2str(scans(i).meanHKL(3)), '\right)$, ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'])
    xlabel('A3 (deg.)')
    if scans(i).norm==0 % Normalized to time
        ylabel('Intensity $\left( \frac{\mbox{cts}}{\mbox{sec}} \right)$')
    elseif scans(i).norm==1 % Normalized to monitor
        ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$')
    end
    errorbar(scans(i).a3, scans(i).int, scans(i).intErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    pbaspect([1 1.618 1])
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\BT7Reduction\fpx', num2str(scans(i).fileNum), 'Curve.png'])
    
    EuPd3S4PlotBGGaussianErr(scans, i);
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\BT7Reduction\fpx', num2str(scans(i).fileNum), 'Error.png'])
end

% Calculate structure factors
for i = 1:length(scans)
    [scans(i).R0, scans(i).M] = ResMat(scans(i).Q, scans(i).W, scans(i).Exp);
    [scans(i).f2, scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);

    if ismember(fileNum(i), expt1a)
        scans(i).f2 = factor1*scans(i).f2;
        scans(i).f2Err = factor1*scans(i).f2Err;
    elseif ismember(fileNum(i), expt1b)
        scans(i).f2 = factor1*scans(i).f2;
        scans(i).f2Err = factor1*scans(i).f2Err;
    elseif ismember(fileNum(i), expt2a)
        scans(i).f2 = factor2*scans(i).f2;
        scans(i).f2Err = factor2*scans(i).f2Err;
    elseif ismember(fileNum(i), expt2b)
        scans(i).f2 = factor2*scans(i).f2;
        scans(i).f2Err = factor2*scans(i).f2Err;
    elseif ismember(fileNum(i), expt2c)
        scans(i).f2 = factor2*scans(i).f2;
        scans(i).f2Err = factor2*scans(i).f2Err;
    elseif ismember(fileNum(i), expt2d)
        scans(i).f2 = factor2*scans(i).f2;
        scans(i).f2Err = factor2*scans(i).f2Err;
    end
end
save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\StrFact.mat') % Save the Workspace

%% Save structure factors to .int files.

clear
close all

% Load peak-fitted workspace
load('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\StrFact.mat')

fileID5K=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S45K.int','w');
formatSpec='EuPd3S4 BT7 (January 2020) INT File: 5K, 35meV\n';
fprintf(fileID5K, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID5K, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID5K, formatSpec);
fileID1p5K=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S41p5K.int','w');
formatSpec='EuPd3S4 BT7 (January 2020) INT File: 1.5 K, 35 meV\n';
fprintf(fileID1p5K, formatSpec);
formatSpec='(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5K, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID1p5K, formatSpec);
formatSpec='1\n';
fprintf(fileID1p5K, formatSpec);
formatSpec='1 0.5 0.5 0.5\n';
fprintf(fileID1p5K, formatSpec);
fileID4p5K=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S44p5K.int','w');
formatSpec='EuPd3S4 BT7 (March 2020) INT File: 4.5 K, 35 meV\n';
fprintf(fileID4p5K, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID4p5K, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID4p5K, formatSpec);
fileID0p5TAFM=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S40p5TAFM.int','w');
formatSpec='EuPd3S4 BT7 (March 2020) INT File: 0.3 K, 0.5 T, 35 meV, AFM\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec='(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec='1\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec='1 0.5 0.5 0.5\n';
fprintf(fileID0p5TAFM, formatSpec);
fileID0p5TFM=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S40p5TFM.int','w');
formatSpec='EuPd3S4 BT7 (March 2020) INT File: 0.3 K, 0.5 T, 35 meV, FM\n';
fprintf(fileID0p5TFM, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID0p5TFM, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID0p5TFM, formatSpec);
fileID3T=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S43T.int','w');
formatSpec='EuPd3S4 BT7 (March 2020) INT File: 0.4 K, 3 T, 35 meV\n';
fprintf(fileID3T, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID3T, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID3T, formatSpec);
for i=1:length(scans)
    if scans(i).scanType==1 && scans(i).include
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID5K, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType==2 && scans(i).include
        formatSpec='%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID1p5K, formatSpec, scans(i).meanH-0.5, scans(i).meanK-0.5, scans(i).meanL-0.5, 1, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType==3 && scans(i).include
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID4p5K, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType==4 && scans(i).include
        formatSpec='%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID0p5TAFM, formatSpec, scans(i).meanH-0.5, scans(i).meanK-0.5, scans(i).meanL-0.5, 1, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType==5 && scans(i).include
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID0p5TFM, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType==6 && scans(i).include
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID3T, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
    end
end
fclose(fileID5K);
fclose(fileID1p5K);
fclose(fileID4p5K);
fclose(fileID0p5TAFM);
fclose(fileID0p5TFM);
fclose(fileID3T);

%% Save spherical absorption corrected structure factors to .int files. Only including Eu for now.

clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')

% Load peak-fitted workspace
load('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\StrFact.mat')

% Determine the absorption factors
volSamp = 5; % Volume of sample in mm3. See EuPd3S4 033022.
atoms = 2; % Per unit cell
vol = 297.502; % Of unit cell in angstroms, our XRD from SI 213 K.
energy = 35; % Energy in meV

sigAbs = 4530*sqrt(25.3/energy); % barns. Times 10^-22 for mm^2
dens = atoms/vol; % Inverse cubic Ang. Times 10^21 for mm^-3
muAbs = dens*sigAbs/10; % Absorption coefficient in inverse mm. 
%radEff = (volSamp*3/4/pi)^(1/3); % Effective radius for spherical absorption correction in mm
radEff = 0.5; % Effective radius for spherical absorption correction in mm

fun = @(R, theta, mu, r, alpha, phi) sin(alpha).*(3./(4.*pi.*R.^3)).*exp(-mu.*(sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)+phi).^2)+sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)-phi).^2)-2.*r.*sin(theta.*(pi/180)).*sin(alpha).*sin(phi))).*r.^2; % Picked up the sin(alpha) because want to integrate over alpha rather than cos(alpha), then flipped integral bounds.
funA = @(R, theta, mu) integral3(@(r, alpha, phi) fun(R, theta, mu, r, alpha, phi), 0, R, 0, pi, 0, 2*pi); % Transmission coefficient. integral3 seems to plug in an array of values so need dots in front of operations in function. Mentioned in documentation.

AStar = zeros(length(scans), 1);
for i = 1:length(scans)
    AStar(i) = 1./funA(radEff, scans(i).meanA4/2, muAbs);
end

j = 1;
k = 1;
for i = 1:length(scans)
    if scans(i).scanType == 1 && scans(i).include
        a45p1K(j) = scans(i).meanA4;
        AStar5p1K(j) = AStar(i);
        j = j+1;
    elseif scans(i).scanType == 2 && scans(i).include
        a41p5K(k) = scans(i).meanA4;
        AStar1p5K(k) = AStar(i);
        k = k+1;
    end
end

figure(1)
hold on
s1 = scatter(a45p1K, AStar5p1K, 'r');
s2 = scatter(a41p5K, AStar1p5K, 'b');
legend([s1, s2], {'5.1 K', '1.5 K'}, 'Location', 'best')
ylabel('\it{A*}')
xlabel('2\theta (deg.)')
box on
axis square
hold off

% Save the files
fileID5K = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S45KSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (January 2020) INT File: 5K, 35meV, Spherical Correction\n';
fprintf(fileID5K, formatSpec);
formatSpec = '(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID5K, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID5K, formatSpec);
fileID1p5K = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S41p5KSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (January 2020) INT File: 1.5 K, 35 meV, Spherical Correction\n';
fprintf(fileID1p5K, formatSpec);
formatSpec = '(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5K, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID1p5K, formatSpec);
formatSpec = '1\n';
fprintf(fileID1p5K, formatSpec);
formatSpec = '1 0.5 0.5 0.5\n';
fprintf(fileID1p5K, formatSpec);
fileID4p5K = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S44p5KSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (March 2020) INT File: 4.5 K, 35 meV, Spherical Correction\n';
fprintf(fileID4p5K, formatSpec);
formatSpec = '(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID4p5K, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID4p5K, formatSpec);
fileID0p5TAFM = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S40p5TAFMSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (March 2020) INT File: 0.3 K, 0.5 T, 35 meV, AFM, Spherical Correction\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec='(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec = '1\n';
fprintf(fileID0p5TAFM, formatSpec);
formatSpec = '1 0.5 0.5 0.5\n';
fprintf(fileID0p5TAFM, formatSpec);
fileID0p5TFM = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S40p5TFMSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (March 2020) INT File: 0.3 K, 0.5 T, 35 meV, FM, Spherical Correction\n';
fprintf(fileID0p5TFM, formatSpec);
formatSpec = '(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID0p5TFM, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID0p5TFM, formatSpec);
fileID3T = fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\EuPd3S4\EuPd3S43TSphAbs.int','w');
formatSpec = 'EuPd3S4 BT7 (March 2020) INT File: 0.4 K, 3 T, 35 meV, Spherical Correction\n';
fprintf(fileID3T, formatSpec);
formatSpec = '(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID3T, formatSpec);
formatSpec = '1.5289 0 0\n';
fprintf(fileID3T, formatSpec);
for i = 1:length(scans)
    if scans(i).scanType == 1 && scans(i).include
        formatSpec = '%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID5K, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType == 2 && scans(i).include
        formatSpec = '%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID1p5K, formatSpec, scans(i).meanH-0.5, scans(i).meanK-0.5, scans(i).meanL-0.5, 1, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType == 3 && scans(i).include
        formatSpec = '%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID4p5K, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType == 4 && scans(i).include
        formatSpec = '%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID0p5TAFM, formatSpec, scans(i).meanH-0.5, scans(i).meanK-0.5, scans(i).meanL-0.5, 1, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType == 5 && scans(i).include
        formatSpec = '%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID0p5TFM, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    elseif scans(i).scanType == 6 && scans(i).include
        formatSpec = '%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        fprintf(fileID3T, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
    end
end
fclose(fileID5K);
fclose(fileID1p5K);
fclose(fileID4p5K);
fclose(fileID0p5TAFM);
fclose(fileID0p5TFM);
fclose(fileID3T);
