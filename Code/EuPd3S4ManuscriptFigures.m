%% Figures for EuPd3S4. Vincent Morano, 01/29/2021
%% HC and Neutron Figure

% Field scan of (1/2,1/2,1/2)
clear
close all

set(0, 'defaultaxesfontsize', 7)

scans=struct;
scans.filename=strings(5,1);
scans.formatSpec=strings(5,1);
scans.filename(1)='FieldScan0to1_847390.bt7';
scans.filename(2)='FieldScan0p5to3_847488.bt7';
scans.filename(3)='FieldScan3to0p35_847522.bt7';
scans.filename(4)='FieldScan0p35to0p15_847523.bt7';
scans.filename(5)='FieldScan0p15to0_847524.bt7';
scans.formatSpec(1:5)=repmat([repmat('%s ',[1,4]), repmat('%f ',[1,36]),'%s %s','%f %f %f','%s ',repmat('%f ',[1,17]),repmat('%s ',[1,4]),repmat('%f ',[1,33])], 5, 1);
scans.detector=cell(5,1);
scans.detectorError=cell(5,1);
scans.monitor=cell(5,1);
scans.monitorError=cell(5,1);
scans.field=cell(5,1);
scans.intensity=cell(5,1);
scans.intensityError=cell(5,1);
scans.temperature=cell(5,1);

for i=1:length(scans.filename)
    fileID=fopen(strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\EuPd3S4\All\', scans.filename(i)));
    fileData=textscan(fileID, scans.formatSpec(i), 'CommentStyle', '#');
    fclose(fileID);
    scans.detector{i}=fileData{9};
    scans.detectorError{i}=sqrt(fileData{9});
    scans.time{i}=fileData{6};
    scans.monitor{i}=fileData{8};
    scans.monitorError{i}=sqrt(fileData{8});
    scans.field{i}=fileData{49};
    scans.intensity{i}=scans.detector{i}./scans.time{i};
    scans.intensityError{i}=sqrt(scans.detector{i})./scans.time{i};
    scans.temperature{i}=fileData{7};
end

decreasingIntensity=[scans.intensity{3}; scans.intensity{4}; scans.intensity{5}];
decreasingIntensityError=[scans.intensityError{3}; scans.intensityError{4}; scans.intensityError{5}];
decreasingField=[scans.field{3}; scans.field{4}; scans.field{5}];

figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.7]) % Left, bottom, width, height
tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile(4)
hold on
text(0.9, 0.9, '\bf{d}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.7, 0.7, '$$\left( \frac{1}{2} \frac{1}{2} \frac{1}{2} \right)$$', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'latex')
text(0.2, 0.5, '\it{T} \rm{= 0.3 K}', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('\mu_0\it{H} \rm{(T)}')
ylabel('\it{I} \rm{(cts / sec.)}')
errorbar(decreasingField, decreasingIntensity, decreasingIntensityError, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
xline(0.275, '--k')
xline(2.6, '--k')
%pbaspect([16 9 1])
axis square
xlim([0 4])
xticks([0 1 2 3 4])
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
box on
hold off


% Add heat capacity field scans. 11/14 data.
% First import the data

factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole

fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
fileName1 = '0p19mg_drhc_EuPd2S4_11132020.dat';
fileLoc2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileName2 = '0p19mg_drhc_EuPd2S4_09182020_fieldsweeps.dat';
file1 = readcell(strcat(fileLoc1, fileName1), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
indF1 = [103, 3130];
file2 = readcell(strcat(fileLoc2, fileName2), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
indF2 = [753, 1677]; % Only going up to 1 K because the field scans become much shorter after that, these are in the supplement.

[temp1, specHeat1, specHeatErr1, field1] = importHC(file1, mass, molarMass, indF1);
[temp2, specHeat2, specHeatErr2, field2] = importHC(file2, mass, molarMass, indF2);

% Then build up the figure

figure(1)
nexttile(2)
hold on
text(0.9, 0.9, '\bf{c}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
clim([0.3, 1])
cmap = jet;
c = colorbar('northoutside', 'Ticks', 0.4:0.2:1, 'TickLength', 0.02);
colormap(jet)
ylabel(c, '\it{T} \rm{(K)}')

unTemps1 = unique(round(temp1, 1));
unTemps2 = unique(round(temp2, 1));
unTemps = [unTemps1; unTemps2];

% Map every field to its corresponding jet cover where jet has 256 rows of rgb values.
slope1 = (length(cmap) - 1)/(max(unTemps)-min(unTemps));
indCmap1 = round(slope1*(unTemps1 - min(unTemps)) + 1);
rgb1 = cmap(indCmap1, :);
slope2 = (length(cmap) - 1)/(max(unTemps)-min(unTemps));
indCmap2 = round(slope2*(unTemps2 - min(unTemps)) + 1);
rgb2 = cmap(indCmap2, :);

for i = 1:length(unTemps1)
    ind = round(temp1, 1)==round(unTemps1(i), 1); % Temperature value
    
    if any(ind)
        tempInd = temp1(ind);
        fieldIndPre = field1(ind);
        specHeatIndPre = specHeat1(ind);
        specHeatErrIndPre = specHeatErr1(ind);
        
        % Sort data from low-B to high-B
        [fieldInd, bInd] = sort(fieldIndPre);
        specHeatInd = specHeatIndPre(bInd);
        specHeatErrInd = specHeatErrIndPre(bInd);

        % Average sets of hc points
        pos = 1;
        k = 1;
        fieldAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(fieldInd)
            indAvg = abs(fieldInd - fieldInd(k)) < 0.00001;
            fieldAvg(pos) = mean(fieldInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        plot(fieldAvg, specHeatAvg, 'LineWidth', 1.5, 'Color', rgb1(i,:))
    end
end
for i = 1:length(unTemps2)
    ind = round(temp2, 1)==round(unTemps2(i), 1); % Temperature value
    
    if any(ind)
        tempInd = temp2(ind);
        fieldIndPre = field2(ind);
        specHeatIndPre = specHeat2(ind);
        specHeatErrIndPre = specHeatErr2(ind);
        
        % Sort data from low-B to high-B
        [fieldInd, bInd] = sort(fieldIndPre);
        specHeatInd = specHeatIndPre(bInd);
        specHeatErrInd = specHeatErrIndPre(bInd);

        % Average sets of hc points
        pos = 1;
        k = 1;
        fieldAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(fieldInd)
            indAvg = abs(fieldInd - fieldInd(k)) < 0.00001;
            fieldAvg(pos) = mean(fieldInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        plot(fieldAvg, specHeatAvg, 'LineWidth', 1.5, 'Color', rgb2(i,:))
    end
end

ylabel('\it{C} \rm{(J / mol-K)}')
%pbaspect([16 9 1])
axis square
xlim([0 4])
ylim([0 3.5])
xticks([0 1 2 3 4])
yticks([0 1 2 3])
xticklabels([])
xline(0.275, '--k')
xline(2.6, '--k')
box on
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
hold off


% Order Parameter Measurement
% Import data

clear

betaMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\betaColorPlot.mat';
load(betaMat)

nexttile(3)
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.71, 0.7, '$$\left( \frac{1}{2} \frac{1}{2} \frac{1}{2} \right)$$', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'latex')
text(0.1, 0.5, '\mu_0\it{H} \rm{= 0 T}', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.1, 0.1, '\it{T}\rm{_N = 2.89(1) K}', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.1, 0.2, '\beta \rm{= 0.45(4)}', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('\it{T} \rm{(K)}')
ylabel('\it{I} \rm{(cts / sec.)}')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
tempPlotPrelim = linspace(temp(1), temp(end), 400);
yPlotPrelim = bgrBest.*ones(size(tempPlotPrelim)) + (I0Best.*(1-tempPlotPrelim./TcBest).^(2.*betaBest)).*(tempPlotPrelim<=TcBest);
tempPlot = linspace(tempFit(1), tempFit(end), 400);
ind3 = tempPlotPrelim>=min(tempFit);
plot(tempPlotPrelim(ind3), yPlotPrelim(ind3), 'r', 'LineWidth', 1.0)
plot(tempPlotPrelim(~ind3), yPlotPrelim(~ind3), '--r', 'LineWidth', 1.0)
%xline(min(tempFit), '--', {'Lower Fitting Limit'})
xline(min(TcBest), '--')
xlim([0 4])
xticks([0 1 2 3 4])
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
box on
%pbaspect([16 9 1])
axis square
hold off

disp(' ')
disp('Optimal Parameters')
disp(['Background: ', num2str(bgrBest)])
disp(['Low-T Intensity: ', num2str(I0Best)])
disp(['T_c: ', num2str(TcBest,5), '+-', num2str((max(TcL) - min(TcL))/2, 2)])
disp(['Beta: ', num2str(betaBest,3), '+-', num2str((max(betaL) - min(betaL))/2, 3)])
disp(['Chi_r^2: ', num2str(chiBest)])


% Add heat capacity above the corresponding order parameter
% First import the data

factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole

% file1 is 9/15 set, file2 9/15 0T, file3 11/14 set
fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileLoc2 = fileLoc1;
fileLoc3 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
fileName1 = '0p19mg_drhc_EuPd2S4_09152020_0T_fieldscan.dat';
fileName2 = '0p19mg_drhc_EuPd2S4_09152020_0T1T3T5T.dat';
fileName3 = '0p19mg_drhc_EuPd2S4_11132020.dat';
file1 = readcell(strcat(fileLoc1, fileName1), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
file2 = readcell(strcat(fileLoc2, fileName2), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
file3 = readcell(strcat(fileLoc3, fileName3), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
ind1 = [1, length(file1)];
ind2 = [1, length(file2)];
ind3 = [3206, 4912]; % 11/14/2020 indices. 5-79 is quick B = 0 T scan. Field scans start at 103. Manual low-T scans end at 3130. High-T T-scans begin at 3206.

[temp1, specHeat1, specHeatErr1, field1] = importHC(file1, mass, molarMass, ind1);
[temp2, specHeat2, specHeatErr2, field2] = importHC(file2, mass, molarMass, ind2);
[temp3, specHeat3, specHeatErr3, field3] = importHC(file3, mass, molarMass, ind3);

temp = [temp1; temp2; temp3];
specHeat = [specHeat1; specHeat2; specHeat3];
specHeatErr = [specHeatErr1; specHeatErr2; specHeatErr3];
field = [field1; field2; field3];

% Then build up the figure

figure(1)
nexttile(1)
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
clim([0, 3])
cmap = jet;
c = colorbar('northoutside', 'Ticks', 0:0.5:3, 'TickLength', 0.02);
colormap(jet)
ylabel(c, '\mu_0\it{H} \rm{(T)}')

unFields = unique(round(field, 3));

% Map every field to its corresponding jet cover where jet has 256 rows of rgb values.
slope = (length(cmap) - 1)/max(unFields);
indCmap = round(slope*unFields + 1);
rgb = cmap(indCmap, :);

for i = 1:length(unFields)
    ind = round(field1, 3)==round(unFields(i), 3); % Field value
    
    if any(ind)
        fieldInd = field1(ind);
        tempIndPre = temp1(ind);
        specHeatIndPre = specHeat1(ind);
        specHeatErrIndPre = specHeatErr1(ind);
        
        % Sort data from low-T to high-T
        [tempInd, tInd] = sort(tempIndPre);
        specHeatInd = specHeatIndPre(tInd);
        specHeatErrInd = specHeatErrIndPre(tInd);

        % Average sets of hc points
        pos = 1;
        k = 1;
        tempAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(tempInd)
            indAvg = abs(tempInd - tempInd(k)) < 0.01;
            tempAvg(pos) = mean(tempInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        plot(tempAvg, specHeatAvg, 'LineWidth', 1.5, 'Color', rgb(i,:))
    end
end
for i = 1:length(unFields)
    ind = round(field2, 3)==round(unFields(i), 3); % Field value
    
    if any(ind)
        fieldInd = field2(ind);
        tempIndPre = temp2(ind);
        specHeatIndPre = specHeat2(ind);
        specHeatErrIndPre = specHeatErr2(ind);
        
        % Sort data from low-T to high-T
        [tempInd, tInd] = sort(tempIndPre);
        specHeatInd = specHeatIndPre(tInd);
        specHeatErrInd = specHeatErrIndPre(tInd);

        % Average sets of hc points
        pos = 1;
        k = 1;
        tempAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(tempInd)
            indAvg = abs(tempInd - tempInd(k)) < 0.01;
            tempAvg(pos) = mean(tempInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        plot(tempAvg, specHeatAvg, 'LineWidth', 1.5, 'Color', rgb(i,:))
    end
end
for i = 1:length(unFields)
    ind = round(field3, 3)==round(unFields(i), 3); %% Field value
    
    if any(ind)
        fieldInd = field3(ind);
        tempIndPre = temp3(ind);
        specHeatIndPre = specHeat3(ind);
        specHeatErrIndPre = specHeatErr3(ind);

        % Sort data from low-T to high-T
        [tempInd, tInd] = sort(tempIndPre);
        specHeatInd = specHeatIndPre(tInd);
        specHeatErrInd = specHeatErrIndPre(tInd);
        
        % Average sets of hc points
        pos = 1;
        k = 1;
        tempAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(tempInd)
            indAvg = abs(tempInd - tempInd(k)) < 0.01;
            tempAvg(pos) = mean(tempInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        plot(tempAvg, specHeatAvg, 'LineWidth', 1.5, 'Color', rgb(i,:))
    end
end

ylabel('\it{C} \rm{(J / mol-K)}')
%pbaspect([16 9 1])
axis square
xlim([0 4])
ylim([0 11])
xticks([0 1 2 3 4])
yticks([0 2 4 6 8 10])
xticklabels([])
xline(min(TcBest), '--')
box on
set(gca, 'TickLength', [0.02, 0.01])
hold off

% Export figure. This is to avoid the fontsize changing when saving manually
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\HCFig.png';
exportgraphics(gcf, fdir, 'resolution', 450)

%% Peaks and phase diagram

clear
close all

set(0, 'defaultaxesfontsize', 7)

expt1Factor = 1.600; % Determined by the ratio of the paramagnetic 112 peak intensities

fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\EuPd3S4\All\';
fileNum = [845721, 845614, 847555, 847424, 847493, 845662]; % 5 K 0 T; 1.5 K 0 T; 4.7 K 0 T; 0.3 K 0.5 T; 0.4 K 3 T; 1.6 K 0 T
fileNames = strcat('fpx', string(fileNum), '.bt7');
headerLines = 41;

scans = importDataBT7(fileLoc, fileNames, headerLines);

% Convert from a3(deg) do delta Q
a = 6.67570;
for i=1:length(fileNames)
    scans(i).Q = 2*pi*sqrt(scans(i).meanH^2 + scans(i).meanK^2 + scans(i).meanL^2)/a; % aStar=2pi/a then Q=aStar*sqrt(H^2+..)
    scans(i).a3Rad = scans(i).a3*(pi/180);
    scans(i).meanA3Rad = scans(i).meanA3*(pi/180);
    scans(i).deltaQ = scans(i).Q.*(scans(i).a3Rad - scans(i).meanA3Rad);
end

% Fit to scans
errPts = 1e2; % For a3/a4 sigma, center fits
fact = [0.2, 0.3, 0.6, 0]; % bg, area, sigma, center
offset = [0.01, 0.6, 1e-3, 0.1];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a3 - x(4))./x(3)).^2./2); % bg, area, sigma, center
for i = 1:length(fileNames)
    modelInput = @(x) model(x, scans(i).a3);
    tmp = sort(scans(i).intTime); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intTime - bg0))./sum(scans(i).intTime - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(scans(i).intTime) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0 = [bg0, area0, sig0, cent0];
    [scans(i).xFit, scans(i).redChi2Fit, scans(i).xErr, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(scans(i).intTime, scans(i).intTimeErr, modelInput, x0, errPts, fact, offset);

    % Plot a3 fit
    scans(i).a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 5e2)';
    scans(i).intCal = model(scans(i).xFit, scans(i).a3Cal);
    scans(i).a3RadCal = scans(i).a3Cal*(pi/180);
    scans(i).deltaQCal = scans(i).Q.*(scans(i).a3RadCal - scans(i).meanA3Rad);
end

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile(1)
hold on
e4p7K = errorbar(scans(3).deltaQ*100, scans(3).intTime, scans(3).intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'MarkerSize', 4);
e0p3K0p5T = errorbar(scans(4).deltaQ*100, scans(4).intTime, scans(4).intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
e0p4K3T = errorbar(scans(5).deltaQ*100, scans(5).intTime, scans(5).intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
e1p6K = errorbar(scans(6).deltaQ*100, scans(6).intTime*expt1Factor, scans(6).intTimeErr*expt1Factor, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'g', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'g', 'MarkerSize', 4);
p4p7K = plot(scans(3).deltaQCal*100, scans(3).intCal, 'LineWidth', 1, 'Color', 'r');
p0p3K0p5T = plot(scans(4).deltaQCal*100, scans(4).intCal, 'LineWidth', 1, 'Color', 'b');
p0p4K3T = plot(scans(5).deltaQCal*100, scans(5).intCal, 'LineWidth', 1, 'Color', 'k');
p1p6K = plot(scans(6).deltaQCal*100, scans(6).intCal*expt1Factor, 'LineWidth', 1, 'Color', 'g');
ylabel('\it{I} \rm{(cts / sec.)}')
text(0.1, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.1, 0.75, '(001)', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'latex')
text(0.58, 0.7, '4.7 K, 0 T', 'Units', 'normalized', 'fontsize', 6, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'r')
text(0.58, 0.8, '0.3 K, 0.5 T', 'Units', 'normalized', 'fontsize', 6, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'b')
text(0.58, 0.9, '0.4 K, 3 T', 'Units', 'normalized', 'fontsize', 6, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
text(0.58, 0.6, '1.6 K, 0 T', 'Units', 'normalized', 'fontsize', 6, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'g')
xlim([-2 2])
ylim([0 50])
xticks((-2:1:2))
yticks([0 10 20 30 40 50])
xticklabels([])
set(gca, 'TickLength', [0.04, 0.01])
set(gca, 'FontName', 'Arial')
%pbaspect([16 9 1])
axis square
box on
hold off

nexttile(5)
hold on
e5K = errorbar(scans(1).deltaQ*100, scans(1).intTime, scans(1).intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'MarkerSize', 4);
e1p5K = errorbar(scans(2).deltaQ*100, scans(2).intTime, scans(2).intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
p5K = plot(scans(1).deltaQCal*100, scans(1).intCal, 'LineWidth', 1, 'Color', 'r');
p1p5K = plot(scans(2).deltaQCal*100, scans(2).intCal, 'LineWidth', 1, 'Color', 'b');
xlabel('$\delta \textsf{Q}_{\perp} ~(10^{-2} ~$\AA$^{-1})$', 'interpreter', 'latex')
ylabel('\it{I} \rm{(cts / sec.)}')
text(0.1, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.1, 0.75, '$\left( \frac{1}{2} \frac{1}{2} \frac{1}{2} \right)$', 'Units', 'normalized', 'fontsize', 7, 'interpreter', 'latex')
text(0.60, 0.22, '5.1 K, 0 T', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'r')
text(0.58, 0.70, '1.5 K, 0 T', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'b')
xlim([-2 2])
ylim([0 70])
xticks((-2:1:2))
yticks([0 20 40 60])
set(gca, 'TickLength', [0.04, 0.01])
set(gca, 'FontName', 'Arial')
%pbaspect([16 9 1])
axis square
box on
hold off

% Unit cell
nexttile([2,1])
text(0.0, 0.9, '\bf{c}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.0, 0.8, '\bf{d}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
text(0.0, 0.7, '\bf{e}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
set(gca, 'visible', 'off')

% Phase diagram
% Record temperature and field coordinates for phase transitions
clear

fieldHCTS=[0; 0.1; 0.233; 0.367; 0.5; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8; 2; 2.1; 2.2];
temperatureHCTS=[2.8641; 2.8272; 2.8344; 2.8284; 2.7651; 2.6984; 2.5631; 2.4297; 2.3; 2.1044; 1.9733; 1.7737; 1.5751; 1.4435; 1.377];
fieldHCFS=[0.2872; 0.33846; 0.33334; 0.35; 0.38889; 0.4; 0.38889; 0.4; 0.38889; 0.37975; 2.30001; 2.4; 2.45001; 2.5; 2.50001; 2.55; 2.60001; 2.55];
temperatureHCFS=[1.2; 2.8; 1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3];
%fieldHCFS=[0.33846; 0.07692; 0.10769; 0.12308; 0.13846; 0.15385; 0.18461; 0.20000; 0.22758; 0.26667; 0.28718; 0.30769; 0.33334; 0.35; 0.38889; 0.4; 0.38889; 0.4; 0.38889; 0.37975; 2.30001; 2.4; 2.45001; 2.5; 2.50001; 2.55; 2.60001; 2.55];
%temperatureHCFS=[2.8; 2.6; 2.5; 2.4; 2.3; 2.2; 2.1; 2.0; 1.8; 1.6; 1.5; 1.2; 1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3];
fieldNFS=[0; 0.27494; 2.6; 0.75];
temperatureNFS=[2.886; 0.3147; 0.3777; 2.5618];
%fieldNFS=[0];
%temperatureNFS=[2.886];

% Arrays for drawing phase boundaries
tempPhBLow = [0.3; 0.5; 0.75; 1.2]; % Estimate of points on the lower phase boundary
fieldPhBLow = [0.38; 0.4; 0.4; 0.29]; % Estimate of points on the lower phase boundary
tempPhBHigh = [0.3; 1; 2; 2.25; 2.5; 2.75; 2.88]; % Estimate of points on the upper phase boundary
fieldPhBHigh = [2.6; 2.38; 1.6; 1.3; 0.93; 0.50; 0.0]; % Estimate of points on the upper phase boundary
tempPhBLowInt = linspace(min(tempPhBLow), max(tempPhBLow), 2e2);
tempPhBHighInt = linspace(min(tempPhBHigh), max(tempPhBHigh), 2e2);
fieldPhBLowInt = interp1(tempPhBLow, fieldPhBLow, tempPhBLowInt, 'spline'); % Change estimate of discrete points to a curve by interpolation
fieldPhBHighInt = interp1(tempPhBHigh, fieldPhBHigh, tempPhBHighInt, 'spline'); % Change estimate of discrete points to a curve by interpolation

% Plot points over color plot for phase diagram
factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole

fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileLoc2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
fileName1 = '0p19mg_drhc_EuPd2S4_09182020_fieldsweeps.dat';
fileName2 = '0p19mg_drhc_EuPd2S4_09152020_0T_fieldscan.dat';
fileName3 = '0p19mg_drhc_EuPd2S4_09152020_0T1T3T5T.dat';
fileName4 = '0p19mg_drhc_EuPd2S4_11132020.dat';
file1 = readcell(strcat(fileLoc1, fileName1), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
file2 = readcell(strcat(fileLoc1, fileName2), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
file3 = readcell(strcat(fileLoc1, fileName3), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
file4 = readcell(strcat(fileLoc2, fileName4), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');

tmp = file1(:,8);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
temp1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file1(:,10);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeat1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from J/g-K to J/mol-K. Convert from old to new mass.
tmp = file1(:,6);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
field1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

tmp = file2(:,8);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
temp2 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file2(:,10);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeat2 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from muJ/K to J/mol-K
tmp = file2(:,6);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
field2 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

tmp = file3(:,8);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
temp3 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file3(:,10);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeat3 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from muJ/K to J/mol-K
tmp = file3(:,6);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
field3 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

tmp = file4(:,8);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
temp4 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file4(:,10);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeat4 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from muJ/K to J/mol-K
tmp = file4(:,6);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
field4 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

temp = [temp1; temp2; temp3; temp4];
specHeat = [specHeat1; specHeat2; specHeat3; specHeat4];
field = [field1; field2; field3; field4];

% Arrange the arrays on a grid with interpolation
[TGrid, BGrid] = meshgrid(0:0.05:4, 0:0.05:4);
CoTGrid = griddata(temp, field, specHeat./temp, TGrid, BGrid); % This interpolates to TGrid and BGrid values, passing through all points given by temperature, field, specific heat.
CoTGrid(TGrid<0.29) = NaN; % Don't interpolate lower than 0.3 K because I don't have much data there.

% Calculate the magnitude of a generalized gradient
% [dCdT, dCdB] = gradient(CGrid, TGrid(1,:), BGrid(:,1));
% grad = sqrt(dCdT.^2 + dCdB.^2);
%grad = dCdB;

% Plot the magnitude of the heat capacity
figure(1)
%gca.Layout.Tile = 3;
nexttile([2,2])
hold on
p1 = pcolor(TGrid, BGrid, CoTGrid);
set(p1, 'EdgeColor', 'none');
xlabel('\it{T} \rm{(K)}')
ylabel('\mu_0\it{H} \rm{(T)}')
c = colorbar;
colormap(jet)
clim([0, 5])
ylabel(c, '\it{C/T} \rm{(J / mol-K^{2})}')
s1 = scatter([temperatureHCTS; temperatureHCFS], [fieldHCTS; fieldHCFS], 20, 'm', '^', 'filled');
s2 = scatter(temperatureNFS, fieldNFS, 20, 'k', 'square', 'filled');
line(tempPhBLowInt, fieldPhBLowInt, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1) % Sketch lower phase boundary
line(tempPhBHighInt, fieldPhBHighInt, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1) % Sketch upper phase boundary
legend([s1, s2], {'Heat Capacity', 'Neutron Diffraction'}, 'Position', [0.705, 0.84, 0.16, 0.14])
legend('boxoff')
axis square
xlim([0 4.05])
ylim([0 4.05])
xticks([0 1 2 3 4])
yticks([0 1 2 3 4])
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
set(gca,'Layer','top')
box on
text(0.08, 0.95, '\bf{f}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'color', 'w', 'interpreter', 'tex')
text(0.2, 0.04, 'AFM', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
text(0.2, 0.25, 'Canted AFM', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
text(0.8, 0.2, 'PM', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'w')
hold off

% hold on
% p2 = pcolor(TGrid, BGrid, grad);
% set(p2, 'EdgeColor', 'none')
% xlabel('$T ~(K)$', 'fontsize', 10)
% ylabel('$B ~(T)$', 'fontsize', 10)
% c = colorbar;
% colormap(jet)
% clim([0, 10])
% ylabel(c, '$| \nabla_{T,B} C |$', 'fontsize', 10, 'interpreter', 'latex') % Note that the units are ill-defined
% s1 = scatter(temperatureHCTS, fieldHCTS, 60, 'm', 'x', 'LineWidth', 2);
% s2 = scatter(temperatureHCFS, fieldHCFS, 60, 'k', 'x', 'LineWidth', 2);
% s3 = scatter(temperatureNFS, fieldNFS, 60, 'g', '+', 'LineWidth', 2);
% legend([s1, s2, s3], {'HC T-Scans', 'HC B-Scans', 'Neutron Scans'}, 'fontsize', 10)
% axis square
% xlim([0 4])
% ylim([0 4])
% xticks([0 1 2 3 4])
% yticks([0 1 2 3 4])
% set(gca,'Layer','top')
% box on
% text(0.2, 0.025, 'AFM', 'FontSize', 12, 'Units', 'normalized', 'Color', 'k')
% text(0.2, 0.25, 'Canted AFM', 'FontSize', 12, 'Units', 'normalized', 'Color', 'k')
% text(0.15, 0.75, 'Field Polarized PM', 'FontSize', 12, 'Units', 'normalized', 'Color', 'k')
% text(0.8, 0.2, 'PM', 'FontSize', 12, 'Units', 'normalized', 'Color', 'k')
% hold off

% Export figure. This is to avoid the fontsize changing when saving manually
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\PhaseDiagramFigMATLAB.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% Plot heat capacity field scans at low-fields with high-T data. 9/15 data.
clear
close all

set(0, 'defaultaxesfontsize', 7)

% First import the data

factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole
axLim = 0.8; % xlim upper bound

fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileName = '0p19mg_drhc_EuPd2S4_09182020_fieldsweeps.dat';
%fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
%fileName = '0p19mg_drhc_EuPd2S4_11132020.dat';
file = readcell(strcat(fileLoc, fileName), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
indF = [1, length(file)];
%indF = [103, 3130]; % Field scans

[temp, specHeat, specHeatErr, field] = importHC(file, mass, molarMass, indF);

% Then build up the figure

unTemps = unique(round(temp, 1));
unTemps = unTemps(unTemps < 2.75);

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 4.2])
tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
hold on
clim([min(unTemps), max(unTemps)])
cmap = jet;
c = colorbar('northoutside', 'Ticks', linspace(0.5, 2.5, 6), 'TickLength', 0.02);
colormap(jet)
ylabel(c, '\it{T} \rm{(K)}')
hold off

% Map every field to its corresponding jet cover where jet has 256 rows of rgb values.
slope = (length(cmap) - 1)/(max(unTemps)-min(unTemps));
indCmap = round(slope*(unTemps - min(unTemps)) + 1);
rgb = cmap(indCmap, :);

for i = 1:length(unTemps)
    ind = round(temp, 1)==round(unTemps(i), 1); % Temperature value
    
    if any(ind)
        tempInd = temp(ind);
        fieldIndPre = field(ind);
        specHeatIndPre = specHeat(ind);
        specHeatErrIndPre = specHeatErr(ind);
        
        % Sort data from low-B to high-B
        [fieldInd, bInd] = sort(fieldIndPre);
        specHeatInd = specHeatIndPre(bInd);
        specHeatErrInd = specHeatErrIndPre(bInd);

        % Average sets of hc points
        pos = 1;
        k = 1;
        fieldAvg = [];
        specHeatAvg = [];
        specHeatErrAvg = [];
        while k <= length(fieldInd)
            indAvg = abs(fieldInd - fieldInd(k)) < 0.001;
            fieldAvg(pos) = mean(fieldInd(indAvg));
            specHeatAvg(pos) = mean(specHeatInd(indAvg));
            specHeatErrAvg(pos) = sqrt(sum(specHeatErrInd(indAvg).^2))/nnz(indAvg);
            k = k + nnz(indAvg);
            pos = pos + 1;
        end
        
        nexttile(1)
        hold on
        plot(fieldAvg, specHeatAvg, 'LineWidth', 2, 'Color', rgb(i,:))
        hold off

        nexttile(2)
        hold on
        plot(fieldAvg, (specHeatAvg-min(specHeatAvg(fieldAvg < axLim)))/max(specHeatAvg(fieldAvg < axLim)-min(specHeatAvg(fieldAvg < axLim))), 'LineWidth', 2, 'Color', rgb(i,:))
        hold off
    end
end

nexttile(1)
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
ylabel('\it{C} \rm{(J / mol-K)}')
pbaspect([16 9 1])
xlim([0 axLim])
xticks(0:0.2:axLim)
xticklabels([])
box on
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
hold off

nexttile(2)
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('\mu_0\it{H} \rm{(T)}')
ylabel('\it{C} \rm{(arb.)}')
pbaspect([16 9 1])
xlim([0 axLim])
ylim([0 1])
xticks(0:0.2:axLim)
box on
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\MetamagneticTransition.png';
%exportgraphics(gcf, fdir, 'resolution', 450)

%% Plot some representative hc pulses from the semi-adiabatic scans from our beta fits
clear
close all

set(0, 'defaultaxesfontsize', 7)

factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole

% First import the data: file1 is 9/15 T-scan at 0 T, file2 is 11/14 T-scans at a few fields
fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileLoc2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
fileName1 = '0p19mg_drhc_EuPd2S4_09152020_0T1T3T5T.raw';
fileName2 = '0p19mg_drhc_EuPd2S4_11132020.raw';
dir1 = strcat(fileLoc1, fileName1);
dir2 = strcat(fileLoc2, fileName2);
indPulseParams1 = [30104, 30152];
indData1 = [30154, 30409];
indPulseParams2 = [606702, 606720];
indData2 = [606722, 606977];

optsPulseParams1 = delimitedTextImportOptions('DataLines', indPulseParams1, 'VariableNamesLine', 10, 'Delimiter', ',');
optsData1 = delimitedTextImportOptions('NumVariables', 7, 'DataLines', indData1, 'VariableNamesLine', 10, 'Delimiter', ',', 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'});
pulseParams1 = readtable(dir1, optsPulseParams1);
data1 = readtable(dir1, optsData1);
optsPulseParams2 = delimitedTextImportOptions('DataLines', indPulseParams2, 'VariableNamesLine', 10, 'Delimiter', ',');
optsData2 = delimitedTextImportOptions('NumVariables', 7, 'DataLines', indData2, 'VariableNamesLine', 10, 'Delimiter', ',', 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'});
pulseParams2 = readtable(dir2, optsPulseParams2);
data2 = readtable(dir2, optsData2);

time1 = data1.TimeStamp_Seconds_ - data1.TimeStamp_Seconds_(1);
platformTemp1 = data1.PlatformTemp_Kelvin_;
platformTempFit1 = data1.PlatformTempFit_Kelvin_;
sampleTempFit1 = data1.SampleTempFit_Kelvin_;
time2 = data2.TimeStamp_Seconds_ - data2.TimeStamp_Seconds_(1);
platformTemp2 = data2.PlatformTemp_Kelvin_;
platformTempFit2 = data2.PlatformTempFit_Kelvin_;
sampleTempFit2 = data2.SampleTempFit_Kelvin_;

% Plot the hc pulses
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.7])
%tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
%nexttile
hold on
xlabel('Time (s)')
ylabel('\it{T} \rm{(K)}')
s1 = scatter(time1, platformTemp1, 15, 'b');
%p1 = plot(time1, platformTempFit1, 'Color', 'b', 'LineWidth', 2);
p2 = plot(time1, sampleTempFit1, 'Color', 'r', 'LineWidth', 2);
%text(0.1, 0.9, 'Experiment 1', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'Color', 'k')
%legend([p1, p2], {'Platform', 'Sample'})
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
xlim([0, max([time1; time2])])
% xticklabels([])
pbaspect([16 9 1])
box on
hold off

% nexttile
% hold on
% xlabel('Time (sec.)')
% ylabel('\it{T} \rm{(K)}')
% s2 = scatter(time2, platformTemp2, 15, 'b');
%p3 = plot(time2, platformTempFit2, 'Color', 'b', 'LineWidth', 2);
% p4 = plot(time2, sampleTempFit2, 'Color', 'r', 'LineWidth', 2);
%text(0.1, 0.9, 'Experiment 2', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'Color', 'k')
%legend([p3, p4], {'Platform', 'Sample'})
% set(gca, 'TickLength', [0.02, 0.01])
% set(gca, 'FontName', 'Arial')
% xlim([0, max([time1; time2])])
% pbaspect([16 9 1])
% box on
% hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\HeatPulse.png';
exportgraphics(gcf, fdir, 'resolution', 450)

%% Supplementary beta figure

clear
close all

betaMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\betaColorPlot.mat';
load(betaMat)

figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 5.9])
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
hold on
xlabel('\it{T}\rm{_{min} (K)}')
ylabel('\beta')
text(0.05, 0.1, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(temp(temp>fitMin*estTN & temp<fitMax*estTN), xFit(4,:), xErr(4,:), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
box on
pbaspect([16 9 1])
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
hold off

nexttile
hold on
s = surf(TcGrid, betaGrid, redChi2Grid);
s.EdgeColor = 'none'; % Don't show individual pixels
c = colorbar('north', 'Ticks', 1:0.1:1.5, 'TickLength', 0.02, 'Color', 'w');
clim([1, 1.5])
%ylabel(c, '\chi_r^2', 'fontsize', 9)
text(0.5, 0.65, '\chi_r^2', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'Color', 'w', 'HorizontalAlignment', 'center')
xlabel('\it{T}\rm{_N (K)}')
ylabel('\beta')
text(0.05, 0.1, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex', 'color', 'w')
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
set(gca, 'FontName', 'Arial')
hold off

nexttile
hold on
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{I}\rm{ (cts / sec.)}')
text(0.05, 0.1, '\bf{c}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
tempPlotPrelim = linspace(temp(1), temp(end), 400);
yPlotPrelim = bgrBest.*ones(size(tempPlotPrelim)) + (I0Best.*(1-tempPlotPrelim./TcBest).^(2.*betaBest)).*(tempPlotPrelim<=TcBest);
tempPlot = linspace(tempFit(1), tempFit(end), 400);
yPlot = bgrBest.*ones(size(tempPlot)) + (I0Best.*(1 - tempPlot./TcBest).^(2.*betaBest)).*(tempPlot<=TcBest);
ind3 = tempPlotPrelim>=min(tempFit);
plot(tempPlotPrelim(ind3), yPlotPrelim(ind3), 'r', 'LineWidth', 1.0)
plot(tempPlotPrelim(~ind3), yPlotPrelim(~ind3), '--r', 'LineWidth', 1.0)
%xline(min(tempFit), '--', {'Lower Fitting Limit'})
xline(min(TcBest), '--', {'\it{T}\rm{_N}'}, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle', 'fontsize', 7)
pbaspect([16 9 1])
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\betaColorPlot.png';
exportgraphics(gcf, fdir, 'resolution', 450)

%% Supplementary EIGER figure: Plot 001 and 002 integrated intensities as a function of temperature at 8 T

clear
close all

set(0, 'defaultaxesfontsize', 7)

% Fit for 001 a3 scan integrated intensities, combining identical scans
errPts = 1e2; % For a3/a4 sigma, center fits
fact = [0.2, 0.3, 0.6, 0]; % bg, area, sigma, center
offset = [0.01, 0.6, 1e-3, 0.1];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a3 - x(4))./x(3)).^2./2); % bg, area, sigma, center

% Import scans
% 001 T-scan filenames: 4108, 4109; 4112, 4113; 4118, 4119; 4123, 4124, 4125; 4130,
% 4131, 4132
% 002 T-scan filenames: 4110; 4114, 4115; 4120, 4121; 4126, 4127, 4128;
% 4133
factor=1e3; % To get reasonable intensities after monitor normalization
directory='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\PSI\EIGER\EuPd3S4\091621\';
filenames={{'eiger2021n004108.scn', 'eiger2021n004109.scn'}; {'eiger2021n004112.scn', 'eiger2021n004113.scn'}; {'eiger2021n004118.scn', 'eiger2021n004119.scn'}; {'eiger2021n004123.scn', 'eiger2021n004124.scn', 'eiger2021n004125.scn'}; {'eiger2021n004130.scn', 'eiger2021n004131.scn', 'eiger2021n004132.scn'}}; % Each cell is a temperature, each element of the cell is a scan at that temperature
for i=1:length(filenames) % Loop through every temperature
    for j=1:length(filenames{i}) % Loop through every scan at a given temperature
        data{i}{j}=readtable([directory, filenames{i}{j}], 'FileType', 'text', 'NumHeaderLines', 52);
    end
end

% Fit combined a3 scans to Gaussian peak
figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, length(data), 'TileSpacing', 'tight', 'Padding', 'tight')
for i=1:length(data)
    if length(filenames{i})==2 % Two scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS);
    elseif length(filenames{i})==3 % Three scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1+data{i}{3}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1 + data{i}{3}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS+data{i}{3}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS + data{i}{3}.CNTS);
    end
    a3{i}=data{i}{1}.A3_1; % Use first scan's a3 values since identical
    temperature(i) = mean(data{i}{1}.TS); % Use first scan's temperature values since identical
    counts{i}=rawCounts{i}./mon{i}.*factor; % Normalize to monitor with overall scale factor
    countsErr{i}=sqrt(rawCounts{i})./mon{i}.*factor;
    %countsErr{i}=(rawCounts{i}./mon{i}).*factor.*sqrt((rawCountsErr{i}./rawCounts{i}).^2 + (monErr{i}./mon{i}).^2);

    % Fit to scans    
    modelInput = @(x) model(x, a3{i});
    tmp = sort(counts{i}); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(a3{i}.*(counts{i} - bg0))./sum(counts{i} - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(counts{i}) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0 = [bg0, area0, sig0, cent0];
    [xFit{i}, redChi2Fit{i}, xErr{i}, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(counts{i}, countsErr{i}, modelInput, x0, errPts, fact, offset);
    area(i) = xFit{i}(2);
    areaErr(i) = xErr{i}(2);

    % Plot a3 fit
    a3Cal{i} = linspace(min(a3{i}), max(a3{i}), 5e2)';
    intCal{i} = model(xFit{i}, a3Cal{i});
    nexttile
    hold on
    errorbar(a3{i}, counts{i}, countsErr{i}, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
    plot(a3Cal{i}, intCal{i}, 'LineWidth', 1, 'Color', 'b');
    xlabel('\it{a3} \rm{(deg.)}')
    ylabel('\it{I} \rm{(arb.)}')
    set(gca, 'TickLength', [0.04, 0.01])
    set(gca, 'FontName', 'Arial')
    axis square
    box on
    hold off
end
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERTScanRocking001.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

% Fit to Brillouin curve
N=1/(297.6e-30); % Number of Eu atoms per unit volume
g=2; % g-factor
J=7/2; % Total angular momentum
brillouin=@(arg) (2.*J+1)./(2.*J).*coth((2.*J+1).*arg./(2.*J))-coth(arg./(2.*J))/(2.*J); % Brillouin function
muB=9.274e-24; % Bohr magneton
kB=1.38e-23; % Boltzmann constant
arg=@(fieldB, tempB) J.*g.*muB.*fieldB./(kB.*tempB); % Ratio of Zeeman energy to thermal energy
mag=@(ratioE) N.*g.*muB.*J.*brillouin(ratioE); % Magnetization
modelB = @(factB, fieldB, tempB) factB.*mag(arg(fieldB, tempB)).^2; % I'm fitting to the integrated intensity, which is proportional to square of magnetization.

field=8; % Field in T
modelBInput = @(factB) modelB(factB, field, temperature);
x0B = [1e-10]; % Using brackets since fitting function expects array.
errPtsB = 1e3;
factB = [0];
offsetB = [x0B*.95];
[xFitB, redChi2FitB, xErrB, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(area, areaErr, modelBInput, x0B, errPtsB, factB, offsetB);

% Plot 001 integrated intensities as a function of temperature at 8 T
temperatureB=linspace(0,300,1e3); % Temperatures to plot
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.25]) % Left, bottom, width, height
hold on
errorbar(temperature, area, areaErr, 'o', 'Color', 'b', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1.0);
plot(temperatureB, modelB(xFitB, field, temperatureB), 'Color', 'b', 'LineWidth', 1.0)
hold off

% Plot 002 integrated intensities as a function of temperature at 8 T
clear

% Fit for 002 a3 scan integrated intensities, combining identical scans
errPts = 1e2; % For a3/a4 sigma, center fits
fact = [0.2, 0.3, 0.6, 0]; % bg, area, sigma, center
offset = [0.01, 0.6, 1e-3, 0.1];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a3 - x(4))./x(3)).^2./2); % bg, area, sigma, center

% Import scans
factor=1e3; % To get reasonable intensities after monitor normalization
directory='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\PSI\EIGER\EuPd3S4\091621\';
filenames={{'eiger2021n004110.scn'}; {'eiger2021n004114.scn', 'eiger2021n004115.scn'}; {'eiger2021n004120.scn', 'eiger2021n004121.scn'}; {'eiger2021n004126.scn', 'eiger2021n004127.scn', 'eiger2021n004128.scn'}; {'eiger2021n004133.scn'}}; % Each cell is a temperature, each element of the cell is a scan at that temperature
for i=1:length(filenames) % Loop through every temperature
    for j=1:length(filenames{i}) % Loop through every scan at a given temperature
        data{i}{j}=readtable([directory, filenames{i}{j}], 'FileType', 'text', 'NumHeaderLines', 52);
    end
end

% Fit combined a3 scans to Gaussian peak
figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8])
tiledlayout(1, length(data), 'TileSpacing', 'tight', 'Padding', 'tight')
for i=1:length(data)
    if length(filenames{i})==1 % One scan
        mon{i}=data{i}{1}.M1;
        monErr{i} = sqrt(data{i}{1}.M1);
        rawCounts{i}=data{i}{1}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS.^2);
    elseif length(filenames{i})==2 % Two scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS);
    elseif length(filenames{i})==3 % Three scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1+data{i}{3}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1 + data{i}{3}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS+data{i}{3}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS + data{i}{3}.CNTS);
    end
    a3{i}=data{i}{1}.A3_1; % Use first scan's a3 values since identical
    temperature(i) = mean(data{i}{1}.TS); % Use first scan's temperature values since identical
    counts{i}=rawCounts{i}./mon{i}.*factor; % Normalize to monitor with overall scale factor
    countsErr{i}=sqrt(rawCounts{i})./mon{i}.*factor;
    %countsErr{i}=(rawCounts{i}./mon{i}).*factor.*sqrt((rawCountsErr{i}./rawCounts{i}).^2 + (monErr{i}./mon{i}).^2);

    % Fit to scans    
    modelInput = @(x) model(x, a3{i});
    tmp = sort(counts{i}); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(a3{i}.*(counts{i} - bg0))./sum(counts{i} - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(counts{i}) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0 = [bg0, area0, sig0, cent0];
    [xFit{i}, redChi2Fit{i}, xErr{i}, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(counts{i}, countsErr{i}, modelInput, x0, errPts, fact, offset);
    area(i) = xFit{i}(2);
    areaErr(i) = xErr{i}(2);

    % Plot a3 fit
    a3Cal{i} = linspace(min(a3{i}), max(a3{i}), 5e2)';
    intCal{i} = model(xFit{i}, a3Cal{i});
    nexttile
    hold on
    errorbar(a3{i}, counts{i}, countsErr{i}, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
    plot(a3Cal{i}, intCal{i}, 'LineWidth', 1, 'Color', 'b');
    xlabel('\it{a3} \rm{(deg.)}')
    ylabel('\it{I} \rm{(arb.)}')
    set(gca, 'TickLength', [0.04, 0.01])
    set(gca, 'FontName', 'Arial')
    axis square
    box on
    hold off
end
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERTScanRocking002.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

% Fit to Brillouin curve
N=1/(297.6e-30); % Number of Eu atoms per unit volume
g=2; % g-factor
J=7/2; % Total angular momentum
brillouin=@(arg) (2.*J+1)./(2.*J).*coth((2.*J+1).*arg./(2.*J))-coth(arg./(2.*J))/(2.*J); % Brillouin function
muB=9.274e-24; % Bohr magneton
kB=1.38e-23; % Boltzmann constant
arg=@(fieldB, tempB) J.*g.*muB.*fieldB./(kB.*tempB); % Ratio of Zeeman energy to thermal energy
mag=@(ratioE) N.*g.*muB.*J.*brillouin(ratioE); % Magnetization
modelB = @(factB, fieldB, tempB) factB.*mag(arg(fieldB, tempB)).^2; % I'm fitting to the integrated intensity, which is proportional to square of magnetization.

field=8; % Field in T
modelBInput = @(factB) modelB(factB, field, temperature);
x0B = [1e-10]; % Using brackets since fitting function expects array.
errPtsB = 1e3;
factB = [0];
offsetB = [x0B*.95];
[xFitB, redChi2FitB, xErrB, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(area, areaErr, modelBInput, x0B, errPtsB, factB, offsetB);

% Plot 002 integrated intensities as a function of temperature at 8 T
temperatureB=linspace(0,300,1e3); % Temperatures to plot
figure(2)
hold on
xlabel('\it{T} \rm{(K)}')
ylabel('\it{I} \rm{(det \times deg. / 10^3 mon)}')
errorbar(temperature, area, areaErr, 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 1.0);
plot(temperatureB, modelB(xFitB, field, temperatureB), 'Color', 'r', 'LineWidth', 1.0)
text(0.2, 0.5, '(001)', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex', 'color', 'b')
text(0.05, 0.15, '(002)', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex', 'color', 'r')
text(0.5, 0.8, '\it{\mu}\rm{_0}\it{H}\rm{ = 8 T}', 'Units', 'normalized', 'fontsize', 12, 'fontname', 'Arial')
%legend('', '(001)', '', '(002)')
xlim([-inf, 100])
ylim([0,inf])
box on
hold off
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERTScanBrillouin.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% Supplementary EIGER figure: Plot 001 and 002 integrated intensities as a function of field at 0.1 K

clear
close all

set(0, 'defaultaxesfontsize', 7)

% Fit for 001 a3 scan integrated intensities, combining identical scans
fieldPM = 2.6; % Field in T beyond which the compound is in the field-polarized paramagnetic state and I can reasonably fit to the Brillouin function. Based on 0.3 K BT7 field scans and DRHC data.
errPts = 1e2; % For a3/a4 sigma, center fits
fact = [0.2, 0.3, 0.6, 0]; % bg, area, sigma, center
offset = [0.01, 0.6, 1e-3, 0.1];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a3 - x(4))./x(3)).^2./2); % bg, area, sigma, center

% Import scans
% 001 B-scan filenames: 4058; 4071, 4072, 4073; 4075, 4076, 4077; 4079,
% 4080, 4081; 4083, 4084; 4086, 4087; 4089, 4090; 4104, 4105; 4108, 4109
% 002 B-scan filenames: 4059; 4095, 4096; 4099, 4100; 4101, 4102; 4106;
% 4110
factor=1e3; % To get reasonable intensities after monitor normalization
directory='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\PSI\EIGER\EuPd3S4\091621\';
filenames={{'eiger2021n004058.scn'}; {'eiger2021n004071.scn', 'eiger2021n004072.scn', 'eiger2021n004073.scn'}; {'eiger2021n004075.scn', 'eiger2021n004076.scn', 'eiger2021n004077.scn'}; {'eiger2021n004079.scn', 'eiger2021n004080.scn', 'eiger2021n004081.scn'}; {'eiger2021n004083.scn', 'eiger2021n004084.scn'}; {'eiger2021n004086.scn', 'eiger2021n004087.scn'}; {'eiger2021n004089.scn', 'eiger2021n004090.scn'}; {'eiger2021n004104.scn', 'eiger2021n004105.scn'}; {'eiger2021n004108.scn', 'eiger2021n004109.scn'}}; % Each cell is a field, each element of the cell is a scan at that field
for i=1:length(filenames) % Loop through every temperature
    for j=1:length(filenames{i}) % Loop through every scan at a given temperature
        data{i}{j}=readtable([directory, filenames{i}{j}], 'FileType', 'text', 'NumHeaderLines', 52);
    end
end
data{1}{1}.MF=0.05.*ones(height(data{1}{1}),1); % Datafile doesn't have field for 4058, 4059 so this is what I recorded in the logbook

% Fit combined a3 scans to Gaussian peak
figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, length(data), 'TileSpacing', 'tight', 'Padding', 'tight')
for i=1:length(data)
    if length(filenames{i})==1 % One scan
        mon{i}=data{i}{1}.M1;
        monErr{i} = sqrt(data{i}{1}.M1);
        rawCounts{i}=data{i}{1}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS.^2);
    elseif length(filenames{i})==2 % Two scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS);
    elseif length(filenames{i})==3 % Three scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1+data{i}{3}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1 + data{i}{3}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS+data{i}{3}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS + data{i}{3}.CNTS);
    end
    a3{i}=data{i}{1}.A3_1; % Use first scan's a3 values since identical
    field(i) = mean(data{i}{1}.MF); % Use first scan's field values since identical
    counts{i}=rawCounts{i}./mon{i}.*factor; % Normalize to monitor with overall scale factor
    countsErr{i}=sqrt(rawCounts{i})./mon{i}.*factor;
    %countsErr{i}=(rawCounts{i}./mon{i}).*factor.*sqrt((rawCountsErr{i}./rawCounts{i}).^2 + (monErr{i}./mon{i}).^2);

    % Fit to scans    
    modelInput = @(x) model(x, a3{i});
    tmp = sort(counts{i}); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(a3{i}.*(counts{i} - bg0))./sum(counts{i} - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(counts{i}) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0 = [bg0, area0, sig0, cent0];
    [xFit{i}, redChi2Fit{i}, xErr{i}, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(counts{i}, countsErr{i}, modelInput, x0, errPts, fact, offset);
    area(i) = xFit{i}(2);
    areaErr(i) = xErr{i}(2);

    % Plot a3 fit
    a3Cal{i} = linspace(min(a3{i}), max(a3{i}), 5e2)';
    intCal{i} = model(xFit{i}, a3Cal{i});
    nexttile
    hold on
    errorbar(a3{i}, counts{i}, countsErr{i}, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
    plot(a3Cal{i}, intCal{i}, 'LineWidth', 1, 'Color', 'b');
    xlabel('\it{a3} \rm{(deg.)}')
    ylabel('\it{I} \rm{(arb.)}')
    set(gca, 'TickLength', [0.04, 0.01])
    set(gca, 'FontName', 'Arial')
    axis square
    box on
    hold off
end
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERBScanRocking001.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

% Fit to Brillouin curve
N=1/(297.6e-30); % Number of Eu atoms per unit volume
g=2; % g-factor
J=7/2; % Total angular momentum
brillouin=@(arg) (2.*J+1)./(2.*J).*coth((2.*J+1).*arg./(2.*J))-coth(arg./(2.*J))/(2.*J); % Brillouin function
muB=9.274e-24; % Bohr magneton
kB=1.38e-23; % Boltzmann constant
arg=@(fieldB, tempB) J.*g.*muB.*fieldB./(kB.*tempB); % Ratio of Zeeman energy to thermal energy
mag=@(ratioE) N.*g.*muB.*J.*brillouin(ratioE); % Magnetization
modelB = @(factB, fieldB, tempB) factB.*mag(arg(fieldB, tempB)).^2; % I'm fitting to the integrated intensity, which is proportional to square of magnetization.

temperature=0.1; % Temperature in K
modelBInput = @(factB) modelB(factB, field(field>fieldPM), temperature);
x0B = [1e-10]; % Using brackets since fitting function expects array.
errPtsB = 1e3;
factB = [0];
offsetB = [x0B*.95];
[xFitB, redChi2FitB, xErrB, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(area(field>fieldPM), areaErr(field>fieldPM), modelBInput, x0B, errPtsB, factB, offsetB);

% Plot 001 integrated intensities as a function of field at 0.1 K
fieldBDash=linspace(0,fieldPM,1e3); % Fields to plot beyond fitting range
fieldBSol=linspace(fieldPM,20,1e3); % Fields to plot within fitting range
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.25]) % Left, bottom, width, height
hold on
errorbar(field, area, areaErr, 'o', 'Color', 'b', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1.0);
plot(fieldBDash, modelB(xFitB, fieldBDash, temperature), 'Color', 'b', 'LineStyle', '--','LineWidth', 1.0)
plot(fieldBSol, modelB(xFitB, fieldBSol, temperature), 'Color', 'b', 'LineWidth', 1.0)
hold off

% Plot 002 integrated intensities as a function of field at 0.1 K
clear

% Fit for 002 a3 scan integrated intensities, combining identical scans
fieldPM = 2.7; % Field in T beyond which the compound is in the field-polarized paramagnetic state and I can reasonably fit to the Brillouin function. Based on 0.3 K BT7 field scans and DRHC data.
errPts = 1e2; % For a3/a4 sigma, center fits
fact = [0.2, 0.3, 0.6, 0]; % bg, area, sigma, center
offset = [0.01, 0.6, 1e-3, 0.1];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a3 - x(4))./x(3)).^2./2); % bg, area, sigma, center

% Import scans
factor=1e3; % To get reasonable intensities after monitor normalization
directory='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\PSI\EIGER\EuPd3S4\091621\';
filenames={{'eiger2021n004059.scn'}; {'eiger2021n004095.scn', 'eiger2021n004096.scn'}; {'eiger2021n004099.scn', 'eiger2021n004100.scn'}; {'eiger2021n004101.scn', 'eiger2021n004102.scn'}; {'eiger2021n004106.scn'}; {'eiger2021n004110.scn'}}; % Each cell is a field, each element of the cell is a scan at that field
for i=1:length(filenames) % Loop through every temperature
    for j=1:length(filenames{i}) % Loop through every scan at a given temperature
        data{i}{j}=readtable([directory, filenames{i}{j}], 'FileType', 'text', 'NumHeaderLines', 52);
    end
end
data{1}{1}.MF=0.05.*ones(height(data{1}{1}),1); % Datafile doesn't have field for 4058, 4059 so this is what I recorded in the logbook

% Fit combined a3 scans to Gaussian peak
figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8])
tiledlayout(1, length(data), 'TileSpacing', 'tight', 'Padding', 'tight')
for i=1:length(data)
    if length(filenames{i})==1 % One scan
        mon{i}=data{i}{1}.M1;
        monErr{i} = sqrt(data{i}{1}.M1);
        rawCounts{i}=data{i}{1}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS.^2);
    elseif length(filenames{i})==2 % Two scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS);
    elseif length(filenames{i})==3 % Three scans
        mon{i}=data{i}{1}.M1+data{i}{2}.M1+data{i}{3}.M1;
        monErr{i} = sqrt(data{i}{1}.M1 + data{i}{2}.M1 + data{i}{3}.M1);
        rawCounts{i}=data{i}{1}.CNTS+data{i}{2}.CNTS+data{i}{3}.CNTS;
        rawCountsErr{i} = sqrt(data{i}{1}.CNTS + data{i}{2}.CNTS + data{i}{3}.CNTS);
    end
    a3{i}=data{i}{1}.A3_1; % Use first scan's a3 values since identical
    field(i) = mean(data{i}{1}.MF); % Use first scan's field values since identical
    counts{i}=rawCounts{i}./mon{i}.*factor; % Normalize to monitor with overall scale factor
    countsErr{i}=sqrt(rawCounts{i})./mon{i}.*factor;
    %countsErr{i}=(rawCounts{i}./mon{i}).*factor.*sqrt((rawCountsErr{i}./rawCounts{i}).^2 + (monErr{i}./mon{i}).^2);

    % Fit to scans    
    modelInput = @(x) model(x, a3{i});
    tmp = sort(counts{i}); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(a3{i}.*(counts{i} - bg0))./sum(counts{i} - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(counts{i}) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0 = [bg0, area0, sig0, cent0];
    [xFit{i}, redChi2Fit{i}, xErr{i}, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(counts{i}, countsErr{i}, modelInput, x0, errPts, fact, offset);
    area(i) = xFit{i}(2);
    areaErr(i) = xErr{i}(2);

    % Plot a3 fit
    a3Cal{i} = linspace(min(a3{i}), max(a3{i}), 5e2)';
    intCal{i} = model(xFit{i}, a3Cal{i});
    nexttile
    hold on
    errorbar(a3{i}, counts{i}, countsErr{i}, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
    plot(a3Cal{i}, intCal{i}, 'LineWidth', 1, 'Color', 'b');
    xlabel('\it{a3} \rm{(deg.)}')
    ylabel('\it{I} \rm{(arb.)}')
    set(gca, 'TickLength', [0.04, 0.01])
    set(gca, 'FontName', 'Arial')
    axis square
    box on
    hold off
end
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERBScanRocking002.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

% Fit to Brillouin curve
N=1/(297.6e-30); % Number of Eu atoms per unit volume
g=2; % g-factor
J=7/2; % Total angular momentum
brillouin=@(arg) (2.*J+1)./(2.*J).*coth((2.*J+1).*arg./(2.*J))-coth(arg./(2.*J))/(2.*J); % Brillouin function
muB=9.274e-24; % Bohr magneton
kB=1.38e-23; % Boltzmann constant
arg=@(fieldB, tempB) J.*g.*muB.*fieldB./(kB.*tempB); % Ratio of Zeeman energy to thermal energy
mag=@(ratioE) N.*g.*muB.*J.*brillouin(ratioE); % Magnetization
modelB = @(factB, fieldB, tempB) factB.*mag(arg(fieldB, tempB)).^2; % I'm fitting to the integrated intensity, which is proportional to square of magnetization.

temperature=0.1; % Temperature in K
modelBInput = @(factB) modelB(factB, field(field>fieldPM), temperature);
x0B = [1e-10]; % Using brackets since fitting function expects array.
errPtsB = 1e3;
factB = [0];
offsetB = [x0B*.95];
[xFitB, redChi2FitB, xErrB, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(area(field>fieldPM), areaErr(field>fieldPM), modelBInput, x0B, errPtsB, factB, offsetB);

% Plot 002 integrated intensities as a function of field at 0.1 K
fieldBDash=linspace(0,fieldPM,1e3); % Fields to plot beyond fitting range
fieldBSol=linspace(fieldPM,20,1e3); % Fields to plot within fitting range
figure(2)
hold on
xlabel('\mu_0\it{H} \rm{(T)}')
ylabel('\it{I} \rm{(det \times deg. / 10^3 mon)}')
errorbar(field, area, areaErr, 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 1.0);
plot(fieldBDash, modelB(xFitB, fieldBDash, temperature), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0)
plot(fieldBSol, modelB(xFitB, fieldBSol, temperature), 'Color', 'r', 'LineWidth', 1.0)
xline(fieldPM, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.0)
text(0.8, 0.8, '(001)', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex', 'color', 'b')
text(0.8, 0.4, '(002)', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex', 'color', 'r')
text(0.5, 0.2, '\it{T}\rm{ = 0.1 K}', 'Units', 'normalized', 'fontsize', 12, 'fontname', 'Arial')
%legend('', '', '(001)', '', '', '(002)', 'Location', 'best')
%legend('(001)', '(002)', 'Location', 'best')
xlim([-inf, 10])
ylim([0,inf])
box on
hold off
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\EIGERBScanBrillouin.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% SI: (001) correlation length

clear
close all

set(0, 'defaultaxesfontsize', 7)

% Fit to the (001) reflection with Voigt whose Gaussian width is fixed to
% the average value for prominent nuclear peaks, then pick out Lorentzian
% hwhm for correlation length. Normalize to time.

% Import prominent 4.7 K, 0 T, experiment with magnet peaks. 112, 220, 222
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\EuPd3S4\All\';
fileNum = [847558, 847559, 847560];
fileNames = strcat('fpx', string(fileNum), '.bt7');
headerLines = 41;

scans = importDataBT7(fileLoc, fileNames, headerLines);

a = 6.67570;
aStar = 2*pi/a;
for i=1:length(fileNames)
    scans(i).Q = 2*pi*sqrt(scans(i).meanH^2 + scans(i).meanK^2 + scans(i).meanL^2)/a; % aStar=2pi/a then Q=aStar*sqrt(H^2+..)
end

% Fit nuclear peaks to a Gaussian and flat background, solving for the
% Gaussian sigma
errPts = 1e2; % For a3 sigma fits
fact = [0.2, 0.3, 0, 0.3]; % bg, area, center, sigma
offset = [0.01, 0.6, 0.05, 1e-3];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(4).*exp(-((a3 - x(3))./x(4)).^2./2); % bg, area, center, sigma
for i = 1:length(fileNames)
    modelInput = @(x) model(x, scans(i).a3);
    tmp = sort(scans(i).intTime); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intTime - bg0))./sum(scans(i).intTime - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(scans(i).intTime) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0, area0, cent0, sig0];
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2Err(scans(i).intTime, scans(i).intTimeErr, modelInput, scans(i).x0, errPts, fact, offset);

    % Plot a3 fit
    scans(i).a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 5e2)';
    scans(i).intCal = model(scans(i).xFit, scans(i).a3Cal);

    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: (', num2str(scans(i).meanHKL), '), ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '(', num2str(i), '/', num2str(length(scans)), ')'])
    xlabel('A_3 (deg.)')
    ylabel('I (det. / sec.)')
    errorbar(scans(i).a3, scans(i).intTime, scans(i).intTimeErr, 'o')
    plot(scans(i).a3Cal, scans(i).intCal)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    tiledlayout(2, 2, 'TileSpacing', 'compact')
    for j = 1:length(scans(i).x0)
        nexttile(j)
        hold on
        ylabel('\chi^2_r')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        box on
        hold off
    end
    nexttile(1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '\pm', num2str(round(scans(i).xErr(1), 4))])
    xlabel('Background (det. / sec.)')
    hold off
    nexttile(2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '\pm', num2str(round(scans(i).xErr(2), 4))])
    xlabel('Integrated Intensity (det. \times deg. / sec.)')
    hold off
    nexttile(3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '\pm', num2str(round(scans(i).xErr(3), 4))])
    xlabel('Gaussian Center (deg.)')
    hold off
    nexttile(4)
    hold on
    title(['Gaussian \sigma: ', num2str(round(scans(i).xFit(4), 4)), '\pm', num2str(round(scans(i).xErr(4), 4))])
    xlabel('Gaussian \sigma (deg.)')
    hold off
    pause(0.1)
end

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('Nuclear \sigma Versus Q')
ylabel('\sigma')
xlabel('\textsf{Q} $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scans), arrayfun(@(x) x.xFit(4), scans), arrayfun(@(x) x.xErr(4), scans), 'o')
yline(mean(arrayfun(@(x) x.xFit(4), scans)), '-', ['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])
box on
hold off
sigGau = mean(arrayfun(@(x) x.xFit(4), scans));
sigGauErr = sqrt(sum(arrayfun(@(x) x.xErr(4), scans).^2)/length(arrayfun(@(x) x.xErr(4), scans))); % 01/12/23 134 notes
disp(['Mean sigma: ', num2str(sigGau)])

% Import (001) 0.4 K, 3 T
fileNum001 = [847493];
fileName001 = strcat('fpx', string(fileNum001), '.bt7');
headerLines = 41;

scan001 = importDataBT7(fileLoc, fileName001, headerLines);
scan001.Q = 2*pi*sqrt(scan001.meanH^2 + scan001.meanK^2 + scan001.meanL^2)/a; % aStar=2pi/a then Q=aStar*sqrt(H^2+..)
scan001.a3Rad = scan001.a3*(pi/180);
scan001.meanA3Rad = scan001.meanA3*(pi/180);
scan001.deltaQ = scan001.Q.*(scan001.a3Rad - scan001.meanA3Rad);
scanGau = scan001;

% Fit (001) to Voigt and solve for correlation length
errPtsV = 0; % Number of parameter iterations when calculating the errorbar
factV = [0.4, 0.3, 0.002, 0.0]; % Determines range of iterated values for errorbar. [bg, peak, center, gamma]
offsetV = [0.0, 0.2, 0.002, 0.1]; % Determines range of iterated values for errorbar
modelV = @(x, a3, sig) x(1) + x(2).*voigt(a3, x(3), sigGau, x(4)); % Voigt function with flat background.

tmpV = sort(scan001.intTime); % Sorted intensities
bg0V = mean(tmpV(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
cent0V = sum(scan001.a3.*(scan001.intTime - bg0V))./sum(scan001.intTime - bg0V); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
gam0 = 0.01;
peak0 = (max(scan001.intTime) - bg0V)/voigt(cent0V, cent0V, sigGau, gam0); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
scan001.x0 = [bg0V, peak0, cent0V, gam0];
modelInputV = @(x) modelV(x, scan001.a3, sigGau);
[scan001.xFit, scan001.redChi2Fit, scan001.xErr, scan001.chiUpper, scan001.chiTrial, scan001.paramTrial, scan001.interpPts, ~, ~, scan001.paramLower, scan001.paramUpper] = fitRedChi2Err(scan001.intTime, scan001.intTimeErr, modelInputV, scan001.x0, errPtsV, factV, offsetV);

fwhm = 2*scan001.Q*scan001.xFit(4)*pi/180; % Fwhm from Lorentzian component of Voigt function
fwhmErr = 2*scan001.Q*scan001.xErr(4)*pi/180;
corr = 2/fwhm; % Factor of 2 is important to remember
corrErr = 2*fwhmErr/fwhm^2; % Square is from the quadrature sum partial derivatives times errors
fwhmMin = 2*scan001.Q*scan001.paramUpper(4)*pi/180; % Using fwhm max b/c refers to corr min
corrMin = 2/fwhmMin; % Factor of 2 is important to remember
fwhmMax = 2*scan001.Q*scan001.paramLower(4)*pi/180; % Fwhm from Lorentzian component of Voigt function
corrMax = 2/fwhmMax; % Factor of 2 is important to remember

% Values for plotting the fit
scan001.a3Cal = linspace(min(scan001.a3), max(scan001.a3), 500)';
scan001.intCal = modelV(scan001.xFit, scan001.a3Cal, sigGau);

% Convert to rlu, in HHL looking at 001 with a3 scan so taking HH0 projection
qH = scan001.Q.*sind(mean(scan001.a3, 1).*ones(size(scan001.a3)) - scan001.a3); % Here the a3 angle in real space corresponds to the rotated angle in reciprocal space because the lattice is cubic.
H = qH./aStar;
qHCal = scan001.Q.*sind(mean(scan001.a3Cal, 1).*ones(size(scan001.a3Cal))-scan001.a3Cal);
HCal = qHCal./aStar;

% Plot scans and fit, display correlation length, display errorbars in a
% separate figure for reference
figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
hold on
e = errorbar(H, scan001.intTime, scan001.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
p = plot(HCal, scan001.intCal, 'Color', 'b', 'LineWidth', 1);
xlabel('(HH0) (r.l.u.)')
ylabel('I (det. / sec.)')
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
box on
hold off

figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
tiledlayout(2, 2, 'TileSpacing', 'compact')
for i = 1:length(scan001.x0)
    nexttile(i)
    hold on
    ylabel('\chi^2_r')
    if ~any(isnan(scan001.interpPts(:, i))) % Check to see if an errorbar was determined before plotting
        scatter(scan001.paramTrial(:, i), scan001.chiTrial(:, i))
        plot([scan001.paramTrial(scan001.interpPts(1, i), i), scan001.paramTrial(scan001.interpPts(2, i), i)], [scan001.chiTrial(scan001.interpPts(1, i), i), scan001.chiTrial(scan001.interpPts(2, i), i)], 'b')
        plot([scan001.paramTrial(scan001.interpPts(3, i), i), scan001.paramTrial(scan001.interpPts(4, i), i)], [scan001.chiTrial(scan001.interpPts(3, i), i), scan001.chiTrial(scan001.interpPts(4, i), i)], 'b')
        yline(scan001.chiUpper, 'Color', 'r', 'LineWidth', 3.0)
        plot([scan001.paramLower(i), scan001.paramUpper(i)], [scan001.chiUpper, scan001.chiUpper], 'k-.o', 'LineWidth', 2.0)
        xlim([-inf inf])
        ylim([scan001.redChi2Fit-0.1 scan001.chiUpper+(scan001.chiUpper-scan001.redChi2Fit)])
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['Background: ', num2str(round(scan001.xFit(1), 4)), '\pm', num2str(round(scan001.xErr(1), 4))])
xlabel('Background (det. / sec.)')
hold off
nexttile(2)
hold on
title(['Voigt Peak: ', num2str(round(scan001.xFit(2), 4)), '\pm', num2str(round(scan001.xErr(2), 4))])
xlabel('Voigt Peak (det. / sec.)')
hold off
nexttile(3)
hold on
title(['Voigt Center: ', num2str(round(scan001.xFit(3), 4)), '\pm', num2str(round(scan001.xErr(3), 4))])
xlabel('Voigt Center (deg.)')
hold off
nexttile(4)
hold on
title(['Lorentzian HWHM: ', num2str(round(scan001.xFit(4), 5)), '\pm', num2str(round(scan001.xErr(4), 4))])
xlabel('Lorentzian HWHM (deg.)')
hold off
pause(0.1)

% Plot gof as a function of correlation length
gam = linspace(0.005, 0.05, 2e2); % Range of Lorentzian hwhm to plot
errPtsG = 0; % Number of parameter iterations when calculating the errorbar
factG = [0.4, 0.3, 0.005]; % Determines range of iterated values for errorbar. [bg, peak, center]
offsetG = [0.0, 0.2, 0.005]; % Determines range of iterated values for errorbar
modelG = @(x, a3, sig, gam) x(1) + x(2).*voigt(a3, x(3), sig, gam); % Voigt function with flat background. center, x0, sigma, and gamma

% Fit for correlation length
for i = 1:length(gam)
    tmp = sort(scan001.intTime); % Sorted intensities
    bg0G = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0G = sum(scan001.a3.*(scan001.intTime-bg0G))./sum(scan001.intTime-bg0G); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    peak0G = (max(scan001.intTime)-bg0G)/voigt(cent0G, cent0G, sigGau, gam(i)); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
    x0G = [bg0G, peak0G, cent0G];
    %dof = length(scan001.intTime)-length(x0G); % Number of datapoints minus number of free parameters
    modelInputG = @(x) modelG(x, scan001.a3, sigGau, gam(i));
    [xFitG(:,i), redChi2FitG(i), xErrG(:,i), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(scan001.intTime, scan001.intTimeErr, modelInputG, x0G, errPtsG, factG, offsetG);
    
    fwhmG(i) = 2*scan001.Q*gam(i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    corrG(i) = 2/fwhmG(i); % Factor of 2 is important to remember

    % Values for plotting the fit
    motorCalG(:,i) = linspace(min(scan001.a3), max(scan001.a3), 500)';
    intCalG(:,i) = modelG(xFitG(:,i), motorCalG, sigGau, gam(i));
end

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, 2, 'TileSpacing', 'Compact')
nexttile
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('$\xi (\mathrm{\AA})$', 'Interpreter', 'latex')
ylabel('\chi^2_r')
%yline(min(redChi2FitG)*(1+1/dof), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
yline(scan001.redChi2Fit*1.2, '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corrG, redChi2FitG, 'o')
pbaspect([16 9 1])
box on
hold off

% Plot the best fits and the fits right next to the gof threshold
chiUpper = scan001.redChi2Fit*1.2;
indMin = find(redChi2FitG<chiUpper & corrG<corr, 1, 'last');
indMax = find(redChi2FitG<chiUpper & corrG>corr, 1, 'first');

figure(6)
nexttile
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(H, scan001.intTime, scan001.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%p1 = plot(HCal, intCalG(:,indMin), 'Color', 'b', 'LineWidth', 1);
p2 = plot(HCal, scan001.intCal, 'Color', 'b', 'LineWidth', 1);
%p3 = plot(HCal, intCalG(:,indMax), 'Color', 'r', 'LineWidth', 1);
xlabel('(HH0) (r.l.u)')
ylabel('I (det. / sec.)')
%legend([p1, p2, p3], {['Min \xi: ', num2str(round(corrG(indMin), -2))], ['Best \xi: ', num2str(round(corr, -2))], ['Max \xi: ', num2str(round(corrG(indMax), -2))]}, 'Location', 'Best')
%legend('boxoff')
pbaspect([16 9 1])
box on
hold off

% Plot the Gaussian fit with width fixed to paramagnetic average next to the gof threshold
errPts = 1e2; % For a3 sigma fits
fact = [0.2, 0.3, 0]; % bg, area, center
offset = [0.01, 0.6, 0.05];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./sigGau.*exp(-((a3 - x(3))./sigGau).^2./2); % bg, area, center, sigma
modelInput = @(x) model(x, scanGau.a3);
tmp = sort(scanGau.intTime); % Sorted intensities
bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
cent0 = sum(scanGau.a3.*(scanGau.intTime - bg0))./sum(scanGau.intTime - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
area0 = (max(scanGau.intTime) - bg0)*sqrt(2*pi)*sigGau; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
scanGau.x0 = [bg0, area0, cent0];
[scanGau.xFit,scanGau.redChi2Fit,scanGau.xErr,scanGau.chiUpper,scanGau.chiTrial,scanGau.paramTrial,scanGau.interpPts,scanGau.slopes,scanGau.intercepts,scanGau.paramLower,scanGau.paramUpper] = fitRedChi2Err(scanGau.intTime, scanGau.intTimeErr, modelInput, scanGau.x0, errPts, fact, offset);
scanGau.a3Cal = linspace(min(scanGau.a3), max(scanGau.a3), 5e2)';
scanGau.intCal = model(scanGau.xFit, scanGau.a3Cal);
qHGau = scanGau.Q.*sind(mean(scanGau.a3, 1).*ones(size(scanGau.a3)) - scanGau.a3); % Here the a3 angle in real space corresponds to the rotated angle in reciprocal space because the lattice is cubic.
HGau = qHGau./aStar;
qHCalGau = scanGau.Q.*sind(mean(scanGau.a3Cal, 1).*ones(size(scanGau.a3Cal))-scanGau.a3Cal);
HCalGau = qHCalGau./aStar;

figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
tiledlayout(3, 1, 'TileSpacing', 'compact')
for j = 1:length(scanGau.x0)
    nexttile(j)
    hold on
    ylabel('\chi^2_r')
    if ~any(isnan(scanGau.interpPts(:, j))) % Check to see if an errorbar was determined before plotting
        scatter(scanGau.paramTrial(:, j), scanGau.chiTrial(:, j))
        plot([scanGau.paramTrial(scanGau.interpPts(1, j), j), scanGau.paramTrial(scanGau.interpPts(2, j), j)], [scanGau.chiTrial(scanGau.interpPts(1, j), j), scanGau.chiTrial(scanGau.interpPts(2, j), j)], 'b')
        plot([scanGau.paramTrial(scanGau.interpPts(3, j), j), scanGau.paramTrial(scanGau.interpPts(4, j), j)], [scanGau.chiTrial(scanGau.interpPts(3, j), j), scanGau.chiTrial(scanGau.interpPts(4, j), j)], 'b')
        yline(scanGau.chiUpper, 'Color', 'r', 'LineWidth', 3.0)
        plot([scanGau.paramLower(j), scanGau.paramUpper(j)], [scanGau.chiUpper, scanGau.chiUpper], 'k-.o', 'LineWidth', 2.0)
        xlim([-inf inf])
        ylim([scanGau.redChi2Fit-0.1 scanGau.chiUpper+(scanGau.chiUpper-scanGau.redChi2Fit)])
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['Background: ', num2str(round(scanGau.xFit(1), 4)), '\pm', num2str(round(scanGau.xErr(1), 4))])
xlabel('Background (det. / sec.)')
hold off
nexttile(2)
hold on
title(['Integrated Intensity: ', num2str(round(scanGau.xFit(2), 4)), '\pm', num2str(round(scanGau.xErr(2), 4))])
xlabel('Integrated Intensity (det. \times deg. / sec.)')
hold off
nexttile(3)
hold on
title(['Gaussian Center: ', num2str(round(scanGau.xFit(3), 4)), '\pm', num2str(round(scanGau.xErr(3), 4))])
xlabel('Gaussian Center (deg.)')
hold off

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, 2, 'TileSpacing', 'Compact')
nexttile
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('$\xi (\mathrm{\AA})$', 'Interpreter', 'latex')
ylabel('\chi^2_r')
%yline(min(redChi2FitG)*(1+1/dof), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
yline(scan001.redChi2Fit*1.2, '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corrG, redChi2FitG, 'o')
pbaspect([16 9 1])
box on
hold off
nexttile
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(HGau, scanGau.intTime, scanGau.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(HCalGau, scanGau.intCal, 'Color', 'b', 'LineWidth', 1);
xlabel('(HH0) (r.l.u.)')
ylabel('I (det. / sec.)')
pbaspect([16 9 1])
box on
hold off

disp(['0.4 K, 3 T (001) Correlation Length: ', num2str(corr), '+-', num2str(corrErr), ' Angstroms'])
disp(['0.4 K, 3 T (001) Correlation Length Min: ', num2str(corrMin), ' Angstroms'])
disp(['0.4 K, 3 T (001) Correlation Length Max: ', num2str(corrMax), ' Angstroms'])

corrBnd = 2/(scan001.Q*(sigGau*sqrt(8*log(2))*pi/180)); % Because the lower bound above is nonsensically large since the peak is basically fit by a Gaussian, we take the lower bound to be given by 1/hwhwm=2/delQ=2/(delOmeg*Q)
corrBndErr = 2/(scan001.Q*(sqrt(8*log(2))*pi/180))*sigGauErr/sigGau^2; % 01/12/23 134 notes
disp(['001 Correlation Length Lower Bound: ', num2str(corrBnd), '+-', num2str(corrBndErr)])

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\001Corr.png';
exportgraphics(gcf, fdir, 'resolution', 450)

%% SI: (1/2,1/2,1/2) correlation length

clear
close all

set(0, 'defaultaxesfontsize', 7)

% Fit to the (1/2,1/2,1/2) reflection with Voigt whose Gaussian width is fixed to
% the average value for prominent nuclear peaks, then pick out Lorentzian
% hwhm for correlation length. Normalize to time.

% Import prominent 5 K, experiment without magnet peaks.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\EuPd3S4\All\';
fileNum = [845748, 845749, 845750, 845772, 845774, 845786, 845795, 845800, 845803];
fileNames = strcat('fpx', string(fileNum), '.bt7');
headerLines = 41;

scans = importDataBT7(fileLoc, fileNames, headerLines);

a = 6.67570;
aStar = 2*pi/a;
for i=1:length(fileNames)
    scans(i).Q = 2*pi*sqrt(scans(i).meanH^2 + scans(i).meanK^2 + scans(i).meanL^2)/a; % aStar=2pi/a then Q=aStar*sqrt(H^2+..)
end

% Fit nuclear peaks to a Gaussian and flat background, solving for the
% Gaussian sigma
errPts = 1e2; % For a3 sigma fits
fact = [0.2, 0.3, 0, 0.3]; % bg, area, center, sigma
offset = [0.01, 0.6, 0.05, 1e-3];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./x(4).*exp(-((a3 - x(3))./x(4)).^2./2); % bg, area, center, sigma
for i = 1:length(fileNames)
    modelInput = @(x) model(x, scans(i).a3);
    tmp = sort(scans(i).intTime); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intTime - bg0))./sum(scans(i).intTime - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.1;
    area0 = (max(scans(i).intTime) - bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0, area0, cent0, sig0];
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2Err(scans(i).intTime, scans(i).intTimeErr, modelInput, scans(i).x0, errPts, fact, offset);

    % Plot a3 fit
    scans(i).a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 5e2)';
    scans(i).intCal = model(scans(i).xFit, scans(i).a3Cal);

    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: (', num2str(scans(i).meanHKL), '), ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '(', num2str(i), '/', num2str(length(scans)), ')'])
    xlabel('A_3 (deg.)')
    ylabel('I (det. / sec.)')
    errorbar(scans(i).a3, scans(i).intTime, scans(i).intTimeErr, 'o')
    plot(scans(i).a3Cal, scans(i).intCal)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    tiledlayout(2, 2, 'TileSpacing', 'compact')
    for j = 1:length(scans(i).x0)
        nexttile(j)
        hold on
        ylabel('\chi^2_r')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        box on
        hold off
    end
    nexttile(1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '\pm', num2str(round(scans(i).xErr(1), 4))])
    xlabel('Background (det. / sec.)')
    hold off
    nexttile(2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '\pm', num2str(round(scans(i).xErr(2), 4))])
    xlabel('Integrated Intensity (det. \times deg. / sec.)')
    hold off
    nexttile(3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '\pm', num2str(round(scans(i).xErr(3), 4))])
    xlabel('Gaussian Center (deg.)')
    hold off
    nexttile(4)
    hold on
    title(['Gaussian \sigma: ', num2str(round(scans(i).xFit(4), 4)), '\pm', num2str(round(scans(i).xErr(4), 4))])
    xlabel('Gaussian \sigma (deg.)')
    hold off
    pause(0.1)
end

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('Nuclear \sigma Versus Q')
ylabel('\sigma')
xlabel('\textsf{Q} $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scans), arrayfun(@(x) x.xFit(4), scans), arrayfun(@(x) x.xErr(4), scans), 'o')
yline(mean(arrayfun(@(x) x.xFit(4), scans)), '-', ['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])
box on
hold off
sigGau = mean(arrayfun(@(x) x.xFit(4), scans));
sigGauErr = sqrt(sum(arrayfun(@(x) x.xErr(4), scans).^2)/length(arrayfun(@(x) x.xErr(4), scans))); % 01/12/23 134 notes
disp(['Mean sigma: ', num2str(sigGau)])

% Import (1/2,1/2,1/2) 1.5 K, 0 T
fileNum0p5 = [845614];
fileName0p5 = strcat('fpx', string(fileNum0p5), '.bt7');
headerLines = 41;

scan0p5 = importDataBT7(fileLoc, fileName0p5, headerLines);
scan0p5.Q = 2*pi*sqrt(scan0p5.meanH^2 + scan0p5.meanK^2 + scan0p5.meanL^2)/a; % aStar=2pi/a then Q=aStar*sqrt(H^2+..)
scan0p5.a3Rad = scan0p5.a3*(pi/180);
scan0p5.meanA3Rad = scan0p5.meanA3*(pi/180);
scan0p5.deltaQ = scan0p5.Q.*(scan0p5.a3Rad - scan0p5.meanA3Rad);
scanGau = scan0p5;

% Fit (1/2,1/2,1/2) to Voigt and solve for correlation length
errPtsV = 0; % Number of parameter iterations when calculating the errorbar
factV = [0.4, 0.3, 0.002, 0.0]; % Determines range of iterated values for errorbar. [bg, peak, center, gamma]
offsetV = [0.0, 0.2, 0.002, 0.1]; % Determines range of iterated values for errorbar
modelV = @(x, a3, sig) x(1) + x(2).*voigt(a3, x(3), sigGau, x(4)); % Voigt function with flat background.

tmpV = sort(scan0p5.intTime); % Sorted intensities
bg0V = mean(tmpV(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
cent0V = sum(scan0p5.a3.*(scan0p5.intTime - bg0V))./sum(scan0p5.intTime - bg0V); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
gam0 = 0.01;
peak0 = (max(scan0p5.intTime) - bg0V)/voigt(cent0V, cent0V, sigGau, gam0); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
scan0p5.x0 = [bg0V, peak0, cent0V, gam0];
modelInputV = @(x) modelV(x, scan0p5.a3, sigGau);
[scan0p5.xFit, scan0p5.redChi2Fit, scan0p5.xErr, scan0p5.chiUpper, scan0p5.chiTrial, scan0p5.paramTrial, scan0p5.interpPts, ~, ~, scan0p5.paramLower, scan0p5.paramUpper] = fitRedChi2Err(scan0p5.intTime, scan0p5.intTimeErr, modelInputV, scan0p5.x0, errPtsV, factV, offsetV);

fwhm = 2*scan0p5.Q*scan0p5.xFit(4)*pi/180; % Fwhm from Lorentzian component of Voigt function
fwhmErr = 2*scan0p5.Q*scan0p5.xErr(4)*pi/180;
corr = 2/fwhm; % Factor of 2 is important to remember
corrErr = 2*fwhmErr/fwhm^2; % Square is from the quadrature sum partial derivatives times errors
fwhmMin = 2*scan0p5.Q*scan0p5.paramUpper(4)*pi/180; % Using fwhm max b/c refers to corr min
corrMin = 2/fwhmMin; % Factor of 2 is important to remember
fwhmMax = 2*scan0p5.Q*scan0p5.paramLower(4)*pi/180; % Fwhm from Lorentzian component of Voigt function
corrMax = 2/fwhmMax; % Factor of 2 is important to remember

% Values for plotting the fit
scan0p5.a3Cal = linspace(min(scan0p5.a3), max(scan0p5.a3), 500)';
scan0p5.intCal = modelV(scan0p5.xFit, scan0p5.a3Cal, sigGau);

% Convert to rlu, in HHL looking at (1/2,1/2,1/2) with a3 scan so taking
% HH-H projection (see 1/18/2023 notes)
qH = scan0p5.Q.*sind(mean(scan0p5.a3, 1).*ones(size(scan0p5.a3)) - scan0p5.a3); % Here the a3 angle in real space corresponds to the rotated angle in reciprocal space because the lattice is cubic.
H = qH./aStar;
qHCal = scan0p5.Q.*sind(mean(scan0p5.a3Cal, 1).*ones(size(scan0p5.a3Cal))-scan0p5.a3Cal);
HCal = qHCal./aStar;

% Plot scans and fit, display correlation length, display errorbars in a
% separate figure for reference
figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
hold on
e = errorbar(H, scan0p5.intTime, scan0p5.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
p = plot(HCal, scan0p5.intCal, 'Color', 'b', 'LineWidth', 1);
xlabel('(\textsf{HH}$\overline{\textsf{H}}$) (r.l.u.)', 'Interpreter', 'latex')
ylabel('I (det. / sec.)')
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
box on
hold off

figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
tiledlayout(2, 2, 'TileSpacing', 'compact')
for i = 1:length(scan0p5.x0)
    nexttile(i)
    hold on
    ylabel('\chi^2_r')
    if ~any(isnan(scan0p5.interpPts(:, i))) % Check to see if an errorbar was determined before plotting
        scatter(scan0p5.paramTrial(:, i), scan0p5.chiTrial(:, i))
        plot([scan0p5.paramTrial(scan0p5.interpPts(1, i), i), scan0p5.paramTrial(scan0p5.interpPts(2, i), i)], [scan0p5.chiTrial(scan0p5.interpPts(1, i), i), scan0p5.chiTrial(scan0p5.interpPts(2, i), i)], 'b')
        plot([scan0p5.paramTrial(scan0p5.interpPts(3, i), i), scan0p5.paramTrial(scan0p5.interpPts(4, i), i)], [scan0p5.chiTrial(scan0p5.interpPts(3, i), i), scan0p5.chiTrial(scan0p5.interpPts(4, i), i)], 'b')
        yline(scan0p5.chiUpper, 'Color', 'r', 'LineWidth', 3.0)
        plot([scan0p5.paramLower(i), scan0p5.paramUpper(i)], [scan0p5.chiUpper, scan0p5.chiUpper], 'k-.o', 'LineWidth', 2.0)
        xlim([-inf inf])
        ylim([scan0p5.redChi2Fit-0.1 scan0p5.chiUpper+(scan0p5.chiUpper-scan0p5.redChi2Fit)])
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['Background: ', num2str(round(scan0p5.xFit(1), 4)), '\pm', num2str(round(scan0p5.xErr(1), 4))])
xlabel('Background (det. / sec.)')
hold off
nexttile(2)
hold on
title(['Voigt Peak: ', num2str(round(scan0p5.xFit(2), 4)), '\pm', num2str(round(scan0p5.xErr(2), 4))])
xlabel('Voigt Peak (det. / sec.)')
hold off
nexttile(3)
hold on
title(['Voigt Center: ', num2str(round(scan0p5.xFit(3), 4)), '\pm', num2str(round(scan0p5.xErr(3), 4))])
xlabel('Voigt Center (deg.)')
hold off
nexttile(4)
hold on
title(['Lorentzian HWHM: ', num2str(round(scan0p5.xFit(4), 5)), '\pm', num2str(round(scan0p5.xErr(4), 4))])
xlabel('Lorentzian HWHM (deg.)')
hold off
pause(0.1)

% Plot gof as a function of correlation length
gam = linspace(0.002, 0.02, 2e2); % Range of Lorentzian hwhm to plot
errPtsG = 0; % Number of parameter iterations when calculating the errorbar
factG = [0.4, 0.3, 0.005]; % Determines range of iterated values for errorbar. [bg, peak, center]
offsetG = [0.0, 0.2, 0.005]; % Determines range of iterated values for errorbar
modelG = @(x, a3, sig, gam) x(1) + x(2).*voigt(a3, x(3), sig, gam); % Voigt function with flat background. center, x0, sigma, and gamma

% Fit for correlation length
for i = 1:length(gam)
    tmp = sort(scan0p5.intTime); % Sorted intensities
    bg0G = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0G = sum(scan0p5.a3.*(scan0p5.intTime-bg0G))./sum(scan0p5.intTime-bg0G); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    peak0G = (max(scan0p5.intTime)-bg0G)/voigt(cent0G, cent0G, sigGau, gam(i)); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
    x0G = [bg0G, peak0G, cent0G];
    %dof = length(scan0p5.intTime)-length(x0G); % Number of datapoints minus number of free parameters
    modelInputG = @(x) modelG(x, scan0p5.a3, sigGau, gam(i));
    [xFitG(:,i), redChi2FitG(i), xErrG(:,i), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(scan0p5.intTime, scan0p5.intTimeErr, modelInputG, x0G, errPtsG, factG, offsetG);
    
    fwhmG(i) = 2*scan0p5.Q*gam(i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    corrG(i) = 2/fwhmG(i); % Factor of 2 is important to remember

    % Values for plotting the fit
    motorCalG(:,i) = linspace(min(scan0p5.a3), max(scan0p5.a3), 500)';
    intCalG(:,i) = modelG(xFitG(:,i), motorCalG, sigGau, gam(i));
end

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, 2, 'TileSpacing', 'Compact')
nexttile
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('$\xi (\mathrm{\AA})$', 'Interpreter', 'latex')
ylabel('\chi^2_r')
%yline(min(redChi2FitG)*(1+1/dof), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
yline(scan0p5.redChi2Fit*1.2, '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corrG, redChi2FitG, 'o')
pbaspect([16 9 1])
box on
hold off

% Plot the best fits and the fits right next to the gof threshold
chiUpper = scan0p5.redChi2Fit*1.2;
indMin = find(redChi2FitG<chiUpper & corrG<corr, 1, 'last');
indMax = find(redChi2FitG<chiUpper & corrG>corr, 1, 'first');

figure(6)
nexttile
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(H, scan0p5.intTime, scan0p5.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%p1 = plot(HCal, intCalG(:,indMin), 'Color', 'b', 'LineWidth', 1);
p2 = plot(HCal, scan0p5.intCal, 'Color', 'b', 'LineWidth', 1);
%p3 = plot(HCal, intCalG(:,indMax), 'Color', 'r', 'LineWidth', 1);
xlabel('(\textsf{HH}$\overline{\textsf{H}}$) (r.l.u.)', 'Interpreter', 'latex')
ylabel('I (det. / sec.)')
%legend([p1, p2, p3], {['Min \xi: ', num2str(round(corrG(indMin), -2))], ['Best \xi: ', num2str(round(corr, -2))], ['Max \xi: ', num2str(round(corrG(indMax), -2))]}, 'Location', 'Best')
%legend('boxoff')
pbaspect([16 9 1])
box on
hold off

% Plot the Gaussian fit with width fixed to paramagnetic average next to the gof threshold
errPts = 1e2; % For a3 sigma fits
fact = [0.2, 0.3, 0]; % bg, area, center
offset = [0.01, 0.6, 0.05];
model = @(x, a3) x(1) + x(2)./sqrt(2.*pi)./sigGau.*exp(-((a3 - x(3))./sigGau).^2./2); % bg, area, center, sigma
modelInput = @(x) model(x, scanGau.a3);
tmp = sort(scanGau.intTime); % Sorted intensities
bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
cent0 = sum(scanGau.a3.*(scanGau.intTime - bg0))./sum(scanGau.intTime - bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
area0 = (max(scanGau.intTime) - bg0)*sqrt(2*pi)*sigGau; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
scanGau.x0 = [bg0, area0, cent0];
[scanGau.xFit,scanGau.redChi2Fit,scanGau.xErr,scanGau.chiUpper,scanGau.chiTrial,scanGau.paramTrial,scanGau.interpPts,scanGau.slopes,scanGau.intercepts,scanGau.paramLower,scanGau.paramUpper] = fitRedChi2Err(scanGau.intTime, scanGau.intTimeErr, modelInput, scanGau.x0, errPts, fact, offset);
scanGau.a3Cal = linspace(min(scanGau.a3), max(scanGau.a3), 5e2)';
scanGau.intCal = model(scanGau.xFit, scanGau.a3Cal);
qHGau = scanGau.Q.*sind(mean(scanGau.a3, 1).*ones(size(scanGau.a3)) - scanGau.a3); % Here the a3 angle in real space corresponds to the rotated angle in reciprocal space because the lattice is cubic.
HGau = qHGau./aStar;
qHCalGau = scanGau.Q.*sind(mean(scanGau.a3Cal, 1).*ones(size(scanGau.a3Cal))-scanGau.a3Cal);
HCalGau = qHCalGau./aStar;

figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
tiledlayout(3, 1, 'TileSpacing', 'compact')
for j = 1:length(scanGau.x0)
    nexttile(j)
    hold on
    ylabel('\chi^2_r')
    if ~any(isnan(scanGau.interpPts(:, j))) % Check to see if an errorbar was determined before plotting
        scatter(scanGau.paramTrial(:, j), scanGau.chiTrial(:, j))
        plot([scanGau.paramTrial(scanGau.interpPts(1, j), j), scanGau.paramTrial(scanGau.interpPts(2, j), j)], [scanGau.chiTrial(scanGau.interpPts(1, j), j), scanGau.chiTrial(scanGau.interpPts(2, j), j)], 'b')
        plot([scanGau.paramTrial(scanGau.interpPts(3, j), j), scanGau.paramTrial(scanGau.interpPts(4, j), j)], [scanGau.chiTrial(scanGau.interpPts(3, j), j), scanGau.chiTrial(scanGau.interpPts(4, j), j)], 'b')
        yline(scanGau.chiUpper, 'Color', 'r', 'LineWidth', 3.0)
        plot([scanGau.paramLower(j), scanGau.paramUpper(j)], [scanGau.chiUpper, scanGau.chiUpper], 'k-.o', 'LineWidth', 2.0)
        xlim([-inf inf])
        ylim([scanGau.redChi2Fit-0.1 scanGau.chiUpper+(scanGau.chiUpper-scanGau.redChi2Fit)])
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['Background: ', num2str(round(scanGau.xFit(1), 4)), '\pm', num2str(round(scanGau.xErr(1), 4))])
xlabel('Background (det. / sec.)')
hold off
nexttile(2)
hold on
title(['Integrated Intensity: ', num2str(round(scanGau.xFit(2), 4)), '\pm', num2str(round(scanGau.xErr(2), 4))])
xlabel('Integrated Intensity (det. \times deg. / sec.)')
hold off
nexttile(3)
hold on
title(['Gaussian Center: ', num2str(round(scanGau.xFit(3), 4)), '\pm', num2str(round(scanGau.xErr(3), 4))])
xlabel('Gaussian Center (deg.)')
hold off

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.8]) % Left, bottom, width, height
tiledlayout(1, 2, 'TileSpacing', 'Compact')
nexttile
hold on
text(0.9, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
xlabel('$\xi (\mathrm{\AA})$', 'Interpreter', 'latex')
ylabel('\chi^2_r')
%yline(min(redChi2FitG)*(1+1/dof), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
yline(scan0p5.redChi2Fit*1.2, '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corrG, redChi2FitG, 'o')
pbaspect([16 9 1])
box on
hold off
nexttile
hold on
text(0.9, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(HGau, scanGau.intTime, scanGau.intTimeErr, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(HCalGau, scanGau.intCal, 'Color', 'b', 'LineWidth', 1);
xlabel('(\textsf{HH}$\overline{\textsf{H}}$) (r.l.u.)', 'Interpreter', 'latex')
ylabel('I (det. / sec.)')
pbaspect([16 9 1])
box on
hold off

disp(['1.5 K, 0 T (1/2,1/2,1/2) Correlation Length: ', num2str(corr), '+-', num2str(corrErr), ' Angstroms'])
disp(['1.5 K, 0 T (1/2,1/2,1/2) Correlation Length Min: ', num2str(corrMin), ' Angstroms'])
disp(['1.5 K, 0 T (1/2,1/2,1/2) Correlation Length Max: ', num2str(corrMax), ' Angstroms'])

corrBnd = 2/(scan0p5.Q*(sigGau*sqrt(8*log(2))*pi/180)); % Because the lower bound above is nonsensically large since the peak is basically fit by a Gaussian, we take the lower bound to be given by 1/hwhwm=2/delQ=2/(delOmeg*Q)
corrBndErr = 2/(scan0p5.Q*(sqrt(8*log(2))*pi/180))*sigGauErr/sigGau^2; % 01/12/23 134 notes
disp(['(1/2,1/2,1/2) Correlation Length Lower Bound: ', num2str(corrBnd), '+-', num2str(corrBndErr)])

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\0p5Corr.png';
exportgraphics(gcf, fdir, 'resolution', 450)

%% SI: BT7 Refinement GOF figure. r=0.5 mm spherical absorption correction.

clear
close all

set(0, 'defaultaxesfontsize', 7)

% 5.1 K & 0 T Nuclear, 1.5 K & 0 T AFM, 4.7 K & 0 T Nuclear, 0.3 K & 0.5 T AFM, 0.3 K & 0.5 T FM, 0.4 K & 3 T FM
scale = [2.232, 2.232, 1.463, 1.463, 1.463, 1.463]; % Refined scale factor from FullProf. For converting units back to barns.
frame = ["\bf{a}"; "\bf{b}"; "\bf{c}"; "\bf{d}"; "\bf{e}"; "\bf{f}"];
label = ["5.1 K, 0 T"; "1.5 K, 0 T"; "4.7 K, 0 T"; "0.3 K, 0.5 T"; "0.3 K, 0.5 T"; "0.4 K, 3 T"];
rad = ["\it{r}\rm{ = 0.5 mm}"; "\it{r}\rm{ = 0.5 mm}"; "\it{r}\rm{ = 0.5 mm}"; "\it{r}\rm{ = 0.5 mm}"; "\it{r}\rm{ = 0.5 mm}"; "\it{r}\rm{ = 0.5 mm}"];
labelStr = ["Nuclear"; "AFM"; "Nuclear"; "AFM"; "FM"; "FM"];

file(1) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\5K\EuPd3S45K.prf');
file(2) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\1p5K\sc\EuPd3S41p5K.prf');
file(3) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\4p5K\EuPd3S44p5K.prf');
file(4) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\0p5TAFM\ab\sc\EuPd3S40p5TAFM.prf');
file(5) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\0p5TFM\sc\EuPd3S40p5TFM.prf');
file(6) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\0p5mm\3T\sc\EuPd3S43T.prf');

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 8.25])
tiledlayout(4,3)
for i=1:length(file)
    data = file(i).data;
    
    f2Obs = data(:, 2)/scale(i);
    sigObs = data(:, 4)/scale(i); % For the observed structure factors
    f2Cal = data(:, 3)/scale(i);
    
    lin = linspace(min([0; f2Cal; f2Obs - sigObs])*1.5, max([f2Cal; f2Obs + sigObs])*1.5);

    nexttile
    hold on
    text(0.05, 0.95, frame(i), 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.85, label(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.75, rad(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.65, labelStr(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    xlabel('\rm{\mid}\it{F}\rm{\mid^2_{cal} (b)}')
    ylabel('\rm{\mid}\it{F}\rm{\mid^2_{obs} (b)}')
    errorbar(f2Cal, f2Obs, sigObs, 'o', 'MarkerFaceColor', 'w')
    plot(lin, lin, 'k', 'LineWidth', 1)
    axis([min([0; f2Cal; f2Obs - sigObs])*1.1 max([f2Cal; f2Obs + sigObs])*1.1 min([0; f2Cal; f2Obs - sigObs])*1.1 max([f2Cal; f2Obs + sigObs])*1.1])
    axis square
    set(gca, 'XTick', get(gca, 'YTick'));
    set(gca, 'TickLength', [0.02, 0.01])
    box on
    hold off
end

% SI: BT7 Refinement GOF figure. r=1.5 mm spherical absorption correction.

clear

% 5.1 K & 0 T Nuclear, 1.5 K & 0 T AFM, 4.7 K & 0 T Nuclear, 0.3 K & 0.5 T AFM, 0.3 K & 0.5 T FM, 0.4 K & 3 T FM
scale = [11.45, 11.45, 10.55, 10.55, 10.55, 10.55]; % Refined scale factor from FullProf. For converting units back to barns.
frame = ["\bf{g}"; "\bf{h}"; "\bf{i}"; "\bf{j}"; "\bf{k}"; "\bf{l}"];
label = ["5.1 K, 0 T"; "1.5 K, 0 T"; "4.7 K, 0 T"; "0.3 K, 0.5 T"; "0.3 K, 0.5 T"; "0.4 K, 3 T"];
rad = ["\it{r}\rm{ = 1.5 mm}"; "\it{r}\rm{ = 1.5 mm}"; "\it{r}\rm{ = 1.5 mm}"; "\it{r}\rm{ = 1.5 mm}"; "\it{r}\rm{ = 1.5 mm}"; "\it{r}\rm{ = 1.5 mm}"];
labelStr = ["Nuclear"; "AFM"; "Nuclear"; "AFM"; "FM"; "FM"];

file(1) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\5K\EuPd3S45K.prf');
file(2) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\1p5K\sc\EuPd3S41p5K.prf');
file(3) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\4p5K\EuPd3S44p5K.prf');
file(4) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\0p5TAFM\ab\sc\EuPd3S40p5TAFM.prf');
file(5) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\0p5TFM\sc\EuPd3S40p5TFM.prf');
file(6) = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\EuPd3S4\Finished\SphericalCorrection\1p5mm\3T\sc\EuPd3S43T.prf');

figure(1)
for i=1:length(file)
    data = file(i).data;
    
    f2Obs = data(:, 2)/scale(i);
    sigObs = data(:, 4)/scale(i); % For the observed structure factors
    f2Cal = data(:, 3)/scale(i);
    
    lin = linspace(min([0; f2Cal; f2Obs - sigObs])*1.5, max([f2Cal; f2Obs + sigObs])*1.5);

    nexttile
    hold on
    text(0.05, 0.95, frame(i), 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.85, label(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.75, rad(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    text(0.05, 0.65, labelStr(i), 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'interpreter', 'tex')
    xlabel('\rm{\mid}\it{F}\rm{\mid^2_{cal} (b)}')
    ylabel('\rm{\mid}\it{F}\rm{\mid^2_{obs} (b)}')
    errorbar(f2Cal, f2Obs, sigObs, 'o', 'MarkerFaceColor', 'w')
    plot(lin, lin, 'k', 'LineWidth', 1)
    axis([min([0; f2Cal; f2Obs - sigObs])*1.1 max([f2Cal; f2Obs + sigObs])*1.1 min([0; f2Cal; f2Obs - sigObs])*1.1 max([f2Cal; f2Obs + sigObs])*1.1])
    axis square
    set(gca, 'XTick', get(gca, 'YTick'));
    set(gca, 'TickLength', [0.02, 0.01])
    box on
    hold off
end

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\BT7GOF.emf';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% Compare some November to September 2020 hc pulses from field scans at higher fields

clear
close all

set(0, 'defaultaxesfontsize', 7)

factorM = 1.07; % Calibration factor for thin plate hc determined by overplotting with 2.34 mg data
massNom = 0.194; % Nominal sample mass in mg
mass = massNom*factorM;
molarMass = 599.46; % Sample grams per mole

% First import the data: file1 is 9/15 T-scan at 0 T, file2 is 11/14 T-scans at a few fields
fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileLoc2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity111420\';
fileName1 = '0p19mg_drhc_EuPd2S4_09182020_fieldsweeps.raw';
fileName2 = '0p19mg_drhc_EuPd2S4_11132020.raw';
dir1 = strcat(fileLoc1, fileName1);
dir2 = strcat(fileLoc2, fileName2);
indPulseParams1 = [274601, 274649];
indData1 = [274651, 274906];
indPulseParams2 = [73804, 73852];
indData2 = [73854, 74109];

optsPulseParams1 = delimitedTextImportOptions('DataLines', indPulseParams1, 'VariableNamesLine', 10, 'Delimiter', ',');
optsData1 = delimitedTextImportOptions('NumVariables', 7, 'DataLines', indData1, 'VariableNamesLine', 10, 'Delimiter', ',', 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'});
pulseParams1 = readtable(dir1, optsPulseParams1);
data1 = readtable(dir1, optsData1);
optsPulseParams2 = delimitedTextImportOptions('DataLines', indPulseParams2, 'VariableNamesLine', 10, 'Delimiter', ',');
optsData2 = delimitedTextImportOptions('NumVariables', 7, 'DataLines', indData2, 'VariableNamesLine', 10, 'Delimiter', ',', 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'});
pulseParams2 = readtable(dir2, optsPulseParams2);
data2 = readtable(dir2, optsData2);

time1 = data1.TimeStamp_Seconds_ - data1.TimeStamp_Seconds_(1);
platformTemp1 = data1.PlatformTemp_Kelvin_;
platformTempFit1 = data1.PlatformTempFit_Kelvin_;
sampleTempFit1 = data1.SampleTempFit_Kelvin_;
time2 = data2.TimeStamp_Seconds_ - data2.TimeStamp_Seconds_(1);
platformTemp2 = data2.PlatformTemp_Kelvin_;
platformTempFit2 = data2.PlatformTempFit_Kelvin_;
sampleTempFit2 = data2.SampleTempFit_Kelvin_;

% Plot the hc pulses
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.7])
tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
hold on
ylabel('\it{T} \rm{(K)}')
s1 = scatter(time1, platformTemp1, 15, 'b');
%p1 = plot(time1, platformTempFit1, 'Color', 'b', 'LineWidth', 2);
p2 = plot(time1, sampleTempFit1, 'Color', 'r', 'LineWidth', 2);
text(0.1, 0.9, 'Sept. 2020', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'Color', 'k')
%legend([p1, p2], {'Platform', 'Sample'})
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
xlim([0, max([time1; time2])])
xticklabels([])
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
xlabel('Time (sec.)')
ylabel('\it{T} \rm{(K)}')
s2 = scatter(time2, platformTemp2, 15, 'b');
%p3 = plot(time2, platformTempFit2, 'Color', 'b', 'LineWidth', 2);
p4 = plot(time2, sampleTempFit2, 'Color', 'r', 'LineWidth', 2);
text(0.1, 0.9, 'Nov. 2020', 'Units', 'normalized', 'fontsize', 7, 'fontname', 'Arial', 'Color', 'k')
%legend([p3, p4], {'Platform', 'Sample'})
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'FontName', 'Arial')
xlim([0, max([time1; time2])])
pbaspect([16 9 1])
box on
hold off

%% TF Dipolar Heisenberg Plots. Plot susc*T/C.

clear
close all

set(0, 'defaultaxesfontsize', 7)

% First import the data
molarMass = 599.46; % Grams per mole
mass1 = 12.6;
% mass2 = 1; %?
fileLoc1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\ACMS\EuPd3S4\';
% fileLoc2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\EuPd3S4\';
file(1).filename = 'ML_Namond_02092019_1_TBe220-ACMS-magvT-1T-12point6mg.dat';
% file(2).filename = 'ML_TBe350_09192021_EuPd3S4_singlecrystal_test sample_huge sample.dat';
file(1).data = readtable(strcat(fileLoc1, file(1).filename), 'FileType', 'text', 'NumHeaderLines', 30, 'Delimiter', ',');
% file(2).data = readtable(strcat(fileLoc2, file(2).filename), 'FileType', 'text', 'NumHeaderLines', 44, 'Delimiter', ',');
temperature1 = file(1).data.Temperature_K_;
field1 = file(1).data.MagneticField_Oe_; % Oersted
moment1 = file(1).data.M_DC_emu_/(mass1/1e3)*molarMass; % Convert to emu/mol
momentErr1 = file(1).data.M_Std_Dev__emu_/(mass1/1e3)*molarMass; % Convert to emu/mol
% temperature2 = file(2).data.Temperature_K_;
% field2 = file(2).data.MagneticField_Oe_; % Oersted
% moment2 = file(2).data.DCMomentFixedCtr_emu_/(mass2/1e3)*molarMass; % Convert to emu/mol
% momentErr2 = file(2).data.DCMomentErrFixedCtr_emu_/(mass2/1e3)*molarMass; % Convert to emu/mol

susc1 = moment1./field1;
suscErr1 = momentErr1./field1;
suscRecip1 = 1./susc1;
suscRecipErr1 = suscErr1./susc1.^2; % 01/12/23 notes
% susc2 = moment2./field2;
% suscErr2 = momentErr2./field2;

% Plot the data
% figure('Units', 'inches', 'Position', [0.0, 1.0, 10, 6])
% tiledlayout(1,2)
% nexttile
% hold on
% title(file(1).filename, 'Interpreter', 'none')
% xlabel('\it{T}\rm{ (K)}')
% ylabel('\rm{H/M (mol/emu)}')
% errorbar(temperature1, suscRecip1, suscRecipErr1, 'o', 'MarkerFaceColor', 'w')
% box on
% hold off
% nexttile
% hold on
% title(file(2).filename, 'Interpreter', 'none')
% xlabel('\it{T}\rm{ (K)}')
% ylabel('\rm{M/H (emu/mol)}')
% errorbar(temperature2, susc2, suscErr2, 'o', 'MarkerFaceColor', 'w')
% box on
% hold off

% Plot chi*T/C
C = 4.62; % emu*K/mol
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.7])
hold on
%title(file(1).filename, 'Interpreter', 'none')
xlabel('1/\it{T}\rm{ (K^{-1})}')
ylabel('\rm{\chi*}\it{T}\rm{/}\it{C}')
errorbar(1./temperature1, susc1.*temperature1./C, suscErr1.*temperature1./C, 'o', 'MarkerFaceColor', 'w')
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\ChiTC.png';
%exportgraphics(gcf, fdir, 'resolution', 450)

%% Global fit using GlobalSearch to power law assuming alpha=alpha' and B=B' as in Ahlers 1974 EuO and 1975 paper

clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 7)

% First try fitting to 0 T data
mass = 2.34; % Sample mass in mg
molarMass = 599.46; % Sample grams per mole

fileLoc='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\Alpha\';
fileName1='IQM_IQMPPMS_20200218_1_TB_0-DRHCtry2.dat'; % J/g-K
fileName2='IQM_IQMPPMS_20200223_1_TB_0-DRHCextra0T.dat';
file1=readcell(strcat(fileLoc, fileName1), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', '\t');
file2=readcell(strcat(fileLoc, fileName2), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
ind1=[536, 662]; % 0 T
ind2=[1, length(file2)];

tmp = file1(ind1(1):ind1(2), 8);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
tempPre1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file1(ind1(1):ind1(2), 10);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
hcPre1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass; % Converting from J/g-K to J/mol-K
tmp = file1(ind1(1):ind1(2), 11);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
hcErrPre1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass; % Converting from J/g-K to J/mol-K
tmp = file1(ind1(1):ind1(2), 6);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
fieldPre1 = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

[tempPre2, hcPre2, hcErrPre2, fieldPre2] = importHC(file2, mass, molarMass, ind2);

tempPre = [tempPre1; tempPre2];
hcPre = [hcPre1; hcPre2];
hcErrPre = [hcErrPre1; hcErrPre2];
fieldPre = [fieldPre1; fieldPre2];

% Sort data from low-T to high-T
[temp, ind] = sort(tempPre);
hc = hcPre(ind);
hcErr = hcErrPre(ind);

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.95])
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
errorbar(temp, hc, hcErr, 'o', 'MarkerFaceColor', 'w')
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
box on
hold off

%tC0 = 2.83; % Estimated transition temperature in K for initial fitting value
tC0 = 2.85; % Rerunning fit to adjust plot
%x0 = [1, 1, 1, 0.05, 1, tC0]; % A, B, E, alpha, AP, tCrit
x0 = [2.0, 18, 0.005, -0.11, 1.38, tC0]; % A, B, E, alpha, AP, tCrit
fitBnds = [0.01, 0.05];
tR=@(tc, t) t./tc-1;
fitInd = @(tC, t) abs(tR(tC, t)) < fitBnds(2) & abs(tR(tC, t)) > fitBnds(1); % Indices for fitted range of data given some Tc and temperatures. E.g. use for temperatures within fitBnds of the critical temperature
model = @(x, t) ((x(1)/x(4)).*abs(tR(x(6), t(fitInd(x(6), t)))).^(-x(4))+x(2)+x(3).*tR(x(6), t(fitInd(x(6), t)))).*(t(fitInd(x(6), t)) > x(6))+((x(5)/x(4)).*abs(tR(x(6), t(fitInd(x(6), t)))).^(-x(4))+x(2)+x(3).*tR(x(6), t(fitInd(x(6), t)))).*(t(fitInd(x(6), t)) < x(6)); % A, B, E, alpha, AP, tCrit, t. Eq. (1).
modelFit = @(x, t) ((x(1)/x(4)).*abs(tR(x(6), t)).^(-x(4))+x(2)+x(3).*tR(x(6), t)).*(t > x(6))+((x(5)/x(4)).*abs(tR(x(6), t)).^(-x(4))+x(2)+x(3).*tR(x(6), t)).*(t < x(6)); % Model without fitting range excluded, so can model the hc anywhere after fitting.
tempCalPre = linspace(min(temp(fitInd(tC0, temp))), max(temp(fitInd(tC0, temp))), 5e2)';
hcCalPre = modelFit(x0, tempCalPre);

% figure(1)
% hold on
% plot(tempCalPre, hcCalPre, 'Color', 'g', 'LineWidth', 1)
% ylim([0 max(hc+hcErr)*1.1])
% hold off

% Comment out section to manually fit with identical alpha's
modelInput = @(x) model(x, temp);
gs = GlobalSearch('PlotFcn', @gsplotbestf);
gs.NumTrialPoints = 6e5;
gs.MaxTime = 390;
lb = [0, 0, 0, -1, 0, 2.82]; % A, B, E, alpha, AP, tCrit
ub = [5, 100, 10, 1, 5, 2.87]; % A, B, E, alpha, AP, tCrit
redChi2 = @(fnObs,fnErr,fnCal,nParam) sum((fnObs-fnCal).^2./fnErr.^2, 'all')/(length(fnObs)-nParam);
obj = @(x) redChi2(hc(fitInd(x(6), temp)),hcErr(fitInd(x(6), temp)),modelInput(x),length(x0));
opts = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','x0',x0,'objective',obj,'lb',lb,'ub',ub);
[xFit,fval] = run(gs,problem);
tempCal = linspace(min(temp(fitInd(xFit(6), temp))), max(temp(fitInd(xFit(6), temp))), 5e2)';
hcCal = modelFit(xFit, tempCal);

figure(1)
hold on
text(0.05, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
plot(tempCal, hcCal, 'Color', 'r', 'LineWidth', 1)
xline(min(temp(fitInd(xFit(6), temp))), '--k')
xline(max(temp(fitInd(xFit(6), temp))), '--k')
xline((fitBnds(1)+1)*xFit(6), '--k')
xline((-fitBnds(1)+1)*xFit(6), '--k')
ylim([0 max(hc+hcErr)*1.1])
axis square
hold off
nexttile
hold on
text(0.05, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(temp, hc, hcErr, 'o', 'MarkerFaceColor', 'w')
plot(tempCal, hcCal, 'Color', 'r', 'LineWidth', 1)
xline(min(temp(fitInd(xFit(6), temp))), '--k')
xline(max(temp(fitInd(xFit(6), temp))), '--k')
xline((fitBnds(1)+1)*xFit(6), '--k')
xline((-fitBnds(1)+1)*xFit(6), '--k')
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
ylim([min(hcCal)*0.9, max(hcCal)*1.1])
xlim([min(tempCal)*0.97, max(tempCal)*1.03])
errorbar(xFit(6), 13.5, 0.06, 'ok', 'horizontal', 'LineWidth', 1) % Plot 0.06 K temperature change
box on
axis square
hold off
drawnow

%% Globally fit thin plate 9/15 0 T data using GlobalSearch to power law assuming alpha=alpha' and B=B' as in Ahlers 1974 EuO and 1975 paper

clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 7)

factor = 1.07; % Factor to be applied to the plate mass
massP = 0.194; % Nominal sample mass in mg
massPCal = massP*factor;
molarMass = 599.46; % Sample grams per mole

fileLocP = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileNameP = '0p19mg_drhc_EuPd2S4_09152020_0T1T3T5T.dat';
fileP = readcell(strcat(fileLocP, fileNameP), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
indP = [1, length(fileP)];

[temp, hc, hcErr, ~] = importHC(fileP, massPCal, molarMass, indP);

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.95])
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
errorbar(temp, hc, hcErr, 'o', 'MarkerFaceColor', 'w')
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
box on
hold off

tC0 = 2.87; % Estimated transition temperature in K for initial fitting value
x0 = [1, 1, 1, 0.05, 1, tC0]; % A, B, E, alpha, AP, tCrit
fitBnds = [0.015, 0.04];
tR=@(tc, t) t./tc-1;
fitInd = @(tC, t) abs(tR(tC, t)) < fitBnds(2) & abs(tR(tC, t)) > fitBnds(1); % Indices for fitted range of data given some Tc and temperatures. E.g. use for temperatures within fitBnds of the critical temperature
model = @(x, t) ((x(1)/x(4)).*abs(tR(x(6), t(fitInd(x(6), t)))).^(-x(4))+x(2)+x(3).*tR(x(6), t(fitInd(x(6), t)))).*(t(fitInd(x(6), t)) > x(6))+((x(5)/x(4)).*abs(tR(x(6), t(fitInd(x(6), t)))).^(-x(4))+x(2)+x(3).*tR(x(6), t(fitInd(x(6), t)))).*(t(fitInd(x(6), t)) < x(6)); % A, B, E, alpha, AP, tCrit, t. Eq. (1).
modelFit = @(x, t) ((x(1)/x(4)).*abs(tR(x(6), t)).^(-x(4))+x(2)+x(3).*tR(x(6), t)).*(t > x(6))+((x(5)/x(4)).*abs(tR(x(6), t)).^(-x(4))+x(2)+x(3).*tR(x(6), t)).*(t < x(6)); % Model without fitting range excluded, so can model the hc anywhere after fitting.
tempCalPre = linspace(min(temp(fitInd(tC0, temp))), max(temp(fitInd(tC0, temp))), 5e2)';
hcCalPre = modelFit(x0, tempCalPre);

% figure(1)
% hold on
% plot(tempCalPre, hcCalPre, 'Color', 'g', 'LineWidth', 1)
% ylim([0 max(hc+hcErr)*1.1])
% hold off

% Comment out section to manually fit with identical alpha's
modelInput = @(x) model(x, temp);
gs = GlobalSearch('PlotFcn', @gsplotbestf);
gs.NumTrialPoints = 1e5;
gs.MaxTime = 60;
lb = [0, 0, 0, -1, 0, 2.82]; % A, B, E, alpha, AP, tCrit
ub = [5, 100, 10, 1, 5, 2.91]; % A, B, E, alpha, AP, tCrit
redChi2 = @(fnObs,fnErr,fnCal,nParam) sum((fnObs-fnCal).^2./fnErr.^2, 'all')/(length(fnObs)-nParam);
obj = @(x) redChi2(hc(fitInd(x(6), temp)),hcErr(fitInd(x(6), temp)),modelInput(x),length(x0));
opts = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','x0',x0,'objective',obj,'lb',lb,'ub',ub);
[xFit,fval] = run(gs,problem);
tempCal = linspace(min(temp(fitInd(xFit(6), temp))), max(temp(fitInd(xFit(6), temp))), 5e2)';
hcCal = modelFit(xFit, tempCal);

figure(1)
hold on
text(0.05, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
plot(tempCal, hcCal, 'Color', 'r', 'LineWidth', 1)
xline(min(temp(fitInd(xFit(6), temp))), '--k')
xline(max(temp(fitInd(xFit(6), temp))), '--k')
xline((fitBnds(1)+1)*xFit(6), '--k')
xline((-fitBnds(1)+1)*xFit(6), '--k')
ylim([0 max(hc+hcErr)*1.1])
axis square
hold off
nexttile
hold on
text(0.05, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(temp, hc, hcErr, 'o', 'MarkerFaceColor', 'w')
plot(tempCal, hcCal, 'Color', 'r', 'LineWidth', 1)
xline(min(temp(fitInd(xFit(6), temp))), '--k')
xline(max(temp(fitInd(xFit(6), temp))), '--k')
xline((fitBnds(1)+1)*xFit(6), '--k')
xline((-fitBnds(1)+1)*xFit(6), '--k')
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
ylim([min(hcCal)*0.9, max(hcCal)*1.1])
xlim([min(tempCal)*0.97, max(tempCal)*1.03])
box on
axis square
hold off
drawnow

%% Now run code from fitRedChi2Err.m to get an errorbar using the best-fit parameters. But use fmincon so bounds consistent with previous fitting, and can't get unphysical solution that weren't allowed in GlobalSearch with fmincon.

data = @(x) hc(fitInd(x, temp));
dataErr = @(x) hcErr(fitInd(x, temp));
modelInput = @(x) model(x, temp);
x0 = xFit;
xFitLoc = xFit';
errPts = 1e2;
fact = [0.2, 0.3, 0.0, 0.5, 0.15, 0.005]; % A, B, E, alpha, AP, tCrit
offset = [0.0, 20.0, 10.0, 0.2, 0.0, 0.0]; % A, B, E, alpha, AP, tCrit
redChi2Fit = fval;

paramTrial = nan(errPts, length(xFitLoc));
for i = 1:length(xFitLoc)
    paramTrial(:, i) = linspace(xFitLoc(i)*(1-fact(i))-offset(i), xFitLoc(i)*(1+fact(i))+offset(i), errPts);
    for j = 1:length(paramTrial(:, i)) % Make sure the iterated parameters obey the fitting bounds
        if paramTrial(j, i)<lb(i)
            paramTrial(j, i) = lb(i);
        elseif paramTrial(j, i)>ub(i)
            paramTrial(j, i) = ub(i);
        end
    end
end
xTrialFit = nan(errPts, length(xFitLoc)-1, length(xFitLoc)); % Each column is the parameter value after fitting, depth is going from one iterated parameter to the next while holding other constant.
chiTrial = nan(errPts, length(xFitLoc)); % Chi2 for every fit for the errorbar
xErr = nan(size(xFitLoc'));
slopes = nan(2, length(xFitLoc)); % First element is for the lower-parameter, second upper
intercepts = nan(2, length(xFitLoc));
interpPts = nan(4, length(xFitLoc));
paramLower = nan(1, length(xFitLoc));
paramUpper = nan(1, length(xFitLoc));
A = [];
b = [];
Aeq = [];
beq = [];
for j = 1:length(xFitLoc)
    initParam = xFitLoc';  % Fitting to one fewer parameter, only send the initial values for the fitted parameters. Can make this more efficient.
    lbLoc = lb;
    ubLoc = ub;
    initParam(j) = [];
    lbLoc(j) = [];
    ubLoc(j) = [];
    if j==1
        modelCalcErr = @(x, iterParam) modelInput([iterParam, x(j:length(xFitLoc)-1)]); % Notation is fun(var, x(1), x(2), x(3), x(4)), fun(x(1), var, x(2), x(3), x(4)), etc.
    elseif j<length(xFitLoc)
        modelCalcErr = @(x, iterParam) modelInput([x(1:j-1), iterParam, x(j:length(xFitLoc)-1)]);
    else
        modelCalcErr = @(x, iterParam) modelInput([x(1:j-1), iterParam]);
    end
    
    % Fit to the peak
    for k = 1:errPts
        modelInputErr = @(x) modelCalcErr(x, paramTrial(k, j));
        if j ~= 6
            obj = @(x) redChi2(data(x(5)), dataErr(x(5)), modelInputErr(x), length(x0)); % Using nParam=length(x0) rather than length(initParam) so gof's match when iterating parameter
        else % Plug in the fixed value of Tc I'm iterating over for the errorbar.
            obj = @(x) redChi2(data(paramTrial(k,j)), dataErr(paramTrial(k,j)), modelInputErr(x), length(x0)); % Using nParam=length(x0) rather than length(initParam) so gof's match when iterating parameter
        end
        [xTrialFit(k,:,j), chiTrial(k,j)] = fmincon(obj, initParam, A, b, Aeq, beq, lbLoc, ubLoc);
    end
    chiUpper = redChi2Fit*(1+1/(length(data(xFit(6)))-length(xFitLoc))); % Upper bound of the gof for determining the errorbar. Check this since dof changes and may be undefined.
    int1 = find(chiTrial(:, j)>chiUpper & paramTrial(:, j)<xFitLoc(j), 1, 'last');
    int2 = find(chiTrial(:, j)<chiUpper & paramTrial(:, j)<xFitLoc(j), 1, 'first');
    int3 = find(chiTrial(:, j)<chiUpper & paramTrial(:, j)>xFitLoc(j), 1, 'last');
    int4 = find(chiTrial(:, j)>chiUpper & paramTrial(:, j)>xFitLoc(j), 1, 'first');
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

figure('Units', 'normalized', 'Position', [0.0, 0.3, 1.0, 0.6])
tiledlayout(1, 6, 'TileSpacing', 'compact')
for j = 1:length(x0)
    nexttile(j)
    hold on
    ylabel('\chi^2_r')
    scatter(paramTrial(:, j), chiTrial(:, j))
    yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
    xlim([-inf inf])
    ylim([redChi2Fit-0.1 chiUpper+(chiUpper-redChi2Fit)])
    if ~any(isnan(interpPts(:, j))) % Check to see if an errorbar was determined before plotting
        plot([paramTrial(interpPts(1, j), j), paramTrial(interpPts(2, j), j)], [chiTrial(interpPts(1, j), j), chiTrial(interpPts(2, j), j)], 'b')
        plot([paramTrial(interpPts(3, j), j), paramTrial(interpPts(4, j), j)], [chiTrial(interpPts(3, j), j), chiTrial(interpPts(4, j), j)], 'b')
        plot([paramLower(j), paramUpper(j)], [chiUpper, chiUpper], 'k-.o', 'LineWidth', 2.0)
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['A: ', num2str(round(xFitLoc(1), 4)), '\pm', num2str(round(xErr(1), 4))])
xlabel('A')
hold off
nexttile(2)
hold on
title(['B: ', num2str(round(xFitLoc(2), 4)), '\pm', num2str(round(xErr(2), 4))])
xlabel('B')
hold off
nexttile(3)
hold on
title(['E: ', num2str(round(xFitLoc(3), 4)), '\pm', num2str(round(xErr(3), 4))])
xlabel('E')
nexttile(4)
hold on
title(['\alpha: ', num2str(round(xFitLoc(4), 4)), '\pm', num2str(round(xErr(4), 4))])
xlabel('\alpha')
nexttile(5)
hold on
title(['AP: ', num2str(round(xFitLoc(5), 4)), '\pm', num2str(round(xErr(5), 4))])
xlabel('AP')
hold off
nexttile(6)
hold on
title(['T_C: ', num2str(round(xFitLoc(6), 4)), '\pm', num2str(round(xErr(6), 4))])
xlabel('T_C')
hold off

p = (1-xFit(1)/xFit(5))/xFit(4);
pErr = sqrt(xErr(1)^2/xFit(5)^2/+xFit(1)^2*xErr(5)^2/xFit(5)^4+(1-xFit(1)/xFit(5))^2*xErr(4)^2/xFit(4)^2)/abs(xFit(4));

disp(['Universal Parameter: ', num2str(p)])
disp(['Universal Parameter Uncertainty: ', num2str(pErr)])
disp(['Reduced Chi-Squared: ', num2str(fval)])

% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\AlphaFit.png';
% figure(1)
% exportgraphics(gcf, fdir, 'resolution', 450)
% 
% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\AlphaErr.png';
% figure(3)
% exportgraphics(gcf, fdir, 'resolution', 450)
% 
% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\AlphaData.mat';
% save(fdir)

%% log(C-B) vs log(t) plot for comparison with log(C) vs log(t) that had found two different positive alpha's

% clear
% close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 7)

% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\AlphaData.mat';
% load(fdir)

% First plot above the transition, then plot below.

redTAbove = tR(xFit(6),temp(temp>xFit(6)));
redTBelow = abs(tR(xFit(6),temp(temp<xFit(6))));
hcAbove = hc(temp>xFit(6));
hcBelow = hc(temp<xFit(6));
hcErrAbove = hcErr(temp>xFit(6));
hcErrBelow = hcErr(temp<xFit(6));
redTCalAbove = tR(xFit(6),tempCal((tempCal>xFit(6)) & fitInd(xFit(6), tempCal)~=0)); % Get the reduced temperature of those model temperatures that only lie in the fitted range. This is where we expect a linear heat capacity.
redTCalBelow = abs(tR(xFit(6),tempCal(tempCal<xFit(6) & fitInd(xFit(6), tempCal)~=0))); % Get the reduced temperature of those model temperatures that only lie in the fitted range. This is where we expect a linear heat capacity.
hcCalAbove = hcCal(tempCal>xFit(6) & fitInd(xFit(6), tempCal)~=0); % Get the hc of those model temperatures that only lie in the fitted range. This is where we expect a linear heat capacity.
hcCalBelow = hcCal(tempCal<xFit(6) & fitInd(xFit(6), tempCal)~=0); % Get the hc of those model temperatures that only lie in the fitted range. This is where we expect a linear heat capacity.

% First show log(C)
figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.95])
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
text(0.05, 0.9, '\bf{a}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(redTAbove, hcAbove, hcErrAbove, 'o', 'MarkerFaceColor', 'w')
plot(redTCalAbove, hcCalAbove, 'Color', 'r', 'LineWidth', 1)
xlabel('\it{\tau}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
axis square
set(gca,'XScale','log','YScale','log')
box on
xlim([-inf inf])
ylim([-inf inf])
hold off
nexttile
hold on
text(0.05, 0.9, '\bf{b}', 'Units', 'normalized', 'fontsize', 9, 'fontname', 'Arial', 'interpreter', 'tex')
errorbar(redTBelow, hcBelow, hcErrBelow, 'o', 'MarkerFaceColor', 'w')
plot(redTCalBelow, hcCalBelow, 'Color', 'r', 'LineWidth', 1)
xlabel('\rm{|}\it{\tau}\rm{| (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
axis square
set(gca,'XScale','log','YScale','log')
box on
xlim([-inf inf])
ylim([-inf inf])
hold off

% Now plot log(C-B)
figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 3.11])
hold on
errorbar(redTAbove, xFit(2)-hcAbove, hcErrAbove, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'Color', 'b')
plot(redTCalAbove, xFit(2)-hcCalAbove, 'Color', 'r', 'LineWidth', 1)
errorbar(redTBelow, xFit(2)-hcBelow, hcErrBelow, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k')
plot(redTCalBelow, xFit(2)-hcCalBelow, 'Color', 'r', 'LineWidth', 1)
legend('\it{T>T_N}', '', '\it{T<T_N}', '', 'Location', 'northwest')
legend('boxoff')
xlabel('\rm{|}\it{\tau}\rm{| (K)}')
ylabel('\it{C-B}\rm{ (J / mol-K)}')
axis square
set(gca,'XScale','log','YScale','log')
box on
xlim([-inf inf])
ylim([-inf inf])
hold off

% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\AlphaLogLogC.png';
% figure(1)
% exportgraphics(gcf, fdir, 'resolution', 450)
% 
% fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\AlphaLogLog.png';
% figure(2)
% exportgraphics(gcf, fdir, 'resolution', 450)

%% Plot absorption correction as a function of 2theta for nuclear and magnetic peaks

clear
close all

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
radEff = (volSamp*3/4/pi)^(1/3); % Effective radius for spherical absorption correction in mm

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

%% Compare small plate to larger hc to get scale factor

clear
close all

% First try fitting to 0 T data
massL = 2.34; % Sample mass in mg
molarMass = 599.46; % Sample grams per mole

fileLocL='C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\Alpha\';
fileName1L='IQM_IQMPPMS_20200218_1_TB_0-DRHCtry2.dat'; % J/g-K
fileName2L='IQM_IQMPPMS_20200223_1_TB_0-DRHCextra0T.dat';
file1L=readcell(strcat(fileLocL, fileName1L), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', '\t');
file2L=readcell(strcat(fileLocL, fileName2L), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
ind1L=[536, 662]; % 0 T
ind2L=[1, length(file2L)];

tmp = file1L(ind1L(1):ind1L(2), 8);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
tempPre1L = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file1L(ind1L(1):ind1L(2), 10);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
hcPre1L = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass; % Converting from J/g-K to J/mol-K
tmp = file1L(ind1L(1):ind1L(2), 11);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
hcErrPre1L = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass; % Converting from J/g-K to J/mol-K
tmp = file1L(ind1L(1):ind1L(2), 6);
mask = cellfun(@ismissing, tmp);
tmp(mask) = {[]};
fieldPre1L = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T

[tempPre2L, hcPre2L, hcErrPre2L, fieldPre2L] = importHC(file2L, massL, molarMass, ind2L);

tempPreL = [tempPre1L; tempPre2L];
hcPreL = [hcPre1L; hcPre2L];
hcErrPreL = [hcErrPre1L; hcErrPre2L];
fieldPreL = [fieldPre1L; fieldPre2L];

% Sort data from low-T to high-T
[tempL, indL] = sort(tempPreL);
hcL = hcPreL(indL);
hcErrL = hcErrPreL(indL);

% Repeat for the small plate using 9/15 data since the 11/14 data had a
% quick 0 T T-scan

factor = 1.07; % Factor to be applied to the plate mass
massP = 0.194; % Nominal sample mass in mg
massPCal = massP*factor;

fileLocP = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\EuPd3S4\HeatCapacity091520\';
fileNameP = '0p19mg_drhc_EuPd2S4_09152020_0T1T3T5T.dat';
fileP = readcell(strcat(fileLocP, fileNameP), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');
indP = [1, length(fileP)];

[tempP, hcP, hcErrP, fieldP] = importHC(fileP, massP, molarMass, indP);
[tempPCal, hcPCal, hcErrPCal, fieldPCal] = importHC(fileP, massPCal, molarMass, indP);

% Plot

figure('Units', 'inches', 'Position', [0.0, 1.0, 6.5, 2.95])
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
errorbar(tempL, hcL, hcErrL, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'Color', 'b')
errorbar(tempP, hcP, hcErrP, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'Color', 'r')
legend(["Large","Small"])
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
title('Uncalibrated')
box on
hold off
nexttile
hold on
errorbar(tempL, hcL, hcErrL, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'Color', 'b')
errorbar(tempPCal, hcPCal, hcErrPCal, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'Color', 'r')
legend(["Large","Small"])
xlabel('\it{T}\rm{ (K)}')
ylabel('\it{C}\rm{ (J / mol-K)}')
title('Calibrated')
box on
hold off

%% Figure 2c second order susceptibility fitting. 09/06/2023 MPMS, 05/10/2019 powder sample, 74.93 mg. Tcw, muEff. No EuS.

clear
close all

set(0, 'defaultaxesfontsize', 18)

tempLB = 100; % Lower fitting limit in K
fdirMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\Susc.mat'; % Directory to save data
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\Susc.eps';

% Import data
molarMass = 599.46; % Grams per mole f.u.
mass = 74.93; % Sample mass in mg
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\EuPd3S4\';
file.filename = 'EuPd3S4_MvT_051019Powder_74p93mg_090623.dat';
file.data = readtable(strcat(fileLoc, file.filename), 'FileType', 'text', 'NumHeaderLines', 44, 'Delimiter', ','); % Sample
fileLocLa = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Susceptibility\LaPd3S4\';
fileLa.filename = 'LaPd3S4 mvT.xlsx';
fileLa.data = readtable(strcat(fileLocLa, fileLa.filename), 'FileType', 'spreadsheet'); % Offset to isolate Eu spin susceptibility

indEu = 1:length(file.data.Temperature_K_); % 500 Oe data
temp = file.data.Temperature_K_(indEu); % K
field = file.data.MagneticField_Oe_(indEu); % Oersted
moment = file.data.DCMomentFixedCtr_emu_(indEu)/(mass/1e3)*molarMass; % Convert to emu/mol
momentErr = file.data.DCMomentErrFixedCtr_emu_(indEu)/(mass/1e3)*molarMass; % Convert to emu/mol
[tempLaPre, indLa] = sort(fileLa.data.Temperature_K_); % K
fieldLaPre = fileLa.data.MagneticField_Oe_(indLa); % Oersted
suscPre = moment./field; % emu/mol*Oe
suscErrPre = momentErr./field; % emu/mol*Oe
suscLaPre = fileLa.data.moment_emu_mol_Oe_(indLa); % emu/mol*Oe
suscErrLaPre = fileLa.data.M_Std_Err__emu_(indLa).*fileLa.data.moment_emu_mol_Oe_(indLa)./fileLa.data.Moment_emu_(indLa); % Convert error from emu to emu/mol*Oe

% Interpolate and subtract the La data from Eu
[B, tempLa, ~] = groupcounts(tempLaPre); % Average susceptibilty at identical temperatures so interp1 works. B is number of occurance for each unique temperature tempLa.
suscLa = nan(size(tempLa));
suscErrLa = nan(size(tempLa));
for i = 1:length(tempLa)
    if B(i)>1 % Temperatures that are repeated
        suscLa(i) = mean(suscLaPre(tempLaPre==tempLa(i))); % Average susceptibilities at identical temperature
        suscErrLa(i) = sqrt(1/length(suscLaPre(tempLaPre==tempLa(i)))*sum(suscErrLaPre(tempLaPre==tempLa(i)).^2)); % Propagate errors of average. 01/12/2023 EuPd3S4 notes.
    else
        suscLa(i) = suscLaPre(tempLaPre==tempLa(i)); % Extract the susceptibility at a unique temperature
        suscErrLa(i) = suscErrLaPre(tempLaPre==tempLa(i));
    end
end
suscLaInterp = interp1(tempLa, suscLa, temp, 'linear', 'extrap');
suscLaInterp(end) = suscLa(end); % Extrapolation of last point is problematic so taking last La value. See 04/28/2023 EuPd3S4 notes.
suscErrLaInterp = nan(size(temp));
for i = 1:length(temp)
    dif = temp(i)-tempLa;
    indLow = length(dif(dif>0)); % Find the neighboring temperatures to the interpolated point. This is the index of the La temperature immediately below the interpolated point. indHigh would be indLow+1.
    if indLow>0 && indLow<length(suscErrLa) % Elseif's are in case there's only one neighboring temperature in the La data. Update for extrapolation if applied to new data.
        suscErrLaInterp(i) = sqrt( ((tempLa(indLow+1)-temp(i))/(tempLa(indLow+1)-tempLa(indLow)))^2*suscErrLa(indLow)^2 + ((temp(i)-tempLa(indLow))/(tempLa(indLow+1)-tempLa(indLow)))^2*suscErrLa(indLow+1)^2 ); % See 05/01/2023 EuPd3S4 notes.
    elseif indLow==0
        suscErrLaInterp(i) = suscErrLa(indLow+1);
    elseif indLow==length(suscErrLa)
        suscErrLaInterp(i) = suscErrLa(indLow);
    end
end
susc = suscPre-suscLaInterp;
suscErr = sqrt(suscErrPre.^2+suscErrLaInterp.^2);

indEuZFC = 1:length(file.data.Temperature_K_); % 500 Oe ZFC data
%indEuFC = 6918:7777; % 500 Oe FC data
suscZFC = susc(1:length(indEuZFC));
%suscFC = susc(length(indEuZFC)+1:end);
suscErrZFC = suscErr(1:length(indEuZFC));
%suscErrFC = suscErr(length(indEuZFC)+1:end);
tempZFC = temp(1:length(indEuZFC));
%tempFC = temp(length(indEuZFC)+1:end);

% Fit to ZFC data above tempLB
tempFit = tempZFC(tempZFC>tempLB);
suscFit = suscZFC(tempZFC>tempLB);
suscErrFit = suscErrZFC(tempZFC>tempLB);

% Wakshima 2001 fitting. 05/10/2023 EuPd3S4 notes.
NA = 6.022e23; % Avogadro's number
kB = 1.38e-16; % Boltzmann constant in cgs
muB = 9.274e-21; % Bohr magneton in cgs
aEq = @(lambda, T) lambda./T;
a = @(T) aEq(480, T); % Hinatsu 2001 is 519 K, Van Vleck 1932 is 365 K or 418 K. Van Vleck 1968 is 480 K. Takikawa 2010 is 471 K.
euThreeEq = @(T, x) 0.5.*NA.*muB.^2./(3.*kB.*a(T).*T).*((24 + (13.5.*a(T)-1.5).*exp(-a(T)) + (67.5.*a(T)-2.5).*exp(-3.*a(T)) + (189.*a(T)-3.5).*exp(-6.*a(T)))./(1 + 3.*exp(-a(T)) + 5.*exp(-3.*a(T)) + 7.*exp(-6.*a(T))));
euTwoEq = @(T, x) 0.5.*NA.*x(2).^2./(3*kB*T).*(1 + x(1)./T);
suscEq = @(T, x) euThreeEq(T, x) + euTwoEq(T, x); % x(1) is TCW, x(2) is muEff. First term is Van Vleck, second is high-T expansion of susceptibility to second order in 1/T
suscEqFit = @(x) suscEq(tempFit, x);
x0 = [-2.6, 7.94*muB];% muEff: Hinatsu 2001 is 7.58, theory 7.94
errPts = 1e3;
fact = [0, 0];
offset = [0.3, 0.03*muB];
[xFit, redChi2Fit, xErr, chiUpper, chiTrial, paramTrial, interpPts, slopes, intercepts, paramLower, paramUpper] = fitRedChi2Err(suscFit, suscErrFit, suscEqFit, x0, errPts, fact, offset);
tempModel = linspace(min(tempFit), max(tempFit), 5e2);
suscModel = suscEq(tempModel, xFit);
disp(['Tcw: ', num2str(xFit(1))])
disp(['Tcw Uncertainty: ', num2str(xErr(1))])
disp(['muEff: ', num2str(xFit(2)/muB)])
disp(['muEff Uncertainty: ', num2str(xErr(2)/muB)])
disp(['Reduced chi-squared: ', num2str(redChi2Fit)])
disp(' ')

% Model the 3+ and 2+ contributions
euThreeFit = euThreeEq(tempFit, xFit); % Fitted data
euThreeModel = euThreeEq(tempModel, xFit); % Model
euTwoFit = euTwoEq(tempFit, xFit); % Fitted data
euTwoModel = euTwoEq(tempModel, xFit); % Model

disp(['Tcw (K): ', num2str(xFit(1))])
disp(['Tcw Uncertainty (K): ', num2str(xErr(1))])
disp(['muEff (muB): ', num2str(xFit(2)/muB)])
disp(['muEff Uncertainty (muB): ', num2str(xErr(2)/muB)])
disp(['Reduced chi-squared: ', num2str(redChi2Fit)])
disp(' ')

% Plot the fit
figure(1)
ax1=gca;
hold on
xlabel('\it{T}\rm{ (K)}')
ylabel('\rm{1/\chi}\rm{ (mol\cdotemu^{-1})}')
%errorbar(temp, 1./susc, suscErr./susc.^2, 'o', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'CapSize', 8) % Uncertainty from 7/18/23 notes
scatter(temp, 1./susc, 'o', 'MarkerFaceColor', 'w')
p1=plot(tempModel, 1./suscModel, 'LineWidth', 2.0);
text(0.5, 0.1, '\it{H}\rm{ = 500 Oe}', 'Units', 'normalized', 'fontsize', 18, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
set(gca, 'LineWidth', 2.0, 'TickLength', [0.02, 0.01], 'XMinorTick', 'on', 'YMinorTick', 'on')
pbaspect([12 10 1])
ylim([-inf, 160])
box on
hold off
ax2=axes('Position', [.31 .52 .33 .33]);
hold on
xlabel('\rm{1/}\it{T}\rm{ (K^{-1})}')
ylabel('\rm{\chi\cdot}\it{T}\rm{ (emu\cdotK\cdotmol^{-1})}')
%errorbar(1./temp, susc.*temp, suscErr.*temp, 'o', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'CapSize', 8)
scatter(1./temp, susc.*temp, 'o', 'MarkerFaceColor', 'w')
plot(1./tempModel, suscModel.*tempModel, 'LineWidth', 2.0);
set(gca, 'LineWidth', 1.5, 'TickLength', [0.04, 0.01], 'XMinorTick', 'on', 'YMinorTick', 'on')
pbaspect([12 10 1])
xlim([0, 0.014])
ylim([3.5, 5])
box on
hold off
exportgraphics(gcf, fdir, 'ContentType', 'vector')

figure(2)
hold on
xlabel('\it{T}\rm{ (K)}')
ylabel('\rm{\chi}\rm{ (mol\cdotemu^{-1})}')
errorbar(temp, susc, suscErr, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'CapSize', 8) % Uncertainty from 7/18/23 notes
plot(tempModel, suscModel, 'g', 'LineWidth', 1.5)
%plot(tempModel, euThreeModel, 'm', 'LineWidth', 1.5)
text(0.54, 0.8, '\it{H}\rm{ = 500 Oe}', 'Units', 'normalized', 'fontsize', 14, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
set(gca, 'LineWidth', 2.0)
pbaspect([40 31 1])
box on
hold off

figure('Units', 'normalized', 'Position', [0.0, 0.3, 1.0, 0.6])
tiledlayout(1, 2, 'TileSpacing', 'compact')
for j = 1:length(x0)
    nexttile(j)
    hold on
    ylabel('\chi^2_r')
    scatter(paramTrial(:, j), chiTrial(:, j))
    yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
    xlim([-inf inf])
    ylim([0.99*redChi2Fit chiUpper+(chiUpper-redChi2Fit)])
    if ~any(isnan(interpPts(:, j))) % Check to see if an errorbar was determined before plotting
        plot([paramTrial(interpPts(1, j), j), paramTrial(interpPts(2, j), j)], [chiTrial(interpPts(1, j), j), chiTrial(interpPts(2, j), j)], 'b')
        plot([paramTrial(interpPts(3, j), j), paramTrial(interpPts(4, j), j)], [chiTrial(interpPts(3, j), j), chiTrial(interpPts(4, j), j)], 'b')
        plot([paramLower(j), paramUpper(j)], [chiUpper, chiUpper], 'k-.o', 'LineWidth', 2.0)
    end
    box on
    hold off
end
nexttile(1)
hold on
title(['\it{\Theta}\rm{_{CW}: ', num2str(round(xFit(1), 2)), '\pm', num2str(round(xErr(1), 2)), ' K}'])
xlabel('\it{\Theta}\rm{_{CW} (K)}')
hold off
nexttile(2)
hold on
title(['\mu\rm{_{Eff} (\mu_B): ', num2str(round(xFit(2)/muB, 2)), '\pm', num2str(round(xErr(2)/muB, 2)), '}'])
xlabel('\mu\rm{_{Eff} (\mu_B)}')
hold off

%save(fdirMat)

%% MPMS SI figure

clear
close all

set(0, 'defaultaxesfontsize', 14)

tempLB = 75; % Lower fitting limit in K
fdirMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\500OeMPMSSusc.mat'; % Directory to save data
fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\500OeMPMSSusc.eps';

% Import and plot 500 Oe ZFC, FC MPMS data
molarMass = 599.46; % Grams per mole f.u.
mass = 15.54; % Sample mass in mg
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\EuPd3S4\';
file.filename = 'EuPd3S4_TBe350_15p54mg_MPMS_070823.dat';
file.data = readtable(strcat(fileLoc, file.filename), 'FileType', 'text', 'NumHeaderLines', 44, 'Delimiter', ','); % Sample
fileLocLa = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Susceptibility\LaPd3S4\';
fileLa.filename = 'LaPd3S4 mvT.xlsx';
fileLa.data = readtable(strcat(fileLocLa, fileLa.filename), 'FileType', 'spreadsheet'); % Offset to isolate Eu spin susceptibility

indEu = 6059:7777; % 500 Oe data
temp = file.data.Temperature_K_(indEu); % K
field = file.data.MagneticField_Oe_(indEu); % Oersted
moment = file.data.DCMomentFixedCtr_emu_(indEu)/(mass/1e3)*molarMass; % Convert to emu/mol
momentErr = file.data.DCMomentErrFixedCtr_emu_(indEu)/(mass/1e3)*molarMass; % Convert to emu/mol
[tempLaPre, indLa] = sort(fileLa.data.Temperature_K_); % K
fieldLaPre = fileLa.data.MagneticField_Oe_(indLa); % Oersted
suscPre = moment./field; % emu/mol*Oe
suscErrPre = momentErr./field; % emu/mol*Oe
suscLaPre = fileLa.data.moment_emu_mol_Oe_(indLa); % emu/mol*Oe
suscErrLaPre = fileLa.data.M_Std_Err__emu_(indLa).*fileLa.data.moment_emu_mol_Oe_(indLa)./fileLa.data.Moment_emu_(indLa); % Convert error from emu to emu/mol*Oe

% Interpolate and subtract the La data from Eu
[B, tempLa, ~] = groupcounts(tempLaPre); % Average susceptibilty at identical temperatures so interp1 works. B is number of occurance for each unique temperature tempLa.
suscLa = nan(size(tempLa));
suscErrLa = nan(size(tempLa));
for i = 1:length(tempLa)
    if B(i)>1 % Temperatures that are repeated
        suscLa(i) = mean(suscLaPre(tempLaPre==tempLa(i))); % Average susceptibilities at identical temperature
        suscErrLa(i) = sqrt(1/length(suscLaPre(tempLaPre==tempLa(i)))*sum(suscErrLaPre(tempLaPre==tempLa(i)).^2)); % Propagate errors of average. 01/12/2023 EuPd3S4 notes.
    else
        suscLa(i) = suscLaPre(tempLaPre==tempLa(i)); % Extract the susceptibility at a unique temperature
        suscErrLa(i) = suscErrLaPre(tempLaPre==tempLa(i));
    end
end
suscLaInterp = interp1(tempLa, suscLa, temp, 'linear', 'extrap');
suscLaInterp(end) = suscLa(end); % Extrapolation of last point is problematic so taking last La value. See 04/28/2023 EuPd3S4 notes.
suscErrLaInterp = nan(size(temp));
for i = 1:length(temp)
    dif = temp(i)-tempLa;
    indLow = length(dif(dif>0)); % Find the neighboring temperatures to the interpolated point. This is the index of the La temperature immediately below the interpolated point. indHigh would be indLow+1.
    if indLow>0 && indLow<length(suscErrLa) % Elseif's are in case there's only one neighboring temperature in the La data. Update for extrapolation if applied to new data.
        suscErrLaInterp(i) = sqrt( ((tempLa(indLow+1)-temp(i))/(tempLa(indLow+1)-tempLa(indLow)))^2*suscErrLa(indLow)^2 + ((temp(i)-tempLa(indLow))/(tempLa(indLow+1)-tempLa(indLow)))^2*suscErrLa(indLow+1)^2 ); % See 05/01/2023 EuPd3S4 notes.
    elseif indLow==0
        suscErrLaInterp(i) = suscErrLa(indLow+1);
    elseif indLow==length(suscErrLa)
        suscErrLaInterp(i) = suscErrLa(indLow);
    end
end
susc = suscPre-suscLaInterp;
suscErr = sqrt(suscErrPre.^2+suscErrLaInterp.^2);

indEuZFC = 6059:6917; % 500 Oe ZFC data
indEuFC = 6918:7777; % 500 Oe FC data
suscZFC = susc(1:length(indEuZFC));
suscFC = susc(length(indEuZFC)+1:end);
suscErrZFC = suscErr(1:length(indEuZFC));
suscErrFC = suscErr(length(indEuZFC)+1:end);
tempZFC = temp(1:length(indEuZFC));
tempFC = temp(length(indEuZFC)+1:end);

figure(1)
hold on
xlabel('\it{T}\rm{ (K)}')
ylabel('\rm{\chi (emu\cdotmol^{-1})}')
%e1 = errorbar(tempZFC, suscZFC, suscErrZFC, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'CapSize', 6, 'Color', 'b');
%e2 = errorbar(tempFC, suscFC, suscErrFC, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'CapSize', 6, 'Color', 'r');
s1 = scatter(tempZFC, suscZFC, 'ob', 'MarkerFaceColor', 'w');
s2 = scatter(tempFC, suscFC, 'or', 'MarkerFaceColor', 'w');
text(0.6, 0.6, '\it{H}\rm{ = 500 Oe}', 'Units', 'normalized', 'fontsize', 14, 'fontname', 'Arial', 'interpreter', 'tex', 'Color', 'k')
%legend([e1, e2], {'ZFC', 'FC'}, 'Location', 'best')
legend([s1, s2], {'ZFC', 'FC'}, 'Location', 'best')
%legend('box', 'off')
pbaspect([40 31 1])
xlim([0, 40])
%ylim([4, 6.5])
box on
set(gca, 'TickLength', [0.02, 0.01])
hold off
exportgraphics(gcf, fdir, 'ContentType', 'vector')
save(fdirMat)

%% TOC graphic

% Order parameter
clear
close all

set(0, 'defaultaxesfontsize', 10)

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\EuPd3S4\Manuscript\TOCMATLAB.eps';

betaMat = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\EuPd3S4\betaColorPlot.mat';
load(betaMat)

figure('Units', 'inches', 'Position', [0.0, 1.0, 3.25, 1.75]) % Left, bottom, width, height
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
text(0.4, 0.1, '$\left( \frac{1}{2} \frac{1}{2} \frac{1}{2} \right)$', 'Units', 'normalized', 'fontsize', 8, 'fontname', 'Arial', 'interpreter', 'latex')
xlabel('\it{T} \rm{(K)}')
ylabel('\it{I} \rm{(Arbs)}')
errorbar(temp, (int-bgrBest)./max(int-bgrBest), intErr./max(int-bgrBest), 'o', 'MarkerFaceColor', 'w', 'MarkerSize', 3, 'CapSize', 3)
tempPlotPrelim = linspace(temp(1), temp(end), 400);
yPlotPrelim = bgrBest.*ones(size(tempPlotPrelim)) + (I0Best.*(1-tempPlotPrelim./TcBest).^(2.*betaBest)).*(tempPlotPrelim<=TcBest);
tempPlot = linspace(tempFit(1), tempFit(end), 400);
ind3 = tempPlotPrelim>=min(tempFit);
%plot(tempPlotPrelim(ind3), yPlotPrelim(ind3), 'r', 'LineWidth', 1.0)
%plot(tempPlotPrelim(~ind3), yPlotPrelim(~ind3), '--r', 'LineWidth', 1.0)
%xline(min(tempFit), '--', {'Lower Fitting Limit'})
xline(min(TcBest), '--', '\fontsize{8}\it{T}\rm{_N}', 'LabelHorizontalAlignment', 'center', 'LabelOrientation', 'horizontal')
xlim([0 4])
ylim([-0.1 1.2])
xticks([0 1 2 3 4])
set(gca, 'TickLength', [0.03, 0.01])
set(gca, 'FontName', 'Arial')
box on
%pbaspect([16 9 1])
axis square
hold off

% Add XRD plot

TXRD = [305, 315, 325, 335, 340, 350, 370];
IXRD = [2.75, 2.5, 1.81, 1.19, 0.00, 0.00, 0.00];
errXRD = [0.1, 0.1, 0.08, 0.03, 0.001, 0.001, 0.001];

nexttile
hold on
text(0.23, 0.1, '$\left( 100 \right)$', 'Units', 'normalized', 'fontsize', 8, 'fontname', 'Arial', 'interpreter', 'latex')
xlabel('\it{T} \rm{(K)}')
%ylabel('\it{I} \rm{(Arbs)}')
errorbar(TXRD, IXRD./max(IXRD), errXRD./max(IXRD), 'o', 'MarkerFaceColor', 'w')
xline(340, '--', '\fontsize{8}\it{T}\rm{_{CO}}', 'LabelHorizontalAlignment', 'center', 'LabelOrientation', 'horizontal')
xlim([300 380])
ylim([-0.1 1.2])
set(gca, 'YTickLabel', [])
set(gca, 'TickLength', [0.03, 0.01])
set(gca, 'FontName', 'Arial')
box on
%pbaspect([16 9 1])
axis square
hold off
exportgraphics(gcf, fdir, 'ContentType', 'vector')
%% Functions

function [temp, specHeat, specHeatErr, field] = importHC(file, mass, molarMass, ind)
    tmp = file(ind(1):ind(2), 8);
    mask = cellfun(@ismissing, tmp);
    tmp(mask) = {[]};
    temp = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
    tmp = file(ind(1):ind(2), 10);
    mask = cellfun(@ismissing, tmp);
    tmp(mask) = {[]};
    specHeat = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass/(mass/1e3)/(1e6); % Converting from muJ/K to J/mol-K
    tmp = file(ind(1):ind(2), 11);
    mask = cellfun(@ismissing, tmp);
    tmp(mask) = {[]};
    specHeatErr = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))*molarMass/(mass/1e3)/(1e6); % Converting from muJ/K to J/mol-K
    tmp = file(ind(1):ind(2), 6);
    mask = cellfun(@ismissing, tmp);
    tmp(mask) = {[]};
    field = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false))./1e4; % Converting from Oe to T
end
