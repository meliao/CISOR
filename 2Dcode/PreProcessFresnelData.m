% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

%%% This script loads Fresnel dataset, orgenizes it into a convenient form,
%%% and saves the result as a MAT file.
%%%
%%% U. S. Kamilov, MERL, 2017.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; home;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Speed of light [m/s]
c = 299792458;

%%% Set of measured frequencies [1/s]
frequencySet = [2, 3, 4, 5, 6, 7, 8, 9, 10];
numFrequencies = length(frequencySet);

%%% Radius of a rig where sensors are located
sensorRadius = 1.7840625; % (1.67) [m]

lambdaSet = c./frequencySet/1e9; % wavelength [m]
kbSet = 2*pi./lambdaSet; % wavenumber [1/m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderName = ['.' filesep...
    'data' filesep];

fileName = 'FoamTwinDielTM';
extension = '.exp';

fullFilePath = [folderName fileName extension];

startRow = 10;
startCol = 0;

fullData = dlmread(fullFilePath, '', startRow, startCol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling grid for the object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size of the reconstruction domain
Lx = 0.1536; % [m]
Ly = 0.1536; % [m]

% Number of pixels
Nx = 256;
Ny = 256;

% Sampling step
dx = Lx/Nx;
dy = Ly/Ny;

% Locations of the pixels
xPix = (-Nx/2:Nx/2-1)*dx;
yPix = (-Ny/2:Ny/2-1)*dy;

% Meshgrid the pixel locations
[XPix, YPix] = meshgrid(xPix, yPix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute locations of transmitters and receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Location of the transmissions

switch fileName
    case {'FoamDielIntTM', 'FoamDielIntTE','FoamDielExtTM', 'FoamDielExtTE'}
        transmitterAngles = (0:45:315)*pi/180; % [rad]
    case {'FoamTwinDielTM', 'FoamTwinDielTE'}
        transmitterAngles = (0:20:340)*pi/180; % [rad]
    otherwise
end

x_transmit = sensorRadius * cos(transmitterAngles); % [m]
y_transmit = sensorRadius * sin(transmitterAngles); % [m]
numTransmitters = numel(transmitterAngles);

% Location of the receivers
receiverAngles = (0:1:359)*pi/180; % [rad]
x_receive = sensorRadius * cos(receiverAngles);
y_receive = sensorRadius * sin(receiverAngles);
numReceivers = numel(receiverAngles);

% Compute receiver locations
numActiveReceivers = 241;
receiverIndicesSet = zeros(numTransmitters, numActiveReceivers, numFrequencies);
receiverMaskSet = false(numTransmitters, numReceivers, numFrequencies);

for indFreq = 1:numFrequencies

    frequency = frequencySet(indFreq);

    for ind = 1:numTransmitters
        I = fullData(:, 1) == ind & fullData(:, 3) == frequency;
        receiverIndicesSet(ind, :, indFreq) = fullData(I, 2);
        receiverMaskSet(ind, fullData(I, 2), indFreq) = true;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot measurement set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('Color', 'w', 'Name', 'Acquisition Set-up');
%
% for ind = 1:numTransmitters
%     I = receiverIndicesSet(ind, :);
%
%     subplot(3, 3, ind);
%     plot(x_transmit(ind), y_transmit(ind), '+', 'MarkerSize', 12);
%     hold on;
%     plot(x_receive(I), y_receive(I), '.');
%     plot(XPix(:), YPix(:), '.');
%     axis equal tight;
%     grid on;
%     xlim([-2, 2]);
%     ylim([-2, 2]);
%     title(sprintf('Tx %d', ind));
%     xlabel('x [m]');
%     ylabel('y [m]');
% end
%
% drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance data between sensors and pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transmitter to receiver
diff_x_tr = repmat(x_transmit(:), [1, numReceivers])-...
    repmat(x_receive, [numTransmitters, 1]);
diff_y_tr = repmat(y_transmit(:), [1, numReceivers])-...
    repmat(y_receive, [numTransmitters, 1]);
distanceTransToRec = sqrt(diff_x_tr.^2 + diff_y_tr.^2);

% Receiver to pixel
diff_x_rp = repmat(XPix, [1, 1, numReceivers])-...
    repmat(reshape(x_receive, [1, 1, numReceivers]), [Ny, Nx]);
diff_y_rp = repmat(YPix, [1, 1, numReceivers])-...
    repmat(reshape(y_receive, [1, 1, numReceivers]), [Ny, Nx]);
distanceRecToPix = sqrt(diff_x_rp.^2 + diff_y_rp.^2);

% Transmitter to pixel
diff_x_tp = repmat(XPix, [1, 1, numTransmitters])-...
    repmat(reshape(x_transmit, [1, 1, numTransmitters]), [Ny, Nx]);
diff_y_tp = repmat(YPix, [1, 1, numTransmitters])-...
    repmat(reshape(y_transmit, [1, 1, numTransmitters]), [Ny, Nx]);
distanceTransToPix = sqrt(diff_x_tp.^2 + diff_y_tp.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load the measured data
utotMeasSet = zeros(numTransmitters, numReceivers, numFrequencies);
uincMeasSet = zeros(numTransmitters, numReceivers, numFrequencies);

for indFreq = 1:numFrequencies

    frequency = frequencySet(indFreq);

    for ind = 1:numTransmitters

        % Extract data for current transmission and frequency
        I = fullData(:, 1) == ind & fullData(:, 3) == frequency;
        activeData = fullData(I, :);

        indRec = receiverIndicesSet(ind, :, indFreq);

        % Fill in data
        utotMeasSet(ind, indRec, indFreq)=...
            conj(activeData(:,4)+1j*activeData(:,5));
        uincMeasSet(ind, indRec, indFreq)=...
            conj(activeData(:,6)+1j*activeData(:,7));
    end
end

%%% Scattered data
uscatMeasSet = utotMeasSet-uincMeasSet;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibrate the input fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uincPredSet = zeros(size(uincMeasSet));
ampSet = zeros(numFrequencies);

for indFreq = 1:numFrequencies

    %%% Extract measured incident field
    uincMeas = uincMeasSet(:,:,indFreq);

    %%% Extract wavenumber in air
    kb = kbSet(indFreq);

    %%% Predicted input field
    uincPred = 1j*0.25*besselh(0, 1, kb*distanceTransToRec);

    %%% Find ratio in the middleswitch fileName
    amp = zeros(numTransmitters,1);
    switch fileName
        case {'FoamDielIntTM', 'FoamDielIntTE','FoamDielExtTM', 'FoamDielExtTE'}
            for indTrans = 1:numTransmitters
                amp(indTrans) = uincMeas(indTrans,mod(181+(indTrans-1)*45,360))...
                    /uincPred(indTrans,mod(181+(indTrans-1)*45,360));
            end
        case {'FoamTwinDielTM', 'FoamTwinDielTE'}
            for indTrans = 1:numTransmitters
                amp(indTrans) = uincMeas(indTrans,mod(181+(indTrans-1)*20,360))...
                    /uincPred(indTrans,mod(181+(indTrans-1)*20,360));
            end
        otherwise
    end

    %%% Find average
    amp = mean(amp);

    %%% Store
    ampSet(indFreq) = amp;

    %%% Scale the predicted input field
    uincPred = uincPred .* amp;

    %%% Store
    uincPredSet(:,:,indFreq) = uincPred;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure phase error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorField = uincMeasSet(receiverMaskSet)./uincPredSet(receiverMaskSet);
errorField = errorField./abs(errorField);

error = norm(1-errorField)^2;

fprintf('Total Error: %.2f\n', error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot calibrated fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w', 'Name', 'Input fields [Tx = 1]');

for indFreq = 1:numFrequencies

    %%% Active receivers
    I = receiverIndicesSet(1, :, indFreq);

    %%% Extract fields
    uincMeas = uincMeasSet(1, I, indFreq);
    uincPred = uincPredSet(1, I, indFreq);

    subplot(numFrequencies, 3, 3*indFreq-2);
    plot(I, abs(uincMeas), 'b.',...
        I, abs(uincPred), 'r.',...
        'LineWidth', 1.2);
    grid on;
    xlim([0, 359]);
    title(sprintf('Amplitude (%d GHz)', frequencySet(indFreq)));

    subplot(numFrequencies, 3, 3*indFreq-1);
    plot(I, angle(uincMeas), 'b.',...
        I, angle(uincPred), 'r.',...
        'LineWidth', 1.2);
    grid on;
    xlim([0, 359]);
    title(sprintf('Phase (%d GHz)', frequencySet(indFreq)));

    subplot(numFrequencies, 3, 3*indFreq);
    plot(I, abs(angle(uincMeas./uincPred)), 'k.',...
        'LineWidth', 1.2);
    grid on;
    xlim([0, 359]);
    title(sprintf('Phase Difference (%d GHz)', frequencySet(indFreq)));
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([fileName '.mat'], 'uincMeasSet', 'utotMeasSet', 'ampSet',...
    'receiverMaskSet', 'receiverIndicesSet');
