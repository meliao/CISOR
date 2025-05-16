% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

%%% This script compares CISOR with first Born and iterative linearization
%%% for different contrast level using simulated data.
%%% Yanting Ma, MERL, 2017 July


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; home;

addpath functions

fprintf('Current directory: %s\n', pwd);

%%% Select reconstruction algorithm from

% 'FB_TV' : first Born approximation with TV regularization
% 'IL_TV' : iterative linearization with TV regularization
% 'CSI_TV': contrast source inversion with TV regularization
% 'CISOR_TV' : CISOR (our algorithm) with TV regularization
% 'SEAGLE_TV' : SEAGLE with TV regularization

algo = {'FB_TV','IL_TV','CSI_TV','CISOR_TV','SEAGLE_TV'};
algo = algo{5};

%%% Handle to the Hankel function
hankFun = @(x) 1j*0.25*besselh(0, 1, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Speed of light [m/s]
c = 299792458;

%%% Set of measured frequencies [1/s]
frequencySet = [10]; % in GHz
numFrequencies = length(frequencySet);

lambdaSet = c./frequencySet/1e9; % wavelength [m]
kbSet = 2*pi./lambdaSet; % wavenumber [1/m]

% contrast = 0.002:0.04:0.4;
contrast = 2.0;
numIter = 500;
alpha = 0.98;


%%% parameters for simulation

% stepSize = 1;

% if strcmp(algo,'CSI_TV')
%     a = 10;
%     stepSize = 0.5;
% end

% if strcmp(algo,'CISOR_TV')
%     alpha = 0.98;
% end

% if strcmp(algo,'SEAGLE_TV')
%     numIter = 10;
%     stepSize = 2;10;100;
% end


recSNRFinal = zeros(length(contrast),length(lamScale));
recrelL2Final = zeros(length(contrast),length(lamScale));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling grid for the object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size of the reconstruction domain
% Size of the reconstruction domain
Lx = 0.16; % [m]
Ly = 0.16; % [m]

% Number of pixels
Nx = 128;
Ny = 128;
% Sampling step
dx = Lx/Nx;
dy = Ly/Ny;


% Locations of the pixels
x = (-Nx/2:Nx/2-1)*dx;
y = (-Ny/2:Ny/2-1)*dy;

% Meshgrid the pixel locations
[XPix, YPix] = meshgrid(x, y);

% Save the final reconstructions here
ohats = zeros(length(contrast),length(lamScale),Ny,Nx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute locations of transmitters and receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Location of the transmissions
transmitterAngles = 0:45:315; % [degree]
transmitter_rad = 1.67;
numTransmitters = numel(transmitterAngles);
x_transmit = transmitter_rad * cosd(transmitterAngles); % [m]
y_transmit = transmitter_rad * sind(transmitterAngles); % [m]



% Location of the receivers
numReceivers = 360;
receiver_rad = 1.67;

x_receive = receiver_rad * cosd(1:360); % [m]
y_receive = receiver_rad * sind(1:360); % [m]

% y_receive(1:169) = (-84:84)*0.0384;
% y_receive = repmat(y_receive,[1,2]);


receiverMaskSet = ones(numTransmitters,numReceivers,numFrequencies);
% Print shape of receiverMaskSet
fprintf('Shape of receiverMaskSet: %d x %d x %d\n', size(receiverMaskSet, 1), size(receiverMaskSet, 2), size(receiverMaskSet, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance data between sensors and pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%% Define input fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uincDomSet = zeros(Ny, Nx, numTransmitters, numFrequencies);

for indFreq = 1:numFrequencies
    uincDom =  hankFun(kbSet(indFreq)*distanceTransToPix);
    uincDomSet(:,:,:,indFreq) = uincDom;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define sensor Green's functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sensorGreensFunctionSet = zeros(Ny, Nx, numReceivers, numFrequencies);

for indFreq = 1:numFrequencies
    sensorGreensFunction = hankFun(kbSet(indFreq)*distanceRecToPix);
    sensorGreensFunction = (kbSet(indFreq)^2)*sensorGreensFunction;
    sensorGreensFunctionSet(:,:,:,indFreq) = sensorGreensFunction;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define domain Green's functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locations of the Green's function's pixels
xGreen = (-Nx:Nx-1)*dx;
yGreen = (-Ny:Ny-1)*dy;

% Meshgrid the pixel locations
[XGreen, YGreen] = meshgrid(xGreen, yGreen);
R = sqrt(XGreen.^2+YGreen.^2);

domainGreensFunctionSet = zeros([size(R), numFrequencies]);

for indFreq = 1:numFrequencies
    
    %%% Extract wave-number
    kb = kbSet(indFreq);
    
    %%% generate Hankel function and remove singularity (singular at R=0)
    domainGreensFunction = hankFun(kb*R);
    domainGreensFunction(Ny+1,Nx+1) = quad2d(@(x,y) hankFun(kb*sqrt(x.^2+y.^2)),...
        -dx/2, dx/2, -dy/2, dy/2)/(dx*dy);
    domainGreensFunction = (kb^2)*domainGreensFunction;
    
    %%% Store
    domainGreensFunctionSet(:,:,indFreq) = domainGreensFunction;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate ground truth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for indContr = 1:length(contrast)
    
    o = gen_FoamDeilExt_phantom(XPix, YPix);
    
    % Show a plot of the object
    % figure(1); clf;
    % imagesc(x, y, o); axis image; colormap(parula);
    % colorbar;
    % title('Ground Truth');
    % xlabel('x [m]'); ylabel('y [m]');
    % set(gca, 'YDir', 'normal');
    % set(gca, 'FontSize', 16);
    % drawnow;
    
    % if strcmp(algo,'CISOR_TV')
    %     if indContr ==  1
    %         stepSize = 1.5*sqrt(1/contrast(indContr));
    %     else
    %         stepSize = 2*sqrt(1/contrast(indContr));
    %     end
    % else
    %     stepSize = sqrt(1/contrast(indContr));
    % end

    % if strcmp(algo, 'SEAGLE_TV')
    %     stepSize = 1.0;
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load Measurements from File
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This loads variables 'uincMeasSet', 'utotMeasSet', 'ampSet', 'receiverMaskSet', 'recieverIndicesSet'
    load('./FoamDielExtTM.mat');
    
    % For utotMeasSet, uincMeasSet, and receiverMaskSet, the last dimension is the frequency. We only want the frequency corresponding to index frequencySet(1).
    uincMeasSet = uincMeasSet(:,:,frequencySet(1)-1);
    utotMeasSet = utotMeasSet(:,:,frequencySet(1)-1);
    receiverMaskSet = receiverMaskSet(:,:,frequencySet(1)-1);
    
    uscatPredSet = utotMeasSet-uincMeasSet;
    
    % Print shape of receiverMaskSet
    fprintf('Shape of receiverMaskSet: %d x %d x %d\n', size(receiverMaskSet, 1), size(receiverMaskSet, 2), size(receiverMaskSet, 3));
    
    % utotDomSet = forwardProp(uincDomSet, o, domainGreensFunctionSet, uincDomSet, dx, dy);
    % uscatPredSet = fullPropagateToSensor(o, utotDomSet,...
    % sensorGreensFunctionSet, receiverMaskSet, dx, dy);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Reconstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = uscatPredSet;
    tol = 1e-4;
    lam = (5e-5)*0.5*norm(data(:))^2;
    lamAll = lam*lamScale;
    
    % parameters needed for plots
    plotRec.flag = 0; % 1 -- plot result at every iteration; 0 -- don't plot
    plotRec.Lx = Lx; plotRec.Ly = Ly;
    plotRec.emax = max(o(:));
    
    
    for indLam = 1:length(lamScale)
        lam = lamAll(indLam);
        filename = sprintf('../data/result_SEAGLE_TV_lam_%d_stepSize_%d.mat', lam, stepSize);
        switch algo
            case 'FB_TV'
                [ohat, outs] = firstBornTV(data,uincDomSet,...
                    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,...
                    dx,dy,numIter,plotRec,o,tol,lam,stepSize);
                recSNRFinal(indContr,indLam) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));
                save ContrastFB.mat  contrast numIter recSNRFinal
            case 'IL_TV'
                [ohat, outs] = iterativeLinearizationTV(data,uincDomSet,...
                    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,...
                    dx,dy,numIter,plotRec,o,tol,lam,stepSize);
                recSNRFinal(indContr,indLam) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));
                save ContrastIL.mat  contrast numIter recSNRFinal
            case 'CSI_TV'
                [ohat, outs] = contrastSourceInverseTV(data,uincDomSet,...
                    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,...
                    dx,dy,numIter,plotRec,o,tol,lam,stepSize,a);
                save ContrastCSI.mat  contrast numIter recSNRFinal
            case 'CISOR_TV'
                [ohat, outs] = cisorTV(data,uincDomSet,...
                    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,...
                    dx,dy,numIter,plotRec,alpha,o,tol,lam,stepSize);
                recSNRFinal(indContr,indLam) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));
                save ContrastCISOR.mat  contrast numIter recSNRFinal
            case 'SEAGLE_TV'
                [ohat, outs, relCost, tvCost, signalCost, times] = seagleTV(data,uincDomSet,...
                    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,...
                    dx,dy,numIter,plotRec,o,tol,lam,stepSize);
                recSNRFinal(indContr,indLam) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));
                recrelL2Final(indContr,indLam) = norm(o(:)-ohat(:))/norm(o(:));
                ohats(indContr,indLam,:,:) = ohat;
                save(filename,"ohat", "relCost", "tvCost", "signalCost", "times");
                disp(['Saving to ' filename])
            otherwise
                error('No such algorithm found!')
        end
        
    end
    
end
