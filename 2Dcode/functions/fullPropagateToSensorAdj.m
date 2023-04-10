% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function f = fullPropagateToSensorAdj(uscatPredSet, utotDomSet,...
    sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Propagate all the scattered fields to computational domain
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numTrans = size(uscatPredSet, 1);
numFreq = size(uscatPredSet, 3);
Ny = size(sensorGreensFunctionSet, 1);
Nx = size(sensorGreensFunctionSet, 2);

%%% Initialize
f = zeros(Ny, Nx);

for indFreq = 1:numFreq

    sensorGreensFunction = sensorGreensFunctionSet(:,:,:,indFreq);

    for indTrans = 1:numTrans

        %%% Extract
        receiverMask = receiverMaskSet(indTrans,:,indFreq);
        utotDom = utotDomSet(:,:,indTrans,indFreq);
        uscatPred = uscatPredSet(indTrans, :, indFreq);

        %%% Predicted scattered field
        uscatPred = uscatPred .* receiverMask;
        contSrc = propagateToSensorAdj(uscatPred, utotDom, sensorGreensFunction, dx, dy);

        %%% Add contribution
        f = f + contSrc;
    end
end
