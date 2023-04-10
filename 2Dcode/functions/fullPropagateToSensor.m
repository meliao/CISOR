% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later


function uscatPredSet = fullPropagateToSensor(f, utotDomSet,...
    sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Propagate all the total fields to the sensors
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numTrans = size(utotDomSet, 3);
numFreq = size(utotDomSet, 4);
numRec = size(sensorGreensFunctionSet, 3);

%%% Initialize
uscatPredSet = zeros(numTrans, numRec, numFreq);

for indFreq = 1:numFreq

    sensorGreensFunction = sensorGreensFunctionSet(:,:,:,indFreq);

    for indTrans = 1:numTrans

        %%% Extract
        receiverMask = receiverMaskSet(indTrans,:,indFreq);
        utotDom = utotDomSet(:,:,indTrans,indFreq);

        %%% Predicted scattered field
        uscatPred = propagateToSensor(f, utotDom, sensorGreensFunction, dx, dy);
        uscatPred = receiverMask(:) .* uscatPred;

        %%% Store
        uscatPredSet(indTrans,:,indFreq) = uscatPred;

    end
end
