% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later


function uscat = propagateToSensor(f, uin, g, dx, dy)
%%% Computes the scattered field at the sensor locations specified by the
%%% set of Green's functions.
%%%
%%% Input:
%%% - f: (Ny x Nx) scattering potential
%%% - uin: (Ny x Nx) input field
%%% - g: (Ny x Nx x Nr) Green's functions
%%% - dx, dy: sampling steps
%%%
%%% Ouput:
%%% - z: (Nr x 1) scattered field
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% number of transmissions
Nr = size(g, 3);

%%% contrast-source
contSrc = f .* uin;
contSrc = repmat(contSrc,[1,1,Nr]);

%%% compute propagated field
uscat = dx*dy*sum(sum(g.*contSrc,1),2); % measure
uscat = uscat(:);
