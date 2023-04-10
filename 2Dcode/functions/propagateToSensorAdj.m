% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later


function fhat = propagateToSensorAdj(uscat, uin, g, dx, dy)
%%% Adjoint of the propagation function
%%%
%%% Inputs:
%%% - uscat: (Nr x 1) scattered field
%%% - uin: (Ny x Nx) input field inside computational domain
%%% - g: (Ny x Nx x Nr) Green's function
%%% - dx,dy: sampling steps
%%%
%%% Outputs:
%%% - fhat: (Ny x Nx) scattering potential
%%%
%%% U. S. Kamilov, MERL, 2016

%%% dimensions of the problem
Nr = size(g,3);
Ny = size(uin,1);
Nx = size(uin,2);

%%% reorganize scattered field
uscat = reshape(uscat,[1,1,Nr]);
uscat = repmat(uscat,[Ny,Nx]);

%%% compute inner product
contSrc = dx*dy*sum(uscat.*conj(g),3);

%%% Multiply to adjoint field
fhat = conj(uin).*contSrc;
