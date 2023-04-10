% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function TVnorm = tv_cost(f)

[Ny,Nx] = size(f);
X = reshape(f, Ny, Nx);
fx = [X(:,2:end), X(:,1)] - X;
fy = [X(2:end,:); X(1,:)] - X;

TVf=sqrt(abs(fx).^2+abs(fy).^2);% Amplitude of the gradient vector

TVnorm=sum(TVf(:));
