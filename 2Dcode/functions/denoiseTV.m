% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function [x, s, q, iter]=denoiseTV(y,lambda,opts)

%This method uses the TV regularizer and can be applied only to 2D data.
%x: denoised image
%iter: number of iterations for getting to the solution.


if ~isfield(opts, 'maxiter'), maxiter = 100; else, maxiter = opts.maxiter; end
if ~isfield(opts, 'tol'), tol = 1e-4; else, tol = opts.tol; end
if ~isfield(opts, 'bounds'), bounds = [-Inf, Inf]; else, bounds = opts.bounds; end

sigSize = size(y);
y = y(:);
Nx = sigSize(2);
Ny = sigSize(1);
N = 2*numel(y);

% q0 = zeros(N,1);
% s0 = zeros(N,1);
if ~isfield(opts, 'P0'), s0 = zeros(N,1); else, s0 = opts.P0(:); end
if ~isfield(opts, 'Q0'), q0 = zeros(N,1); else, q0 = opts.Q0(:); end

D = @(x) forwardTV(x,Nx,Ny);
Dt = @(x) adjointTV(x,Nx,Ny);

wx = 2*pi*(0:Nx-1)/Nx;
wy = 2*pi*(0:Ny-1)/Ny;

[Wx, Wy] = meshgrid(wx, wy);

freqResponseDTD =...
    (4 - exp(-1j*Wx) - exp(1j*Wx) - exp(-1j*Wy) - exp(1j*Wy));

% project onto bounds interval
y(y < bounds(1)) = bounds(1);
y(y > bounds(2)) = bounds(2);

% TV denoising
x = y;

q = q0;
s = s0;
mu = 1;
iter = 0;

Dx = D(x);
while norm(s - Dx)/norm(Dx) > tol && lambda/mu > 0
    iter = iter+1;
    % update z
    x = y + mu*Dt(s + q/mu);
    if isreal(x)
        x = real(ifft2(fft2(reshape(x,sigSize(1),sigSize(2))) ./ (1 + mu*freqResponseDTD)));
    else
        x = (ifft2(fft2(reshape(x,sigSize(1),sigSize(2))) ./ (1 + mu*freqResponseDTD)));
    end
    x = x(:);

    Dx = D(x);

    % update s
    s = TV_L21_shrink(Dx - q/mu, lambda/mu);

    % update q
    q = q + mu*(s - Dx);

    if iter>=maxiter
        break;
    end
end

x = reshape(x,sigSize);

end


function y = forwardTV(x,Nx,Ny)
    X = reshape(x, Ny, Nx);
    Dx = [X(:,2:end), X(:,1)] - X;
    Dy = [X(2:end,:); X(1,:)] - X;
    y = [Dx(:); Dy(:)];

end

function y = adjointTV(x,Nx,Ny)
    Dx = reshape(x(1:Nx*Ny),Ny,Nx);
    Dy = reshape(x(Nx*Ny+1:2*Nx*Ny), Ny, Nx);

    y = [Dy(end,:); Dy(1:end-1,:)]  - Dy + [Dx(:,end), Dx(:,1:end-1)] - Dx;

    y = y(:);
end

function [s] = TV_L21_shrink(a, gamma)
N = length(a);
% rearrange data
a = [a(1:N/2), a(N/2+1:N)];
s = bsxfun(@times, a, max(sqrt(sum(abs(a).^2,2)) - gamma, 0)./(sqrt(sum(abs(a).^2,2))+1e-6));
s = s(:);

end
