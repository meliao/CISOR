% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function [ohat, outs, relCost, tvCost, signalCost, times] = cisorTV(data,uincDomSet,...
    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,dx,dy,...
    numIter,plotRec,alpha,o,tol,lam,stepSize,ohat)

%%% Implementation of CISOR

%%% At each iteration, gradient of the following cost function w.r.t. f is computed
%%% exactly (without linearization)


%%% Solving the following optimization problem using FISTA:
% min_{f} 0.5*||y-Hdiag{utot}f)||^2+lam*TV(f)
% f: contrast
% utot: total field


%%% Input argument:
% data: y
% uincDomSet: uinc
% domainGreensFunctionSet: G
% sensorGreensFunctionSet: H

% receiverMaskSet: 1 -- active transmitter-receiver pair
%                  0 -- inactive transmitter-receiver pair
% dx, dy: sampling step
% numIter: maximum number of iterations
% plotRec: parameters needed for plotting reconstruction results
% alpha: parameter for relaxed FISTA
% o: groundtruth (only used for computing reconstruction SNR)
% tol: stopping criteria
% lam: regularization parameter
% stepSize: stepSize for relaxed FISTA

%%% Output argument:
% ohat: estimate of constrast f
% outs: store cost functions

% Yanting Ma @ MERL June 2017


[Ny,Nx,numTransmitters,numFrequencies] = size(uincDomSet);

% stepSize =  1;
stopCounter = 0;


relCost = zeros(numIter, 1);
signalCost = zeros(numIter, 1);
times = zeros(numIter, 1);
tvCost = zeros(numIter, 1);
totalCost = zeros(numIter, 1);
gradientNorm = zeros(numIter,1);
recSNR = zeros(numIter,1);

% ohat = zeros([Ny Nx]);
s = ohat;
q = 1;
u = zeros(size(uincDomSet));
v = zeros(size(uincDomSet));
P0 = zeros(Ny*Nx*2, 1);
Q0 = zeros(Ny*Nx*2, 1);

for indIter = 1:numIter

    tstart = tic;

    u = forwardProp(uincDomSet, s, domainGreensFunctionSet, u, dx, dy);
    dataPred = fullPropagateToSensor(s, u,...
        sensorGreensFunctionSet, receiverMaskSet, dx, dy);

    % Print shape of data and dataPred
    

    res = dataPred-data;
    HTres = H_adjoint_full(res, sensorGreensFunctionSet, receiverMaskSet, dx, dy);

    utotRes = repmat(s, [1,1,numTransmitters,numFrequencies]) .* HTres;

    v = backwardProp(utotRes, s, domainGreensFunctionSet, v, dx, dy);
    GTv = G_adjoint_full(v, domainGreensFunctionSet, dx, dy);

    grad = sum(sum(real(conj(u) .* (HTres+GTv)), 3), 4);

    ohatnext = s - stepSize*grad;

    optsTV.maxiter = 100;
    optsTV.bounds = [0, Inf];
    optsTV.P0 = P0;
    optsTV.Q0 = Q0;
    [ohatnext, P0, Q0] = denoiseTV(ohatnext, lam*stepSize, optsTV);

    qnext = 0.5*(1+sqrt(1+4*q*q));
    gradientNorm(indIter) = norm(s(:)-ohatnext(:));
    s = ohatnext + alpha*((q-1)/qnext)*(ohatnext-ohat);

    reldiff = norm(ohatnext(:)-ohat(:))/norm(ohat(:));

    q = qnext;

    ohat = ohatnext;

    times(indIter) = toc(tstart);

    %%% for computing cost only
    u_ohat = forwardProp(uincDomSet, ohat, domainGreensFunctionSet, u, dx, dy);
    dataPred_ohat = fullPropagateToSensor(ohat, u_ohat,...
        sensorGreensFunctionSet, receiverMaskSet, dx, dy);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relCost(indIter) = norm(dataPred_ohat(:)-data(:))/norm(data(:));
    tvCost(indIter) = tv_cost(ohat);
    signalCost(indIter) = norm(ohat-o)/norm(o);
    totalCost(indIter) = relCost(indIter) + tvCost(indIter);
    recSNR(indIter) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));



    outs.relCost = relCost;
    outs.tvCost = tvCost;
    outs.totalCost = totalCost;
    outs.gradientNorm = gradientNorm;
    outs.recSNR = recSNR;

    if indIter > 1
        if reldiff<tol
            stopCounter = stopCounter + 1;
        else
            stopCounter = 0;
        end

        if stopCounter>=10
            fprintf('Stopping criteria met. Stop iteration.\n');
            break;
        end
    end


    fprintf('indIter = %d, signalCost = %e, relCost = %e,  time=%e \n',...
        indIter, signalCost(indIter), relCost(indIter),  times(indIter));

    if plotRec.flag
        Lx = plotRec.Lx;
        Ly = plotRec.Ly;
        emax = plotRec.emax;

        figure(110); plot(recSNR(1:indIter)); drawnow;
        recSNR(indIter)

        figure(105);
        set(gcf, 'Name', sprintf('[%d/%d]', indIter, numIter));

        subplot(3, 4, 1:8);
%         imagesc([-Lx/2, Lx/2], [-Ly/2 Ly/2], ohat, [0, 1.1*emax]);
        imagesc([-Lx/2, Lx/2], [-Ly/2 Ly/2], ohat);
        hold off;
        axis equal tight xy;
        set(gca, 'FontSize', 16);
        colorbar;

        subplot(3, 4, 9);
        semilogy(1:indIter, relCost(1:indIter), 'b-',...
            indIter, relCost(indIter), 'bo',...
            'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Data Cost: %.1e', relCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 10);
        semilogy(1:indIter, tvCost(1:indIter), 'k-',...
            indIter, tvCost(indIter), 'ko', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('TV Cost: %.1e', tvCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 11);
        semilogy(1:indIter, totalCost(1:indIter), 'r-',...
            indIter, totalCost(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Total Cost: %.1e', totalCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 12);
        semilogy(1:indIter, gradientNorm(1:indIter), 'r-',...
            indIter, gradientNorm(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Gradient Norm: %.1e', gradientNorm(indIter)));
        set(gca, 'FontSize', 16);

        drawnow;
    end
end
