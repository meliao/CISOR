% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later


function [ohat, outs] = iterativeLinearizationTV(data,uincDomSet,...
    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,dx,dy,...
    numIter,plotRec,o,tol,lam,stepSize)

%%% At each iteration, in the forward problem, total field is computed using the current
%%% estimate of contrast f, then in the inverse problem, total field is
%%% assumed to be fixed, hence resulting in a linear inverse problem


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

% o: ground-truth (only used for computing reconstructed SNR)
% tol: stopping criteria
% lam: regularization parameter
% stepSize: stepSize for FISTA

%%% Output argument:
% ohat: estimate of constrast f
% outs: store cost functions

% Yanting Ma @ MERL June 2017

[Ny,Nx,~,~] = size(uincDomSet);
% stepSize =  1;
stopCounter = 0;


relCost = zeros(numIter, 1);
tvCost = zeros(numIter, 1);
totalCost = zeros(numIter, 1);
recSNR = zeros(numIter,1);

numInnerIter = 10;
ohat = zeros([Ny Nx]);
u = zeros(size(uincDomSet));
s = ohat;
P = zeros([Ny Nx 2]);

P0 = zeros(Ny*Nx*2,1);
Q0 = zeros(Ny*Nx*2,1);
for indIter = 1:numIter

    u = forwardProp(uincDomSet, ohat, domainGreensFunctionSet, u, dx, dy);

    Hdiagu = @(x)  fullPropagateToSensor(x, u,sensorGreensFunctionSet, receiverMaskSet, dx, dy);
    HdiaguT = @(z) fullPropagateToSensorAdj(z, u, sensorGreensFunctionSet, receiverMaskSet, dx, dy);

    q = 1;
    for innerIter = 1:numInnerIter
        grad = real(HdiaguT(Hdiagu(s) - data));

        ohatnext = s - stepSize*grad;
        optsTV.maxiter = 100;
        optsTV.bounds = [0, Inf];
        optsTV.P0 = P0;
        optsTV.Q0 = Q0;
        [ohatnext, P0, Q0] = denoiseTV(ohatnext, lam*stepSize, optsTV);

        qnext = 0.5*(1+sqrt(1+4*q*q));
        s = ohatnext + ((q-1)/qnext)*(ohatnext-ohat);

        q = qnext;

        if norm(ohatnext - ohat)/norm(ohatnext)<1e-6
            break;
        end

        ohat = ohatnext;
    end

    uscatPredSet = fullPropagateToSensor(ohat, u ,sensorGreensFunctionSet, receiverMaskSet, dx, dy);

    relCost(indIter) = 0.5*norm(uscatPredSet(:)-data(:))^2/norm(data(:))^2;
    tvCost(indIter) = lam*tv_cost(ohat);
    totalCost(indIter) = relCost(indIter) + tvCost(indIter);
    recSNR(indIter) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));

    outs.relCost = relCost;
    outs.tvCost = tvCost;
    outs.totalCost = totalCost;
    outs.recSNR = recSNR;


    if indIter > 1
        if abs(totalCost(indIter-1)-totalCost(indIter))/totalCost(indIter-1)<tol ||...
                totalCost(indIter-1) - totalCost(indIter)<0
            stopCounter = stopCounter + 1;
        else
            stopCounter = 0;
        end

        if stopCounter>=10
            break;
        end
    end

%     fprintf('indIter = %d, totalCost = %e\n',indIter, totalCost(indIter));

    if plotRec.flag
        Lx = plotRec.Lx;
        Ly = plotRec.Ly;
        emax = plotRec.emax;
        figure(105);
        set(gcf, 'Name', sprintf('[%d/%d]', indIter, numIter));

        subplot(3, 3, 1:6);
        imagesc([-Lx/2, Lx/2], [-Ly/2 Ly/2], ohat, [0, 1.1*emax]);
        axis equal tight xy;
        set(gca, 'FontSize', 16);
        colorbar;

        subplot(3, 3, 7);
        semilogy(1:indIter, relCost(1:indIter), 'b-',...
            indIter, relCost(indIter), 'bo',...
            'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Data Cost: %.1e', relCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 3, 8);
        semilogy(1:indIter, tvCost(1:indIter), 'k-',...
            indIter, tvCost(indIter), 'ko', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('TV Cost: %.1e', tvCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 3, 9);
        semilogy(1:indIter, totalCost(1:indIter), 'r-',...
            indIter, totalCost(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Total Cost: %.1e', totalCost(indIter)));
        set(gca, 'FontSize', 16);

        drawnow;
    end
end
