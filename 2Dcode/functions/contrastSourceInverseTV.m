% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function [ohat, outs] = contrastSourceInverseTV(data,uincDomSet,...
    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,dx,dy,...
    numIter,plotRec,o,tol,lam,stepSize,a)

%%% Solving the following optimization problem:
% min_{f,w} 0.5*lam1*||y-Hw||^2+0.5*lam2*||f*uinc-w+diag{f}Gw||^2+lam*TV(f)
% f: contrast
% w: (=f*utot) contrst source
% Alternating between the estimation of f and w.

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

[Ny,Nx,numTransmitters,numFrequencies] = size(uincDomSet);
lam = lam/(norm(data(:))^2);
lam1 = 1/(norm(data(:))^2);
lam2 = a/(norm(uincDomSet(:))^2);
stopCounter = 0;

ostepSize = stepSize;
wStepSize = stepSize;


numIter1 = 100;
numIter2 = 100;


relCost1 = zeros(numIter, 1);
relCost2 = zeros(numIter, 1);
relCost = zeros(numIter, 1);
tvCost = zeros(numIter, 1);
totalCost = zeros(numIter, 1);
recSNR = zeros(numIter, 1);


ohat = zeros([Ny Nx]);
P = zeros([Ny Nx 2]);

w = zeros(Ny,Nx,numTransmitters,numFrequencies);

P0 = zeros(Ny*Nx*2,1);
Q0 = zeros(Ny*Nx*2,1);
for indIter = 1:numIter

    a = a*1.1;

    % sub-problem 1: estimate contrast source w=utot*f
    qw = 1;
    sw = w;
    for indIter1 = 1:numIter1

        % compute gradient of 1st term in cost: H*(Hw - y)
        Hw = H_forward_full(sw, sensorGreensFunctionSet, receiverMaskSet, dx, dy);
        resw1 = Hw - data;
        sensorGradientSet = H_adjoint_full(resw1, sensorGreensFunctionSet, receiverMaskSet, dx, dy);

        % compute gradient of 2nd term in cost: B*(Bw-uinc*f), B=I-diag{f}G
        Bw = B_forward_full(sw, ohat, domainGreensFunctionSet, dx,dy);
        resw2 = Bw - uincDomSet.*repmat(ohat,[1,1,numTransmitters,numFrequencies]);
        domainGradientSet = B_adjoint_full(resw2, ohat, domainGreensFunctionSet,dx,dy);


        % full gradient for sub-problem 1:
        gradsw = lam1*sensorGradientSet + lam2*domainGradientSet;
        if norm(gradsw(:))<1e-6
            break
        end

        wnext = sw - wStepSize * gradsw;

        qwnext = 0.5*(1+sqrt(1+4*qw*qw));
        sw = wnext + ((qw-1)/qwnext)*(wnext-w);

        qw = qwnext;
        w = wnext;
    end

    % sub-problem 2: estimate contrast f
    q = 1;
    s = ohat;
    for indIter2 = 1:numIter2
        s_full = repmat(s,[1,1,numTransmitters,numFrequencies]);
        reso = s_full.*G_forward_full(w,domainGreensFunctionSet, dx, dy)...
            + uincDomSet.*s_full - w;
        grads = conj(uincDomSet).*reso + conj(w).*G_adjoint_full(reso,domainGreensFunctionSet,dx,dy);
        grads = real(grads);
        grads = sum(sum(grads,4),3);

        ohatnext = s - ostepSize*grads;
        optsTV.maxiter = 100;
        optsTV.bounds = [0, Inf];
        optsTV.P0 = P0;
        optsTV.Q0 = Q0;
        [ohatnext, P0, Q0] = denoiseTV(ohatnext, lam*stepSize, optsTV);

        if norm(ohatnext(:)-ohat(:))/norm(ohat(:))<1e-4
            break
        end

        qnext = 0.5*(1+sqrt(1+4*q*q));
        s = ohatnext + ((q-1)/qnext)*(ohatnext-ohat);

        q = qnext;
        ohat = ohatnext;
    end

    relCost1(indIter) = 0.5*lam1*norm(resw1(:))^2;
    relCost2(indIter) = 0.5*lam2*norm(resw2(:))^2;

    relCost(indIter) = relCost1(indIter) + relCost2(indIter);
    tvCost(indIter) = lam*tv_cost(ohat);
    totalCost(indIter) = relCost(indIter) + tvCost(indIter);
    recSNR(indIter) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));

    outs.relCost1 = relCost1;
    outs.relCost2 = relCost2;
    outs.relCost = relCost;
    outs.tvCost = tvCost;
    outs.totalCost = totalCost;
    outs.recSNR = recSNR;

    if indIter > 1
        if abs(totalCost(indIter-1)-totalCost(indIter))/totalCost(indIter-1)<tol
            stopCounter = stopCounter + 1;
        else
            stopCounter = 0;
        end

        if stopCounter>=10
            break;
        end
    end

    fprintf('indIter = %d, totalCost = %e\n',indIter, totalCost(indIter));


    if plotRec.flag
        Lx = plotRec.Lx;
        Ly = plotRec.Ly;
        emax = plotRec.emax;

        figure(105);
        set(gcf, 'Name', sprintf('[%d/%d]', indIter, numIter));

        subplot(3, 4, 1:8);
        imagesc([-Lx/2, Lx/2], [-Ly/2 Ly/2], ohat, [0, 1.1*emax]);
        axis equal tight xy;
        set(gca, 'FontSize', 16);
%         colormap gray;
        colorbar;

        subplot(3, 4, 9);
        semilogy(1:indIter, relCost1(1:indIter), 'b-',...
            indIter, relCost1(indIter), 'bo',...
            'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Data Cost1: %.1e', relCost1(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 10);
        semilogy(1:indIter, relCost2(1:indIter), 'b-',...
            indIter, relCost2(indIter), 'bo',...
            'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Data Cost2: %.1e', relCost2(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 11);
        semilogy(1:indIter, tvCost(1:indIter), 'k-',...
            indIter, tvCost(indIter), 'ko', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('TV Cost: %.1e', tvCost(indIter)));
        set(gca, 'FontSize', 16);

        subplot(3, 4, 12);
        semilogy(1:indIter, totalCost(1:indIter), 'r-',...
            indIter, totalCost(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('Total Cost: %.1e', totalCost(indIter)));
        set(gca, 'FontSize', 16);
        drawnow;
    end
end
end

function u = B_adjoint(z, f, g, dx, dy)
%%% Adjoint of operator B = I - diag(f)*G
%%%
%%% Input:
%%% - z: (Ny x Nx) field
%%% - f: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) Green's function
%%%
%%% Output:
%%% - u: (Ny x Nx) field

assert(all(size(g)==2*size(z)), 'all(size(g)==2*size(z))');
assert(all(size(f)==size(z)), 'all(size(f)==size(z))');

u = z - dx*dy*conv2DAdj(z.*conj(f), g);

end

function wSet = B_adjoint_full(zSet, f, domainGreensFunctionSet, dx, dy)


[~,~,numTransmitters,numFrequencies] = size(zSet);
wSet = zeros(size(zSet));

for indFreq = 1:numFrequencies

    domainGreensFunction = domainGreensFunctionSet(:,:,indFreq);

    for indTrans = 1:numTransmitters

        z = zSet(:,:,indTrans,indFreq);
        w = B_adjoint(z, f, domainGreensFunction, dx, dy);
        wSet(:,:,indTrans,indFreq) = w;
    end
end

end


function z = B_forward(w, f, g, dx, dy)
%%% Operator B = I - diag(f)*G
%%%
%%% Input:
%%% - w: (Ny x Nx) contrast source
%%% - f: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) Green's function
%%%
%%% Output:
%%% - z: (Ny x Nx) field

assert(all(size(g)==2*size(w)), 'all(size(g)==2*size(w))');
assert(all(size(f)==size(w)), 'all(size(f)==size(w))');

z = w - dx*dy*f.*conv2D(w, g);

end

function zSet = B_forward_full(wSet, f, domainGreensFunctionSet, dx, dy)


[~,~,numTransmitters,numFrequencies] = size(wSet);
zSet = zeros(size(wSet));

for indFreq = 1:numFrequencies

    domainGreensFunction = domainGreensFunctionSet(:,:,indFreq);

    for indTrans = 1:numTransmitters

        w = wSet(:,:,indTrans,indFreq);
        z = B_forward(w, f, domainGreensFunction, dx, dy);
        zSet(:,:,indTrans,indFreq) = z;
    end
end
end
