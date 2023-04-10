% Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later


function [ohat, outs] = seagleTV(data,uincDomSet,...
    domainGreensFunctionSet,sensorGreensFunctionSet,receiverMaskSet,dx,dy,...
    numIter,plotRec,o,tol,lam,stepSize,fileName)

%%% This implementation is only for single frequency case, need to modify
%%% the code for multiple frequencies.

%%% Implementation of SEAGLE
%%% Yanting Ma, MERL, 2017 August

[Ny,Nx,numTransmitters,numFrequencies] = size(uincDomSet);
K = 200;
tol_seagle = norm(uincDomSet(:))^2*5*1e-7;
stopCounter = 0;

relCost = zeros(numIter, 1);
tvCost = zeros(numIter, 1);
totalCost = zeros(numIter, 1);
gradientNorm = zeros(numIter,1);
recSNR = zeros(numIter,1);

ohat = zeros([Ny Nx]);
s = ohat;
ohatnext = ohat;
q = 1;
P = zeros([Ny Nx 2]);
reductionRate = 0.99;

P0 = zeros(Ny*Nx*2,1);
Q0 = zeros(Ny*Nx*2,1);
for indIter = 1:numIter

    outsForward = seagle_forward(s,domainGreensFunctionSet,sensorGreensFunctionSet,...
        receiverMaskSet,K,tol_seagle,uincDomSet,dx,dy);
    grad = seagle_backprop(s,outsForward,domainGreensFunctionSet,sensorGreensFunctionSet,...
        receiverMaskSet,outsForward.K,data,uincDomSet,dx,dy);

    stepSize = stepSize/reductionRate;
    if stepSize < 1e-6
        break;
    end
    z = outsForward.z;
    dataCost_s = 0.5*norm(data(:)-z(:))^2;
    keepSearch = 1;
    while keepSearch
        stepSize = stepSize*reductionRate
        if stepSize < 1e-2
            break
        end
        Q =@(f) dataCost_s + sum(sum(grad.*(f-s))) + 1/(2*stepSize)*sum(sum((s-f).^2));
        outsForward = seagle_forward(ohatnext,domainGreensFunctionSet,sensorGreensFunctionSet,...
        receiverMaskSet,K,tol_seagle,uincDomSet,dx,dy);
        z = outsForward.z;
        dataCost_ohat = 0.5*norm(data(:)-z(:))^2;

        ohatnext = s - stepSize*grad;
        optsTV.maxiter = 100;
        optsTV.bounds = [0, Inf];
        optsTV.P0 = P0;
        optsTV.Q0 = Q0;
        [ohatnext, P0, Q0] = denoiseTV(ohatnext, lam*stepSize, optsTV);

        if dataCost_ohat < Q(ohatnext)
            keepSearch = 0;
        end

    end

    qnext = 0.5*(1+sqrt(1+4*q*q));
    gradientNorm(indIter) = norm(s(:)-ohatnext(:));
    s = ohatnext + ((q-1)/qnext)*(ohatnext-ohat);

    q = qnext;

    ohat = ohatnext;



    %%% for computing cost only
    outsForward = seagle_forward(ohat,domainGreensFunctionSet,sensorGreensFunctionSet,...
        receiverMaskSet,K,tol_seagle,uincDomSet,dx,dy);
    dataPred_ohat = outsForward.z;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relCost(indIter) = norm(dataPred_ohat(:)-data(:))/norm(data(:));
    tvCost(indIter) = lam*tv_cost(ohat)/(0.5*norm(data(:)));
    totalCost(indIter) = relCost(indIter) + tvCost(indIter);
    recSNR(indIter) = 20*log10(norm(o(:))/norm(ohat(:)-o(:)));

    outs.relCost = relCost;
    outs.tvCost = tvCost;
    outs.totalCost = totalCost;
    outs.gradientNorm = gradientNorm;
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


    fprintf('indIter = %d, totalCost = %e, K=%d, stepSize = %4.4f \n',...
        indIter, totalCost(indIter),outsForward.K,stepSize);

    if plotRec.flag
        Lx = plotRec.Lx;
        Ly = plotRec.Ly;
        emax = plotRec.emax;

        figure(110); plot(recSNR(1:indIter)); drawnow;
        recSNR(indIter)

        figure(105);
        set(gcf, 'Name', sprintf('[%d/%d]', indIter, numIter));

        subplot(3, 4, 1:8);
        imagesc([-Lx/2, Lx/2], [-Ly/2 Ly/2], ohat, [0, 1.1*emax]);
%         hold on;
%         switch fileName
%             case {'FoamDielIntTM', 'FoamDielIntTE'}
%                 viscircles([0, 0], 0.08/2,'LineStyle','-','EdgeColor','b');
%                 viscircles([-0.005, 0], 0.031/2,'LineStyle','-','EdgeColor','r');
%             case {'FoamDielExtTM', 'FoamDielExtTE'}
%                 viscircles([0, 0], 0.08/2,'LineStyle','-','EdgeColor','b');
%                 viscircles([-0.0555, 0], 0.031/2,'LineStyle','-','EdgeColor','r');
%         end
        hold off;
        axis equal tight xy;
        set(gca, 'FontSize', 16);
%         colormap gray;
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
end

function outs = seagle_forward(f,domainGreensFunctionSet,sensorGreensFunctionSet,...
    receiverMaskSet,K,tol,uinc,dx,dy)

[Ny,Nx,numTransmitters] = size(uinc);
gammar = zeros(K,1);
mu = zeros(K,1);
s = zeros(Ny,Nx,numTransmitters,K);


%%% initialization
t = 0;
u = uinc;
uprev = uinc;
unext = uinc;

%%% Scattering operators
A = @(u) A_forward(u, f, domainGreensFunctionSet, dx, dy);
AT = @(z) A_adjoint(z, f, domainGreensFunctionSet, dx, dy);


for k = 1:K
    tnext = (1+sqrt(1 + 4*t^2))/2;
    mu(k) = (1-t)/tnext;
    s(:,:,:,k) = (1-mu(k))*u + mu(k)*uprev;
    As = zeros(size(uinc));
    for indTrans = 1:numTransmitters
        As(:,:,indTrans) = A(s(:,:,indTrans,k));
    end
    res = As - uinc;
    g = zeros(size(uinc));
    for indTrans = 1:numTransmitters
        g(:,:,indTrans) = AT(res(:,:,indTrans));
    end

    if norm(g(:)) < tol
        K = k-1;
        gammar = gammar(1:K);
        s = s(:,:,:,1:K);
        mu = mu(1:K);
        break
    end

    Ag = zeros(size(uinc));
    for indTrans = 1:numTransmitters
        Ag(:,:,indTrans) = A(g(:,:,indTrans));
    end
    gammar(k) = norm(g(:))^2/norm(Ag(:))^2;
    unext = s(:,:,:,k) - gammar(k)*g;

    uprev = u;
    u = unext;


end
uhat = unext;

uf = uhat.*repmat(f,[1,1,numTransmitters]);
z = H_forward_full(uf,sensorGreensFunctionSet,receiverMaskSet,dx,dy);

outs.z = z;
outs.uhat = uhat;
outs.s = s;
outs.gammar = gammar;
outs.mu = mu;
outs.K = K;
end

function grad = seagle_backprop(f,outs,domainGreensFunctionSet,...
    sensorGreensFunctionSet,receiverMaskSet,K,y,uinc,dx,dy)

%%% Scattering operators
A = @(u) A_forward(u, f, domainGreensFunctionSet, dx, dy);
AT = @(z) A_adjoint(z, f, domainGreensFunctionSet, dx, dy);
GT = @(z) G_adjoint_full(z,domainGreensFunctionSet, dx,dy);

z = outs.z;
uhat = outs.uhat;
gammar = outs.gammar;
s = outs.s;
mu = outs.mu;

[Ny,Nx,numTransmitters] = size(uinc);
HTres = H_adjoint_full(z-y,sensorGreensFunctionSet,receiverMaskSet,dx,dy);

qnext = zeros(Ny,Nx,numTransmitters); % q^{K+1}
q = repmat(conj(f),[1,1,numTransmitters]).*HTres; % q^K
r = conj(uhat).*HTres;

Sq = zeros(size(q));
Sqnext = zeros(size(qnext));

for k = K:-1:1
    if k == K
        for indTrans = 1:numTransmitters
            Sq(:,:,indTrans) = q(:,:,indTrans) - gammar(k)*AT(A(q(:,:,indTrans)));
        end
        qprev = (1-mu(k))*Sq;
    else
        for indTrans = 1:numTransmitters
            Sqnext(:,:,indTrans) = qnext(:,:,indTrans) - gammar(k+1)*AT(A(qnext(:,:,indTrans)));
        end
        qprev = (1-mu(k))*Sq + mu(k+1)*Sqnext;
    end

    As = zeros(Ny,Nx,numTransmitters);
    Aq = zeros(Ny,Nx,numTransmitters);
    for indTrans = 1:numTransmitters
        As(:,:,indTrans) = A(s(:,:,indTrans,k));
        Aq(:,:,indTrans) = A(q(:,:,indTrans));
    end
    r = r + gammar(k)*(q.*conj(GT(As - uinc)) + conj(s(:,:,:,k)).*GT(Aq));

    qnext = q;
    q = qprev;
    Sq = Sqnext;
end
grad = real(sum(r,3));


end
