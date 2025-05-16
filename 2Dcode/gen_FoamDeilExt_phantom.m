function phantom = gen_FoamDeilExt_phantom(XPix, YPix)
% gen_FoamDeilExt_phantom generates the FoamDeilExt phantom
% with a circle of radius 0.04 centered at the origin and contrast 2.0.
% A second circle of radius 0.0031 /2 is added with contrast 0.45 with center (-0.0555, 0.0)
%
% Inputs:
%   XPix - x-coordinates of the pixel grid. Has shape (n, n)
%   YPix - y-coordinates of the pixel grid. Has shape (n, n)
%
% Outputs:
%   phantom - n x n matrix representing the FoamDeilExt phantom

% Define the circle
radius = 0.04;
contrast = 0.45;

n = size(XPix, 1); % Assuming XPix and YPix are square matrices

% Create the phantom
phantom = zeros(n, n);
circleMask = (XPix.^2 + YPix.^2) <= radius^2;
phantom(circleMask) = contrast;


% Add the second circle
radius2 = 0.0155;
contrast2 = 2.0;
center2 = [-0.0555, 0.0];
circleMask2 = ((XPix - center2(1)).^2 + (YPix - center2(2)).^2) <= radius2^2;
phantom(circleMask2) = contrast2;
end