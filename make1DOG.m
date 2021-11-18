
function RF = make1DOG(pWorld,pRF)

%
%pWorld.n([1,2]))       size of image in pixels
%pWorld.visAngle        width of image in degrees of visual angle
%
%pRF.center             [x,y] center of DOG
%pRF.sig                1x2 vector of standard deviations of the two Gaussians
%pRF.contrast           1x2 vector of two amplitudes

%
%Creates a 3D receptive field (filter) that is a 2D difference of Gaussian (DOG)
%in space and a difference of Gamma functions in time.  
%
%Parameters:
%   pWorld.t            Time vector (see initWorld.m)
%   pWorld.n            1x2 vector containing image size
%   
%   pRF.center          Center of RF (0,0) is middle of image
%   pRF.sf              Spatial frequency of Gabor
%   pRF.sig             Standard deviation of Gabor
%   pRF.contrast        Amplitude of Gabor (max of 1)

%   pRF.n               1x2 vector contaning # of cascades for the positive
%                       (1st) and negative (2nd) Gamma function
%   pRF.theta           1x2 vector containing time constants for two Gammas
%   pRF.a               1x2 vector containing the amplitudes of two Gammas
%   pRF.nt              Number of time frames for the RF
%
%SEE ALSO:  Gamma.m  makeGaborRF.m

if ~isfield(pRF,'delay')
    pRF.delay = 0;
end

t = pWorld.t(1:(min(pRF.nt,length(pWorld.t))))-pRF.delay;

g1 = Gamma(pRF.n(1),pRF.theta(1),t);
g2 = Gamma(pRF.n(2),pRF.theta(2),t);

g = pRF.a(1)*g1 - pRF.a(2)*g2;

g(t<=0) = 0;
g = -g./norm(g);


RF = fliplr(g);
RF = reshape(RF,[pWorld.n(1),pWorld.n(2),pRF.nt]);
