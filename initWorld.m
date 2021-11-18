function pWorld = initWorld(pWorld)
%pWorld = initWorld(pWorld)
%
%Adds parameters to the 'pWorld' structure for use in the vision demo
%programs.  Specifically it adds matrices x and y from a meshgrid
%containing visual angle coordinates, and a time vector t based on the
%parameters:
%
%   pWorld.n            1x3 vector containing image size and number of frames(h w t)
%   pWorld.visAngle     width of image in degrees of visual angle
%   pWorld.dur          duration of stimulus (seconds)

%%10/7/08 Written by G.M. Boynton at the University of Washington

ratio = pWorld.n(1)/pWorld.n(2);

[x,y] = meshgrid(linspace(-pWorld.visAngle/2,pWorld.visAngle/2,pWorld.n(2)),...
                (linspace(-ratio*pWorld.visAngle/2,ratio*pWorld.visAngle/2,pWorld.n(1))));
            
pWorld.x = x;
pWorld.y = y;
pWorld.t = linspace(0,pWorld.dur,pWorld.n(3));



