function [rOverR,r,preTwist,chord,percT,t,pitchAxis] =...
    interpolateBladeData(dataFile,numSects,interpType)
% This program takes a tabluated blade design
% and interpolates between it to create more airfoil sections
This is useful for creating 3d blades using cad packages.

% import the data
% this code assumes that the data is formatted in columns of
% r/R, r, preTwist, chord, % thickness, thickness, pitch axis location
data = importdata(dataFile);

% name the data
harpOpt.r = data(:,2);
harpOpt.rOverR = data(:,1);
raidus = harpOpt.r(1)/harpOpt.rOverR(1);
harpOpt.preTwist = data(:,3);
harpOpt.chord = data(:,4);
harpOpt.percT = data(:,5);
harpOpt.t = data(:,6);
harpOpt.pitchAxis = data(:,7);
harpOpt.length = raidus - harpOpt.r(1);
harpOptEnd = harpOpt.r(end);
r = harpOpt.r(1):harpOpt.length/numSects:harpOpt.r(end);

% interpolated between the airofil profiles data
rOverR = r/raidus;
preTwist    = interp1(harpOpt.r,harpOpt.preTwist,r,interpType);
chord       = interp1(harpOpt.r,harpOpt.chord,r,interpType);
percT       = interp1(harpOpt.r,harpOpt.percT,r,interpType);
t           = interp1(harpOpt.r,harpOpt.t,r,interpType);
pitchAxis   = interp1(harpOpt.r,harpOpt.pitchAxis,r,interpType);
length = harpOpt.length;
