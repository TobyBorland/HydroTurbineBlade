function [blade] =...
    readHarpOptData(dataFile)
% This program takes a harpOptdesign developed by Harp _Opt and interpolates
% between it to create more airfoil sections than the origional Harp_Opt
% design. This is useful for creating 3d blades using cad packages.
% This code assumes that the data is formatted in columns of
% r/R, r, preTwist, chrod, % thickness, thickness, pitch axis

% import the data
dataFile
data = importdata(dataFile);

% name the data
blade.r = data(:,2);
blade.rOverR = data(:,1);
raidus = blade.r(1)/blade.rOverR(1);
blade.preTwist = data(:,3);
blade.chord = data(:,4);
blade.percT = data(:,5);
blade.t = data(:,6);
blade.pitchAxis = data(:,7);
blade.length = raidus - blade.r(1);

