function [Output] = Struct_Analysis(RElm,TWIST,CHORD,PERCENT_THICKNESS,AoA,AFang,LocVel,Cl,Cd,Thrust)

global Normalized_AF_Coordinates HubRad RotorRad...
       NumCases rho NumSeg NumAFPoints...
       Einput allowableStrain matDensity DecreasingST STmin STdel;
%=========================================================================%   
% DUMMY VALUES USED FOR DEBUGGING ONLY
% load Coordinates.mat %Load the normalized airfoil coordinates FFA-W3-XXX family
% global Type
% Type = 2;
% % %GUI User Inputs
% Einput = 27.6; %GPa bulk material elasticity
% allowableStrain = 0.003; %(m/m) i.e. 3000 microstrain
% matDensity = 1800; %(kg/m^3) bulk material density
% DecreasingST = 0; %if 1, shell thickness will be monotonicaly decreasing
% STmin = 0.001; %minimum shell thickness (m)
% STdel = 0.00025; %shell thickness increment (m)
% SFstruct = 1.5; %safety factor multiplied to bending moments
% 
% %-------------------------------------------------------------------------%
% %data from circle root blade
% CircularRoot = 1;
% RootTranSt = 0.09;
% RootDir = pwd;
% HubRad = 0.19;
% RotorRad = 2.5;
% SpdSt = 0.2;
% SpdEnd = 3;
% SpdDel = 0.1;
% rho = 1025; 
% NumCases = round((SpdEnd - SpdSt)/SpdDel + 1);
% if isnan(NumCases); NumCases = 1; end
% NumSeg = 40;
% NumAFPoints = 40;
% RElm = [0.192000000000000;0.199000000000000;0.213000000000000;0.234000000000000;0.262000000000000;0.297000000000000;0.338000000000000;0.385000000000000;0.439000000000000;0.498000000000000;0.562000000000000;0.630000000000000;0.704000000000000;0.781000000000000;0.862000000000000;0.946000000000000;1.03200000000000;1.12000000000000;1.20900000000000;1.30000000000000;1.39000000000000;1.48100000000000;1.57000000000000;1.65800000000000;1.74400000000000;1.82800000000000;1.90900000000000;1.98600000000000;2.06000000000000;2.12800000000000;2.19200000000000;2.25100000000000;2.30500000000000;2.35200000000000;2.39300000000000;2.42800000000000;2.45600000000000;2.47700000000000;2.49100000000000;2.49800000000000;];
% TWIST = [22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;22.5600000000000;19.4600000000000;17.3200000000000;15.7000000000000;14.3800000000000;13.2800000000000;12.3200000000000;11.4500000000000;10.6600000000000;9.93000000000000;9.24000000000000;8.59000000000000;7.98000000000000;7.40000000000000;6.85000000000000;6.34000000000000;5.86000000000000;5.42000000000000;5.02000000000000;4.65000000000000;4.33000000000000;4.04000000000000;3.80000000000000;3.59000000000000;3.43000000000000;3.31000000000000;3.23000000000000;3.19000000000000;];
% CHORD = [0.270000000000000;0.270000000000000;0.270000000000000;0.273000000000000;0.286000000000000;0.318000000000000;0.370000000000000;0.439000000000000;0.518000000000000;0.598000000000000;0.669000000000000;0.719000000000000;0.743000000000000;0.743000000000000;0.723000000000000;0.692000000000000;0.654000000000000;0.615000000000000;0.577000000000000;0.540000000000000;0.506000000000000;0.475000000000000;0.447000000000000;0.420000000000000;0.396000000000000;0.374000000000000;0.352000000000000;0.333000000000000;0.314000000000000;0.296000000000000;0.280000000000000;0.265000000000000;0.251000000000000;0.238000000000000;0.227000000000000;0.217000000000000;0.210000000000000;0.204000000000000;0.200000000000000;0.198000000000000;];
% PERCENT_THICKNESS = [100;100;100;98.9000000000000;94.5000000000000;84.9000000000000;72.9000000000000;61.3000000000000;51.7000000000000;44.4000000000000;39.1000000000000;35.6000000000000;33.2000000000000;31.5000000000000;30;28.4000000000000;26.5000000000000;24.5000000000000;22.2000000000000;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;21;];
% fid = fopen([RootDir '\Output_Files\Circle_Root2\Circle_Root2_1.bed'],'rt');
% [BED_Error,BED,AoA,AFang,LocVel,Cl,Cd,CpMin,Thrust] = Read_BED(fid,NumCases,NumSeg,1);
% fclose(fid);
%=========================================================================%   

%Only need the profiles for the current percent thickness values
perT = cell2mat(Normalized_AF_Coordinates(:,2));
Scaled_AF_Coordinates = cell(NumSeg,2);
for n = 1:length(PERCENT_THICKNESS)
    Ind = FindInd(perT,PERCENT_THICKNESS(n)); 
    Scaled_AF_Coordinates(n,1) = {CHORD(n).*Normalized_AF_Coordinates{Ind,1}};
    Scaled_AF_Coordinates(n,2) = {perT(Ind)};
end


%Calculate the span (width) of each blade segment
EndPts(1,1) = HubRad;
EndPts(NumSeg+1,1) = RotorRad;
LocSpan = zeros(NumSeg,1);
for n = 2:NumSeg+1
EndPts(n,1) = RElm(n-1) + RElm(n-1)-EndPts(n-1);
LocSpan(n-1,1) = EndPts(n) - EndPts(n-1); %(m)
end

%From the BEM data, only use the Lift and Drag from the flow speed that corresponds to the
%maximum root flap bending moment
RootFlap = zeros(NumCases,1);
for n = 1:NumCases
T = Thrust(NumSeg*(n-1)+1:NumSeg*n);
RootFlap(n,1) = sum((RElm-HubRad).*LocSpan.*T./1000);
end
[MaxFlapMom N] = max(RootFlap);
N = N(1); %just incase duplicate maximums exist
AoA_n = AoA(NumSeg*(N-1)+1:NumSeg*N);
AFang_n = AFang(NumSeg*(N-1)+1:NumSeg*N);
LocVel_n = LocVel(NumSeg*(N-1)+1:NumSeg*N);
Cl_n = Cl(NumSeg*(N-1)+1:NumSeg*N);
Cd_n = Cd(NumSeg*(N-1)+1:NumSeg*N);

Area = LocSpan.*CHORD; %area of each blade segment (m^2)
Lift = 0.5.*rho.*(LocVel_n.^2).*Area.*Cl_n; %(N)
Drag = 0.5.*rho.*(LocVel_n.^2).*Area.*Cd_n; %(N), positive towards TE
Fn = Lift.*cosd(AoA_n) + Drag.*sind(AoA_n); %force perpendicular to chordline (N), positive towards upper surface
Ft = Drag.*cosd(AoA_n) - Lift.*sind(AoA_n); %force parallel to chordline (N), positive towards TE

%Local bending moments at the center of the elements, normal and tangent to the chordline,
Mnorm = zeros(NumSeg,1); %using right hand rule, positive when thumb points towards LE along chordline (upper in compression, lower in tension)
Mtan = zeros(NumSeg,1); %using right hand rule, positive when thumb points perpendicular to chordline and towards lower surface (LE in compression, TE in tension)
%Local bending moments at the center of the elements, normal and tangent to the plane of rotation
Mflap = zeros(NumSeg,1); %using right hand rule, positive when thumb points in direction of blade rotation along plane of rotation
Medge = zeros(NumSeg,1); %using right hand rule, positive when thumb points perpendicular to plane of rotation, finger curl in rotation of blade
for n = 2:NumSeg
    Twist_diff = abs((TWIST(n:end) - TWIST(n-1))); %units: radians
    r_arm = RElm(n:end) - RElm(n-1); %moment arm (m)
    Mnorm(n-1,1) = sum((Fn(n:end).*cosd(Twist_diff)...
                      + Ft(n:end).*sind(Twist_diff)).*r_arm)./1000; %(kN-m)
    Mtan(n-1,1)  = sum((Fn(n:end).*sind(Twist_diff)...
                      - Ft(n:end).*cosd(Twist_diff)).*r_arm)./1000; %(kN-m)
    Mflap(n-1,1) = sum((Lift(n:end).*sind(AFang_n(n:end))...
                      + Drag(n:end).*cosd(AFang_n(n:end))).*r_arm)./1000; %(kN-m)
    Medge(n-1,1) = sum((Lift(n:end).*cosd(AFang_n(n:end))...
                      - Drag(n:end).*sind(AFang_n(n:end))).*r_arm)./1000; %(kN-m)
end

% [Lift Drag sqrt(Lift.^2 + Drag.^2) Fn Ft sqrt(Fn.^2 + Ft.^2)]
% [Mnorm Mtan sqrt(Mnorm.^2 + Mtan.^2) Mflap Medge sqrt(Mflap.^2 + Medge.^2)]

%Preallocate array
ShellThickness(1:NumSeg,1) = STmin;

%pause;%__________________________________________________________________________________________

%determine the shell thickness at each radial station
for N = NumSeg:-1:1
    
    %determine STmin if we want the shell thickness to be monotonically decreasing
    if DecreasingST == 1 && N < NumSeg
       varSTmin = ShellThickness(N+1,1);
    else
       varSTmin = STmin;
    end
    
    %Maximum possible shell thickness, depends on how thick the current section is
    %either assume the max thickness based on the "nameplate" % thickness,
    %or actually measure it based on the x-y coordinates
    %     STmax = CHORD(N)*PERCENT_THICKNESS(N)/200;

    %NumAFPoints = 40; Changed to a global variable in
    %HARP_Opt.m __________TB_FEB13
    STmax = max(Scaled_AF_Coordinates{N}(1:NumAFPoints/2,2) - ...
                flipud(Scaled_AF_Coordinates{N}(NumAFPoints/2:end,2)))/2;
    
    % Case of scaled chord thickness < STmin
    if varSTmin > STmax
        varSTmin = STmax;
        disp('WARNING: chord thickness < minimum specified structural thickness');
    end
    
    ST = [varSTmin:STdel:STmax]';
    %%
    for M = 1:length(ST)
        
        %calculate the strain and shell properties
        [Strain I A] = getStrains(ST(M),STmax,CHORD(N),Mtan(N),Mnorm(N),Scaled_AF_Coordinates{N});
        maxStrain = max(abs([Strain.LE Strain.TE Strain.Upper Strain.Lower]));
        
        ShellThickness(N,1) = ST(M);
        StrainLE(N,1) = Strain.LE;
        StrainTE(N,1) = Strain.TE;
        StrainUpper(N,1) = Strain.Upper;
        StrainLower(N,1) = Strain.Lower;
        Iuu(N,1) = I.uu;
        Ivv(N,1) = I.vv;
        ShellArea(N,1) = A;
           
        %check if strain meets requirements
        if maxStrain(1) < allowableStrain;
           %this shell thickness is sufficient to meet strain requirements
           break %exits this for loop, and continues to next blade station
        end       
    end    
end
%%
%pause on;
%[StrainLE StrainTE StrainUpper StrainLower]

%Calculate stresses
StressLE = Einput.*StrainLE.*1000; %(MPa)
StressTE = Einput.*StrainTE.*1000; %(MPa)
StressUpper = Einput.*StrainUpper.*1000; %(MPa)
StressLower = Einput.*StrainLower.*1000; %(MPa)

%Calculate the blade mass distribution
%disp('ShellArea size');disp(size(ShellArea));
MassDist = ShellArea.*LocSpan.*matDensity;
BladeMass = sum(MassDist);

%store output variables in a structure
Output.TotalMass = BladeMass;
Output.MassUnitLen = MassDist./LocSpan;
Output.ShellThickness = ShellThickness.*1000; %converted to (mm)
Output.Iuu = Iuu;
Output.Ivv = Ivv;
Output.StrainLE = StrainLE;
Output.StrainTE = StrainTE;
Output.StrainUpper = StrainUpper;
Output.StrainLower = StrainLower;
Output.StressLE = StressLE;
Output.StressTE = StressTE;
Output.StressUpper = StressUpper;
Output.StressLower = StressLower;
Output.Mnorm = Mnorm;
Output.Mtan = Mtan;
